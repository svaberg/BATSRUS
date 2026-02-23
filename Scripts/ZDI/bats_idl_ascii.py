#!/usr/bin/env python3
"""Read BATSRUS/SWMF shell #SAVEPLOT outputs in idl_ascii (.out) format.

This module currently targets regular lon-lat shell plots (e.g. `shl VAR idl_ascii`).
It also provides a simple bilinear resampler so movie tools can handle cases where
frames (or runs) are saved on different shell grids.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math
import re
from typing import Dict, Iterable, List, Sequence

import numpy as np


STEP_RE = re.compile(r"_n(\d+)\.")


@dataclass
class ShellFrame:
    path: Path
    units: List[str]
    step: int
    time: float
    ndim: int
    nparam: int
    nvar: int
    dims: Sequence[int]
    var_names: List[str]
    lon: np.ndarray  # shape (nLon,)
    lat: np.ndarray  # shape (nLat,)
    data: Dict[str, np.ndarray]  # each shape (nLat, nLon)
    resampled: bool = False

    @property
    def nlon(self) -> int:
        return int(self.dims[0])

    @property
    def nlat(self) -> int:
        return int(self.dims[1])


def _parse_float_rows(lines: Sequence[str], ncol: int) -> np.ndarray:
    rows = []
    for line in lines:
        s = line.strip()
        if not s:
            continue
        parts = s.split()
        if len(parts) != ncol:
            continue
        try:
            rows.append([float(x) for x in parts])
        except ValueError:
            continue
    if not rows:
        raise ValueError("No numeric rows found in idl_ascii file")
    return np.asarray(rows, dtype=float)


def _extract_step(path: Path, fallback_step: int) -> int:
    match = STEP_RE.search(path.name)
    if match:
        try:
            return int(match.group(1))
        except ValueError:
            pass
    return fallback_step


def read_shell_idl_ascii(path: str | Path) -> ShellFrame:
    """Read a BATSRUS shell `idl_ascii` file into structured arrays."""

    file_path = Path(path)
    lines = file_path.read_text().splitlines()
    if len(lines) < 5:
        raise ValueError(f"File too short for idl_ascii shell plot: {file_path}")

    units = lines[0].split()
    meta_parts = lines[1].split()
    if len(meta_parts) < 5:
        raise ValueError(f"Malformed metadata line in {file_path}: {lines[1]!r}")

    step_header = int(meta_parts[0])
    time = float(meta_parts[1])
    ndim = int(meta_parts[2])
    nparam = int(meta_parts[3])
    nvar = int(meta_parts[4])
    if ndim != 2:
        raise ValueError(f"Only 2D shell plots are supported, got ndim={ndim} in {file_path}")

    dims = tuple(int(x) for x in lines[2].split())
    if len(dims) != ndim:
        raise ValueError(f"Expected {ndim} dims, got {len(dims)} in {file_path}")

    var_names = lines[3].split()
    expected_cols = ndim + nvar
    if len(var_names) != expected_cols:
        raise ValueError(
            f"Expected {expected_cols} column names (coords + vars), got {len(var_names)} in {file_path}"
        )

    if nparam != 0:
        # Shell VAR plots used here are nParam=0. Raise loudly if a different flavor appears.
        raise ValueError(
            f"Reader currently expects nParam=0 for shell plots, got nParam={nparam} in {file_path}"
        )

    arr = _parse_float_rows(lines[4:], expected_cols)
    n_expected = int(np.prod(dims))
    if arr.shape[0] != n_expected:
        raise ValueError(
            f"Data row count mismatch in {file_path}: got {arr.shape[0]}, expected {n_expected}"
        )

    lon_vals = np.unique(arr[:, 0])
    lat_vals = np.unique(arr[:, 1])
    if lon_vals.size != dims[0] or lat_vals.size != dims[1]:
        raise ValueError(
            f"Unexpected grid uniqueness in {file_path}: unique lon/lat = "
            f"{lon_vals.size}/{lat_vals.size}, dims={dims}"
        )

    lon_vals = np.asarray(lon_vals, dtype=float)
    lat_vals = np.asarray(lat_vals, dtype=float)

    lon_index = {float(val): i for i, val in enumerate(lon_vals)}
    lat_index = {float(val): j for j, val in enumerate(lat_vals)}

    data = {
        name: np.full((lat_vals.size, lon_vals.size), np.nan, dtype=float)
        for name in var_names
    }

    for row in arr:
        i_lon = lon_index[float(row[0])]
        j_lat = lat_index[float(row[1])]
        for i_col, name in enumerate(var_names):
            data[name][j_lat, i_lon] = row[i_col]

    for name, grid in data.items():
        if np.isnan(grid).any():
            raise ValueError(f"Incomplete grid for variable {name} in {file_path}")

    # Keep only physical columns in data map? No: retaining lon/lat grids is useful for debugging.
    step = _extract_step(file_path, step_header)
    return ShellFrame(
        path=file_path,
        units=units,
        step=step,
        time=time,
        ndim=ndim,
        nparam=nparam,
        nvar=nvar,
        dims=dims,
        var_names=var_names,
        lon=lon_vals,
        lat=lat_vals,
        data=data,
    )


def discover_shell_files(io2_dir: str | Path, pattern: str = "shl_var_*.out") -> List[Path]:
    paths = sorted(Path(io2_dir).glob(pattern))
    def sort_key(p: Path):
        m = STEP_RE.search(p.name)
        return (int(m.group(1)) if m else 10**18, p.name)
    return sorted(paths, key=sort_key)


def _drop_duplicate_lon_endpoint(lon: np.ndarray, field: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if lon.size < 2:
        return lon, field
    span = lon[-1] - lon[0]
    if math.isclose(abs(span), 360.0, rel_tol=0.0, abs_tol=1e-6):
        step = lon[1] - lon[0]
        if not math.isclose(step, 0.0, abs_tol=1e-12):
            return lon[:-1], field[:, :-1]
    return lon, field


def _interp_periodic_lon(src_lon: np.ndarray, src_vals: np.ndarray, dst_lon: np.ndarray) -> np.ndarray:
    """Interpolate along longitude with periodic wrap (degrees)."""
    src_lon, src_vals = _drop_duplicate_lon_endpoint(src_lon, src_vals)
    if src_lon.ndim != 1:
        raise ValueError("src_lon must be 1D")
    if src_vals.shape[-1] != src_lon.size:
        raise ValueError("Longitude axis mismatch")
    if src_lon.size < 2:
        return np.repeat(src_vals[..., :1], dst_lon.size, axis=-1)

    period = 360.0
    base = src_lon[0]
    x = src_lon.astype(float)
    y = src_vals.astype(float)
    # Extend one point for wrap-around interpolation.
    x_ext = np.concatenate([x, [x[0] + period]])
    y_ext = np.concatenate([y, y[..., :1]], axis=-1)

    dst_mod = ((dst_lon - base) % period) + base
    out = np.empty((y.shape[0], dst_mod.size), dtype=float)
    for j in range(y.shape[0]):
        out[j, :] = np.interp(dst_mod, x_ext, y_ext[j, :])
    return out


def _interp_lat(src_lat: np.ndarray, src_vals: np.ndarray, dst_lat: np.ndarray) -> np.ndarray:
    if src_vals.shape[0] != src_lat.size:
        raise ValueError("Latitude axis mismatch")
    out = np.empty((dst_lat.size, src_vals.shape[1]), dtype=float)
    for i in range(src_vals.shape[1]):
        out[:, i] = np.interp(dst_lat, src_lat, src_vals[:, i])
    return out


def resample_field_to_grid(
    src_lon: np.ndarray,
    src_lat: np.ndarray,
    src_field: np.ndarray,
    dst_lon: np.ndarray,
    dst_lat: np.ndarray,
) -> np.ndarray:
    """Bilinearly resample a regular lon-lat field to another regular lon-lat grid."""
    if src_field.shape != (src_lat.size, src_lon.size):
        raise ValueError("src_field shape does not match src lat/lon")

    lon_interp = _interp_periodic_lon(src_lon, src_field, dst_lon)
    return _interp_lat(src_lat, lon_interp, dst_lat)


def resample_frame(frame: ShellFrame, ref_lon: np.ndarray, ref_lat: np.ndarray) -> ShellFrame:
    if np.array_equal(frame.lon, ref_lon) and np.array_equal(frame.lat, ref_lat):
        return frame
    if np.allclose(frame.lon, ref_lon) and np.allclose(frame.lat, ref_lat):
        return frame

    resampled_data: Dict[str, np.ndarray] = {}
    for name, field in frame.data.items():
        if field.shape == (frame.lat.size, frame.lon.size):
            resampled_data[name] = resample_field_to_grid(frame.lon, frame.lat, field, ref_lon, ref_lat)
        else:
            # Non-gridded extras are copied as-is (not expected for shell VAR plots).
            resampled_data[name] = field.copy()

    return ShellFrame(
        path=frame.path,
        units=list(frame.units),
        step=frame.step,
        time=frame.time,
        ndim=frame.ndim,
        nparam=frame.nparam,
        nvar=frame.nvar,
        dims=(len(ref_lon), len(ref_lat)),
        var_names=list(frame.var_names),
        lon=np.asarray(ref_lon, dtype=float),
        lat=np.asarray(ref_lat, dtype=float),
        data=resampled_data,
        resampled=True,
    )


def read_shell_series(
    io2_dir: str | Path,
    pattern: str = "shl_var_*.out",
    step_min: int | None = None,
    step_max: int | None = None,
    stride: int = 1,
    resample_to_first: bool = True,
) -> List[ShellFrame]:
    files = discover_shell_files(io2_dir, pattern=pattern)
    frames: List[ShellFrame] = []
    for path in files:
        frame = read_shell_idl_ascii(path)
        if step_min is not None and frame.step < step_min:
            continue
        if step_max is not None and frame.step > step_max:
            continue
        frames.append(frame)

    if stride > 1:
        frames = frames[::stride]

    if not frames:
        return []

    if resample_to_first:
        ref_lon = frames[0].lon
        ref_lat = frames[0].lat
        frames = [resample_frame(f, ref_lon, ref_lat) for f in frames]

    return frames


def _cli_summary(paths: Iterable[Path]) -> int:
    for path in paths:
        frame = read_shell_idl_ascii(path)
        print(
            f"{path}: step={frame.step} time={frame.time:g} "
            f"grid={frame.nlon}x{frame.nlat} vars={len(frame.var_names)}"
        )
        print(f"  vars: {' '.join(frame.var_names)}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", help="Shell idl_ascii file(s) to summarize")
    args = parser.parse_args(argv)
    return _cli_summary([Path(p) for p in args.paths])


if __name__ == "__main__":
    raise SystemExit(main())
