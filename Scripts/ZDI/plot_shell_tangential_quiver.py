#!/usr/bin/env python3
"""Plot shell tangential magnetic field arrows on a field-strength background.

Background:
- scalar field strength |B| = sqrt(Br^2 + Bphi^2 + Btheta^2)
- grayscale colormap to avoid clashing with arrow colors

Arrows:
- tangential field on the shell using (Bphi, Bmer), where Bmer = -Btheta
- arrow length scales with tangential field strength
- arrow color is split by sign of Br (two complementary colors)

This targets BATSRUS/SWMF shell `idl_ascii` outputs (`shl VAR idl_ascii`).
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np

from bats_idl_ascii import read_shell_idl_ascii


POS_COLOR = "#d95f02"  # orange (Br > 0)
NEG_COLOR = "#1b9e77"  # teal   (Br < 0)


def _get_unit(frame, var_name: str) -> str:
    names = [n.lower() for n in frame.var_names]
    try:
        i = names.index(var_name.lower())
    except ValueError:
        return ""
    if i < len(frame.units):
        return frame.units[i]
    return ""


def _build_quiver_vectors(
    lon2d: np.ndarray,
    lat2d: np.ndarray,
    bphi: np.ndarray,
    bmer: np.ndarray,
    stride_lon: int,
    stride_lat: int,
    lon_step_deg: float,
    lat_step_deg: float,
    arrow_len_frac: float,
    cos_floor: float,
    arrow_vref: float | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    """Return quiver sample points and vectors in plotted lon/lat degree axes.

    Horizontal axis is longitude [deg], vertical is latitude [deg].
    Convert the physical tangential components approximately into this plotted
    coordinate system by scaling the longitudinal component with sec(lat).
    """

    ii = slice(None, None, stride_lat)
    jj = slice(None, None, stride_lon)

    lon_q = lon2d[ii, jj]
    lat_q = lat2d[ii, jj]
    bphi_q = bphi[ii, jj]
    bmer_q = bmer[ii, jj]

    coslat = np.cos(np.deg2rad(lat_q))
    coslat = np.where(np.abs(coslat) < cos_floor, np.sign(coslat) * cos_floor, coslat)

    # Raw vector components in plotted (lon, lat) axes, before global scaling.
    u_raw = bphi_q / coslat
    v_raw = bmer_q

    mag_raw = np.sqrt(u_raw**2 + v_raw**2)
    finite = np.isfinite(mag_raw)
    if not np.any(finite):
        return lon_q, lat_q, np.zeros_like(u_raw), np.zeros_like(v_raw), 0.0

    # Scale so a reference field magnitude maps to a fixed fraction of the arrow spacing.
    spacing_deg = np.sqrt((stride_lon * lon_step_deg) * (stride_lat * lat_step_deg))
    target_len_deg = max(spacing_deg * arrow_len_frac, 1e-6)
    if arrow_vref is None:
        mag_ref = float(np.quantile(mag_raw[finite], 0.95))
    else:
        mag_ref = float(arrow_vref)

    if mag_ref <= 0.0:
        scale = 0.0
    else:
        scale = target_len_deg / mag_ref

    return lon_q, lat_q, u_raw * scale, v_raw * scale, scale


def plot_shell_quiver(args: argparse.Namespace) -> Path:
    frame = read_shell_idl_ascii(args.shell_file)

    # Case-insensitive lookup but preserve original keys
    key_map = {k.lower(): k for k in frame.data.keys()}
    for need in ("br", "bphi", "btheta"):
        if need not in key_map:
            raise SystemExit(f"Missing required variable '{need}' in {frame.path}")

    br = frame.data[key_map["br"]]
    bphi = frame.data[key_map["bphi"]]
    btheta = frame.data[key_map["btheta"]]
    bmer = -btheta
    bmag = np.sqrt(br**2 + bphi**2 + btheta**2)

    lon2d, lat2d = np.meshgrid(frame.lon, frame.lat)
    lon_step = float(np.median(np.diff(frame.lon))) if frame.nlon > 1 else 1.0
    lat_step = float(np.median(np.diff(frame.lat))) if frame.nlat > 1 else 1.0

    vabs = args.vmax if args.vmax is not None else float(np.quantile(bmag[np.isfinite(bmag)], args.quantile))
    vabs = max(vabs, 1e-12)

    fig, ax = plt.subplots(figsize=(13.5, 6.2), constrained_layout=True)
    mesh = ax.pcolormesh(
        lon2d, lat2d, bmag,
        shading="auto",
        cmap=args.cmap,
        norm=Normalize(vmin=0.0, vmax=vabs),
    )

    if args.draw_br_zero and np.nanmin(br) < 0.0 < np.nanmax(br):
        ax.contour(lon2d, lat2d, br, levels=[0.0], colors="k", linewidths=0.8, alpha=0.7)

    lon_q, lat_q, u_q, v_q, field_to_plot_scale = _build_quiver_vectors(
        lon2d, lat2d, bphi, bmer,
        stride_lon=args.stride_lon,
        stride_lat=args.stride_lat,
        lon_step_deg=lon_step,
        lat_step_deg=lat_step,
        arrow_len_frac=args.arrow_len_frac,
        cos_floor=args.cos_floor,
        arrow_vref=args.arrow_vref,
    )

    br_q = br[::args.stride_lat, ::args.stride_lon]
    finite_vec = np.isfinite(u_q) & np.isfinite(v_q) & np.isfinite(br_q)
    pos = finite_vec & (br_q > 0.0)
    neg = finite_vec & (br_q < 0.0)

    qkw = dict(
        angles="xy",
        scale_units="xy",
        scale=1.0,
        width=args.quiver_width,
        headwidth=args.headwidth,
        headlength=args.headlength,
        headaxislength=args.headaxislength,
        pivot="mid",
        alpha=args.arrow_alpha,
        zorder=4,
    )

    q_obj = None
    if np.any(pos):
        q_obj = ax.quiver(lon_q[pos], lat_q[pos], u_q[pos], v_q[pos], color=POS_COLOR, **qkw)
    if np.any(neg):
        q_neg = ax.quiver(lon_q[neg], lat_q[neg], u_q[neg], v_q[neg], color=NEG_COLOR, **qkw)
        if q_obj is None:
            q_obj = q_neg

    # Tiny legend proxies
    ax.plot([], [], color=POS_COLOR, lw=3, label="Arrow color: Br > 0")
    ax.plot([], [], color=NEG_COLOR, lw=3, label="Arrow color: Br < 0")
    ax.legend(loc="upper right", framealpha=0.9)

    if args.arrow_key_value is not None and q_obj is not None and field_to_plot_scale > 0.0:
        unit = _get_unit(frame, key_map["br"])
        label_unit = f" {unit}" if unit else ""
        ax.quiverkey(
            q_obj,
            X=args.arrow_key_x,
            Y=args.arrow_key_y,
            U=float(args.arrow_key_value) * field_to_plot_scale,
            label=f"{args.arrow_key_value:g}{label_unit} tangential",
            labelpos="E",
            coordinates="axes",
            color="k",
            fontproperties={"size": 9},
        )

    ax.set_title(args.title if args.title else "Shell tangential field vectors on |B| background")
    ax.set_xlabel("Longitude [deg]")
    ax.set_ylabel("Latitude [deg]")
    ax.set_xlim(float(frame.lon.min()), float(frame.lon.max()))
    ax.set_ylim(float(frame.lat.min()), float(frame.lat.max()))

    cbar = fig.colorbar(mesh, ax=ax, fraction=0.03, pad=0.02)
    unit = _get_unit(frame, key_map["br"])
    cbar.set_label(f"|B| [{unit}]" if unit else "|B|")

    if args.annotate:
        txt = (
            f"step={frame.step} | time={frame.time:g} | {frame.path.name}\n"
            f"arrows=(Bphi, Bmer), Bmer=-Btheta | stride=({args.stride_lon},{args.stride_lat})"
        )
        ax.text(
            0.01, 0.01, txt,
            transform=ax.transAxes,
            va="bottom", ha="left",
            fontsize=9,
            color="k",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.7, edgecolor="none"),
        )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=args.dpi)
    plt.close(fig)
    return out


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("shell_file", help="Path to a shell idl_ascii file (shl_var_*.out)")
    p.add_argument("--output", required=True, help="Output image path (PNG/PDF)")
    p.add_argument("--title", help="Figure title")
    p.add_argument("--cmap", default="Greys", help="Background colormap (default: Greys)")
    p.add_argument("--quantile", type=float, default=0.995, help="Quantile for |B| vmax if --vmax not set")
    p.add_argument("--vmax", type=float, help="Fixed background vmax for |B| (vmin is 0)")
    p.add_argument("--stride-lon", type=int, default=20, help="Arrow sampling stride in longitude index")
    p.add_argument("--stride-lat", type=int, default=20, help="Arrow sampling stride in latitude index")
    p.add_argument("--arrow-len-frac", type=float, default=0.55, help="Arrow length vs sample spacing")
    p.add_argument(
        "--arrow-vref",
        type=float,
        help="Fixed tangential-field reference magnitude for arrow scaling (same units as B); "
             "if omitted, use per-frame 95th percentile",
    )
    p.add_argument(
        "--arrow-key-value",
        type=float,
        help="Draw quiver key for this tangential-field magnitude (same units as B)",
    )
    p.add_argument("--arrow-key-x", type=float, default=0.78, help="Quiver key X position in axes coordinates")
    p.add_argument("--arrow-key-y", type=float, default=1.02, help="Quiver key Y position in axes coordinates")
    p.add_argument("--cos-floor", type=float, default=0.25, help="Minimum |cos(lat)| for lon-axis conversion")
    p.add_argument("--quiver-width", type=float, default=0.0022)
    p.add_argument("--headwidth", type=float, default=3.4)
    p.add_argument("--headlength", type=float, default=4.6)
    p.add_argument("--headaxislength", type=float, default=4.2)
    p.add_argument("--arrow-alpha", type=float, default=0.95)
    p.add_argument("--dpi", type=int, default=180)
    p.add_argument("--draw-br-zero", action="store_true", help="Draw Br=0 contour (black)")
    p.add_argument("--annotate", action="store_true", help="Annotate file/step info")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    out = plot_shell_quiver(args)
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
