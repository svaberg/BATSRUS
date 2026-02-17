#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Convert BATSRUS .idl to .dat via pIDL and plot with pcolormesh"
    )
    p.add_argument("idl_file", help="Path to BATSRUS .idl file (typically in GM/IO2)")
    p.add_argument("--var", default=None, help="Variable name to plot (default: first physical variable)")
    p.add_argument("--out", default=None, help="Output PNG path (default: <dat>.png)")
    return p.parse_args()


def base_name_from_idl(idl: Path) -> str:
    # 3d__raw_..._pe0000.idl -> 3d__raw_...
    m = re.match(r"(.+)_pe\d+\.idl$", idl.name)
    if m:
        return m.group(1)
    if idl.suffix == ".idl":
        return idl.stem
    raise ValueError(f"Not an .idl file: {idl}")


def run_pidl_to_dat(idl: Path) -> Path:
    io2_dir = idl.parent
    if io2_dir.name != "IO2":
        raise ValueError(f"Expected file under IO2/, got {idl}")

    gm_dir = io2_dir.parent
    pidl = gm_dir / "pIDL"
    if not pidl.exists():
        raise FileNotFoundError(f"Could not find pIDL at {pidl}")

    base = base_name_from_idl(idl)
    cmd = [str(pidl), "-k", "-f=tec", f"IO2/{base}"]
    subprocess.run(cmd, cwd=gm_dir, check=True)

    dat = io2_dir / f"{base}.dat"
    if not dat.exists():
        raise FileNotFoundError(f"Expected converted dat file not found: {dat}")
    return dat


def parse_tecplot_dat(dat: Path):
    lines = dat.read_text().splitlines()

    var_line = next((ln for ln in lines if ln.strip().startswith("VARIABLES=")), None)
    if var_line is None:
        raise ValueError("VARIABLES line not found")
    vars_ = re.findall(r'"([^"]+)"', var_line)
    vars_clean = [v.strip() for v in vars_]

    zone_line = next((ln for ln in lines if ln.strip().startswith("ZONE")), None)
    if zone_line is None:
        raise ValueError("ZONE line not found")

    mi = re.search(r"\bI=\s*(\d+)", zone_line)
    mj = re.search(r"\bJ=\s*(\d+)", zone_line)
    if not mi or not mj:
        raise ValueError(f"Could not parse I/J from ZONE line: {zone_line}")
    ni = int(mi.group(1))
    nj = int(mj.group(1))

    data_start = None
    for idx, ln in enumerate(lines):
        s = ln.lstrip()
        if s and (s[0].isdigit() or s[0] in "+-"):
            data_start = idx
            break
    if data_start is None:
        raise ValueError("No numeric data lines found")

    npts = ni * nj
    arr = np.loadtxt(dat, skiprows=data_start, max_rows=npts)
    if arr.ndim == 1:
        arr = arr[None, :]
    if arr.shape[0] == 0:
        raise ValueError("Data zone is empty (0 points)")

    lower_to_idx = {v.lower(): i for i, v in enumerate(vars_clean)}

    def pick_col(names: list[str]) -> int:
        for n in names:
            if n.lower() in lower_to_idx:
                return lower_to_idx[n.lower()]
        raise KeyError(f"Could not find any of {names} in variables: {vars_clean}")

    ix = pick_col(["x", "X"])
    iy = pick_col(["y", "Y"])

    return arr, vars_clean, ni, nj, ix, iy


def pick_plot_var(vars_clean: list[str], requested: str | None) -> str:
    if requested is not None:
        return requested
    skip = {"j", "i", "x", "y", "z"}
    for v in vars_clean:
        if v.strip().lower() not in skip:
            return v
    raise ValueError("No plottable variable found")


def main() -> None:
    args = parse_args()
    idl = Path(args.idl_file).resolve()

    dat = run_pidl_to_dat(idl)
    arr, vars_clean, ni, nj, ix, iy = parse_tecplot_dat(dat)

    var_name = pick_plot_var(vars_clean, args.var)
    try:
        ivar = [v.lower() for v in vars_clean].index(var_name.lower())
    except ValueError as e:
        raise ValueError(f"Variable '{var_name}' not found in {vars_clean}") from e

    expected = ni * nj
    if arr.shape[0] < expected:
        raise ValueError(
            f"Not enough points for structured grid: got {arr.shape[0]}, expected {expected}."
        )

    # File layout is point format with I varying fastest, then J.
    x2d = arr[:expected, ix].reshape(nj, ni)
    y2d = arr[:expected, iy].reshape(nj, ni)
    v2d = arr[:expected, ivar].reshape(nj, ni)

    fig, ax = plt.subplots(figsize=(9, 3.5))
    m = ax.pcolormesh(x2d, y2d, v2d, shading="auto")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"{var_name} from {idl.name}")
    fig.colorbar(m, ax=ax, label=var_name)
    fig.tight_layout()

    out = Path(args.out).resolve() if args.out else dat.with_suffix(".png")
    fig.savefig(out, dpi=150)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
