#!/usr/bin/env python3
"""Make a shell-plot movie (br/bphi/btheta) from BATSRUS idl_ascii output.

This script reads `shl VAR idl_ascii` files (e.g. `SC/IO2/shl_var_*.out`) and
creates a 3-panel movie. By default it uses a fixed linear color range
(-25..25) for all panels and draws a black zero-contour in each panel. If the
shell grid changes between frames, later frames are resampled to the first
frame's lon-lat grid.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess
from typing import Dict, List, Sequence

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, SymLogNorm
import numpy as np

from bats_idl_ascii import ShellFrame, read_shell_series


DEFAULT_COMPONENTS = ["br", "bphi", "btheta"]  # customary ZDI panel order: radial, azimuthal, meridional
DEFAULT_CMAPS = {
    "br": "RdBu_r",
    "bphi": "RdBu_r",
    "btheta": "RdBu_r",
}

COMPONENT_LABELS = {
    "br": "Radial (Br)",
    "bphi": "Azimuthal (Bphi)",
    # BATSRUS movie output uses Btheta (co-latitude / polar component).
    # ZDI papers/plots often label the third panel as meridional/latitudinal,
    # where Blat = -Btheta.
    "btheta": "Meridional (=-Btheta in BATSRUS)",
    "zdibr": "ZDI target radial (Br)",
    "zdibphi": "ZDI target azimuthal (Bphi)",
    "zdibtheta": "ZDI target meridional (=-Btheta)",
    "zdibmer": "ZDI target meridional (Bmer)",
}


@dataclass
class ColorScale:
    vmin: float
    vmax: float
    linthresh: float
    unit: str


def _find_unit(frame: ShellFrame, name: str, fallback: str = "") -> str:
    if name in frame.var_names and len(frame.units) == len(frame.var_names):
        return frame.units[frame.var_names.index(name)]
    return fallback


def _compute_color_scale(frames: Sequence[ShellFrame], component: str, q: float, linfrac: float) -> ColorScale:
    vals = np.concatenate([np.ravel(f.data[component]) for f in frames])
    avals = np.abs(vals[np.isfinite(vals)])
    if avals.size == 0:
        return ColorScale(vmin=-1.0, vmax=1.0, linthresh=1e-3, unit="")
    vmax = float(np.quantile(avals, q))
    vmax = max(vmax, float(avals.max()) * 1e-6, 1e-12)
    nonzero = avals[avals > 0.0]
    if nonzero.size:
        linthresh = max(float(np.quantile(nonzero, 0.25)) * linfrac, vmax * 1e-4, 1e-12)
    else:
        linthresh = max(vmax * 1e-4, 1e-12)
    return ColorScale(vmin=-vmax, vmax=vmax, linthresh=linthresh, unit=_find_unit(frames[0], component))


def _fixed_linear_scale(frame: ShellFrame, component: str, vabs: float) -> ColorScale:
    vmax = float(abs(vabs))
    return ColorScale(vmin=-vmax, vmax=vmax, linthresh=max(vmax * 1e-3, 1e-12), unit=_find_unit(frame, component))


def _mismatch_text(frame: ShellFrame) -> str:
    pieces = []
    for comp, target in [("br", "zdibr"), ("bphi", "zdibphi"), ("btheta", "zdibtheta")]:
        if comp in frame.data and target in frame.data:
            diff = np.abs(frame.data[comp] - frame.data[target])
            pieces.append(f"mean|d{comp}|={diff.mean():.3g}")
    return ", ".join(pieces)


def _component_label(component: str) -> str:
    return COMPONENT_LABELS.get(component.lower(), component)


def _component_stats(frames: Sequence[ShellFrame], component: str) -> tuple[float, float, float]:
    vals = np.concatenate([np.ravel(f.data[component]) for f in frames])
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return 0.0, 0.0, 0.0
    vmin = float(np.min(vals))
    vmax = float(np.max(vals))
    amax = float(np.max(np.abs(vals)))
    return vmin, vmax, amax


def _plot_frame(
    frame: ShellFrame,
    components: Sequence[str],
    scales: Dict[str, ColorScale],
    out_png: Path,
    layout: str,
    title: str,
    annotate: bool,
    color_scale_mode: str,
    draw_zero_contour: bool,
) -> None:
    n = len(components)
    if layout == "1x3":
        nrows, ncols = 1, n
        figsize = (5.1 * n, 4.4)
    else:
        nrows, ncols = n, 1
        figsize = (6.4, 3.2 * n)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, constrained_layout=True)
    axes_arr = np.atleast_1d(axes).ravel()

    lon2d, lat2d = np.meshgrid(frame.lon, frame.lat)
    for ax, comp in zip(axes_arr, components):
        scale = scales[comp]
        field = frame.data[comp]
        if color_scale_mode == "linear":
            norm = Normalize(vmin=scale.vmin, vmax=scale.vmax)
            cbar_suffix = ""
        else:
            norm = SymLogNorm(linthresh=scale.linthresh, vmin=scale.vmin, vmax=scale.vmax, base=10)
            cbar_suffix = " (symlog)"
        mesh = ax.pcolormesh(lon2d, lat2d, field, shading="auto", cmap=DEFAULT_CMAPS.get(comp, "RdBu_r"), norm=norm)
        if draw_zero_contour and np.nanmin(field) < 0.0 < np.nanmax(field):
            ax.contour(lon2d, lat2d, field, levels=[0.0], colors="k", linewidths=0.7)
        ax.set_title(_component_label(comp))
        ax.set_xlabel("Longitude [deg]")
        ax.set_ylabel("Latitude [deg]")
        ax.set_xlim(float(frame.lon.min()), float(frame.lon.max()))
        ax.set_ylim(float(frame.lat.min()), float(frame.lat.max()))
        ax.set_aspect("auto")
        cbar = fig.colorbar(mesh, ax=ax, fraction=0.046, pad=0.03)
        unit = scale.unit or ""
        cbar.set_label(unit + cbar_suffix if unit else (cbar_suffix.strip() or ""))

    if annotate:
        notes = [f"step={frame.step}", f"time={frame.time:g}", frame.path.name]
        if frame.resampled:
            notes.append("resampled-to-first-grid")
        mismatch = _mismatch_text(frame)
        if mismatch:
            notes.append(mismatch)
        fig.suptitle(f"{title}\n" + " | ".join(notes), fontsize=11)
    elif title:
        fig.suptitle(title, fontsize=11)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def _run_ffmpeg(frame_dir: Path, output: Path, fps: int) -> None:
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg is None:
        raise RuntimeError("ffmpeg not found in PATH")
    cmd = [
        ffmpeg,
        "-y",
        "-framerate",
        str(fps),
        "-i",
        str(frame_dir / "frame_%05d.png"),
        "-vf",
        "pad=ceil(iw/2)*2:ceil(ih/2)*2",
        "-c:v",
        "libx264",
        "-pix_fmt",
        "yuv420p",
        str(output),
    ]
    subprocess.run(cmd, check=True)


def build_movie(args: argparse.Namespace) -> None:
    frames = read_shell_series(
        args.io2_dir,
        pattern=args.pattern,
        step_min=args.step_min,
        step_max=args.step_max,
        stride=args.stride,
        resample_to_first=True,
    )
    if not frames:
        raise SystemExit(f"No shell frames found in {args.io2_dir}")

    components = [c.strip().lower() for c in args.components.split(",") if c.strip()]
    for comp in components:
        if comp not in frames[0].data:
            raise SystemExit(f"Component '{comp}' not found in {frames[0].path}")

    if args.color_scale == "linear":
        scales = {comp: _fixed_linear_scale(frames[0], comp, args.linear_vabs) for comp in components}
    else:
        scales = {
            comp: _compute_color_scale(frames, comp, q=args.quantile, linfrac=args.linthresh_frac)
            for comp in components
        }

    print(f"Color scale mode: {args.color_scale}")
    if args.color_scale == "linear":
        print(f"Plot range for all panels: [{-abs(args.linear_vabs):g}, {abs(args.linear_vabs):g}]")
    for comp in components:
        vmin, vmax, amax = _component_stats(frames, comp)
        print(f"{comp}: min={vmin:.6g} max={vmax:.6g} absmax={amax:.6g}")

    frame_dir = Path(args.frame_dir) if args.frame_dir else Path(args.output).with_suffix("")
    frame_dir.mkdir(parents=True, exist_ok=True)

    for i, frame in enumerate(frames):
        out_png = frame_dir / f"frame_{i:05d}.png"
        _plot_frame(
            frame=frame,
            components=components,
            scales=scales,
            out_png=out_png,
            layout=args.layout,
            title=args.title,
            annotate=not args.no_annotations,
            color_scale_mode=args.color_scale,
            draw_zero_contour=not args.no_zero_contour,
        )

    if not args.frames_only:
        _run_ffmpeg(frame_dir, Path(args.output), args.fps)

    print(f"Wrote {len(frames)} frame(s) to {frame_dir}")
    if not args.frames_only:
        print(f"Wrote movie {args.output}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("io2_dir", help="Path to SC/IO2 directory containing shl_var_*.out")
    parser.add_argument("--output", required=True, help="Output movie path (.mp4)")
    parser.add_argument("--frame-dir", help="Directory for intermediate PNG frames (default: output stem)")
    parser.add_argument("--pattern", default="shl_var_*_n*.out", help="Glob pattern for shell files")
    parser.add_argument(
        "--components",
        default=",".join(DEFAULT_COMPONENTS),
        help=(
            "Comma-separated variables (default: br,bphi,btheta). "
            "Customary ZDI panel order is radial, azimuthal, meridional; "
            "BATSRUS stores the meridional/polar panel as btheta, so ZDI meridional = -btheta."
        ),
    )
    parser.add_argument("--layout", choices=["1x3", "3x1"], default="1x3")
    parser.add_argument("--title", default="BATSRUS shell magnetic field")
    parser.add_argument(
        "--color-scale",
        choices=["linear", "symlog"],
        default="linear",
        help="Color normalization mode (default: linear)",
    )
    parser.add_argument(
        "--linear-vabs",
        type=float,
        default=25.0,
        help="Fixed symmetric range for linear mode: [-V,+V] (default: 25)",
    )
    parser.add_argument("--fps", type=int, default=10)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--step-min", type=int)
    parser.add_argument("--step-max", type=int)
    parser.add_argument("--quantile", type=float, default=0.995, help="Abs-value quantile for symmetric color limits")
    parser.add_argument("--linthresh-frac", type=float, default=0.5, help="Factor applied to Q25(|data|) for symlog linthresh")
    parser.add_argument("--frames-only", action="store_true", help="Only write PNG frames; skip ffmpeg")
    parser.add_argument("--no-zero-contour", action="store_true", help="Disable black B=0 contour overlay")
    parser.add_argument("--no-annotations", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    build_movie(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
