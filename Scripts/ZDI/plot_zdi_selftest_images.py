#!/usr/bin/env python3
"""Render PNG images from ZDI startup self-test idl_ascii 2D dumps.

This uses the same plotting style as `make_shell_component_movie.py`:
- 3 panels
- fixed linear color range by default
- black zero contour
- annotations/title

It reads the BATSRUS/SWMF ASCII plot files produced by the `#ZDISELFTEST`
startup dump (e.g. `zdi_selftest_field_2d.out`, `zdi_selftest_singlecoeff_*.out`).
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from bats_idl_ascii import read_shell_idl_ascii
from make_shell_component_movie import _fixed_linear_scale, _plot_frame  # reuse exact frame style


DEFAULT_COMPONENTS = "zdibr,zdibphi,zdibmer"


def _resolve_component_name(frame, component: str) -> str | None:
    want = component.strip().lower()
    for name in frame.data.keys():
        if name.lower() == want:
            return name
    return None


def _norm_component_names(frame, components: Sequence[str]) -> list[str]:
    requested = [c.strip().lower() for c in components if c.strip()]
    names: list[str] = []
    missing: list[str] = []
    for c in requested:
        resolved = _resolve_component_name(frame, c)
        if resolved is None:
            missing.append(c)
        else:
            names.append(resolved)
    if missing:
        available = ", ".join(frame.data.keys())
        raise SystemExit(f"Missing component(s) {missing} in {frame.path}; available: {available}")
    return names


def _discover_files(inputs: Sequence[str]) -> list[Path]:
    paths: list[Path] = []
    for item in inputs:
        p = Path(item)
        if any(ch in item for ch in "*?[]"):
            paths.extend(sorted(Path().glob(item)))
        elif p.is_dir():
            paths.extend(sorted(p.glob("zdi_selftest*.out")))
        else:
            paths.append(p)
    # Deduplicate while preserving order
    out: list[Path] = []
    seen: set[Path] = set()
    for p in paths:
        rp = p.resolve()
        if rp in seen:
            continue
        seen.add(rp)
        out.append(p)
    return out


def render_images(args: argparse.Namespace) -> None:
    files = _discover_files(args.inputs)
    if not files:
        raise SystemExit("No input files found")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    count = 0
    for path in files:
        if path.suffix.lower() != ".out":
            continue
        try:
            frame = read_shell_idl_ascii(path)
        except Exception as exc:
            if args.skip_bad:
                print(f"skip {path}: {exc}")
                continue
            raise

        components = _norm_component_names(frame, args.components.split(","))
        scales = {c: _fixed_linear_scale(frame, c, args.linear_vabs) for c in components}
        out_png = out_dir / (path.stem + ".png")

        title = args.title if args.title else "ZDI startup self-test field"
        _plot_frame(
            frame=frame,
            components=components,
            scales=scales,
            out_png=out_png,
            layout=args.layout,
            title=title,
            annotate=not args.no_annotations,
            color_scale_mode="linear",
            draw_zero_contour=not args.no_zero_contour,
        )
        print(f"Wrote {out_png}")
        count += 1

    if count == 0:
        raise SystemExit("No images were rendered")
    print(f"Rendered {count} image(s) to {out_dir}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "inputs",
        nargs="+",
        help="Input .out file(s), glob(s), or directory/directories containing zdi_selftest*.out",
    )
    p.add_argument("--output-dir", required=True, help="Directory for rendered PNG files")
    p.add_argument(
        "--components",
        default=DEFAULT_COMPONENTS,
        help="Comma-separated components (default: zdibr,zdibphi,zdibmer)",
    )
    p.add_argument("--layout", choices=["1x3", "3x1"], default="1x3")
    p.add_argument("--linear-vabs", type=float, default=3.0, help="Fixed linear color range +/-vabs")
    p.add_argument("--title", help="Optional figure title prefix")
    p.add_argument("--no-zero-contour", action="store_true")
    p.add_argument("--no-annotations", action="store_true")
    p.add_argument("--skip-bad", action="store_true", help="Skip files that fail parsing")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    render_images(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
