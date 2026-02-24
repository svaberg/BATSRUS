#!/usr/bin/env python3
"""Prepare AWSoM ZDI experiment run directories using a ZDI-derived PFSS harmonics file.

Creates three runs (free/clamp/nudge) from a template run directory by rewriting
PARAM.in to use the provided PFSS harmonics file for initialization while keeping
the direct ZDI boundary target enabled for the tangential experiment.

This is intentionally a lightweight text-rewrite utility for local science runs.
"""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path


MODE_TO_TYPE = {
    "free": "off",
    "clamp": "clamp",
    "nudge": "nudge",
}


def replace_labeled_line(lines: list[str], label: str, new_value: str) -> None:
    for i, line in enumerate(lines):
        if line.lstrip().startswith("!"):
            continue
        if line.rstrip().endswith(label):
            lines[i] = new_value
            return
    raise ValueError(f"Could not find active line with label {label!r}")


def set_first_stop_max_iteration(lines: list[str], max_iter: int) -> None:
    in_first_stop = False
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped == "#STOP":
            in_first_stop = True
            continue
        if not in_first_stop:
            continue
        if line.lstrip().startswith("!"):
            continue
        if "MaxIteration" in line:
            # Preserve trailing comments if present.
            suffix = ""
            if "MaxIteration" in line:
                idx = line.index("MaxIteration")
                suffix = line[idx:]
            lines[i] = f"{max_iter}\t\t\t{suffix}".rstrip()
            return
    raise ValueError("Could not find first active #STOP MaxIteration line")


def comment_line_if_active(lines: list[str], idx: int) -> None:
    if idx < 0 or idx >= len(lines):
        return
    if lines[idx].strip() and not lines[idx].lstrip().startswith("!"):
        lines[idx] = "!" + lines[idx]


def enable_harmonics_block(lines: list[str]) -> None:
    """Uncomment #HARMONICSFILE and HARMONICSGRID blocks if they were commented."""
    i = 0
    while i < len(lines):
        stripped = lines[i].strip()
        raw = lines[i]
        if stripped in ("!#HARMONICSFILE", "#HARMONICSFILE"):
            # command + value line + likely blank line
            for j in range(i, min(i + 3, len(lines))):
                if lines[j].startswith("!"):
                    lines[j] = lines[j][1:]
            i += 3
            continue
        if stripped in ("!HARMONICSGRID", "HARMONICSGRID"):
            # command + 7 lines + likely blank line = 9 lines in this PARAM.in
            for j in range(i, min(i + 9, len(lines))):
                if lines[j].startswith("!"):
                    lines[j] = lines[j][1:]
            i += 9
            continue
        i += 1


def rewrite_param(
    param_path: Path,
    harmonics_file: Path,
    zdi_coeff_file: Path | None,
    mode: str,
    quick: bool,
    zdi_scale: float | None,
) -> None:
    lines = param_path.read_text().splitlines()

    enable_harmonics_block(lines)

    replace_labeled_line(
        lines, "NameHarmonicsFile", f"{harmonics_file}\tNameHarmonicsFile"
    )
    replace_labeled_line(lines, "UseZdiBoundary", "T\t\t\tUseZdiBoundary")
    replace_labeled_line(
        lines, "TypeZdiBoundary", f"{MODE_TO_TYPE[mode]}\t\t\tTypeZdiBoundary"
    )
    if zdi_coeff_file is not None:
        replace_labeled_line(
            lines, "NameZdiCoeffFile", f"{zdi_coeff_file}\tNameZdiCoeffFile"
        )
    if zdi_scale is not None:
        replace_labeled_line(lines, "ZdiFieldScaleIo", f"{zdi_scale}\t\t\tZdiFieldScaleIo")

    if quick:
        replace_labeled_line(lines, "DoSaveInitial", "T\t\t\tDoSaveInitial")
        set_first_stop_max_iteration(lines, 1)

    param_path.write_text("\n".join(lines) + "\n")


def clear_outputs(run_dir: Path) -> None:
    for rel in ("SC/IO2", "GM/IO2"):
        d = run_dir / rel
        if not d.exists():
            continue
        for child in d.iterdir():
            if child.is_file() or child.is_symlink():
                child.unlink()
    for name in ("runlog", "field_2d.out", "new_field_2d.out"):
        f = run_dir / name
        if f.exists():
            f.unlink()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare free/clamp/nudge AWSoM runs with ZDI-derived PFSS initialization."
    )
    parser.add_argument(
        "--template-run",
        default="run_test_awsom_mel25_zdiBr_free",
        help="Template run directory to copy (default: %(default)s)",
    )
    parser.add_argument(
        "--output-prefix",
        default="run_test_awsom_pfsszdi",
        help="Output run prefix; mode suffixes _free/_clamp/_nudge are added",
    )
    parser.add_argument(
        "--harmonics-file",
        required=True,
        type=Path,
        help="PFSS harmonics file (BATSRUS #HARMONICSFILE format) derived from ZDI Br",
    )
    parser.add_argument(
        "--zdi-coeff-file",
        type=Path,
        default=None,
        help="Optional override for #ZDIMAGNETOGRAM NameZdiCoeffFile",
    )
    parser.add_argument(
        "--zdi-scale",
        type=float,
        default=None,
        help="Optional override for ZdiFieldScaleIo (e.g. 0.1 for stability tests)",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Set DoSaveInitial=T and first #STOP MaxIteration=1",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output run directories",
    )
    args = parser.parse_args()

    repo_root = Path.cwd()
    template_dir = (repo_root / args.template_run).resolve()
    harmonics_file = args.harmonics_file.resolve()
    zdi_coeff_file = args.zdi_coeff_file.resolve() if args.zdi_coeff_file else None

    if not template_dir.is_dir():
        raise SystemExit(f"Template run directory not found: {template_dir}")
    if not harmonics_file.exists():
        raise SystemExit(f"Harmonics file not found: {harmonics_file}")
    if zdi_coeff_file and not zdi_coeff_file.exists():
        raise SystemExit(f"ZDI coeff file not found: {zdi_coeff_file}")

    created: list[Path] = []
    for mode in ("free", "clamp", "nudge"):
        out_dir = (repo_root / f"{args.output_prefix}_{mode}").resolve()
        if out_dir.exists():
            if args.force:
                shutil.rmtree(out_dir)
            else:
                raise SystemExit(f"Output directory exists (use --force): {out_dir}")
        shutil.copytree(template_dir, out_dir, symlinks=True)
        rewrite_param(
            out_dir / "PARAM.in",
            harmonics_file=harmonics_file,
            zdi_coeff_file=zdi_coeff_file,
            mode=mode,
            quick=args.quick,
            zdi_scale=args.zdi_scale,
        )
        clear_outputs(out_dir)
        created.append(out_dir)

    print("Prepared runs:")
    for path in created:
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
