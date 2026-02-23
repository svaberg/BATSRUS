#!/usr/bin/env python3
"""Compare BATSRUS ZDI self-test point dumps against ZDIpy.

Supports:
- main startup dump (single case, all coefficients active)
- single-coefficient sweep dump (multiple '# case ...' sections)

Assumes BATSRUS point dump columns:
idx lon_deg lat_deg br_raw bphi_raw btheta_raw bmer_raw br_io bphi_io btheta_io bmer_io
"""
from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


def parse_case_file(path: Path) -> Dict[str, np.ndarray]:
    cases: Dict[str, List[List[float]]] = {}
    current = None
    for line in path.read_text().splitlines():
        if line.startswith('# case'):
            current = line[len('# case') :].strip()
            cases[current] = []
            continue
        if not line.strip() or line.startswith('#'):
            continue
        if current is None:
            continue
        try:
            cases[current].append([float(x) for x in line.split()])
        except ValueError as exc:
            raise ValueError(f"Could not parse data line in {path}: {line!r}") from exc
    return {k: np.array(v, dtype=float) for k, v in cases.items() if v}


def load_zdipy(zdipy_root: Path):
    sys.path.insert(0, str(zdipy_root))
    from core import magneticGeom  # type: ignore

    return magneticGeom


@dataclass
class ErrStats:
    max_abs: float
    mean_abs: float


def comp_stats(a: np.ndarray, b: np.ndarray) -> ErrStats:
    e = np.abs(a - b)
    return ErrStats(float(np.max(e)), float(np.mean(e)))


def compare_main_case(cases: Dict[str, np.ndarray], magneticGeom, coeff_file: Path, lmax: int | None = None):
    # Prefer the explicit loaded-coeff case if present, else first available.
    key = None
    for k in cases:
        if 'Loaded coefficients' in k:
            key = k
            break
    if key is None:
        key = sorted(cases)[0]
    arr = cases[key]
    lon_deg = arr[:, 1]
    lat_deg = arr[:, 2]
    lon = np.deg2rad(lon_deg)
    clat = np.deg2rad(90.0 - lat_deg)

    kwargs = {'verbose': 0}
    if lmax is not None:
        kwargs['lmax'] = lmax
    mg = magneticGeom.magSphHarmoicsFromFile(str(coeff_file), **kwargs)
    mg.initMagGeom(clat, lon)
    B = mg.getAllMagVectors()  # [Br, Bclat, Blon]

    br = B[0]
    btheta = B[1]
    bphi = B[2]
    bmer = -B[1]

    print(f"main_case: {key}")
    for label, col, ref in [
        ('Br_raw', 3, br),
        ('Bphi_raw', 4, bphi),
        ('Btheta_raw', 5, btheta),
        ('Bmer_raw', 6, bmer),
    ]:
        s = comp_stats(arr[:, col], ref)
        print(f"  {label:10s} max={s.max_abs:.6e} mean={s.mean_abs:.6e}")


def compare_single_coeff(cases: Dict[str, np.ndarray], magneticGeom, lmax: int):
    patt = re.compile(r'single_coeff\d+\s+(alpha|beta|gamma)\s+(re|im)\s+l=(\d+)\s+m=(\d+)')
    scases = {k: v for k, v in cases.items() if patt.match(k)}
    if not scases:
        print('single_coeff: no single-coefficient cases found')
        return

    first = scases[sorted(scases)[0]]
    lon = np.deg2rad(first[:, 1])
    clat = np.deg2rad(90.0 - first[:, 2])
    mg = magneticGeom.magSphHarmoics(lmax)
    mg.initMagGeom(clat, lon)
    lm_to_i = {(int(l), int(m)): i for i, (l, m) in enumerate(zip(mg.l, mg.m))}

    worst = {'Br': (None, -1.0), 'Bphi': (None, -1.0), 'Btheta': (None, -1.0), 'Bmer': (None, -1.0)}
    by_set: Dict[str, Dict[str, List[float]]] = {
        s: {c: [] for c in ['Br', 'Bphi', 'Btheta', 'Bmer']} for s in ['alpha', 'beta', 'gamma']
    }

    for key in sorted(scases):
        m = patt.match(key)
        assert m
        set_name, part, l_str, m_str = m.groups()
        l = int(l_str)
        mm = int(m_str)

        mg.alpha[:] = 0j
        mg.beta[:] = 0j
        mg.gamma[:] = 0j
        getattr(mg, set_name)[lm_to_i[(l, mm)]] = 1.0 + 0.0j if part == 're' else 0.0 + 1.0j
        B = mg.getAllMagVectors()
        refs = {'Br': B[0], 'Btheta': B[1], 'Bphi': B[2], 'Bmer': -B[1]}
        cols = {'Br': 3, 'Bphi': 4, 'Btheta': 5, 'Bmer': 6}

        arr = scases[key]
        for c in ['Br', 'Bphi', 'Btheta', 'Bmer']:
            emax = float(np.max(np.abs(arr[:, cols[c]] - refs[c])))
            by_set[set_name][c].append(emax)
            if emax > worst[c][1]:
                worst[c] = (key, emax)

    print(f"single_coeff cases: {len(scases)}")
    for set_name in ['alpha', 'beta', 'gamma']:
        print(f"  {set_name}:")
        for c in ['Br', 'Bphi', 'Btheta', 'Bmer']:
            vals = by_set[set_name][c]
            if vals:
                print(f"    {c:6s} max={max(vals):.6e} mean={sum(vals)/len(vals):.6e}")
    print('  worst cases:')
    for c in ['Br', 'Bphi', 'Btheta', 'Bmer']:
        print(f"    {c:6s} {worst[c][0]} -> {worst[c][1]:.6e}")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument('--zdipy-root', default='/Users/dagfev/Documents/starwinds/ZDIpy')
    p.add_argument('--coeff-file', default='/Users/dagfev/Documents/toupies-magnetograms-colin/coeff-Mel25-5.dat')
    p.add_argument('--points-file', default='zdi_selftest_points.out')
    p.add_argument('--singlecoeff-file', default='zdi_selftest_singlecoeff_points.out')
    p.add_argument('--singlecoeff-lmax', type=int, default=2)
    p.add_argument('--skip-singlecoeff', action='store_true')
    args = p.parse_args()

    magneticGeom = load_zdipy(Path(args.zdipy_root))
    point_cases = parse_case_file(Path(args.points_file))
    compare_main_case(point_cases, magneticGeom, Path(args.coeff_file))

    if not args.skip_singlecoeff:
        scases = parse_case_file(Path(args.singlecoeff_file))
        compare_single_coeff(scases, magneticGeom, args.singlecoeff_lmax)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
