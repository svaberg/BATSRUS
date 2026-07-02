#!/usr/bin/env python3
"""Compare ZDI-derived PFSS coefficients with a BATSRUS harmonics file.

This is a convention check for ModUserAwsomZdiB0. It verifies that the radial
ZDI alpha coefficients map to the same surface-Br spherical harmonics used by
BATSRUS' existing PFSS path.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


def zdi_to_pfss_conversion(l_degree: int, m_order: int) -> float:
    complex_to_real = 1.0 if m_order == 0 else math.sqrt(2.0)
    condon_shortley = 1.0 if m_order % 2 == 0 else -1.0
    return math.sqrt(2 * l_degree + 1) / (
        condon_shortley * complex_to_real * math.sqrt(4.0 * math.pi)
    )


def read_zdi_alpha(path: Path) -> tuple[dict[tuple[int, int], tuple[float, float]], int]:
    with path.open() as stream:
        stream.readline()
        n_coeff, _n_set, flag = [int(x) for x in stream.readline().split()[:3]]
        coeff: dict[tuple[int, int], tuple[float, float]] = {}
        for _ in range(n_coeff):
            l_text, m_text, re_text, im_text = stream.readline().split()
            l_degree = int(l_text)
            m_order = int(m_text)
            re_val = float(re_text)
            im_val = float(im_text)
            if flag == -3:
                im_val = -im_val
            coeff[(l_degree, m_order)] = (re_val, im_val)
    return coeff, flag


def read_batsrus_harmonics(path: Path) -> dict[tuple[int, int], tuple[float, float]]:
    coeff: dict[tuple[int, int], tuple[float, float]] = {}
    with path.open() as stream:
        for line in stream:
            parts = line.split()
            if len(parts) != 4:
                continue
            try:
                l_degree = int(parts[0])
                m_order = int(parts[1])
                g_val = float(parts[2])
                h_val = float(parts[3])
            except ValueError:
                continue
            coeff[(l_degree, m_order)] = (g_val, h_val)
    return coeff


def legendre_arrays(order: int, theta: float) -> list[list[float]]:
    sin_theta = max(math.sin(theta), 1.0e-10)
    cos_theta = math.cos(theta)

    max_int = max(order * order, 5 * order, 10)
    sqrt_i = [0.0] + [math.sqrt(i) for i in range(1, max_int + 1)]
    sqrt_ratio = [0.0] * (order + 1)
    sqrt_ratio[0] = 1.0
    for m_order in range(1, order + 1):
        sqrt_ratio[m_order] = (
            sqrt_ratio[m_order - 1]
            * sqrt_i[2 * m_order - 1]
            / sqrt_i[2 * m_order]
        )

    p_lm = [[0.0] * (order + 1) for _ in range(order + 1)]
    dp_lm = [[0.0] * (order + 1) for _ in range(order + 1)]

    sin_theta_m = 1.0
    sin_theta_m1 = 1.0
    for m_order in range(order + 1):
        if m_order == 0:
            coef1 = sqrt_i[2 * m_order + 1]
        else:
            coef1 = sqrt_i[2 * (2 * m_order + 1)]

        p_lm[m_order][m_order] = sqrt_ratio[m_order] * coef1 * sin_theta_m
        if m_order < order:
            p_lm[m_order + 1][m_order] = (
                p_lm[m_order][m_order] * sqrt_i[2 * m_order + 3] * cos_theta
            )

        dp_lm[m_order][m_order] = (
            sqrt_ratio[m_order] * coef1 * m_order * cos_theta * sin_theta_m1
        )
        if m_order < order:
            dp_lm[m_order + 1][m_order] = sqrt_i[2 * m_order + 3] * (
                cos_theta * dp_lm[m_order][m_order]
                - sin_theta * p_lm[m_order][m_order]
            )

        sin_theta_m1 = sin_theta_m
        sin_theta_m *= sin_theta

    for m_order in range(order - 1):
        for l_degree in range(m_order + 2, order + 1):
            coef1 = sqrt_i[2 * l_degree + 1] / sqrt_i[l_degree**2 - m_order**2]
            coef2 = sqrt_i[2 * l_degree - 1]
            coef3 = sqrt_i[(l_degree - 1) ** 2 - m_order**2] / sqrt_i[
                2 * l_degree - 3
            ]
            p_lm[l_degree][m_order] = coef1 * (
                coef2 * cos_theta * p_lm[l_degree - 1][m_order]
                - coef3 * p_lm[l_degree - 2][m_order]
            )
            dp_lm[l_degree][m_order] = coef1 * (
                coef2
                * (
                    cos_theta * dp_lm[l_degree - 1][m_order]
                    - sin_theta * p_lm[l_degree - 1][m_order]
                )
                - coef3 * dp_lm[l_degree - 2][m_order]
            )

    for m_order in range(order + 1):
        for l_degree in range(m_order, order + 1):
            coef1 = 1.0 / sqrt_i[2 * l_degree + 1]
            p_lm[l_degree][m_order] *= coef1
            dp_lm[l_degree][m_order] *= coef1

    return p_lm


def surface_br(coeff: dict[tuple[int, int], tuple[float, float]], lon: float, lat: float) -> float:
    order = max(l_degree for l_degree, _m_order in coeff)
    theta = math.pi / 2.0 - lat
    p_lm = legendre_arrays(order, theta)
    br_val = 0.0
    for m_order in range(order + 1):
        cos_m_phi = math.cos(m_order * lon)
        sin_m_phi = math.sin(m_order * lon)
        for l_degree in range(max(1, m_order), order + 1):
            g_val, h_val = coeff[(l_degree, m_order)]
            br_val += p_lm[l_degree][m_order] * (
                g_val * cos_m_phi + h_val * sin_m_phi
            )
    return br_val


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("zdi_coeff_file", type=Path)
    parser.add_argument("batsrus_harmonics_file", type=Path)
    args = parser.parse_args()

    zdi_coeff, flag = read_zdi_alpha(args.zdi_coeff_file)
    batsrus_coeff = read_batsrus_harmonics(args.batsrus_harmonics_file)

    zdi_pfss_coeff: dict[tuple[int, int], tuple[float, float]] = {}
    max_coeff_error = 0.0
    worst_coeff = None
    for key, (re_val, im_val) in zdi_coeff.items():
        l_degree, m_order = key
        factor = zdi_to_pfss_conversion(l_degree, m_order)
        g_val = factor * re_val
        h_val = -factor * im_val
        zdi_pfss_coeff[key] = (g_val, h_val)
        ref_g, ref_h = batsrus_coeff[key]
        error = max(abs(g_val - ref_g), abs(h_val - ref_h))
        if error > max_coeff_error:
            max_coeff_error = error
            worst_coeff = (key, g_val, h_val, ref_g, ref_h)

    max_br_error = 0.0
    for lon_deg in (0, 10, 60, 123, 180, 250, 359):
        for lat_deg in (-80, -30, 0, 45, 80):
            lon = math.radians(lon_deg)
            lat = math.radians(lat_deg)
            error = abs(
                surface_br(zdi_pfss_coeff, lon, lat)
                - surface_br(batsrus_coeff, lon, lat)
            )
            max_br_error = max(max_br_error, error)

    print(f"donati_flag {flag}")
    print(f"max_coeff_abs_error {max_coeff_error:.8e}")
    print(f"worst_coeff {worst_coeff}")
    print(f"max_sample_surface_br_abs_error {max_br_error:.8e}")


if __name__ == "__main__":
    main()
