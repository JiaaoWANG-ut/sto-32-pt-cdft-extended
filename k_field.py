#!/usr/bin/env python3
"""Electrostatic force on K+ when Δn electrons are put on Pt6 (fixed geometry).

    python k_field.py <label> [tags...]      (serial; default tags dn0p50..dn2p00)

F_K = q_K * ∇ΔV, ΔV = V_es(dn) - V_es(dn=0) from the saved fixed-geometry CDFT
densities (Poisson of pseudo density + compensation charges, same dipole-layer
solver). This is partition-free: no λw potential, no weight-derivative term.
The field is averaged over spherical shells around K, so K's own polarization
cancels to first order. Positive F_toward_Pt = K pulled toward Pt6.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from ase.units import Bohr, Ha
from scipy.ndimage import map_coordinates

from gpaw.old.calculator import GPAW

ROOT = Path(__file__).resolve().parent
SHELLS = (1.8, 2.2, 2.6)
DELTA = 0.05


def electrostatic(gpw: Path):
    calc = GPAW(str(gpw), txt=None)
    calc.initialize_positions()
    dens, ham = calc.density, calc.hamiltonian
    dens.interpolate_pseudo_density()
    dens.calculate_pseudo_charge()
    v = ham.finegd.zeros()
    ham.poisson.solve(v, dens.rhot_g, charge=-dens.charge)
    v = ham.finegd.zero_pad(ham.finegd.collect(v, broadcast=True))
    return calc.get_atoms(), v * Ha


def sampler(v, cell, pbc):
    inv = np.linalg.inv(cell)
    n = np.array(v.shape)

    def at(points):
        s = np.asarray(points) @ inv * n
        return map_coordinates(v, s.T, order=1, mode="grid-wrap")

    return at


def shell_points(r, n=146):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5**0.5) * i
    return r * np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]


def grad_avg(at, center, r):
    pts = center + shell_points(r) if r > 0 else center[None, :]
    g = np.zeros(3)
    for k in range(3):
        e = np.zeros(3)
        e[k] = DELTA
        g[k] = np.mean((at(pts + e) - at(pts - e)) / (2 * DELTA))
    return g


def main() -> None:
    label = sys.argv[1]
    tags = sys.argv[2:] or ["dn0p50", "dn1p00", "dn1p50", "dn2p00"]
    work = ROOT / "work" / label
    groups = json.loads((ROOT / "analysis/charges" / label / "gs.json").read_text())["groups"]
    k_idx, pt = groups["K"], groups["Pt6"]

    atoms, v0 = electrostatic(work / "cdft_dn0p00.gpw")
    pos = atoms.positions
    pt_c = pos[pt].mean(0)
    out = {"label": label, "K": k_idx, "shells_A": SHELLS, "points": []}
    for tag in tags:
        gpw = work / f"cdft_{tag}.gpw"
        if not gpw.is_file():
            continue
        _, v = electrostatic(gpw)
        at = sampler(v - v0, atoms.cell.array, atoms.pbc)
        chg = ROOT / "analysis/charges" / label / f"cdft_{tag}.json"
        bader = json.loads(chg.read_text())["bader"] if chg.is_file() else None
        near_pt = float(np.mean(at(pt_c + shell_points(3.0))))
        far = float(np.mean(at(pos[np.argmin(pos[:, 2])] + shell_points(2.0))))
        row = {"tag": tag, "dV_near_Pt_eV": near_pt, "dV_bottom_eV": far, "K": []}
        for k in k_idx:
            q = bader[k] if bader else 0.86
            d = pos[pt] - pos[k]
            u = d[np.argmin(np.linalg.norm(d, axis=1))]
            u = u / np.linalg.norm(u)
            uc = (pt_c - pos[k]) / np.linalg.norm(pt_c - pos[k])
            per_shell = []
            for r in SHELLS:
                f = q * grad_avg(at, pos[k], r)
                per_shell.append({"r": r, "F": f.tolist(), "F_abs": float(np.linalg.norm(f)),
                                  "F_toward_nearest_Pt": float(f @ u),
                                  "F_toward_Pt_centroid": float(f @ uc)})
            row["K"].append({"index": k, "q_bader": q,
                             "d_K_Pt_min": float(np.linalg.norm(d, axis=1).min()),
                             "shells": per_shell})
        out["points"].append(row)
        print(f"{label} {tag}: ΔV near Pt {near_pt:+.3f} eV, bottom {far:+.3f} eV", flush=True)
        for kk in row["K"]:
            s = kk["shells"][1]
            print(f"   K{kk['index']} q={kk['q_bader']:+.2f} d(K-Pt)={kk['d_K_Pt_min']:.2f} Å  "
                  f"|F|={s['F_abs']:.3f} eV/Å  toward nearest Pt {s['F_toward_nearest_Pt']:+.3f}  "
                  f"(shells: " + " ".join(f"{x['F_toward_nearest_Pt']:+.3f}" for x in kk["shells"]) + ")",
                  flush=True)
    dest = ROOT / "analysis" / "k_field"
    dest.mkdir(parents=True, exist_ok=True)
    (dest / f"{label}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
