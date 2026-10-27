#!/usr/bin/env python3
"""First-order electrostatic (de)stabilization of the extra Pt6 electrons by each species.

    python k_stab.py <label> [tags...]

δE_i = q_i * (-ΔV_i): q_i = GS Bader charge of atom i, ΔV_i = change of the electron
potential energy at atom i when Δn electrons sit on Pt6 (fixed geometry), averaged
over a sphere around the atom so its own polarization drops out. Summed per group
(K, OH, H2O, slab). Negative = this species makes the extra electrons on Pt cheaper.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from k_field import ROOT, electrostatic, sampler, shell_points

R_SHELL = 1.5


def main() -> None:
    label = sys.argv[1]
    tags = sys.argv[2:] or ["dn0p50", "dn1p00", "dn1p50", "dn2p00"]
    work = ROOT / "work" / label
    gs = json.loads((ROOT / "analysis/charges" / label / "gs.json").read_text())
    q = np.array(gs["bader"])
    groups = {k: v for k, v in gs["groups"].items() if k != "Pt6" and v}

    atoms, v0 = electrostatic(work / "cdft_dn0p00.gpw")
    shell = shell_points(R_SHELL)
    out = {"label": label, "r_shell_A": R_SHELL, "points": []}
    for tag in tags:
        gpw = work / f"cdft_{tag}.gpw"
        if not gpw.is_file():
            continue
        dn = int(tag[2]) + int(tag[4:]) / 100
        _, v = electrostatic(gpw)
        at = sampler(v - v0, atoms.cell.array, atoms.pbc)
        dv = np.array([np.mean(at(p + shell)) for p in atoms.positions])
        de = -q * dv
        row = {"tag": tag, "dn": dn,
               "groups_eV": {k: float(de[idx].sum()) for k, idx in groups.items()},
               "per_atom": {k: [{"i": i, "sym": atoms[i].symbol, "q": float(q[i]),
                                 "dV": float(dv[i]), "dE": float(de[i])} for i in idx]
                            for k, idx in groups.items() if k in ("K", "OH")}}
        out["points"].append(row)
        print(f"{label} {tag}: " + "  ".join(
            f"{k} {e:+.3f} eV ({e / dn:+.3f}/e)" for k, e in row["groups_eV"].items()), flush=True)
        for k in ("K", "OH"):
            for r in row["per_atom"].get(k, []):
                if r["sym"] in ("K", "O"):
                    print(f"    {k} {r['sym']}{r['i']}: q={r['q']:+.2f} ΔV={r['dV']:+.3f} δE={r['dE']:+.3f}",
                          flush=True)
    dest = ROOT / "analysis" / "k_stab"
    dest.mkdir(parents=True, exist_ok=True)
    (dest / f"{label}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
