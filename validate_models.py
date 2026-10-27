#!/usr/bin/env python3
"""Check modeling rationality for all POSCAR.init under work/."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from ase.io import read

from build_models import (
    FAR_PT_CUT,
    MODEL_SPECS,
    _idx,
    _nH,
    metrics,
    oh_on_pt,
)

ROOT = Path(__file__).resolve().parent
SRC31 = ROOT.parent / "31-Pt-DEDL-CDFT"
PARENT = {
    "pris-koh": "pristine-koh",
    "pris-water": "pristine-water",
    "Vo-koh": "Vo-interface-koh",
    "Vo-water": "Vo-interface-water",
}


def main() -> None:
    rows = []
    n_fail = 0
    print(f"{'label':22s} {'formula':28s} nK OH@Pt  dKmin  dKmax  gap  status")
    for spec in MODEL_SPECS:
        label = spec["label"]
        atoms = read(ROOT / "work" / label / "POSCAR.init")
        m = metrics(atoms)
        kd = [x["d_K_Pt_min"] for x in m["K"]]
        dmin = min(kd) if kd else float("nan")
        dmax = max(kd) if kd else float("nan")
        msgs = []
        if spec["K"] == "noK" and m["nK"] != 0:
            msgs.append("expected nK=0")
        if spec["K"] == "Knear" and (not kd or max(kd) > 4.5):
            msgs.append("Knear d>4.5")
        if spec["K"] == "Kfar" and (not kd or min(kd) < FAR_PT_CUT - 0.05):
            msgs.append(f"Kfar d<{FAR_PT_CUT}")
        if spec["OH"] == "noOH" and m["nOH_on_Pt"] != 0:
            msgs.append("still OH@Pt")
        if spec["OH"] == "OHonPt" and m["nOH_on_Pt"] < 1:
            msgs.append("missing OH@Pt")
        parent = read(SRC31 / "work" / PARENT[spec["parent"]] / "relax.traj")
        if atoms.get_chemical_formula() != parent.get_chemical_formula():
            msgs.append("stoichiometry changed")
        if any(_nH(atoms, oi) > 2 for oi in _idx(atoms, "O")):
            msgs.append("O with >2 H")
        if m["n_clash_lt_0.60"]:
            msgs.append("clash<0.60")
        gap = float(atoms.cell[2, 2] - atoms.positions[:, 2].max())
        if gap < 5.0:
            msgs.append(f"vacuum gap {gap:.1f}<5")
        # Vo: one fewer lattice O than pris partner
        if spec["defect"] == "Vo" and "O78" not in m["formula"] and "O80" not in m["formula"]:
            # Vo-koh family O78; Vo-water family O80
            if m["formula"].count("O") == 0:
                msgs.append("Vo O count odd")
        ok = not msgs
        if not ok:
            n_fail += 1
        status = "OK" if ok else "FAIL:" + ";".join(msgs)
        print(
            f"{label:22s} {m['formula']:28s} {m['nK']:2d} {m['nOH_on_Pt']:5d}  "
            f"{dmin:5.2f}  {dmax:5.2f}  {gap:4.1f}  {status}"
        )
        rows.append({"label": label, "ok": ok, "msgs": msgs, "metrics": m, "gap_A": gap})

    out = ROOT / "structures" / "validation.json"
    out.write_text(json.dumps(rows, indent=2))
    print(f"\n{len(rows)} models, failures={n_fail} → {out}")
    if n_fail:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
