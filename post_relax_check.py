#!/usr/bin/env python3
"""Re-validate factor tags on relaxed.vasp (run after GS opt, before CDFT)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

from ase.io import read

from build_models import FAR_PT_CUT, MODEL_SPECS, metrics

ROOT = Path(__file__).resolve().parent
SPECS = {s["label"]: s for s in MODEL_SPECS}


def check_label(label: str) -> dict:
    spec = SPECS[label]
    path = ROOT / "work" / label / "relaxed.vasp"
    if not path.is_file():
        return {"label": label, "ok": False, "msgs": ["missing relaxed.vasp"]}
    atoms = read(path)
    m = metrics(atoms)
    kd = [x["d_K_Pt_min"] for x in m["K"]]
    msgs = []
    if spec["K"] == "noK" and m["nK"] != 0:
        msgs.append("nK!=0")
    if spec["K"] == "Knear" and (not kd or max(kd) > 4.5):
        msgs.append(f"Knear drifted: {kd}")
    if spec["K"] == "Kfar" and (not kd or min(kd) < FAR_PT_CUT - 0.05):
        msgs.append(f"Kfar collapsed toward Pt: {kd}")
    if spec["OH"] == "noOH" and m["nOH_on_Pt"] != 0:
        msgs.append(f"OH returned to Pt: {m['nOH_on_Pt']}")
    if spec["OH"] == "OHonPt" and m["nOH_on_Pt"] < 1:
        msgs.append("lost OH@Pt")
    return {"label": label, "ok": not msgs, "msgs": msgs, "metrics": m}


def main() -> None:
    labels = sys.argv[1:] or list(SPECS)
    rows = [check_label(lab) for lab in labels]
    n_fail = sum(1 for r in rows if not r["ok"])
    for r in rows:
        st = "OK" if r["ok"] else "FAIL:" + ";".join(r["msgs"])
        m = r.get("metrics") or {}
        kd = ",".join(f"{x['d_K_Pt_min']:.2f}" for x in m.get("K", [])) or "-"
        print(f"{r['label']:22s} OH@Pt={m.get('nOH_on_Pt','?'):}  dK=[{kd}]  {st}")
    out = ROOT / "structures" / "post_relax_validation.json"
    out.write_text(json.dumps(rows, indent=2, default=str))
    print(f"failures={n_fail} → {out}")
    if n_fail:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
