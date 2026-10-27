#!/usr/bin/env python3
"""Build the Tier A+B structure matrix for 32-Pt-DEDL-CDFT-extended.

Parents (relaxed GS from 31-Pt-DEDL-CDFT):
  pris / Vo  ×  {koh (2 K+, OH@Pt), water (no K, OH@Pt)}

Transforms (charge- and stoichiometry-conserving):
  noOH   : proton-transfer — take H from a distant H2O onto each OH@Pt → H2O@Pt + far OH-
  Kfar   : swap each K with a far H2O molecule (O+2H), d(K–Pt) ≳ 6 Å
  Knear  : keep parent K positions (~3.1–3.4 Å)
  noK    : start from water parent (no K atoms)
  OHonPt : keep parent OH with d(O–Pt) < 2.5 Å and nH == 1

Writes:
  structures/<label>.vasp
  work/<label>/POSCAR.init   (Selective dynamics: bottom TiO2 fixed, z ≲ 4.04 Å)
  structures/manifest.json
"""
from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.io import read, write

ROOT = Path(__file__).resolve().parent
SRC31 = ROOT.parent / "31-Pt-DEDL-CDFT"

# OH@Pt: solvent O with exactly one H and min distance to any Pt below this.
OH_PT_CUT = 2.50
# Far solvent for proton source / K swap.
FAR_PT_CUT = 6.0
FAR_Z_LO, FAR_Z_HI = 11.5, 17.0
H_BOND = 1.30
FIXED_ZMAX = 4.05  # bottom TiO2 layer from 31 POSCAR.init


def _sym(atoms: Atoms):
    return np.array(atoms.get_chemical_symbols())


def _idx(atoms: Atoms, symbol: str) -> list[int]:
    return [i for i, s in enumerate(_sym(atoms)) if s == symbol]


def _pt_pos(atoms: Atoms) -> np.ndarray:
    return atoms.positions[_idx(atoms, "Pt")]


def _nH(atoms: Atoms, oi: int) -> int:
    """Count covalently bound H. No MIC: cell edges otherwise fake a second H."""
    h = _idx(atoms, "H")
    if not h:
        return 0
    return int((np.linalg.norm(atoms.positions[oi] - atoms.positions[h], axis=1) < H_BOND).sum())


def _H_of(atoms: Atoms, oi: int) -> list[int]:
    h = _idx(atoms, "H")
    d = np.linalg.norm(atoms.positions[oi] - atoms.positions[h], axis=1)
    return [int(h[i]) for i in range(len(h)) if d[i] < H_BOND]


def oh_on_pt(atoms: Atoms) -> list[int]:
    """Oxygen indices of OH directly coordinating Pt."""
    pt = _pt_pos(atoms)
    out = []
    for oi in _idx(atoms, "O"):
        if _nH(atoms, oi) != 1:
            continue
        d = float(np.linalg.norm(atoms.positions[oi] - pt, axis=1).min())
        if d < OH_PT_CUT:
            out.append(oi)
    return out


def far_waters(atoms: Atoms, exclude: set[int] | None = None) -> list[tuple[float, int]]:
    """(dmin_Pt, O_index) for H2O with d(Pt)>FAR and z in solvent window, farthest first."""
    exclude = exclude or set()
    pt = _pt_pos(atoms)
    cands = []
    for oi in _idx(atoms, "O"):
        if oi in exclude:
            continue
        if _nH(atoms, oi) != 2:
            continue
        z = atoms.positions[oi, 2]
        if not (FAR_Z_LO <= z <= FAR_Z_HI):
            continue
        d = float(np.linalg.norm(atoms.positions[oi] - pt, axis=1).min())
        if d >= FAR_PT_CUT:
            cands.append((d, oi))
    cands.sort(reverse=True)
    return cands


def apply_fix_bottom(atoms: Atoms) -> Atoms:
    fixed = [i for i, p in enumerate(atoms.positions) if p[2] < FIXED_ZMAX]
    atoms.set_constraint(FixAtoms(indices=fixed))
    return atoms


def protonate_oh_on_pt(atoms: Atoms) -> dict:
    """Convert every OH@Pt to H2O by taking one H from a distant H2O (→ far OH)."""
    atoms = atoms.copy()
    used_o: set[int] = set()
    transfers = []
    for oi in oh_on_pt(atoms):
        cands = far_waters(atoms, exclude=used_o | {oi})
        if not cands:
            raise RuntimeError(f"no far water to protonate OH O{oi}")
        _, far_o = cands[0]
        hs = _H_of(atoms, far_o)
        if len(hs) != 2:
            raise RuntimeError(f"far O{far_o} is not H2O (nH={len(hs)})")
        # Move the H that is farther from Pt onto the OH oxygen direction.
        pt = _pt_pos(atoms)
        hi = max(hs, key=lambda h: float(np.linalg.norm(atoms.positions[h] - pt, axis=1).min()))
        # Place H along OH…acceptor geometry: 0.98 Å from oi, away from oi's existing H.
        oh_h = _H_of(atoms, oi)
        assert len(oh_h) == 1
        v = atoms.positions[oi] - atoms.positions[oh_h[0]]
        # Prefer a direction roughly tetrahedral / opposite existing H, tilted up.
        direction = v / (np.linalg.norm(v) + 1e-12)
        # Mix with +z to avoid clash into the slab.
        direction = direction + np.array([0.0, 0.0, 0.4])
        direction /= np.linalg.norm(direction)
        atoms.positions[hi] = atoms.positions[oi] + 0.98 * direction
        used_o.add(far_o)
        transfers.append({"oh_O": oi, "from_O": far_o, "H": hi})
    return {"atoms": atoms, "transfers": transfers}


def swap_k_with_far_water(atoms: Atoms, n_swap: int | None = None) -> dict:
    """Exchange each K with a far H2O (O + 2H)."""
    atoms = atoms.copy()
    ks = _idx(atoms, "K")
    if n_swap is None:
        n_swap = len(ks)
    if n_swap > len(ks):
        raise RuntimeError("n_swap > nK")
    used: set[int] = set()
    swaps = []
    for ki in ks[:n_swap]:
        cands = far_waters(atoms, exclude=used)
        # Also require the far water not too close to this K already
        cands = [(d, oi) for d, oi in cands if atoms.get_distance(ki, oi, mic=True) > 2.0]
        if not cands:
            raise RuntimeError(f"no far water to swap with K{ki}")
        _, oi = cands[0]
        hs = _H_of(atoms, oi)
        if len(hs) != 2:
            raise RuntimeError(f"O{oi} not H2O")
        # Rigid translate water ↔ K: move O to K position, K to O position;
        # translate H by the same O→K vector.
        pos = atoms.positions
        kpos = pos[ki].copy()
        opos = pos[oi].copy()
        delta = kpos - opos
        pos[oi] = kpos
        for h in hs:
            pos[h] = pos[h] + delta
        pos[ki] = opos
        atoms.set_positions(pos)
        used.add(oi)
        pt = _pt_pos(atoms)
        d_new = float(np.linalg.norm(atoms.positions[ki] - pt, axis=1).min())
        swaps.append({"K": ki, "water_O": oi, "H": hs, "d_K_Pt_after": d_new})
    return {"atoms": atoms, "swaps": swaps}


def load_parent(kind: str) -> Atoms:
    """kind in {pris-koh, pris-water, Vo-koh, Vo-water}."""
    mapping = {
        "pris-koh": "pristine-koh",
        "pris-water": "pristine-water",
        "Vo-koh": "Vo-interface-koh",
        "Vo-water": "Vo-interface-water",
    }
    label = mapping[kind]
    traj = SRC31 / "work" / label / "relax.traj"
    atoms = read(traj)
    atoms.pbc = (True, True, False)
    return atoms


def metrics(atoms: Atoms) -> dict:
    pt = _pt_pos(atoms)
    ks = _idx(atoms, "K")
    oh = oh_on_pt(atoms)
    k_d = []
    for ki in ks:
        k_d.append({
            "index": ki,
            "z": float(atoms.positions[ki, 2]),
            "d_K_Pt_min": float(np.linalg.norm(atoms.positions[ki] - pt, axis=1).min()),
        })
    oh_info = []
    for oi in oh:
        oh_info.append({
            "O": oi,
            "z": float(atoms.positions[oi, 2]),
            "d_O_Pt_min": float(np.linalg.norm(atoms.positions[oi] - pt, axis=1).min()),
            "H": _H_of(atoms, oi),
        })
    # min interatomic distance (sanity)
    from ase.neighborlist import neighbor_list
    i, j, d = neighbor_list("ijd", atoms, 1.2)
    min_d = float(d.min()) if len(d) else 99.0
    # count pairs closer than 0.6 (clash)
    n_clash = int((d < 0.60).sum()) if len(d) else 0
    return {
        "formula": atoms.get_chemical_formula(),
        "natoms": len(atoms),
        "nK": len(ks),
        "nOH_on_Pt": len(oh),
        "K": k_d,
        "OH_on_Pt": oh_info,
        "min_pair_lt_1.2": min_d,
        "n_clash_lt_0.60": n_clash,
        "n_fixed": len(atoms.constraints[0].index) if atoms.constraints else 0,
    }


# label -> builder
# Tier A (1-8) + Tier B (9-12)
MODEL_SPECS = [
    # Tier A
    {"label": "pris-noK-noOH", "tier": "A", "parent": "pris-water",
     "K": "noK", "OH": "noOH", "defect": "pris",
     "note": "Baseline: no K, no OH@Pt (Pt may still bind H2O)."},
    {"label": "pris-Knear-noOH", "tier": "A", "parent": "pris-koh",
     "K": "Knear", "OH": "noOH", "defect": "pris",
     "note": "K near Pt without OH@Pt cancellation — key positive control."},
    {"label": "pris-Kfar-noOH", "tier": "A", "parent": "pris-koh",
     "K": "Kfar", "OH": "noOH", "defect": "pris",
     "note": "Distance control: same ions, K far from Pt."},
    {"label": "pris-Knear-OHonPt", "tier": "A", "parent": "pris-koh",
     "K": "Knear", "OH": "OHonPt", "defect": "pris",
     "note": "Parent pristine-koh: K near + OH@Pt (cancel case)."},
    {"label": "Vo-noK-noOH", "tier": "A", "parent": "Vo-water",
     "K": "noK", "OH": "noOH", "defect": "Vo",
     "note": "Vo baseline without K / OH@Pt."},
    {"label": "Vo-Knear-noOH", "tier": "A", "parent": "Vo-koh",
     "K": "Knear", "OH": "noOH", "defect": "Vo",
     "note": "Vo + K near, no OH@Pt."},
    {"label": "Vo-Kfar-noOH", "tier": "A", "parent": "Vo-koh",
     "K": "Kfar", "OH": "noOH", "defect": "Vo",
     "note": "Vo + K far, no OH@Pt."},
    {"label": "Vo-Knear-OHonPt", "tier": "A", "parent": "Vo-koh",
     "K": "Knear", "OH": "OHonPt", "defect": "Vo",
     "note": "Parent Vo-interface-koh: K near + OH@Pt."},
    # Tier B
    {"label": "pris-noK-OHonPt", "tier": "B", "parent": "pris-water",
     "K": "noK", "OH": "OHonPt", "defect": "pris",
     "note": "OH@Pt only (no K) — isolates OH penalty."},
    {"label": "Vo-noK-OHonPt", "tier": "B", "parent": "Vo-water",
     "K": "noK", "OH": "OHonPt", "defect": "Vo",
     "note": "Vo + OH@Pt only."},
    {"label": "pris-Kfar-OHonPt", "tier": "B", "parent": "pris-koh",
     "K": "Kfar", "OH": "OHonPt", "defect": "pris",
     "note": "OH@Pt + K far — OH penalty without near-K stabilization."},
    {"label": "Vo-Kfar-OHonPt", "tier": "B", "parent": "Vo-koh",
     "K": "Kfar", "OH": "OHonPt", "defect": "Vo",
     "note": "Vo + OH@Pt + K far."},
]


def build_one(spec: dict) -> tuple[Atoms, dict]:
    atoms = load_parent(spec["parent"])
    log: dict = {"label": spec["label"], "parent": spec["parent"], "ops": []}

    if spec["OH"] == "noOH":
        if oh_on_pt(atoms):
            r = protonate_oh_on_pt(atoms)
            atoms = r["atoms"]
            log["ops"].append({"op": "protonate_OH_on_Pt", "transfers": r["transfers"]})
        else:
            log["ops"].append({"op": "protonate_OH_on_Pt", "transfers": [], "skip": "already no OH@Pt"})

    if spec["K"] == "Kfar":
        if not _idx(atoms, "K"):
            raise RuntimeError(f"{spec['label']}: Kfar requested but no K in parent")
        r = swap_k_with_far_water(atoms)
        atoms = r["atoms"]
        log["ops"].append({"op": "swap_K_far_water", "swaps": r["swaps"]})
    elif spec["K"] == "Knear":
        log["ops"].append({"op": "keep_K_near"})
    elif spec["K"] == "noK":
        if _idx(atoms, "K"):
            raise RuntimeError(f"{spec['label']}: noK parent still has K")
        log["ops"].append({"op": "no_K_from_water_parent"})

    atoms = apply_fix_bottom(atoms)
    m = metrics(atoms)
    log["metrics"] = m

    # Hard checks against intended factors
    errors = []
    if spec["K"] == "noK" and m["nK"] != 0:
        errors.append("expected nK=0")
    if spec["K"] in ("Knear", "Kfar") and m["nK"] < 1:
        errors.append("expected K present")
    if spec["K"] == "Knear":
        for kd in m["K"]:
            if kd["d_K_Pt_min"] > 4.5:
                errors.append(f"Knear but d={kd['d_K_Pt_min']:.2f}>4.5")
    if spec["K"] == "Kfar":
        for kd in m["K"]:
            if kd["d_K_Pt_min"] < FAR_PT_CUT - 0.05:
                errors.append(f"Kfar but d={kd['d_K_Pt_min']:.2f}<{FAR_PT_CUT}")
    if spec["OH"] == "noOH" and m["nOH_on_Pt"] != 0:
        errors.append(f"expected no OH@Pt, found {m['nOH_on_Pt']}")
    if spec["OH"] == "OHonPt" and m["nOH_on_Pt"] < 1:
        errors.append("expected OH@Pt")
    if m["n_clash_lt_0.60"] > 0:
        errors.append(f"atomic clashes <0.60 Å: {m['n_clash_lt_0.60']}")
    log["errors"] = errors
    return atoms, log


def write_poscar(atoms: Atoms, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # ASE vasp writer preserves FixAtoms as Selective dynamics
    write(str(path), atoms, format="vasp", direct=True, vasp5=True, sort=False)


def main() -> None:
    struct = ROOT / "structures"
    work = ROOT / "work"
    struct.mkdir(parents=True, exist_ok=True)
    manifest = []
    n_err = 0
    for spec in MODEL_SPECS:
        atoms, log = build_one(spec)
        label = spec["label"]
        write_poscar(atoms, struct / f"{label}.vasp")
        wdir = work / label
        wdir.mkdir(parents=True, exist_ok=True)
        write_poscar(atoms, wdir / "POSCAR.init")
        # also save xyz for quick viz
        write(str(struct / f"{label}.xyz"), atoms)
        entry = {**spec, **log}
        manifest.append(entry)
        status = "OK" if not log["errors"] else "FAIL"
        if log["errors"]:
            n_err += 1
        kd = ", ".join(f"{x['d_K_Pt_min']:.2f}" for x in log["metrics"]["K"]) or "-"
        print(
            f"[{status}] {label:22s}  parent={spec['parent']:12s}  "
            f"nK={log['metrics']['nK']}  OH@Pt={log['metrics']['nOH_on_Pt']}  "
            f"d(K-Pt)=[{kd}]  fixed={log['metrics']['n_fixed']}",
            flush=True,
        )
        if log["errors"]:
            for e in log["errors"]:
                print(f"         ! {e}", flush=True)

    (struct / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"\nWrote {len(manifest)} models → {struct}/ and work/*/POSCAR.init")
    print(f"Validation errors: {n_err}")
    if n_err:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
