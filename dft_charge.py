#!/usr/bin/env python3
"""Unconstrained DFT charging (compare to CDFT).

After GS relax + gs.gpw exist, for each Δn:
  - set calculator charge = -dn (uniform background compensation)
  - SCF at fixed relaxed geometry
  - write dft_dnXpYY.gpw and append dft_scan.csv

Environment: OPT_LABEL, optional DN_SCAN.
Does NOT constrain charge onto Pt — electrons go where DFT puts them.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
from ase.io import read
from ase.parallel import world
from gpaw import MixerDif
from gpaw.old.calculator import GPAW

from common import DN_SCAN, ROOT, append_csv_row, make_parallel, pt_indices, write_json


def parse_dn_scan() -> list[float]:
    raw = os.environ.get("DN_SCAN", "").strip()
    if not raw:
        return list(DN_SCAN)
    return [float(x) for x in raw.split(",") if x.strip()]


def main() -> None:
    label = os.environ.get("OPT_LABEL")
    if not label:
        raise SystemExit("OPT_LABEL is not set")
    work = ROOT / "work" / label
    os.chdir(work)

    atoms = read("relaxed.vasp")
    atoms.pbc = (True, True, False)
    atoms.set_initial_magnetic_moments([0.0] * len(atoms))
    pt = pt_indices(atoms)
    parallel = make_parallel(atoms)
    gs = work / "gs.gpw"
    if not gs.is_file():
        raise SystemExit("gs.gpw missing — run pipeline static first")

    header = ["dn", "charge", "E_dft", "status"]
    csv = work / "dft_scan.csv"

    for dn in parse_dn_scan():
        tag = f"dn{dn:.2f}".replace(".", "p")
        out = work / f"dft_{tag}.gpw"
        if out.is_file():
            if world.rank == 0:
                print(f"[dft] skip {out.name}", flush=True)
            continue
        chg = -float(dn)  # extra electrons → negative total charge
        if world.rank == 0:
            print(f"[dft] Δn={dn:.2f}  charge={chg:+.2f}", flush=True)
        calc = GPAW(
            str(gs),
            txt=str(work / f"dft_{tag}.txt"),
            parallel=parallel,
            charge=chg,
            maxiter=500,
            mixer=MixerDif(
                beta=0.10, nmaxold=8, weight=10.0,
                beta_m=0.10, nmaxold_m=8, weight_m=10.0,
            ),
            convergence={
                "eigenstates": 1.0e-7,
                "density": 1.0e-4,
                "energy": 1.0e-4,
            },
        )
        atoms_c = atoms.copy()
        atoms_c.pbc = atoms.pbc
        atoms_c.set_initial_magnetic_moments([0.0] * len(atoms_c))
        atoms_c.calc = calc
        try:
            e = float(atoms_c.get_potential_energy())
            calc.write(str(out), mode="all")
            row = {"dn": dn, "charge": chg, "E_dft": e, "status": "ok"}
            if world.rank == 0:
                print(f"[dft] E={e:.6f} eV", flush=True)
        except Exception as exc:  # noqa: BLE001
            if world.rank == 0:
                print(f"[dft] FAILED Δn={dn}: {exc}", flush=True)
            row = {"dn": dn, "charge": chg, "E_dft": "", "status": f"fail:{type(exc).__name__}"}
        append_csv_row(csv, header, row)
        write_json(work / f"dft_{tag}_meta.json", {"label": label, "pt": pt, **row})
        world.barrier()

    if world.rank == 0:
        print(f"[{label}] dft_charge finished", flush=True)


if __name__ == "__main__":
    main()
