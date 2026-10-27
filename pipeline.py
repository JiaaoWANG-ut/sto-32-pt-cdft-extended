#!/usr/bin/env python3
"""Per-system pipeline: static GS → Pt Q0 → fixed-geometry CDFT Δn scan.

Environment:
  OPT_LABEL   one of common.LABELS
  DN_SCAN     optional comma list, default common.DN_SCAN
"""
from __future__ import annotations

import os
import traceback
from pathlib import Path

import numpy as np
from ase.io import read
from ase.parallel import world
from gpaw.cdft.cdft import CDFT
from gpaw.old.calculator import GPAW

from common import (
    DN_SCAN,
    ROOT,
    append_csv_row,
    make_dft,
    make_parallel,
    measure_q0,
    pt_indices,
    write_json,
)


def parse_dn_scan() -> list[float]:
    raw = os.environ.get("DN_SCAN", "").strip()
    if not raw:
        return list(DN_SCAN)
    return [float(x) for x in raw.split(",") if x.strip()]


def run_static(atoms, work: Path, parallel: dict) -> GPAW:
    gpw = work / "gs.gpw"
    if gpw.is_file():
        if world.rank == 0:
            print(f"[static] restart from {gpw}", flush=True)
        calc = GPAW(str(gpw), txt=str(work / "static.txt"))
        atoms.calc = calc
        # Ensure density is available without a full redo if already converged.
        atoms.get_potential_energy()
        return calc

    if world.rank == 0:
        print("[static] fresh ground-state SCF", flush=True)
    calc = make_dft(str(work / "static.txt"), parallel)
    atoms.calc = calc
    e = atoms.get_potential_energy()
    if world.rank == 0:
        print(f"[static] E = {e:.6f} eV", flush=True)
    calc.write(str(gpw), mode="all")
    if world.rank == 0:
        print(f"[static] wrote {gpw}", flush=True)
    return calc


def run_cdft_point(
    atoms,
    work: Path,
    parallel: dict,
    pt: list[int],
    q0: float,
    dn: float,
) -> dict:
    q_target = q0 - dn
    tag = f"dn{dn:.2f}".replace(".", "p")
    txt = work / f"cdft_{tag}.txt"
    log = work / f"cdft_{tag}.cdft"
    out_gpw = work / f"cdft_{tag}.gpw"

    if world.rank == 0:
        print(
            f"[cdft] Δn={dn:.2f}  Q_target={q_target:.4f}  "
            f"(Q0={q0:.4f})",
            flush=True,
        )

    # Fresh calculator from GS density for each constraint.
    # Do NOT reuse Rc/mu from the Q0 probe: CDFT/WeightFunc(new=True) expects
    # Angstrom and would convert Bohr values twice, changing the partitioning.
    gs = work / "gs.gpw"
    from gpaw import MixerDif

    calc = GPAW(
        str(gs),
        txt=str(txt),
        parallel=parallel,
        maxiter=500,
        mixer=MixerDif(
            beta=0.10,
            nmaxold=8,
            weight=10.0,
            beta_m=0.10,
            nmaxold_m=8,
            weight_m=10.0,
        ),
        convergence={
            "eigenstates": 1.0e-7,
            "density": 1.0e-4,
            "energy": 1.0e-4,
        },
    )
    atoms_pt = atoms.copy()
    atoms_pt.pbc = atoms.pbc
    atoms_pt.set_initial_magnetic_moments([0.0] * len(atoms_pt))

    coef0 = 0.0

    cdft = CDFT(
        calc=calc,
        atoms=atoms_pt,
        charge_regions=[list(pt)],
        charges=[q_target],
        charge_coefs=[coef0],
        method="L-BFGS-B",
        txt=str(log),
        minimizer_options={"gtol": 0.01},
        compute_forces=False,
        restart=True,
        maxstep=20,
        tol=1e-3,
    )
    atoms_pt.calc = cdft
    e_cdft = atoms_pt.get_potential_energy()
    e_dft = float(cdft.dft_energy())
    lam = float(np.asarray(cdft.get_lagrangians()).ravel()[0])

    # err = Ne - (Z - Q_target)  →  Q_achieved = Q_target - err
    err = float(np.asarray(cdft.dn_i).ravel()[0])
    z = float(sum(atoms_pt[i].number for i in pt))
    q_achieved = q_target - err
    n_el = z - q_achieved
    cdft.calc.write(str(out_gpw), mode="all")

    row = {
        "dn": dn,
        "Q0": q0,
        "Q_target": q_target,
        "Q_achieved": q_achieved,
        "Ne": n_el,
        "Z": z,
        "E_cdft": e_cdft,
        "E_dft": e_dft,
        "lambda_eV": lam,
        "status": "ok",
    }
    if world.rank == 0:
        print(
            f"[cdft] done Δn={dn:.2f}  Q={q_achieved:.4f}  "
            f"λ={lam:.4f} eV  Ecdft={e_cdft:.6f}",
            flush=True,
        )
    return row


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
    if len(pt) != 6:
        raise SystemExit(f"expected 6 Pt, found {len(pt)}: {pt}")

    if world.rank == 0:
        print(
            f"[{label}] natoms={len(atoms)} Pt={pt} ranks={world.size}",
            flush=True,
        )

    parallel = make_parallel(atoms)
    calc = run_static(atoms, work, parallel)
    # Free the static calculator before the Q0 probe opens gs.gpw again.
    del calc

    qinfo = measure_q0(atoms, work / "gs.gpw", pt, parallel)
    q0 = float(qinfo["Q"])
    write_json(
        work / "q0.json",
        {
            "label": label,
            "pt_indices": pt,
            **{
                k: qinfo[k]
                for k in ("Q", "Ne", "Z", "Rc", "mu", "Error_at_lambda0")
            },
        },
    )
    if world.rank == 0:
        print(
            f"[Q0] Pt6  Z={qinfo['Z']:.1f}  Ne={qinfo['Ne']:.4f}  "
            f"Q0={q0:.4f}",
            flush=True,
        )

    csv = work / "cdft_scan.csv"
    header = [
        "dn",
        "Q0",
        "Q_target",
        "Q_achieved",
        "Ne",
        "Z",
        "E_cdft",
        "E_dft",
        "lambda_eV",
        "status",
    ]
    rc = qinfo["Rc"]
    mu = qinfo["mu"]
    dn_list = parse_dn_scan()

    for dn in dn_list:
        tag = f"dn{dn:.2f}".replace(".", "p")
        done_gpw = work / f"cdft_{tag}.gpw"
        if done_gpw.is_file():
            if world.rank == 0:
                print(f"[cdft] skip existing {done_gpw.name}", flush=True)
            continue
        try:
            row = run_cdft_point(atoms, work, parallel, pt, q0, dn)
        except Exception as exc:  # noqa: BLE001 — log and continue scan
            if world.rank == 0:
                print(f"[cdft] FAILED Δn={dn}: {exc}", flush=True)
                traceback.print_exc()
            row = {
                "dn": dn,
                "Q0": q0,
                "Q_target": q0 - dn,
                "Q_achieved": "",
                "Ne": "",
                "Z": qinfo["Z"],
                "E_cdft": "",
                "E_dft": "",
                "lambda_eV": "",
                "status": f"fail:{type(exc).__name__}",
            }
        append_csv_row(csv, header, row)
        world.barrier()

    if world.rank == 0:
        print(f"[{label}] pipeline finished", flush=True)


if __name__ == "__main__":
    main()
