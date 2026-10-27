#!/usr/bin/env python3
"""CDFT geometry optimization: extra electrons on Pt6, one Δn per process.

Environment:
  OPT_LABEL   system folder under work/
  DN_OPT      the single Δn to optimize

Reads Q0 from work/<label>/q0.json. Starts from the optimized geometry of the
largest finished smaller Δn (else relaxed.vasp), or resumes from this Δn's
.traj/.pckl if a previous run was interrupted. FixAtoms come from relaxed.vasp.

A second CDFT/GPAW calculation inside the same MPI process hangs (seen on all
four systems that reached Δn=0.5), so only one Δn is allowed per process.
"""
from __future__ import annotations

import json
import os
import traceback
from pathlib import Path

import numpy as np
from ase.io import read, write
from ase.optimize.bfgs import BFGS
from ase.parallel import world
from dftd3.ase import DFTD3
from gpaw import MixerDif
from gpaw.cdft.cdft import CDFT
from gpaw.old.calculator import GPAW

from common import ROOT, append_csv_row, make_parallel, pt_indices

DN_OPT_DEFAULT = [0.1, 0.2, 0.5, 1.0, 1.5, 2.0]
FMAX = 0.05
MAX_STEPS = 300


def parse_dn() -> list[float]:
    """DN_OPT accepts comma, colon, or whitespace (colon preferred for Slurm --export)."""
    raw = os.environ.get("DN_OPT", "").strip()
    if not raw:
        return list(DN_OPT_DEFAULT)
    for sep in (",", ":", " "):
        if sep in raw:
            return [float(x) for x in raw.replace(",", " ").replace(":", " ").split() if x.strip()]
    return [float(raw)]


def tag_dn(dn: float) -> str:
    return f"dn{dn:.2f}".replace(".", "p")


def load_q0(work: Path) -> float:
    path = work / "q0.json"
    if not path.is_file():
        raise SystemExit(f"missing {path}; run charge scan pipeline first")
    return float(json.loads(path.read_text())["Q"])


def make_gpaw(gs: Path, txt: Path, parallel: dict) -> GPAW:
    return GPAW(
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


def previous_dn_structure(work: Path, dn: float, atoms0):
    """Optimized geometry of the largest finished Δn below this one, else relaxed.vasp."""
    best = None
    for p in work.glob("cdft_opt_dn*.vasp"):
        d = float(p.stem.removeprefix("cdft_opt_dn").replace("p", "."))
        if d < dn - 1e-9 and (best is None or d > best[0]):
            best = (d, p)
    if best is None:
        return atoms0, "relaxed.vasp"
    return read(str(best[1])), best[1].name


def run_opt(
    atoms0,
    work: Path,
    parallel: dict,
    pt: list[int],
    q0: float,
    dn: float,
) -> dict:
    q_target = q0 - dn
    tag = tag_dn(dn)
    out_vasp = work / f"cdft_opt_{tag}.vasp"
    out_gpw = work / f"cdft_opt_{tag}.gpw"
    traj = work / f"cdft_opt_{tag}.traj"
    pckl = work / f"cdft_opt_{tag}.pckl"
    log = work / f"cdft_opt_{tag}.log"
    cdft_log = work / f"cdft_opt_{tag}.cdft"
    txt = work / f"cdft_opt_{tag}.txt"

    if out_vasp.is_file():
        if world.rank == 0:
            print(f"[opt] skip existing {out_vasp.name}", flush=True)
        return None

    resuming = traj.is_file() and traj.stat().st_size > 0
    if resuming:
        start = read(str(traj), index=-1)
        source = traj.name
    else:
        start, source = previous_dn_structure(work, dn, atoms0)
        if pckl.is_file():
            pckl.unlink()

    if world.rank == 0:
        print(
            f"[opt] Δn={dn:.2f}  Q_target={q_target:.4f}  "
            f"(Q0={q0:.4f})  fmax={FMAX}  start={source}",
            flush=True,
        )

    atoms = atoms0.copy()
    atoms.positions = start.positions
    atoms.pbc = (True, True, False)
    # FixAtoms always come from relaxed.vasp so every Δn uses the same frozen layer.
    if getattr(atoms0, "constraints", None):
        atoms.set_constraint(atoms0.constraints)
    atoms.set_initial_magnetic_moments([0.0] * len(atoms))

    gs = work / "gs.gpw"
    calc = make_gpaw(gs, txt, parallel)
    cdft = CDFT(
        calc=calc,
        atoms=atoms,
        charge_regions=[list(pt)],
        charges=[q_target],
        charge_coefs=[0.0],
        method="L-BFGS-B",
        txt=str(cdft_log),
        minimizer_options={"gtol": 0.01},
        compute_forces=True,
        restart=True,
        maxstep=20,
        tol=1e-3,
    )
    # Same D3 protocol as the neutral geometry relax.
    atoms.calc = DFTD3(method="PBE", damping="d3zero").add_calculator(cdft)

    opt_log = str(log) if world.rank == 0 else None
    dyn = BFGS(
        atoms,
        trajectory=str(traj),
        append_trajectory=resuming,
        restart=str(pckl),
        logfile=opt_log,
    )
    dyn.run(fmax=FMAX, steps=MAX_STEPS)

    e_cdft = float(atoms.get_potential_energy())
    forces = atoms.get_forces()
    # Mask fixed atoms for reported fmax if constraints present
    fmax = float(np.sqrt((forces**2).sum(axis=1).max()))
    steps = int(getattr(dyn, "nsteps", -1))

    err = float(np.asarray(cdft.dn_i).ravel()[0])
    q_achieved = q_target - err

    if world.rank == 0:
        write(
            str(out_vasp),
            atoms,
            format="vasp",
            direct=True,
            sort=False,
            vasp5=True,
        )
        print(
            f"[opt] done Δn={dn:.2f}  steps={steps}  fmax={fmax:.4f}  "
            f"Q={q_achieved:.4f}  E={e_cdft:.6f}",
            flush=True,
        )
    # Density only; wavefunctions (mode="all") were 7–11 GB per point.
    cdft.calc.write(str(out_gpw))
    world.barrier()

    return {
        "dn": dn,
        "Q0": q0,
        "Q_target": q_target,
        "Q_achieved": q_achieved,
        "E_cdft": e_cdft,
        "fmax": fmax,
        "steps": steps,
        "status": "ok" if fmax <= FMAX + 1e-6 else "not_converged",
    }


def main() -> None:
    label = os.environ.get("OPT_LABEL")
    if not label:
        raise SystemExit("OPT_LABEL is not set")

    work = ROOT / "work" / label
    os.chdir(work)

    atoms0 = read("relaxed.vasp")
    atoms0.pbc = (True, True, False)
    pt = pt_indices(atoms0)
    if len(pt) != 6:
        raise SystemExit(f"expected 6 Pt, found {len(pt)}: {pt}")

    q0 = load_q0(work)
    parallel = make_parallel(atoms0)
    dn_list = parse_dn()
    if len(dn_list) != 1:
        raise SystemExit(f"DN_OPT must hold exactly one Δn per process, got {dn_list}")

    if world.rank == 0:
        print(
            f"[{label}] CDFT opt  natoms={len(atoms0)}  Pt={pt}  "
            f"Q0={q0:.6f}  dn={dn_list}  ranks={world.size}",
            flush=True,
        )

    csv = work / "cdft_opt_scan.csv"
    header = [
        "dn",
        "Q0",
        "Q_target",
        "Q_achieved",
        "E_cdft",
        "fmax",
        "steps",
        "status",
    ]

    for dn in dn_list:
        try:
            row = run_opt(atoms0, work, parallel, pt, q0, dn)
        except Exception as exc:  # noqa: BLE001
            if world.rank == 0:
                print(f"[opt] FAILED Δn={dn}: {exc}", flush=True)
                traceback.print_exc()
            row = {
                "dn": dn,
                "Q0": q0,
                "Q_target": q0 - dn,
                "Q_achieved": "",
                "E_cdft": "",
                "fmax": "",
                "steps": "",
                "status": f"fail:{type(exc).__name__}",
            }
        if row is not None:
            append_csv_row(csv, header, row)
        world.barrier()

    if world.rank == 0:
        print(f"[{label}] CDFT opt pipeline finished", flush=True)


if __name__ == "__main__":
    main()
