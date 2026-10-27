#!/usr/bin/env python3
"""Shared helpers for Pt-region CDFT workflow (static Q0 + charge scan)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from ase.dft.kpoints import monkhorst_pack
from ase.parallel import world
from gpaw import Davidson, FermiDirac, MixerDif
from gpaw.old.calculator import GPAW

ROOT = Path(__file__).resolve().parent
# Tier A+B matrix (see MODELS.md). Built by build_models.py.
LABELS = [
    # Tier A
    "pris-noK-noOH",
    "pris-Knear-noOH",
    "pris-Kfar-noOH",
    "pris-Knear-OHonPt",
    "Vo-noK-noOH",
    "Vo-Knear-noOH",
    "Vo-Kfar-noOH",
    "Vo-Knear-OHonPt",
    # Tier B
    "pris-noK-OHonPt",
    "Vo-noK-OHonPt",
    "pris-Kfar-OHonPt",
    "Vo-Kfar-OHonPt",
]
# Extra electrons on Pt relative to unconstrained Q0 (Q_target = Q0 - dn).
DN_SCAN = [0.0, 0.5, 1.0, 1.5, 2.0]


def grid_points(cell, h: float = 0.15) -> tuple[int, int, int]:
    from gpaw.utilities import h2gpts

    return tuple(int(n) for n in h2gpts(h, cell, idiv=4))


def choose_parallel(size: int, gpts: tuple[int, int, int], lengths):
    def divisors(n: int):
        return [d for d in range(1, n + 1) if n % d == 0]

    best = None
    for nb in (2, 4, 1):
        if size % nb:
            continue
        nd = size // nb
        for nx in divisors(gpts[0]):
            if nd % nx or gpts[0] // nx < 4:
                continue
            rest = nd // nx
            for ny in divisors(gpts[1]):
                if rest % ny or gpts[1] // ny < 4:
                    continue
                nz = rest // ny
                if nz < 1 or gpts[2] % nz or gpts[2] // nz < 4:
                    continue
                edges = (
                    lengths[0] / nx,
                    lengths[1] / ny,
                    lengths[2] / nz,
                )
                aspect = max(edges) / min(edges)
                score = (aspect, 0 if nb == 2 else 1, nb)
                cand = (score, nb, (nx, ny, nz), aspect)
                if best is None or cand[0] < best[0]:
                    best = cand
    if best is None:
        raise RuntimeError(f"no domain/band split for size={size} gpts={gpts}")
    return best[1], best[2], best[3]


def make_parallel(atoms) -> dict:
    gpts = grid_points(atoms.cell, 0.15)
    lengths = np.linalg.norm(np.asarray(atoms.cell), axis=1)
    nb, domain, aspect = choose_parallel(world.size, gpts, lengths)
    if world.rank == 0:
        print(
            f"gpts={gpts} domain={domain} band={nb} aspect={aspect:.2f}",
            flush=True,
        )
    return {
        "domain": domain,
        "band": nb,
        "augment_grids": True,
        "sl_auto": True,
    }


def make_dft(txt: str, parallel: dict) -> GPAW:
    return GPAW(
        h=0.15,
        mode="fd",
        charge=0,
        xc="PBE",
        kpts=monkhorst_pack((1, 1, 1)),
        symmetry="off",
        occupations=FermiDirac(0.1, fixmagmom=False),
        setups={"Ti": ":d,5.8"},
        eigensolver=Davidson(3),
        spinpol=True,
        mixer=MixerDif(
            beta=0.20,
            nmaxold=10,
            weight=10.0,
            beta_m=0.20,
            nmaxold_m=10,
            weight_m=10.0,
        ),
        maxiter=400,
        txt=txt,
        convergence={
            "eigenstates": 1.0e-8,
            "density": 1.0e-4,
            "energy": 1.0e-5,
        },
        poissonsolver={"dipolelayer": "xy"},
        parallel=parallel,
    )


def pt_indices(atoms) -> list[int]:
    return [i for i, s in enumerate(atoms.get_chemical_symbols()) if s == "Pt"]


def measure_q0(atoms, gs_gpw: Path, indices: list[int], parallel: dict, mu: dict | None = None):
    """Unconstrained Pt charge from one CDFT evaluation at λ=0, Q_target=0.

    With Q_target=0 the constraint is Ne=Z, so Error=Ne-Z=-Q and Q0=-Error.
    gtol is set huge so the λ-optimizer stops after the first SCF (λ stays 0).
    This uses GPAW CDFT's own weight + AE path — same as later scan points.
    """
    from gpaw.cdft.cdft import CDFT

    atoms_q = atoms.copy()
    atoms_q.pbc = atoms.pbc
    atoms_q.set_initial_magnetic_moments([0.0] * len(atoms_q))

    calc = GPAW(str(gs_gpw), txt="q0_scf.txt", parallel=parallel, maxiter=400)
    cdft = CDFT(
        calc=calc,
        atoms=atoms_q,
        charge_regions=[list(indices)],
        charges=[0.0],
        charge_coefs=[0.0],
        method="L-BFGS-B",
        txt="q0_probe.cdft",
        # |Error| << 1e12 after first eval → λ-opt declares convergence at λ=0
        minimizer_options={"gtol": 1.0e12},
        compute_forces=False,
        restart=True,
        **({"mu": mu} if mu is not None else {}),
    )
    atoms_q.calc = cdft
    atoms_q.get_potential_energy()

    err = float(np.asarray(cdft.dn_i).ravel()[0])  # Ne - Z
    q = -err
    z = float(sum(atoms_q[i].number for i in indices))
    n_el = z - q
    rc_use = {k: float(v) for k, v in cdft.Rc.items()}
    mu_use = {k: float(v) for k, v in cdft.mu.items()}
    if world.rank == 0:
        print(
            f"[Q0] CDFT λ=0 probe  Error={err:.4f}  →  Q0={q:.4f}  "
            f"Ne={n_el:.4f}  Z={z:.1f}",
            flush=True,
        )
    return {
        "Q": q,
        "Ne": n_el,
        "Z": z,
        "Rc": rc_use,
        "mu": mu_use,
        "indices": list(indices),
        "Error_at_lambda0": err,
    }


def write_json(path: Path, data: dict) -> None:
    if world.rank == 0:
        path.write_text(json.dumps(data, indent=2) + "\n")


def append_csv_row(path: Path, header: list[str], row: dict) -> None:
    if world.rank != 0:
        return
    new = not path.is_file()
    with path.open("a") as f:
        if new:
            f.write(",".join(header) + "\n")
        f.write(",".join(str(row[h]) for h in header) + "\n")
