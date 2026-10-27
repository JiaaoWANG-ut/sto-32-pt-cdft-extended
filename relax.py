#!/usr/bin/env python3
"""GPAW DFT-D3 GS relax — same protocol as 07-STO-Pt6-KOH-water/opt-dft-clean.

Bottom TiO2 fixed (Selective dynamics in POSCAR.init).
PBE+U(Ti=5.8), h=0.15, Gamma, spinpol, D3-zero, dipole xy, BFGS fmax=0.05.

128 ranks on one lr8 node: domain groups that divide the FD grid, plus a
band group of 2 when that makes the real-space domains more cubic.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
from ase.dft.kpoints import monkhorst_pack
from ase.io import read, write
from ase.optimize.bfgs import BFGS
from ase.parallel import world
from dftd3.ase import DFTD3
from gpaw import Davidson, FermiDirac, MixerDif
from gpaw.old.calculator import GPAW

ROOT = Path(__file__).resolve().parent


def grid_points(cell, h: float = 0.15) -> tuple[int, int, int]:
    from gpaw.utilities import h2gpts

    return tuple(int(n) for n in h2gpts(h, cell, idiv=4))


def choose_parallel(size: int, gpts: tuple[int, int, int], lengths):
    """Pick (domain, band) with size = nx*ny*nz*nb and each n dividing the grid."""

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
                # Band parallel of 2 is the useful split for Davidson on this
                # slab; more than that lengthens the z-domains.
                score = (aspect, 0 if nb == 2 else 1, nb)
                cand = (score, nb, (nx, ny, nz), aspect)
                if best is None or cand[0] < best[0]:
                    best = cand
    if best is None:
        raise RuntimeError(f"no domain/band split for size={size} gpts={gpts}")
    return best[1], best[2], best[3]


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
            beta=0.20, nmaxold=10, weight=10.0,
            beta_m=0.20, nmaxold_m=10, weight_m=10.0,
        ),
        maxiter=300,
        txt=txt,
        convergence={"eigenstates": 1.0e-8, "density": 1.0e-4, "energy": 1.0e-5},
        poissonsolver={"dipolelayer": "xy"},
        parallel=parallel,
    )


def main() -> None:
    label = os.environ.get("OPT_LABEL")
    if not label:
        raise SystemExit("OPT_LABEL is not set")
    work = ROOT / "work" / label
    os.chdir(work)

    atoms = read("POSCAR.init")
    atoms.pbc = (True, True, False)
    atoms.set_initial_magnetic_moments([0.0] * len(atoms))

    fixed = []
    for c in atoms.constraints:
        fixed.extend(getattr(c, "index", getattr(c, "indices", [])))
    if world.rank == 0:
        print(
            f"[{label}] natoms={len(atoms)} fixed={len(fixed)} "
            f"free={len(atoms) - len(fixed)} ranks={world.size}",
            flush=True,
        )

    gpts = grid_points(atoms.cell, 0.15)
    lengths = np.linalg.norm(np.asarray(atoms.cell), axis=1)
    nb, domain, aspect = choose_parallel(world.size, gpts, lengths)
    parallel = {
        "domain": domain,
        "band": nb,
        "augment_grids": True,
        "sl_auto": True,
    }
    if world.rank == 0:
        print(
            f"[{label}] gpts={gpts} domain={domain} band={nb} "
            f"aspect={aspect:.2f}",
            flush=True,
        )

    dft = make_dft("relax.txt", parallel)
    atoms.calc = DFTD3(method="PBE", damping="d3zero").add_calculator(dft)

    restarting = Path("relax.pckl").is_file()
    if not restarting:
        atoms.get_potential_energy()
        atoms.get_forces()

    log = "relax.log" if world.rank == 0 else None
    relax = BFGS(atoms, trajectory="relax.traj", restart="relax.pckl", logfile=log)
    relax.run(fmax=0.05)

    if world.rank == 0:
        write("relaxed.vasp", atoms, format="vasp", direct=True, sort=True, vasp5=True)
        print(f"[{label}] wrote relaxed.vasp", flush=True)


if __name__ == "__main__":
    main()
