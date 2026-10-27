#!/usr/bin/env python3
"""Hirshfeld + Bader charges for one saved GPAW density (.gpw).

usage: python charge_analysis.py <label> <gpw-name>
Writes analysis/charges/<label>/<gpw-stem>.json with per-atom charges and
group sums (Pt6, K, OH/H2O solvent, slab). Charge = Z_valence - N (Hirshfeld)
or Z - N_all-electron (Bader); positive = electron-poor.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from ase.io import write
from ase.units import Bohr
from gpaw.analyse.hirshfeld import HirshfeldPartitioning
from gpaw.old.calculator import GPAW

ROOT = Path(__file__).resolve().parent
BADER = Path.home() / "software/bin/bader"


def groups(atoms) -> dict[str, list[int]]:
    sym = np.array(atoms.get_chemical_symbols())
    h = np.where(sym == "H")[0]
    solvent_o = []
    for i in np.where(sym == "O")[0]:
        if len(h) and atoms.get_distances(i, h, mic=True).min() < 1.25:
            solvent_o.append(i)
    oh, water = [], []
    for i in solvent_o:
        nh = int((atoms.get_distances(i, h, mic=True) < 1.25).sum())
        bonded = list(h[atoms.get_distances(i, h, mic=True) < 1.25])
        (oh if nh == 1 else water).extend([i, *bonded])
    solvent = set(solvent_o) | set(h.tolist())
    slab = [i for i in range(len(atoms)) if sym[i] in ("Sr", "Ti", "O") and i not in solvent]
    return {
        "Pt6": np.where(sym == "Pt")[0].tolist(),
        "K": np.where(sym == "K")[0].tolist(),
        "OH": sorted(int(i) for i in oh),
        "H2O": sorted(int(i) for i in water),
        "slab": slab,
    }


def bader_charges(calc, atoms, workdir: Path) -> list[float]:
    rho = calc.get_all_electron_density(gridrefinement=2) * Bohr**3
    workdir.mkdir(parents=True, exist_ok=True)
    cube = workdir / "rho.cube"
    write(str(cube), atoms, data=rho)
    subprocess.run([str(BADER), cube.name], cwd=workdir, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    n = []
    for line in (workdir / "ACF.dat").read_text().splitlines()[2:]:
        parts = line.split()
        if len(parts) >= 5 and parts[0].isdigit():
            n.append(float(parts[4]))
    z = atoms.get_atomic_numbers()
    return [float(zi - ni) for zi, ni in zip(z, n)]


def main() -> None:
    label, gpw = sys.argv[1], sys.argv[2]
    work = ROOT / "work" / label
    out_dir = ROOT / "analysis" / "charges" / label
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = Path(gpw).stem
    t0 = time.time()

    calc = GPAW(str(work / gpw), txt=None)
    atoms = calc.get_atoms()
    hir = [float(x) for x in HirshfeldPartitioning(calc).get_charges()]
    t1 = time.time()
    tmp = out_dir / f"_bader_{stem}"
    bad = bader_charges(calc, atoms, tmp)
    shutil.rmtree(tmp, ignore_errors=True)

    g = groups(atoms)
    res = {
        "label": label,
        "gpw": gpw,
        "symbols": atoms.get_chemical_symbols(),
        "groups": g,
        "hirshfeld": hir,
        "bader": bad,
        "sum_hirshfeld": {k: float(sum(hir[i] for i in v)) for k, v in g.items()},
        "sum_bader": {k: float(sum(bad[i] for i in v)) for k, v in g.items()},
        "total_hirshfeld": float(sum(hir)),
        "total_bader": float(sum(bad)),
        "seconds": {"hirshfeld": t1 - t0, "bader": time.time() - t1},
    }
    (out_dir / f"{stem}.json").write_text(json.dumps(res, indent=1))
    print(f"{label}/{gpw}: Pt6 H={res['sum_hirshfeld']['Pt6']:+.3f} B={res['sum_bader']['Pt6']:+.3f} "
          f"K H={res['sum_hirshfeld']['K']:+.3f} B={res['sum_bader']['K']:+.3f} "
          f"({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
