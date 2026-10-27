# 32 — Pt/STO CDFT extended matrix (K / OH / Vo)

Extension of `31-Pt-DEDL-CDFT` to prove **K⁺ near Pt stabilizes extra electrons on Pt₆**, with controls for **K–Pt distance**, **OH⁻ on Pt**, and **interface Vo**.

Full model table and construction rules: **[MODELS.md](MODELS.md)**.

## Layout

```
32-Pt-DEDL-CDFT-extended/
  build_models.py      # build 12 POSCAR.init from 31 relaxed GS
  validate_models.py   # factor / clash / stoichiometry checks
  structures/          # *.vasp *.xyz manifest.json validation.json
  work/<label>/POSCAR.init
  relax.py             # GS PBE+U+D3 BFGS (same as 31)
  pipeline.py          # static gs.gpw + fixed-geo CDFT Δn scan
  dft_charge.py        # unconstrained DFT charging (CDFT comparison)
  k_stab.py / k_field.py / charge_analysis.py
  submit.sh / submit_pipeline.sh / submit_dft_charge.sh
  launch_all_relax.sh
```

## Workflow per label

```text
POSCAR.init  →  relax.py  →  relaxed.vasp
                         →  pipeline.py   (gs.gpw + cdft_dn*.gpw, cdft_scan.csv)
                         →  dft_charge.py (dft_dn*.gpw, dft_scan.csv)   [after gs.gpw]
                         →  charge_analysis / k_stab / k_field
```

### 1. Build (already done)

```bash
cd /global/home/users/jiaao/workspace/STO/32-Pt-DEDL-CDFT-extended
source env.sh
python build_models.py
python validate_models.py
```

### 2. GS structure optimization

```bash
# all 12
./launch_all_relax.sh

# or one
OPT_LABEL=pris-Knear-noOH sbatch --export=ALL,OPT_LABEL submit.sh
```

Protocol (unchanged from 31): PBE+U(Ti=5.8), DFT-D3 zero, h=0.15, Γ, spinpol, dipole xy, bottom TiO₂ fixed, BFGS fmax=0.05, 128 ranks.

### 3. CDFT (fixed GS geometry)

Needs `work/<label>/relaxed.vasp`.

```bash
OPT_LABEL=pris-Knear-noOH sbatch --export=ALL,OPT_LABEL submit_pipeline.sh
# optional: DN_SCAN=0,0.5,1,1.5,2
```

### 4. Unconstrained DFT charging (compare to CDFT)

Needs `gs.gpw` from the pipeline static step.

```bash
OPT_LABEL=pris-Knear-noOH sbatch --export=ALL,OPT_LABEL submit_dft_charge.sh
```

Compare for each Δn:

| Quantity | CDFT | DFT (`dft_charge`) |
|----------|------|---------------------|
| Where electrons go | Forced onto Pt₆ | Wherever SCF puts them (Bader) |
| Energy | `E_cdft` / `E_dft` in `cdft_scan.csv` | `E_dft` in `dft_scan.csv` |
| Interpretation | Cost to put Δn on Pt | Natural charging response |

If unconstrained DFT already piles charge on Pt when K is near (and less when K is far / absent), that corroborates the CDFT stabilization story without the constraint.

### 5. Analysis (after densities exist)

```bash
python charge_analysis.py <label> gs.gpw
python charge_analysis.py <label> cdft_dn1p00.gpw
python k_stab.py <label>
python k_field.py <label>   # K-containing labels only
```

## Suggested run order (cost)

1. **pris Tier A** (`pris-noK-noOH`, `pris-Knear-noOH`, `pris-Kfar-noOH`, `pris-Knear-OHonPt`) — enough for the main claim  
2. **pris Tier B** (`pris-noK-OHonPt`, `pris-Kfar-OHonPt`)  
3. **Vo Tier A+B** — defect modulation  

## Inputs ready checklist

- [x] 12 named models under `work/*/POSCAR.init`
- [x] Structures + xyz in `structures/`
- [x] `manifest.json` / `validation.json` (all OK)
- [x] `common.LABELS` updated
- [x] `relax` / `pipeline` / `dft_charge` + Slurm submit scripts
- [x] `MODELS.md` documentation
- [ ] GS relax jobs (submit when ready)
- [ ] CDFT + DFT charge scans
