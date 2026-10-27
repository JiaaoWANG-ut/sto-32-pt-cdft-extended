# Model matrix — 32-Pt-DEDL-CDFT-extended

Purpose: isolate whether **K⁺ near Pt stabilizes extra electrons on Pt₆**, and how **OH⁻ on Pt**, **K–Pt distance**, and **interface Vo** change that.

Parents are the **relaxed GS** structures from `31-Pt-DEDL-CDFT` (`pristine-koh/water`, `Vo-interface-koh/water`). New models are **stoichiometry-conserving edits** of those parents, then re-optimized.

## Naming

```
{defect}-{K}-{OH}
```

| Token | Meaning | Operational definition |
|-------|---------|------------------------|
| `pris` | No oxygen vacancy | Parent `pristine-*` |
| `Vo` | Interface Vo | Parent `Vo-interface-*` (one fewer lattice O) |
| `noK` | No potassium | Water parent |
| `Knear` | K⁺ near Pt | Keep parent K; d(K–Pt)_min ≈ 3.1–3.4 Å |
| `Kfar` | K⁺ far from Pt | Swap each K with a far H₂O; d(K–Pt)_min ≳ 6 Å |
| `noOH` | No OH directly on Pt | Proton-transfer: H from far H₂O → each OH@Pt → H₂O@Pt + far OH⁻ |
| `OHonPt` | ≥1 OH bound to Pt | d(O–Pt) < 2.5 Å and nH(O)=1 (no MIC for O–H) |

## Full list (12)

### Tier A — required for the K story

| Label | Defect | K | OH@Pt | Parent | Role |
|-------|--------|---|-------|--------|------|
| `pris-noK-noOH` | pris | none | 0 | pris-water | Baseline charging |
| `pris-Knear-noOH` | pris | near | 0 | pris-koh | **K stabilizes alone** |
| `pris-Kfar-noOH` | pris | far | 0 | pris-koh | **Distance causal control** |
| `pris-Knear-OHonPt` | pris | near | 2 | pris-koh | K + OH cancel (≈ old `pristine-koh`) |
| `Vo-noK-noOH` | Vo | none | 0 | Vo-water | Vo baseline |
| `Vo-Knear-noOH` | Vo | near | 0 | Vo-koh | Vo + K alone |
| `Vo-Kfar-noOH` | Vo | far | 0 | Vo-koh | Vo + distance control |
| `Vo-Knear-OHonPt` | Vo | near | 2 | Vo-koh | Vo + cancel (≈ old `Vo-interface-koh`) |

### Tier B — isolate OH

| Label | Defect | K | OH@Pt | Parent | Role |
|-------|--------|---|-------|--------|------|
| `pris-noK-OHonPt` | pris | none | 1 | pris-water | OH penalty without K |
| `Vo-noK-OHonPt` | Vo | none | 1 | Vo-water | Vo + OH only |
| `pris-Kfar-OHonPt` | pris | far | 2 | pris-koh | OH penalty, K far |
| `Vo-Kfar-OHonPt` | Vo | far | 2 | Vo-koh | Vo + OH, K far |

## Key comparisons

| Question | Pair |
|----------|------|
| Does near K stabilize Pt electrons? | `pris-Knear-noOH` − `pris-noK-noOH` |
| Is it distance / electrostatic? | `pris-Knear-noOH` − `pris-Kfar-noOH` |
| Does OH@Pt cancel K? | `pris-Knear-OHonPt` − `pris-Knear-noOH` |
| OH alone? | `pris-noK-OHonPt` − `pris-noK-noOH` |
| Does Vo change the K effect? | `Vo-Knear-noOH` − `pris-Knear-noOH` (and same for noK baselines) |

Repeat the first four pairs with `Vo-*` labels.

## Built geometry (after `build_models.py`)

| Label | Formula | nK | OH@Pt | d(K–Pt) / Å |
|-------|---------|----|-------|-------------|
| pris-noK-noOH | H72O81Pt6Sr9Ti18 | 0 | 0 | — |
| pris-Knear-noOH | H66K2O79Pt6Sr9Ti18 | 2 | 0 | 3.16, 3.26 |
| pris-Kfar-noOH | H66K2O79Pt6Sr9Ti18 | 2 | 0 | 6.32, 6.49 |
| pris-Knear-OHonPt | H66K2O79Pt6Sr9Ti18 | 2 | 2 | 3.16, 3.26 |
| Vo-noK-noOH | H72O80Pt6Sr9Ti18 | 0 | 0 | — |
| Vo-Knear-noOH | H66K2O78Pt6Sr9Ti18 | 2 | 0 | 3.24, 3.37 |
| Vo-Kfar-noOH | H66K2O78Pt6Sr9Ti18 | 2 | 0 | 6.46, 6.59 |
| Vo-Knear-OHonPt | H66K2O78Pt6Sr9Ti18 | 2 | 2 | 3.24, 3.37 |
| pris-noK-OHonPt | H72O81Pt6Sr9Ti18 | 0 | 1 | — |
| Vo-noK-OHonPt | H72O80Pt6Sr9Ti18 | 0 | 1 | — |
| pris-Kfar-OHonPt | H66K2O79Pt6Sr9Ti18 | 2 | 2 | 8.14, 8.36 |
| Vo-Kfar-OHonPt | H66K2O78Pt6Sr9Ti18 | 2 | 2 | 8.30, 8.78 |

Vo formulas have **one fewer O** than the matching pris series (interface vacancy).

## Construction rules (charge / stoichiometry)

1. **Proton transfer (`noOH`)** — does not add/remove atoms. Each OH@Pt becomes H₂O; a distant H₂O becomes OH⁻. Formal charge unchanged.
2. **K ↔ far H₂O swap (`Kfar`)** — rigid translate water (O+2H) onto the old K site and place K at the old O site. Atom count unchanged.
3. **`noK`** — only from water parents (never delete K from koh, which would break neutrality unless OH is removed together).
4. **Bottom fix** — atoms with z < 4.05 Å fixed (same bottom TiO₂ as 31).

## Rationality notes / caveats

| Item | Assessment |
|------|------------|
| Knear vs Kfar separation | Clear: ~3.2 Å vs ≳6.3 Å (≳8 Å when far waters not used as H donors). |
| OH@Pt count | koh parents: **2**; water parents: **1**. Asymmetric but intentional (parents differ). For clean 1-vs-1 OH counts, compare within koh-derived or water-derived families. |
| H₂O still on Pt after `noOH` | Allowed: factor is **OH** (nH=1), not “no adsorbate”. |
| Lattice O near Pt (nH=0) | STO surface O; not counted as OH@Pt. |
| Kfar after `noOH` | Slightly closer (~6.3–6.6 Å) than Kfar+OHonPt (~8 Å) because protonation consumes the farthest waters first; still above the 6 Å cut. |
| Initial (pre-relax) clashes | Validator requires no pairs < 0.60 Å; GS relax will settle H placements. |
| Neutrality | GPAW `charge=0` GS; K⁺ balanced by solvent OH⁻ as in parents. |
| CDFT vs DFT | Same `relaxed.vasp`: CDFT constrains Δn on Pt₆; `dft_charge.py` adds total charge −Δn without constraint. |

Re-run checks anytime:

```bash
source env.sh
python build_models.py      # regenerate from 31 traj
python validate_models.py
```
