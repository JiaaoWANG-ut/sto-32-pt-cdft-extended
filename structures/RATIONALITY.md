# Modeling rationality check (pre-relax)

Generated from `build_models.py` + `validate_models.py` on the initial `POSCAR.init` set.

## Pass / fail summary

All **12 / 12** models passed automated checks:

- Factor tags match geometry (nK, OH@Pt, d(K–Pt) cuts)
- Stoichiometry identical to 31 parent (edits are atom-conserving)
- No H₃O⁺ / O with >2 H
- No interatomic pairs < 0.60 Å
- Vacuum gap above adsorbates > 5 Å
- Bottom 27 atoms fixed (z < 4.05 Å)

## Design OK

1. **Causal K test is clean** — `Knear-noOH` vs `Kfar-noOH` vs `noK-noOH` share the same formal composition within a defect series; only K placement / presence changes (plus the controlled OH→far-OH proton hops for `noOH`).
2. **OH cancellation is separable** — `Knear-OHonPt` vs `Knear-noOH` keeps K near and toggles only whether OH sits on Pt.
3. **Vo is a true O vacancy** — Vo formulas have one fewer O than pris partners; Sr/Ti/Pt counts match.
4. **Neutrality** — noK from water parents; K₂ from koh parents. No orphaned K⁺ without compensating OH⁻ inventory.

## Acceptable asymmetries (document, do not “fix” by force)

| Issue | Why OK | How to compare fairly |
|-------|--------|------------------------|
| Water parents have **1** OH@Pt; koh parents have **2** | Inherited from 31 GS | Compare OH effects inside the same parent family |
| `Kfar` after `noOH` is ~6.3–6.6 Å, not ~8 Å | Farthest waters used as H donors first | Still ≥ 6 Å cut; both far vs ~3.2 Å near |
| `noOH` may leave **H₂O** on Pt | Factor is OH (nH=1), not bare Pt | State explicitly in paper/talk |
| Proton-transfer H placement is approximate | GS relax will refine | Require post-relax re-validation of OH@Pt / d(K–Pt) |

## Post-relax checklist (required before CDFT)

After `relax.py` finishes, re-run a geometry audit on `relaxed.vasp`:

1. `nOH_on_Pt` still matches the label (`noOH` → 0; `OHonPt` → ≥1)
2. `Knear` still d < 4.5 Å; `Kfar` still d > 6 Å (K must not drift back to Pt during neutral GS opt — if it does, constrain K or note failure)
3. No new clashes; vacuum gap intact
4. Then proceed to `pipeline.py` / `dft_charge.py`

## Residual risk

Neutral GS opt of `Kfar-*` may pull K⁺ toward Pt⁶ if the electrostatic drive is strong even at Q≈0. If post-relax d(K–Pt) collapses below ~5 Å, either (a) freeze K xy during GS, or (b) treat that label as “K wants to be near” evidence and rely on fixed-geometry CDFT from the **pre-relax** `Kfar` POSCAR for the distance control (document which geometry was used).
