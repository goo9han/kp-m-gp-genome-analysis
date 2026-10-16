#!/usr/bin/env python3
"""Define the informative positions from the two isolates of this study alone, so that the
other seven strains stay independent of the definition.

  informative : all three Kp-GP copies carry the same allele, that allele differs from the Kp-M base,
                and both are called.

No other strain and no external genome enters the definition, so Kp-MTX, the five public mouse-derived
strains and Kp-GP2 can all be scored as independent controls.

Contig set follows gp_sites.py so the comparison is like for like: nuclear contigs with allele calls,
with the mitochondrion-derived contig and the short HOG1-AHP1 contig excluded by name.

Writes sites_twoisolate.npz with, per contig, inf (bool), gp (GP allele), ref (Kp-M base), cal.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
KP = os.path.dirname(os.path.dirname(HERE))
B = f"{KP}/het_landscape/b_origin"
DROP = {"scaffold_17", "scaffold_18"}

ab = np.load(f"{B}/ab_kpgp3_kp1.npz")
ref = np.load(f"{B}/asm_kp1.npz")

out = {}
tot = dict(inf=0, cal=0, contigs=0)
for c in ref.files:
    if c.endswith("_MT") or c in DROP or f"A_{c}" not in ab.files:
        continue
    r = ref[c]; n = len(r)
    cut = lambda x: x[1:n + 1] if len(x) == n + 1 else x[:n]
    A, Bc = cut(ab[f"A_{c}"]), cut(ab[f"B_{c}"])
    cal = (A > 0) & (r > 0)
    hom = (A > 0) & (A == Bc)
    inf = hom & (A != r) & (r > 0)
    al = np.zeros(n, np.uint8); al[inf] = A[inf]
    out[f"inf_{c}"] = inf; out[f"gp_{c}"] = al; out[f"ref_{c}"] = r
    out[f"cal_{c}"] = cal.astype(np.uint8)
    tot["inf"] += int(inf.sum()); tot["cal"] += int(cal.sum()); tot["contigs"] += 1

np.savez_compressed(f"{HERE}/sites_twoisolate.npz", **out)
print(f"  contigs {tot['contigs']}   callable {tot['cal']:,} bp")
print(f"  informative positions {tot['inf']:,}")
