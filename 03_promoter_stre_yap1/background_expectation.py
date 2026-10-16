#!/usr/bin/env python3
"""Background expectation for the promoter cis-element figure: how many promoter site changes between the two
clades a randomly chosen gene carries, under the same window, alignment, motif definitions and scanner.

Criteria, fixed before the observed values were read:
  denominator   251 background loci (stre_null_model, 300 random kp1 genes, seed 20260925, none in the stress set,
                survivors of the true-ATG rule) and the 70 stress loci. Same -700..-1 window, same MAFFT alignment.
  event         a clade-level site loss or gain: present in all seven mouse strains and neither guinea-pig strain,
                or the reverse. Counted per locus.
  two levels    site level   = every such site, no read filter. Available for STRE and YAP1, both sets.
                read-filtered = the guinea-pig allele fixed in both guinea-pig strains. Available for STRE in both
                sets (the status column of stre_events_filtered.tsv); not for YAP1 in the background set, which has
                no allele-fraction data, so the two levels are reported separately and never mixed.
  statement     where the eight drawn genes fall in that distribution. Nothing beyond "inside / outside the range".
Output: background_expectation.tsv
"""
import csv, json, os, collections
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
KP = os.path.abspath(os.path.join(HERE, "..", ".."))
NB, NS = 251, 70
DRAWN = ["AHP1", "GPX2", "SOD1", "SKN7", "HSP26", "CRH1", "CWP2", "HAP2"]
rd = lambda p: list(csv.DictReader([l for l in open(p) if not l.startswith("#")], delimiter="\t"))
ev = rd(f"{KP}/stre_null_model/stre_events_filtered.tsv")
ys, yb = rd(f"{HERE}/yap1_sites_stress.tsv"), rd(f"{HERE}/yap1_sites_background.tsv")
meta = json.load(open(f"{KP}/stre_null_model/work/stress_meta.json"))
name2id = {v["label"]: k for k, v in meta.items()}
vec = lambda c, n: np.array(list(c.values()) + [0] * (n - len(c)))

out = []
def level(tag, motifs, fixed_only):
    cs, cb = collections.Counter(), collections.Counter()
    if "STRE" in motifs:
        for r in ev:
            if fixed_only and r["status"] != "PASS":
                continue
            (cs if r["set"] == "STRESS" else cb)[r["kp1_gene"]] += 1
    if "YAP1" in motifs:
        for r in ys:
            if r["kind"] in ("loss", "gain"): cs[r["kp1_gene"]] += 1
        for r in yb:
            if r["kind"] in ("loss", "gain"): cb[r["kp1_gene"]] += 1
    vs, vb = vec(cs, NS), vec(cb, NB)
    out.append(dict(level=tag, motifs="+".join(motifs), read_filtered=str(fixed_only).lower(),
                    stress_events=vs.sum(), stress_loci=NS, stress_per_locus=f"{vs.sum()/NS:.4f}",
                    bg_events=vb.sum(), bg_loci=NB, bg_per_locus=f"{vb.sum()/NB:.4f}",
                    ratio=f"{(vs.sum()/NS)/(vb.sum()/NB):.2f}", gene="", n="", bg_ge_n_pct=""))
    return cs, vb

print("site level, no read filter")
for m in (["STRE"], ["YAP1"], ["STRE", "YAP1"]):
    cs, vb = level("site", m, False)
    print(f"  {'+'.join(m):<10} stress/locus {out[-1]['stress_per_locus']}  background/locus {out[-1]['bg_per_locus']}  ratio {out[-1]['ratio']}")
    if m == ["STRE", "YAP1"]:
        db = collections.Counter(vb.tolist())
        print("    background per-locus: " + ", ".join(f"{k}:{db[k]} ({db[k]/NB*100:.1f}%)" for k in sorted(db)))
        for g in DRAWN:
            n = cs.get(name2id[g], 0); pct = (vb >= n).mean() * 100
            out.append(dict(level="per_gene", motifs="STRE+YAP1", read_filtered="false", stress_events="", stress_loci="",
                            stress_per_locus="", bg_events="", bg_loci=NB, bg_per_locus="", ratio="",
                            gene=g, n=n, bg_ge_n_pct=f"{pct:.1f}"))
            print(f"    {g:<7} {n} change(s); {pct:.1f}% of background loci carry at least that many")
print("read-filtered (fixed in both guinea-pig strains), STRE only")
level("read_filtered", ["STRE"], True)
print(f"  STRE       stress/locus {out[-1]['stress_per_locus']}  background/locus {out[-1]['bg_per_locus']}  ratio {out[-1]['ratio']}")
cols = list(out[0])
with open(f"{HERE}/background_expectation.tsv", "w") as fh:
    fh.write("# Background expectation for the promoter cis-element figure. Criteria in background_expectation.py,\n"
             "# fixed before the observed values were read. ratio = stress rate / background rate.\n")
    w = csv.DictWriter(fh, cols, delimiter="\t"); w.writeheader(); w.writerows(out)
print("\nsaved background_expectation.tsv")
