#!/usr/bin/env python3
"""YAP1 (MA0415.2) site scan of the STRESS and BACKGROUND promoter sets, for the background expectation
of the promoter cis-element figure.

Why this exists: the STRE scan already covers both sets (stre_null_model), but YAP1 was only ever scanned
on the 14 H2O2 loci, so four of the eight drawn genes had no YAP1 call and no YAP1 background existed.

Nothing new is aligned or re-extracted. It reuses, unchanged:
  promoters   stre_null_model/work/stress.fa, work/bg/BG_promoters_trueATG_9strains.fa  (true-ATG, -700..-1)
  alignments  stre_null_model/work/aln_stress, work/aln_bg                              (MAFFT, cached)
  scorer      stre_null_model/stre_core.py                                              (FIMO-equivalent:
              scorer_validation.txt reports 3085/3085 and 7314/7314 hits, 0 differences)
  matrix      stre_null_model/JASPAR2024_CORE_fungi_non-redundant_pfms_meme.txt, MOTIF MA0415.2 YAP1
Site calling and the cross-strain merge are the same code path as stre_turnover_figure/build_data.py:
p < 1e-3, same-strand overlapping hits merged, sites equated across strains by MAFFT column overlap.

Calibration gate: on the stress set the per-strain presence pattern must reproduce every YAP1 row of
h2o2_analysis/true_atg_all/new/sites_all.tsv (the existing FIMO-based table). The script aborts if not.

Output: yap1_sites_stress.tsv, yap1_sites_background.tsv
"""
import csv, json, os, re, sys
KP = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
KP = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, f"{KP}/stre_null_model")
from stre_core import *                                   # noqa  (STRAINS, read_fasta, load_bg, read_meme, PWM, site_list, load_set)

OUT = os.path.dirname(os.path.abspath(__file__))
HN = f"{KP}/stre_null_model"
SEVEN = ["kp1", "kpmtx", "C89-2", "ID1LA", "IDNYC", "QC1", "W26G"]
GP = ["kpgp3", "kpgp6"]
PROM = 700

bg = load_bg()
M = read_meme(f"{HN}/JASPAR2024_CORE_fungi_non-redundant_pfms_meme.txt")
YAP1 = [PWM("YAP1", *M["MA0415.2"], bg)]
print(f"YAP1 MA0415.2: width {YAP1[0].w}, score cutoff {YAP1[0].cut:.4f} at p < 1e-3")


def col2pos(row):
    m, i = {}, -1
    for c, ch in enumerate(row):
        if ch != "-":
            i += 1; m[c] = i
    return m


def scan(alns, labels):
    """Site-level clade classification per locus, on alignment columns."""
    sites = []
    for g, A in alns.items():
        c2p = {s: col2pos(A[s]) for s in STRAINS}
        site_cols = {}
        for s in STRAINS:
            p2c = {v: k for k, v in c2p[s].items()}
            site_cols[s] = [(p2c[a], p2c[b], st) for a, b, st in site_list(A[s].replace("-", ""), YAP1)]
        canon = []
        for s in STRAINS:                                  # identical merge order to build_data.py
            for c0, c1, st in site_cols[s]:
                for cs in canon:
                    if cs["strand"] == st and c0 <= cs["c1"] and c1 >= cs["c0"]:
                        cs["c0"], cs["c1"] = min(cs["c0"], c0), max(cs["c1"], c1); cs["in"].add(s); break
                else:
                    canon.append(dict(c0=c0, c1=c1, strand=st, **{"in": {s}}))
        # a + and a - canonical site with the same span is one palindromic site ("+/-"),
        # the merge rule documented in h2o2_analysis/true_atg_all/h2o2_sites_all.py
        byspan = {}
        for cs in canon:
            byspan.setdefault((cs["c0"], cs["c1"]), []).append(cs)
        canon = []
        for (c0, c1), group in byspan.items():
            strands = {c["strand"] for c in group}
            canon.append(dict(c0=c0, c1=c1, strand="+/-" if strands == {"+", "-"} else group[0]["strand"],
                              **{"in": set().union(*[c["in"] for c in group])}))
        for cs in canon:
            in7 = [s in cs["in"] for s in SEVEN]; ingp = [s in cs["in"] for s in GP]
            kind = ("loss" if all(in7) and not any(ingp) else "gain" if not any(in7) and all(ingp) else
                    "shared" if all(in7) and all(ingp) else "other")
            pos = f"{c2p['kp1'].get(cs['c0'], 0) - PROM}..{c2p['kp1'].get(cs['c1'], 0) - PROM}"
            sites.append(dict(kp1_gene=g, locus=labels.get(g, g), motif="YAP1", kind=kind, pos_from_ATG=pos,
                              strand=cs["strand"], present_in=",".join(s for s in STRAINS if s in cs["in"]),
                              n_present=len(cs["in"]), seq_kp1=A["kp1"][cs["c0"]:cs["c1"] + 1],
                              seq_kpgp3=A["kpgp3"][cs["c0"]:cs["c1"] + 1]))
    return sites


meta = json.load(open(f"{HN}/work/stress_meta.json"))
labels = {k: v["label"] for k, v in meta.items()}
print("loading alignments (cached) ...")
aln_s = load_set(f"{HN}/work/stress.fa", f"{HN}/work/aln_stress")
aln_b = load_set(f"{HN}/work/bg/BG_promoters_trueATG_9strains.fa", f"{HN}/work/aln_bg")
print(f"  stress {len(aln_s)} loci, background {len(aln_b)} loci")
s_sites = scan(aln_s, labels)
b_sites = scan(aln_b, {})
print(f"  YAP1 sites: stress {len(s_sites)}, background {len(b_sites)}")

# ---------------- calibration gate ----------------
name2id = {v: k for k, v in labels.items()}
tgt = [r for r in csv.DictReader(open(f"{KP}/h2o2_analysis/true_atg_all/new/sites_all.tsv"), delimiter="\t")
       if r["family"] == "YAP1"]
byg = {}
for s in s_sites:
    byg.setdefault(s["kp1_gene"], []).append(s)
ok = bad = skip = 0
for r in tgt:
    gid = name2id.get(r["gene"])
    if gid not in byg and gid not in aln_s:
        skip += 1; print(f"  calibration skip {r['gene']}: not in the stress set"); continue
    lo, hi = [int(x) for x in re.match(r"(-\d+)\.\.(-\d+)", r["kp1_from_ATG"]).groups()]
    cand = [s for s in byg.get(gid, []) if not (int(s["pos_from_ATG"].split("..")[1]) < lo - 6
                                                or int(s["pos_from_ATG"].split("..")[0]) > hi + 6)]
    mine = set().union(*[set(s["present_in"].split(",")) for s in cand]) if cand else set()
    theirs = {st for st in STRAINS if r[st] == "1"}
    if mine == theirs:
        ok += 1
    else:
        bad += 1
        print(f"  CALIBRATION DIFF {r['gene']} {r['kp1_from_ATG']}: mine {sorted(mine)} vs sites_all {sorted(theirs)}")
print(f"CALIBRATION: {ok} agree, {bad} differ, {skip} skipped (of {len(tgt)} YAP1 rows in sites_all.tsv)")
if bad:
    sys.exit("ABORT: the YAP1 scan does not reproduce the existing FIMO table; do not use these calls.")

for tag, rows in (("stress", s_sites), ("background", b_sites)):
    cols = ["kp1_gene", "locus", "motif", "kind", "pos_from_ATG", "strand", "present_in", "n_present", "seq_kp1", "seq_kpgp3"]
    with open(f"{OUT}/yap1_sites_{tag}.tsv", "w") as fh:
        fh.write(f"# YAP1 (MA0415.2) sites in the {tag} promoter set, -700..-1 from the true ATG, p<1e-3,\n"
                 "# same-strand overlap merge, sites equated across strains by MAFFT column overlap.\n"
                 "# kind: loss = present in all 7 mouse and neither guinea-pig strain; gain = the reverse;\n"
                 "#       shared = all 9; other = any pattern that is not clade-level. Built by scan_yap1.py.\n")
        w = csv.DictWriter(fh, cols, delimiter="\t"); w.writeheader(); w.writerows(rows)
    import collections
    print(f"  {tag}: " + ", ".join(f"{k} {v}" for k, v in sorted(collections.Counter(r["kind"] for r in rows).items())))
print("saved")
