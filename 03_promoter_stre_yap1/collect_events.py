#!/usr/bin/env python3
"""Collect clade-specific substitutions, lineage-specific substitutions, observed STRE events, and the
read-based filters (kpgp3/kpgp6 Illumina pileup) for the STRESS and BACKGROUND sets.
Writes work/events_<SET>.json. Definitions: analysis_plan.md."""
import json, os, subprocess, sys
from collections import defaultdict, Counter
from stre_core import *

SAMTOOLS = os.path.expanduser("~/miniconda3/envs/buscoenv/bin/samtools")
BAM = {"kpgp3": f"{ROOT}/qpcr_primer_design/validation/pileup/kpgp3_illumina_combined.sorted.bam",
       "kpgp6": f"{ROOT}/qpcr_primer_design/redesign/pileup/kpgp6_illumina_combined.sorted.bam"}
REF = {s: f"{ROOT}/rnaseq_prep/ani/nuclear_clean/{s}_nuclear_clean.fa" for s in GP}
INDEL_ADJ, HOMO_RUN, PROM_LEN = 5, 4, 700
LINEAGES = [(s,) for s in STRAINS] + [tuple(GP)]
bg = load_bg()
M = read_meme(f"{ROOT}/organic_acid_analysis/motifs_organic_acid.meme")
STRE = [PWM(n, *M[n], bg) for n in ("MSN2_STRE", "MSN4_STRE")]
SET = sys.argv[1]
if SET == "STRESS":
    meta = json.load(open(f"{HERE}/work/stress_meta.json"))
    alns = load_set(f"{HERE}/work/stress.fa", f"{HERE}/work/aln_stress")
    ranges = {g: {s: tuple(v) for s, v in m["ranges"].items()} for g, m in meta.items() if m["used"]}
    flags = {g: m.get("flags", []) for g, m in meta.items() if m["used"]}
    labels = {g: m["label"] for g, m in meta.items()}
else:
    alns = load_set(f"{HERE}/work/bg/BG_promoters_trueATG_9strains.fa", f"{HERE}/work/aln_bg")
    ranges = {g: {s: tuple(v) for s, v in json.load(open(f"{HERE}/work/bg/{g}.range.json")).items()} for g in alns}
    flags = {g: [] for g in alns}
    labels = {g: g for g in alns}
print(SET, "loci", len(alns))


def col2pos(row):
    m, i = {}, -1
    for c, ch in enumerate(row):
        if ch != "-":
            i += 1; m[c] = i
    return m


def genome_coord(g, s, i, base):
    contig, a, b, strand = ranges[g][s]
    pad = PROM_LEN - (b - a + 1)
    gi = i - pad
    if gi < 0:
        return None
    return (contig, a + gi, base) if strand == "+" else (contig, b - gi, base.translate(COMP))


def homopolymer(seq, i):
    """base i sits in or next to a run of >= HOMO_RUN identical bases"""
    for j in range(max(0, i - HOMO_RUN), min(len(seq), i + 2)):
        if len(set(seq[j:j + HOMO_RUN])) == 1 and len(seq[j:j + HOMO_RUN]) == HOMO_RUN and j <= i + 1 and j + HOMO_RUN >= i:
            return True
    return False


loci = {}
queries = defaultdict(set)           # strain -> {(contig, pos)}
for g, A in alns.items():
    L = len(A["kp1"])
    cols = [{s: A[s][c] for s in STRAINS} for c in range(L)]
    gapcol = [any(v == "-" for v in col.values()) for col in cols]
    near_gap = [any(gapcol[max(0, c - INDEL_ADJ):c + INDEL_ADJ + 1]) for c in range(L)]
    c2p = {s: col2pos(A[s]) for s in STRAINS}
    full = lambda col: all(v in "ACGT" for v in col.values())
    # ancestral (7-consensus) string for the clade test
    anc_cols = [c for c in range(L) if len({cols[c][s] for s in SEVEN}) == 1 and cols[c]["kp1"] in "ACGT"]
    anc = "".join(cols[c]["kp1"] for c in anc_cols)
    a_idx = {c: i for i, c in enumerate(anc_cols)}
    eligible = defaultdict(list)                    # ancestral base -> ancestral-string indices
    for c in anc_cols:
        if full(cols[c]):
            eligible[cols[c]["kp1"]].append(a_idx[c])
    clade = []
    for c in anc_cols:
        col = cols[c]
        if not full(col):
            continue
        X, Y = col["kp1"], col["kpgp3"]
        if col["kpgp6"] == Y and Y != X:
            ev = dict(col=c, ai=a_idx[c], X=X, Y=Y, kp1_pos=c2p["kp1"][c] - PROM_LEN, indel_adj=near_gap[c],
                      homopolymer=homopolymer(anc, a_idx[c]) or homopolymer(anc[:a_idx[c]] + Y + anc[a_idx[c] + 1:], a_idx[c]),
                      gcoord={s: genome_coord(g, s, c2p[s][c], Y) for s in GP})
            for s in GP:
                if ev["gcoord"][s]:
                    ct, p, _ = ev["gcoord"][s]
                    for d in (range(-5, 6) if ev["homopolymer"] else (0,)):
                        queries[s].add((ct, p + d))
            clade.append(ev)
    # lineage-specific substitutions and their STRE losses/gains (unfiltered, Test 4)
    lin = {}
    for Lg in LINEAGES:
        rest = [s for s in STRAINS if s not in Lg]
        acols = [c for c in range(L) if len({cols[c][s] for s in rest}) == 1 and cols[c][rest[0]] in "ACGT"]
        aseq = "".join(cols[c][rest[0]] for c in acols)
        subs = [(i, cols[c][Lg[0]]) for i, c in enumerate(acols) if full(cols[c]) and len({cols[c][s] for s in Lg}) == 1
                and cols[c][Lg[0]] != cols[c][rest[0]]]
        dseq = list(aseq)
        for i, y in subs:
            dseq[i] = y
        l, gn, _, _ = loss_gain(aseq, "".join(dseq), STRE, [i for i, _ in subs])
        lin["+".join(Lg)] = dict(n_sub=len(subs), loss=l, gain=gn, anc_sites=len(site_list(aseq, STRE)))
    # observed STRE sites on the actual sequences, clade events
    site_cols = {}
    for s in STRAINS:
        p2c = {v: k for k, v in c2p[s].items()}
        site_cols[s] = [(p2c[a], p2c[b], st) for a, b, st in site_list(A[s].replace("-", ""), STRE)]
    canon = []
    for s in STRAINS:
        for c0, c1, st in site_cols[s]:
            for cs in canon:
                if cs["strand"] == st and c0 <= cs["c1"] and c1 >= cs["c0"]:
                    cs["c0"], cs["c1"] = min(cs["c0"], c0), max(cs["c1"], c1); cs["in"].add(s); break
            else:
                canon.append(dict(c0=c0, c1=c1, strand=st, **{"in": {s}}))
    obs = []
    for cs in canon:
        cols_ = range(cs["c0"], cs["c1"] + 1)
        has_n = any("N" in cols[c].values() for c in cols_)
        in7, ingp = [s in cs["in"] for s in SEVEN], [s in cs["in"] for s in GP]
        kind = "loss" if all(in7) and not any(ingp) else "gain" if not any(in7) and all(ingp) else None
        if kind is None:
            continue
        cause_subs = [e for e in clade if cs["c0"] <= e["col"] <= cs["c1"]]
        indel = any(gapcol[c] for c in cols_)
        obs.append(dict(kind=kind, c0=cs["c0"], c1=cs["c1"], strand=cs["strand"], has_N=has_n, indel=indel,
                        kp1_pos=f"{c2p['kp1'].get(cs['c0'], 0) - PROM_LEN}..{c2p['kp1'].get(cs['c1'], 0) - PROM_LEN}",
                        seqs={s: A[s][cs["c0"]:cs["c1"] + 1] for s in STRAINS}, subs=[e["col"] for e in cause_subs],
                        gcoord={s: genome_coord(g, s, c2p[s][next(c for c in cols_ if c in c2p[s])], "N") for s in GP}))
        if indel:
            for s in GP:
                gc = obs[-1]["gcoord"][s]
                if gc:
                    for d in range(-5, cs["c1"] - cs["c0"] + 6):
                        queries[s].add((gc[0], gc[1] + d))
    loci[g] = dict(label=labels[g], flags=flags[g], anc=anc, eligible=eligible, clade=clade, lineages=lin, obs=obs,
                   anc_stre_sites=len(site_list(anc, STRE)))

# ---- pileup
pile = {}
for s in GP:
    bed = f"{HERE}/work/pile_{SET}_{s}.bed"
    with open(bed, "w") as fh:
        for ct, p in sorted(queries[s]):
            fh.write(f"{ct}\t{p - 1}\t{p}\n")
    out = subprocess.run([SAMTOOLS, "mpileup", "-B", "-q", "20", "-Q", "20", "-l", bed, "-f", REF[s], BAM[s]],
                         capture_output=True, text=True, check=True).stdout
    d = {}
    for line in out.splitlines():
        f = line.split("\t")
        ct, p, ref, depth, bases = f[0], int(f[1]), f[2].upper(), int(f[3]), f[4]
        cnt, i, nindel = Counter(), 0, 0
        while i < len(bases):
            ch = bases[i]
            if ch == "^":
                i += 2; continue
            if ch in "+-":
                j = i + 1
                while j < len(bases) and bases[j].isdigit():
                    j += 1
                n = int(bases[i + 1:j]); nindel += 1; i = j + n; continue
            if ch in ".,":
                cnt[ref] += 1
            elif ch.upper() in "ACGT":
                cnt[ch.upper()] += 1
            elif ch == "*":
                cnt["*"] += 1; nindel += 0
            i += 1
        d[(ct, p)] = dict(ref=ref, cnt=dict(cnt), n=sum(v for k, v in cnt.items() if k in "ACGT"), indel=nindel,
                          deleted=cnt.get("*", 0))
    pile[s] = d


def indel_frac(s, ct, p, lo=-5, hi=5):
    fr = []
    for q in range(p + lo, p + hi + 1):
        x = pile[s].get((ct, q))
        if x and x["n"] + x["deleted"] > 0:
            fr.append((x["indel"] + x["deleted"]) / (x["n"] + x["deleted"]))
    return max(fr) if fr else None


for g, D in loci.items():
    for e in D["clade"]:
        st = {}
        for s in GP:
            gc = e["gcoord"][s]
            if not gc:
                st[s] = dict(depth=0, af=None, anc_af=None, indel_frac=None); continue
            x = pile[s].get((gc[0], gc[1]))
            anc_base = e["X"] if ranges[g][s][3] == "+" else e["X"].translate(COMP)
            if not x or x["n"] == 0:
                st[s] = dict(depth=0, af=None, anc_af=None, indel_frac=None); continue
            st[s] = dict(depth=x["n"], af=x["cnt"].get(gc[2], 0) / x["n"], anc_af=x["cnt"].get(anc_base, 0) / x["n"],
                         ref_is_derived=x["ref"] == gc[2],
                         indel_frac=indel_frac(s, gc[0], gc[1]) if e["homopolymer"] else None)
        e["pile"] = st
        dep = [st[s]["depth"] for s in GP]
        af = [st[s]["af"] for s in GP]
        if min(dep) < 10 or any(a is None for a in af) or min(af) < 0.15:
            cls = "UNSUPPORTED"
        elif min(af) >= 0.90:
            cls = "FIXED"
        else:
            cls = "PARTIAL"
        if e["homopolymer"] and any((st[s]["indel_frac"] or 0) > 0.20 for s in GP):
            cls = "UNSUPPORTED"
            e["why"] = "homopolymer: >20% reads with an indel nearby"
        e["af_class"] = cls
    for o in D["obs"]:
        if o["indel"]:
            o["indel_frac"] = {s: indel_frac(s, o["gcoord"][s][0], o["gcoord"][s][1], -5, o["c1"] - o["c0"] + 5)
                               if o["gcoord"][s] else None for s in GP}
json.dump(loci, open(f"{HERE}/work/events_{SET}.json", "w"), default=str)
n = sum(len(D["clade"]) for D in loci.values())
print(SET, "clade substitutions", n, Counter(e["af_class"] for D in loci.values() for e in D["clade"]),
      "obs STRE events", Counter((o["kind"], o["indel"]) for D in loci.values() for o in D["obs"]))
