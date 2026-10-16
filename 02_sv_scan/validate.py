#!/usr/bin/env python3
"""Validate scan intervals with the same evidence types used for HOG1-AHP1 (ahp1_cnv_verdict.md), generalized to any interval:
  depth    : best (lowest-GP for deletions / highest-GP for gains) window inside the interval at its detecting scale:
             Kp-GP, Kp-GP2, the 7 mouse strains and the 4 ONT profiles (normalized to own median), GP percentile rank
  pairs    : Illumina FR read pairs (read MAPQ >= 20, same contig) with one end in the left part and the other in the right part
             of [start-5 kb, end+5 kb] and insert > 1 kb (deletion); RF pairs spanning the interval (tandem gain). Unique (pos, mate pos)
             per strain, 9 strains (illumina_evidence.py counts |TLEN| > 2 kb pairs; same idea, region-generic)
  ONT      : reads (MAPQ >= 20, primary + supplementary) skipping >= 30 % of the interval (min 500 bp) inside [start-5 kb, end+5 kb],
             as one CIGAR deletion or as two same-strand alignments (ont_evidence.py logic without its fixed 5-kb size)
  interchrom: Illumina pairs with the mate on another contig and ONT reads with a supplementary alignment on another contig
             (MAPQ >= 20), per strain, with the most frequent target (rearrangement joining this locus to another contig)
  selfblast: blastn -dust no of the interval against kpgp3_combined.fa; fraction of the interval covered by non-self HSPs
             at >= 98 % identity (a haplotig twin) and >= 90 % (repeat content) (selfblast_table.py logic)
  genes    : kpgp3 funannotate genes overlapping the interval with eggNOG preferred name / description (region_genes.py logic)
  boundary : blastn -dust no of the flanks around the two ends ([start-1 kb, start+h], [end-h, end+1 kb], h = min(1 kb, width/2),
             non-overlapping) against each other (shared repeat at both ends?), and overlap
             with copies of the 316-bp HOG1-AHP1 boundary repeat (repeats/copies_kpgp3.tsv)
Reverse-scan intervals (kp1 coordinates) get: kp1 genes, GP pairs on kp1 spanning the interval, and blastn of the kp1 sequence against
the GP assemblies (kpgp3, kpgp6) and the six other mouse assemblies (nuclear_clean). Output: work/validated_<cat>.json"""
import os, re, sys, json, subprocess, collections
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); R = os.path.dirname(HERE); W = f"{HERE}/work"
B = os.path.expanduser("~/miniconda3/envs/buscoenv/bin")
G = f"{R}/qpcr_primer_design/validation/genomes"
K3, K1 = f"{G}/kpgp3_combined.fa", f"{G}/kp1_combined.fa"
FA3 = "/Users/geongoohan/SynologyDrive/Brown/006_Fungal_WGS/kpgp3/results/funannotate/Kazachstania_pintolopesii"
FA1 = "/Users/geongoohan/SynologyDrive/Brown/006_Fungal_WGS/kp1/scaffolding_medaka_racon_pilon/funannotate/Kazachstania_pintolopesii"
NUC = f"{R}/rnaseq_prep/ani/nuclear_clean"
GP = ["kpgp3_self", "kpgp6_on3"]
MOUSE = ["kp1_on3", "kpmtx_on3", "C89-2_on3", "ID1LA_on3", "IDNYC_on3", "QC1_on3", "W26G_on3"]
ONT = ["kpgp3_ont", "kpgp6_ont", "kp1_ont", "kpmtx_ont"]
BAM = {"kpgp3_self": f"{R}/qpcr_primer_design/validation/pileup/kpgp3_illumina_combined.sorted.bam",
       "kpgp6_on3": f"{R}/ahp1_cnv/bams/kpgp6_on_kpgp3.bam", "kp1_on3": f"{R}/ahp1_cnv/bams/kp1_on_kpgp3.bam",
       "kpmtx_on3": f"{R}/ahp1_cnv/bams/kpmtx_on_kpgp3.bam",
       **{f"{s}_on3": f"{R}/ahp1_cnv/nine_strain/bams/{s}_on_kpgp3.bam" for s in ("C89-2", "ID1LA", "IDNYC", "QC1", "W26G")},
       "kpgp3_ont": f"{R}/ahp1_cnv/bams/kpgp3_ont_on_kpgp3.bam", "kpgp6_ont": f"{R}/ahp1_cnv/bams/kpgp6_ont_on_kpgp3.bam",
       "kp1_ont": f"{R}/ahp1_cnv/nine_strain/bams/kp1_ont_on_kpgp3.bam", "kpmtx_ont": f"{R}/ahp1_cnv/nine_strain/bams/kpmtx_ont_on_kpgp3.bam",
       "kpgp3_on_k1": f"{R}/het_landscape/bams/kpgp3_on_kp1.bam", "kpgp6_on_k1": f"{R}/het_landscape/bams/kpgp6_on_kp1.bam",
       "kp1_self_k1": f"{R}/qpcr_primer_design/redesign/pileup/kp1_illumina_combined.sorted.bam"}
PCT = json.load(open(f"{W}/percentiles.json"))


def sh(*a):
    return subprocess.run(list(a), capture_output=True, text=True).stdout


def load_scale(k, tag="k3"):
    z = np.load(f"{W}/scan_{tag}_{k}kb.npz")
    coords = [x.split(":") for x in z["coords"]]
    return [(c, int(s), int(e)) for c, s, e in coords], {l: z[l] for l in z.files if l != "coords"}


SC = {k: load_scale(k) for k in (1, 5, 9, 20)}


def best_window(c, s, e, scales, cat):
    best = None
    for k in scales:
        coords, V = SC[k]
        idx = [i for i, (cc, ss, ee) in enumerate(coords) if cc == c and ss >= s and ee <= e]
        if not idx: continue
        gpm = np.array([np.mean([V[g][i] for g in GP]) for i in idx])
        i = idx[int(np.argmin(gpm) if cat == "deletion" else np.argmax(gpm))]
        cand = (k, coords[i], {l: float(V[l][i]) for l in V})
        key = np.mean([cand[2][g] for g in GP])
        if best is None or (key < best[0] if cat == "deletion" else key > best[0]): best = (key, cand)
    k, (cc, ss, ee), vals = best[1]
    P = PCT[str(k)]
    # percentile rank of each GP value within its own reference distribution is approximated by the bracket it falls in
    return dict(scale_kb=k, window=f"{cc}:{ss + 1}-{ee}", depth={l: round(v, 3) for l, v in vals.items()},
                gp_p1=[round(P[g]["1"], 3) for g in GP], gp_p99=[round(P[g]["99"], 3) for g in GP],
                mouse_min=round(min(vals[m] for m in MOUSE), 3), mouse_max=round(max(vals[m] for m in MOUSE), 3))


def pairs(bam, c, s, e, kind):
    lo, hi = max(1, s - 5000), e + 5000; mid = (s + e) / 2
    out = sh(f"{B}/samtools", "view", "-F", "0xD04", "-q", "20", bam, f"{c}:{lo}-{int(mid)}")
    seen = set()
    for l in out.splitlines():
        f = l.split("\t"); flag = int(f[1]); pos = int(f[3])
        if f[6] != "=" or not (lo <= pos < mid): continue
        mpos = int(f[7])
        if not (mid < mpos <= hi): continue
        rrev, mrev = bool(flag & 16), bool(flag & 32)
        if kind == "deletion" and not rrev and mrev and mpos - pos > 1000: seen.add((pos, mpos))
        if kind == "gain" and rrev and not mrev and mpos - pos > 1000: seen.add((pos, mpos))
    return len(seen)


def interchrom(bam, c, s, e, ont=False):
    """reads in [start-5 kb, end+5 kb] (primary, MAPQ >= 20) whose mate (Illumina) or supplementary alignment (ONT, MAPQ >= 20)
    lies on another contig; returns (n reads, top target contig:kb, its count)"""
    out = sh(f"{B}/samtools", "view", "-F", "0x904", "-q", "20", bam, f"{c}:{max(1, s - 5000)}-{e + 5000}")
    tg = collections.Counter(); n = 0
    for l in out.splitlines():
        f = l.split("\t")
        if ont:
            sa = [t for t in f[11:] if t.startswith("SA:Z:")]
            hits = [x.split(",") for x in sa[0][5:].strip(";").split(";")] if sa else []
            hits = [(h[0], int(h[1])) for h in hits if h[0] != c and int(h[4]) >= 20]
        else:
            mq = [t for t in f[11:] if t.startswith("MQ:i:")]
            hits = [(f[6], int(f[7]))] if f[6] not in ("=", "*") and (not mq or int(mq[0][5:]) >= 20) else []
        if hits:
            n += 1
            for h in hits: tg[f"{h[0]}:{h[1] // 1000}kb"] += 1
    top = tg.most_common(1)
    return [n, top[0][0] if top else None, top[0][1] if top else 0]


def ont_skips(bam, c, s, e):
    width = e - s; need = max(500, 0.3 * width); lo, hi = max(1, s - 5000), e + 5000
    out = sh(f"{B}/samtools", "view", "-F", "0x104", "-q", "20", bam, f"{c}:{lo}-{hi}")
    al = collections.defaultdict(list); gaps = []
    for l in out.splitlines():
        f = l.split("\t"); flag = int(f[1]); rp = int(f[3]); segs = []; s0 = rp; big = []
        for n, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]):
            n = int(n)
            if op in "M=X": rp += n
            elif op in "DN":
                if n >= need: big.append((rp, rp + n - 1))
                rp += n
        al[f[0]].append((int(f[3]), rp - 1, "-" if flag & 16 else "+", big))
    n = 0; bps = collections.Counter()
    for r, A in al.items():
        hit = None
        for a0, a1, st, big in A:
            for d0, d1 in big:
                if d0 >= lo and d1 <= hi and d0 < e and d1 > s: hit = (d0, d1)
        if not hit:
            A2 = sorted(A)
            for x, y in zip(A2, A2[1:]):
                if x[2] == y[2] and y[0] - x[1] >= need and x[1] >= lo and y[0] <= hi and x[1] < e and y[0] > s: hit = (x[1], y[0])
        if hit: n += 1; bps[(round(hit[0], -2), round(hit[1], -2))] += 1
    return n, (bps.most_common(1)[0][0] if bps else None)


def selfblast(c, s, e):
    q = f"{W}/_q.fa"; open(q, "w").write(sh(f"{B}/samtools", "faidx", K3, f"{c}:{s + 1}-{e}"))
    out = sh(f"{B}/blastn", "-query", q, "-db", f"{G}/kpgp3_db", "-dust", "no", "-evalue", "1e-10", "-max_target_seqs", "500",
             "-max_hsps", "500", "-outfmt", "6 sseqid pident length qstart qend sstart send")
    cov98, cov90 = set(), set(); twins = collections.Counter()
    for l in out.splitlines():
        sid, pid, ln, q0, q1, s0, s1 = l.split("\t"); pid = float(pid); q0, q1, s0, s1 = map(int, (q0, q1, s0, s1))
        if sid == c and min(s0, s1) == s + q0 and max(s0, s1) == s + q1: continue
        rng = range(q0, q1 + 1)
        if pid >= 90: cov90.update(rng)
        if pid >= 98: cov98.update(rng); twins[sid] += q1 - q0 + 1
    L = e - s
    return dict(cov98=round(len(cov98) / L, 3), cov90=round(len(cov90) / L, 3), top98=twins.most_common(2))


EMP_CACHE = {}


def genes(prefix, c, s, e):
    if prefix not in EMP_CACHE:
        emp = {}
        for l in open(prefix + ".emapper.annotations"):
            if l.startswith("#"): continue
            f = l.rstrip("\n").split("\t"); emp[f[0].split("-T")[0]] = (f[8], f[7])
        gl = collections.defaultdict(list)
        for l in open(prefix + ".gff3"):
            if l.startswith("#"): continue
            f = l.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene": continue
            a = dict(x.split("=", 1) for x in f[8].strip(";").split(";") if "=" in x)
            gl[f[0]].append((int(f[3]), int(f[4]), a["ID"], a.get("Name", "")))
        EMP_CACHE[prefix] = (emp, gl)
    emp, gl = EMP_CACHE[prefix]
    out = []
    for gs, ge, gid, nm in gl.get(c, []):
        if ge < s or gs > e: continue
        pref, desc = emp.get(gid, ("", ""))
        lab = nm if nm and not nm.startswith("FUN_") else (pref if pref and pref != "-" else "")
        out.append(dict(id=gid, name=lab, desc=desc[:60], inside=bool(gs >= s and ge <= e)))
    return out


COPIES = [l.rstrip("\n").split("\t") for l in open(f"{HERE}/repeats/copies_kpgp3.tsv")][1:]


def boundary(c, s, e):
    a, b = f"{W}/_l.fa", f"{W}/_r.fa"
    h = min(1000, (e - s) // 2)   # non-overlapping flanks: [s-1 kb, s+h] and [e-h, e+1 kb]
    open(a, "w").write(sh(f"{B}/samtools", "faidx", K3, f"{c}:{max(1, s - 1000)}-{s + h}"))
    open(b, "w").write(sh(f"{B}/samtools", "faidx", K3, f"{c}:{e - h + 1}-{e + 1000}"))
    o = sh(f"{B}/blastn", "-query", a, "-subject", b, "-dust", "no", "-evalue", "1e-10", "-outfmt", "6 pident length")
    hs = sorted(((float(x.split()[0]), int(x.split()[1])) for x in o.splitlines() if int(x.split()[1]) >= 100), key=lambda t: -t[1])
    r316 = [f"{x[0]}:{x[1]}-{x[2]}" for x in COPIES if x[0] == c and (abs(int(x[1]) - s) < 3000 or abs(int(x[2]) - e) < 3000)]
    return dict(shared_flank_hsp=hs[0] if hs else None, rep316_copies_near_ends=r316)


def main(cats):
    for cat in cats:
        if cat == "reverse_gp_absent":
            iv = json.load(open(f"{W}/intervals_reverse_gp_absent.json"))
            z = np.load(f"{W}/scan_k1_1kb.npz"); coords = [x.split(":") for x in z["coords"]]
            res = []
            for r in iv:
                c, s, e = r["contig"], r["start0"], r["end"]
                idx = [i for i, (cc, ss, ee) in enumerate(coords) if cc == c and int(ss) >= s and int(ee) <= e]
                dep = {l: round(float(np.mean(z[l][idx])), 3) for l in ("kpgp3_on_k1", "kpgp6_on_k1", "kp1_self_k1")}
                q = f"{W}/_q1.fa"; open(q, "w").write(sh(f"{B}/samtools", "faidx", K1, f"{c}:{s + 1}-{e}"))
                pres = {}
                for a in ("kpgp3", "kpgp6", "kpmtx", "C89-2", "ID1LA", "IDNYC", "QC1", "W26G"):
                    o = sh(f"{B}/blastn", "-query", q, "-subject", f"{NUC}/{a}_nuclear_clean.fa", "-dust", "no", "-evalue", "1e-10",
                           "-outfmt", "6 qstart qend pident")
                    cov = set()
                    for x in o.splitlines():
                        q0, q1, pid = x.split("\t")
                        if float(pid) >= 90: cov.update(range(int(q0), int(q1) + 1))
                    pres[a] = round(len(cov) / (e - s), 3)
                res.append(dict(**r, depth=dep, assembly_coverage_ge90=pres,
                                gp_pairs_spanning={g: pairs(BAM[g], c, s, e, "deletion") for g in ("kpgp3_on_k1", "kpgp6_on_k1")},
                                genes=genes(FA1, c, s, e)))
                print(cat, c, s, e, dep, pres, flush=True)
        else:
            iv = json.load(open(f"{W}/intervals_{cat}.json"))
            res = []
            for r in iv:
                c, s, e = r["contig"], r["start0"], r["end"]
                rec = dict(**r)
                if cat in ("deletion", "gain"):
                    rec["best"] = best_window(c, s, e, r["scales_kb"], cat)
                    rec["pairs"] = {l: pairs(BAM[l], c, s, e, cat) for l in GP + MOUSE}
                    if cat == "deletion":
                        rec["ont"] = {l: ont_skips(BAM[l], c, s, e) for l in ONT}
                    rec["interchrom_pairs"] = {l: interchrom(BAM[l], c, s, e) for l in GP + MOUSE}
                    rec["interchrom_ont"] = {l: interchrom(BAM[l], c, s, e, ont=True) for l in ONT}
                    rec["selfblast"] = selfblast(c, s, e)
                    rec["boundary"] = boundary(c, s, e)
                else:   # mouse_unmappable: depth summary and mouse ONT reads that skip the interval
                    k = 1; coords, V = SC[k]
                    idx = [i for i, (cc, ss, ee) in enumerate(coords) if cc == c and ss >= s and ee <= e]
                    rec["depth"] = {l: round(float(np.mean(V[l][idx])), 3) for l in GP + MOUSE + ONT}
                    rec["mouse_ont_skip"] = {l: ont_skips(BAM[l], c, s, e)[0] for l in ("kp1_ont", "kpmtx_ont")}
                rec["genes"] = genes(FA3, c, s, e)
                res.append(rec)
                print(cat, c, s, e, flush=True)
        json.dump(res, open(f"{W}/validated_{cat}.json", "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1:] or ["deletion", "gain", "mouse_unmappable", "reverse_gp_absent"])
