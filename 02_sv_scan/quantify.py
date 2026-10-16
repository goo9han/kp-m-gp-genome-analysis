#!/usr/bin/env python3
"""CHK1 deletion, read-level quantification (one BAM per call; results cached in work/q_<label>.json).
usage: quantify.py <label> [<label> ...]      (labels in LAB below)
kp1 coordinates (deleted interval D = scaffold_6:420,988-422,860, 1,873 bp; kp1 CHK1 CDS 420,835-422,547 (-)):
  median     genome median depth as ahp1_cnv/profile.py / sv_scan/depth_bins.py (samtools depth -a -Q 20, contigs >= 20 kb, every
             100th position); reused from the existing median file when there is one
  interior   mean depth over D (MAPQ >= 20 and MAPQ >= 0) / median;  flanks = 4 kb each side (416,988-420,987 and 422,861-426,860)
  reads_in_D reads (primary) with >= 50 % of their aligned bases inside D: count by MAPQ class, list of up to 50 (pos, MAPQ, NM, len)
  noise      genome-wide 1,873-bp windows (contigs >= 100 kb, 20 kb from ends), mean depth (MAPQ >= 20) / median -> percentiles
             (ahp1_cnv/window_noise.py with the window set to the deletion length)
  profile    100-bp bins 405,001-438,000 (MAPQ >= 20) / median, for the supplementary depth panels
  Illumina   junction pairs: unique (pos, mate pos) FR pairs, one mate starting in 417,500-420,988 on the + strand and the mate in
             422,860-426,000 on the - strand, MAPQ >= 20 (a pair that can only exist if D is missing: normal insert ~0.4 kb)
             split reads: primary MAPQ >= 20 with an SA segment, one part ending within 30 bp of 420,988 and the other starting within
             30 bp of 422,860
  ONT        delegated to reconcile_ont.score(), the project's single ONT read definition, so that this file
             and the OG0004482 / positive-control table cannot report different numbers for the same locus. It is the
             definition that was here, with the three window constants made relative to the interval: primary only
             (-F 0x904), MAPQ >= 20, overlapping the interval +/- 1,000 bp; 'junction' = a CIGAR deletion
             >= 0.8 x interval length inside interval +/- 300 or a split (SA) alignment whose parts end/start within
             150 bp of the two boundaries; 'intact' = one block covering interval +/- 100 with no deletion >= 500 bp;
             'other' otherwise (e.g. ends inside). Per-read spans kept for panel C.
kpgp3 coordinates (GP junction scaffold_48:52,670; mouse carries the 1.87-kb sequence here):
  continuous = reads (MAPQ >= 20) aligned without clip/indel >= 20 across 52,640-52,700; clipped = soft clip >= 20 bp within 15 bp of 52,670;
  ONT insertion reads = CIGAR insertion >= 1,500 bp within 300 bp of 52,670 (or split around it)"""
import os, re, sys, json, subprocess
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reconcile_ont            # the single ONT read definition
KP = "/Users/geongoohan/SynologyDrive/Brown/006_Fungal_WGS/kp_wgs_pangenome"; SIBA = "/path/to/ont_alignments"   # ONT alignments of the three isolates, produced
# outside this working tree.
HERE = os.path.dirname(os.path.abspath(__file__)); W = f"{HERE}/work"; os.makedirs(W, exist_ok=True)
SAM = os.path.expanduser("~/miniconda3/envs/buscoenv/bin/samtools")
K1 = f"{KP}/qpcr_primer_design/validation/genomes/kp1_combined.fa"; G3 = f"{KP}/qpcr_primer_design/validation/genomes/kpgp3_combined.fa"
C, D0, D1 = "scaffold_6", 420988, 422860; DL = D1 - D0 + 1
J3C, J3 = "scaffold_48", 52670
SVW = f"{KP}/sv_scan/work"
# label: (strain, tech, coord, bam, median file or None)
LAB = {"kp1_il": ("kp1", "illumina", "kp1", f"{KP}/qpcr_primer_design/redesign/pileup/kp1_illumina_combined.sorted.bam", f"{SVW}/median_kp1_self_k1.txt"),
       "kpgp3_il": ("kpgp3", "illumina", "kp1", f"{KP}/het_landscape/bams/kpgp3_on_kp1.bam", f"{SVW}/median_kpgp3_on_k1.txt"),
       "kpgp6_il": ("kpgp6", "illumina", "kp1", f"{KP}/het_landscape/bams/kpgp6_on_kp1.bam", f"{SVW}/median_kpgp6_on_k1.txt"),
       "kp1_ont": ("kp1", "ont", "kp1", f"{SIBA}/01_ont_boundary/bams/kp1_ont_on_kp1.bam", None),
       "kpgp3_ont": ("kpgp3", "ont", "kp1", f"{SIBA}/01_ont_boundary/bams/kpgp3_ont_on_kp1.bam", None),
       "kpgp6_ont": ("kpgp6", "ont", "kp1", f"{SIBA}/01_ont_boundary/bams/kpgp6_ont_on_kp1.bam", None),
       "kpmtx_ont": ("kpmtx", "ont", "kp1", f"{HERE}/bams/kpmtx_ont_on_kp1.bam", None),
       **{f"{s}_il": (s, "illumina", "kp1", f"{HERE}/bams/{s}_on_kp1.bam", None) for s in ("kpmtx", "C89-2", "ID1LA", "IDNYC", "QC1", "W26G")},
       # kpgp3 coordinates
       "kp1_il_g3": ("kp1", "illumina", "kpgp3", f"{KP}/ahp1_cnv/bams/kp1_on_kpgp3.bam", None),
       "kpmtx_il_g3": ("kpmtx", "illumina", "kpgp3", f"{KP}/ahp1_cnv/bams/kpmtx_on_kpgp3.bam", None),
       **{f"{s}_il_g3": (s, "illumina", "kpgp3", f"{KP}/ahp1_cnv/nine_strain/bams/{s}_on_kpgp3.bam", None) for s in ("C89-2", "ID1LA", "IDNYC", "QC1", "W26G")},
       "kpgp3_il_g3": ("kpgp3", "illumina", "kpgp3", f"{KP}/qpcr_primer_design/validation/pileup/kpgp3_illumina_combined.sorted.bam", None),
       "kpgp6_il_g3": ("kpgp6", "illumina", "kpgp3", f"{KP}/ahp1_cnv/bams/kpgp6_on_kpgp3.bam", None),
       "kp1_ont_g3": ("kp1", "ont", "kpgp3", f"{KP}/ahp1_cnv/nine_strain/bams/kp1_ont_on_kpgp3.bam", None),
       "kpmtx_ont_g3": ("kpmtx", "ont", "kpgp3", f"{KP}/ahp1_cnv/nine_strain/bams/kpmtx_ont_on_kpgp3.bam", None),
       "kpgp3_ont_g3": ("kpgp3", "ont", "kpgp3", f"{KP}/ahp1_cnv/bams/kpgp3_ont_on_kpgp3.bam", None),
       "kpgp6_ont_g3": ("kpgp6", "ont", "kpgp3", f"{KP}/ahp1_cnv/bams/kpgp6_ont_on_kpgp3.bam", None)}


def sh(*a): return subprocess.run(list(a), capture_output=True, text=True, check=True).stdout
def ops(cig): return [(int(n), o) for n, o in re.findall(r"(\d+)([MIDNSHP=X])", cig)]
def refspan(o): return sum(n for n, x in o if x in "MDN=X")


def median(bam, mf, ref):
    if mf and os.path.exists(mf): return float(open(mf).read())
    bed = f"{W}/_med_{os.path.basename(ref)}.bed"
    with open(bed, "w") as fh:
        for l in open(ref + ".fai"):
            c, L = l.split("\t")[:2]
            if int(L) >= 20000 and not c.endswith("_MT"): fh.write(f"{c}\t0\t{L}\n")
    p = subprocess.Popen([SAM, "depth", "-a", "-Q", "20", "-b", bed, bam], stdout=subprocess.PIPE, text=True)
    v = [int(l.split("\t")[2]) for i, l in enumerate(p.stdout) if i % 100 == 0]; p.wait()
    return float(np.median(v))


def mean_depth(bam, c, a, b, q):
    out = sh(SAM, "depth", "-a", "-Q", str(q), "-r", f"{c}:{a}-{b}", bam); v = [int(l.split("\t")[2]) for l in out.splitlines()]
    return float(np.mean(v)) if v else 0.0


def noise(bam, med, ref):
    bed = f"{W}/_noise_{os.path.basename(ref)}.bed"
    if not os.path.exists(bed):
        with open(bed, "w") as fh:
            for l in open(ref + ".fai"):
                c, L = l.split("\t")[0], int(l.split("\t")[1])
                if L < 100000 or c.endswith("_MT"): continue
                for s in range(20000, L - 20000 - DL, DL): fh.write(f"{c}\t{s}\t{s + DL}\n")
    out = sh(SAM, "bedcov", "-Q", "20", bed, bam)
    v = np.array([int(l.split("\t")[3]) / DL / med for l in out.splitlines()])
    ks = [1, 2.5, 5, 25, 50, 75, 95, 97.5, 99]
    return dict(n=int(len(v)), q={str(k): float(x) for k, x in zip(ks, np.percentile(v, ks))}, min=float(v.min()))


def kp1_coords(label, strain, tech, bam, mf):
    med = median(bam, mf, K1); r = dict(label=label, strain=strain, tech=tech, coord="kp1", median=med)
    r["interior_q20"] = mean_depth(bam, C, D0, D1, 20) / med; r["interior_q0"] = mean_depth(bam, C, D0, D1, 0) / med
    r["flank_left_q20"] = mean_depth(bam, C, D0 - 4000, D0 - 1, 20) / med; r["flank_right_q20"] = mean_depth(bam, C, D1 + 1, D1 + 4000, 20) / med
    r["noise"] = noise(bam, med, K1)
    prof = sh(SAM, "depth", "-a", "-Q", "20", "-r", f"{C}:405001-438000", bam)
    v = np.array([int(l.split("\t")[2]) for l in prof.splitlines()]); r["profile_start"] = 405001
    r["profile100"] = [float(x) for x in (v.reshape(-1, 100).mean(1) / med)] if len(v) == 33000 else []
    inside = []; n_in = [0, 0]; pairs = set(); split = 0
    for l in sh(SAM, "view", "-F", "0x904", bam, f"{C}:{D0 - 3600}-{D1 + 3200}").splitlines():
        f = l.split("\t"); flag, pos, mq = int(f[1]), int(f[3]), int(f[4]); o = ops(f[5]); end = pos + refspan(o) - 1
        nm = next((int(t[5:]) for t in f[11:] if t.startswith("NM:i:")), -1)
        # aligned bases inside D
        rp = pos; ins = 0; alg = 0
        for n, x in o:
            if x in "M=X":
                alg += n; ins += max(0, min(rp + n - 1, D1) - max(rp, D0) + 1)
            if x in "MDN=X": rp += n
        if alg and ins / alg >= 0.5:
            n_in[mq >= 20] += 1
            if len(inside) < 50: inside.append((pos, mq, nm, alg))
        sa = next((t[5:] for t in f[11:] if t.startswith("SA:Z:")), None)
        if tech == "illumina":
            if mq >= 20 and flag & 1 and not flag & 12 and f[6] == "=" and not flag & 16 and flag & 32:
                mp = int(f[7])
                if D0 - 3500 <= pos < D0 and D1 < mp <= D1 + 3200: pairs.add((pos, mp))
            if mq >= 20 and sa:
                for s in sa.strip(";").split(";"):
                    c2, p2, st2, cg2, mq2 = s.split(",")[:5]; e2 = int(p2) + refspan(ops(cg2)) - 1
                    if c2 == C and ((abs(end - (D0 - 1)) <= 30 and abs(int(p2) - (D1 + 1)) <= 30) or (abs(pos - (D1 + 1)) <= 30 and abs(e2 - (D0 - 1)) <= 30)): split += 1; break
    r["reads_in_D"] = dict(total=n_in[0] + n_in[1], mapq_ge20=n_in[1], mapq_lt20=n_in[0], examples=inside)
    if tech == "illumina": r["junction_pairs"] = len(pairs); r["pair_list"] = sorted(pairs)[:300]; r["split_reads"] = split
    else:
        sc = reconcile_ont.score(bam, (C, D0, D1))            # one definition, one implementation
        r["ont_reads"] = sc["reads"]
        r["ont_counts"] = {"junction": sc["junction"], "intact": sc["intact"], "other": sc["other"]}
        r["ont_definition"] = "chk1/reconcile_ont.py score(), PAD=1000"
        r["ont_cascade"] = sc["cascade"]; r["ont_mapq_median"] = sc["mapq_med"]
    return r


def g3_coords(label, strain, tech, bam):
    r = dict(label=label, strain=strain, tech=tech, coord="kpgp3", junction=f"{J3C}:{J3}"); cont = clip = insr = 0; spanning = 0
    for l in sh(SAM, "view", "-F", "0x904", "-q", "20", bam, f"{J3C}:{J3 - 30}-{J3 + 30}").splitlines():
        f = l.split("\t"); pos = int(f[3]); o = ops(f[5]); end = pos + refspan(o) - 1
        rp = pos; bad = False; ins_here = False
        if o[0][1] == "S" and o[0][0] >= 20 and abs(pos - J3) <= 15: clip += 1; continue
        if o[-1][1] == "S" and o[-1][0] >= 20 and abs(end - J3) <= 15: clip += 1; continue
        for n, x in o:
            if x in "ID" and n >= 20 and J3 - 300 <= rp <= J3 + 300:
                bad = True
                if x == "I" and n >= 1500: ins_here = True
            if x in "MDN=X": rp += n
        if pos <= J3 - 30 and end >= J3 + 30:
            spanning += 1
            if not bad: cont += 1
        if ins_here: insr += 1
    r.update(spanning=spanning, continuous=cont, clipped=clip, insertion_ge1500=insr)
    return r


if __name__ == "__main__":
    for lab in sys.argv[1:]:
        strain, tech, coord, bam, mf = LAB[lab]
        r = kp1_coords(lab, strain, tech, bam, mf) if coord == "kp1" else g3_coords(lab, strain, tech, bam)
        json.dump(r, open(f"{W}/q_{lab}.json", "w"), default=lambda x: int(x) if isinstance(x, np.integer) else float(x))
        print(lab, {k: v for k, v in r.items() if k not in ("profile100", "pair_list", "ont_reads", "noise", "reads_in_D")},
              "noise q1/q5/q95", [round(r["noise"]["q"][k], 3) for k in ("1", "5", "95")] if "noise" in r else "",
              "reads_in_D", {k: v for k, v in r.get("reads_in_D", {}).items() if k != "examples"}, flush=True)
