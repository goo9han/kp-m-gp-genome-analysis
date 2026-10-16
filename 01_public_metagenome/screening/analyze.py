#!/usr/bin/env python3
"""PRJDB38973 selection against the deposited Kp-GP genome.

The counting set is stated positively rather than by subtraction. The deposited genome ships the nuclear
assembly and the mitochondrion as separate files, so "nuclear" is a file boundary and not a name test,
and only one condition is applied inside it.

  SPACE = contigs of at least 20 kb of the deposited nuclear assembly.

  why 20 kb   detection is a breadth test. It compares the number of occupied 1 kb bins with the number
              expected if the same reads were placed uniformly over one genome. Short contigs are
              enriched for rDNA and repeats, where reads from other organisms pile into few bins, which
              is the failure mode the test exists to catch.

Everything outside SPACE is tabulated, not discarded.

  mitochondrion   kept in the reference so mitochondrial reads are not forced onto nuclear contigs,
                  counted on its own because mitochondrial copy number does not track the nuclear genome.
  short contigs   nuclear sequence below 20 kb, counted on its own.

Kp-confident read = primary, MAPQ >= 20, aligned length >= 100 bp, identity (1 - NM / aligned) >= 0.97.
Detection rule, fixed before the analysis = at least 5 Kp-confident reads in SPACE in at least 5 distinct
1 kb bins, leaving out the recurrent conserved bin scaffold_34:117 kb.

Output: screen_v2.tsv, screen_v2.json
"""
import os, re, json, subprocess, collections, glob

HERE = os.path.dirname(os.path.abspath(__file__))
KP = os.path.dirname(os.path.dirname(HERE))
SAM = os.path.expanduser("~/miniconda3/envs/buscoenv/bin/samtools")
NUCFAI = f"{KP}/public_survey/rigp_ref/nuclear_only.fai"
ALLFAI = f"{KP}/public_survey/rigp_ref/rigp.fa.fai"

NUCLEAR = {c: int(L) for c, L in (l.split("\t")[:2] for l in open(NUCFAI))}
SPACE = {c: L for c, L in NUCLEAR.items() if L >= 20000}
SHORT = {c for c, L in NUCLEAR.items() if L < 20000}
MT = {c for c, _ in (l.split("\t")[:2] for l in open(ALLFAI))} - set(NUCLEAR)
GENOME = sum(SPACE.values())
NBINS = sum((L + 999) // 1000 for L in SPACE.values())
CONSERVED = {("scaffold_34", 117)}


def channel(bam):
    out = subprocess.run([SAM, "view", bam], capture_output=True, text=True).stdout.splitlines()
    bins = collections.Counter(); mt = short = 0
    for l in out:
        f = l.split("\t")
        if int(f[1]) & 0x900 or int(f[4]) < 20:
            continue
        ops = re.findall(r"(\d+)([MIDNSHP=X])", f[5])
        aln = sum(int(n) for n, o in ops if o in "M=X")
        nm = [int(t[5:]) for t in f[11:] if t.startswith("NM:i:")]
        if aln < 100 or not nm or 1 - nm[0] / aln < 0.97:
            continue
        c, pos = f[2], int(f[3])
        if c in SPACE:
            bins[(c, pos // 1000)] += 1
        elif c in MT:
            mt += 1
        elif c in SHORT:
            short += 1
    usable = {k: v for k, v in bins.items() if k not in CONSERVED}
    n = sum(usable.values())
    return dict(space_reads=sum(bins.values()), space_reads_excl_conserved=n,
                bins_excl_conserved=len(usable), conserved_bin_reads=sum(bins.values()) - n,
                mt_reads=mt, short_contig_reads=short, detected=bool(n >= 5 and len(usable) >= 5))


def spots(label):
    if label.startswith("SPIKE"):
        return 1000000 + int(re.search(r"n(\d+)$", label).group(1))
    m = re.search(r"Written (\d+) spots", open(f"{HERE}/logs/fqd_{label}.log").read())
    return int(m.group(1)) if m else None


if __name__ == "__main__":
    labels = sorted(os.path.basename(p)[:-3] for p in glob.glob(f"{HERE}/done/*.ok"))
    res = {}
    for l in labels:
        r = channel(f"{HERE}/bams/{l}_mapped.bam"); r["spots"] = spots(l)
        res[l] = r
    sp = [res[l]["space_reads_excl_conserved"] / int(re.search(r"n(\d+)$", l).group(1))
          for l in labels if l.startswith("SPIKE") and int(re.search(r"n(\d+)$", l).group(1)) >= 100]
    rate = sum(sp) / len(sp) if sp else None
    runs = [l for l in labels if not l.startswith("SPIKE")]
    det = [l for l in runs if res[l]["detected"]]
    json.dump(dict(space_contigs=len(SPACE), genome_bp=GENOME, nbins=NBINS,
                   spike_pass_rate=rate, result=res),
              open(f"{HERE}/screen_v2.json", "w"), indent=1)
    cols = ["label", "spots", "space_reads", "space_reads_excl_conserved", "bins_excl_conserved",
            "conserved_bin_reads", "mt_reads", "short_contig_reads", "detected"]
    with open(f"{HERE}/screen_v2.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for l in labels:
            fh.write("\t".join(str(res[l].get(c, l if c == "label" else "")) for c in cols) + "\n")
    print(f"  analysis space {len(SPACE)} contigs, {GENOME:,} bp, {NBINS:,} 1-kb bins")
    print(f"  {len(det)} of {len(runs)} runs detected")
    print(f"  spike recovery (n>=100) {rate:.3f} reads per spiked pair" if rate else "")
    for l in labels:
        if l.startswith("SPIKE"):
            r = res[l]; print(f"    {l:<18} Kp-confident {r['space_reads_excl_conserved']:>5}  "
                              f"bins {r['bins_excl_conserved']:>5}  detected {r['detected']}")
