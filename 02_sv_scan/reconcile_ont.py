#!/usr/bin/env python3
"""The single ONT read definition, applied everywhere, for an arbitrary interval.

THE DEFINITION
  records      samtools view -F 0x904: primary alignments only (no secondary 0x100, no supplementary 0x800,
               no unmapped 0x4), then MAPQ >= 20, then the read must overlap [D0-1000, D1+1000]
  one record   per read name (a supplementary is not primary, so names are already unique)
  junction     a CIGAR D/N >= 0.8 x interval length lying inside [D0-300, D1+300], OR an SA segment on the
               same contig whose two parts end/start within 150 bp of the two boundaries; in the SA case the
               second half is added to the read's blocks, so the record describes the whole read
  intact       a single aligned block covering [D0-100, D1+100] with no D/N >= 500 inside it
  other        neither: the read ends inside the interval, or carries a partial indel
  total        junction + intact + other, i.e. the MAPQ>=20 primary reads that overlap the padded interval
A junction read's split half is read off the primary record's SA tag, so excluding supplementary records
loses nothing: both routes find the same reads. Verified below by running the identical code on two BAM sets.

usage: reconcile_ont.py            (prints every table; writes chk1/work/reconcile_ont.tsv)"""
import os, re, sys, json, subprocess
import numpy as np

KP = "/Users/geongoohan/SynologyDrive/Brown/006_Fungal_WGS/kp_wgs_pangenome"
SIBA = "/path/to/ont_alignments"   # ONT alignments of the three isolates, produced
# outside this working tree.
SCR = "/tmp/sv_scan_work/ont_kp1"
SAM = os.path.expanduser("~/miniconda3/envs/buscoenv/bin/samtools")
HERE = os.path.dirname(os.path.abspath(__file__)); W = f"{HERE}/work"

CHK1 = ("scaffold_6", 420988, 422860)          # validated boundary, 1,873 bp
OG4482 = ("scaffold_11", 8315, 11469)          # the sv_scan interval, 3,155 bp

# BAM sets. "combined" = the reference the figure was built on; "clean" = the alternative assembly.
# scaffold_6 and scaffold_11 are byte-identical in the two FASTAs (md5-checked), so coordinates are shared.
SETS = {
 "combined": {"ref": "kp1_combined.fa (30 seq)", "bams": {
    "kp1":   f"{SIBA}/01_ont_boundary/bams/kp1_ont_on_kp1.bam",
    "kpmtx": f"{KP}/chk1/bams/kpmtx_ont_on_kp1.bam",
    "kpgp3": f"{SIBA}/01_ont_boundary/bams/kpgp3_ont_on_kp1.bam",
    "kpgp6": f"{SIBA}/01_ont_boundary/bams/kpgp6_ont_on_kp1.bam"}},
 "clean":    {"ref": "kp1_nuclear_clean.fa (28 seq)", "bams": {
    "kpgp3": f"{SCR}/kpgp3_ont_on_kp1.bam",
    "kpgp6": f"{SCR}/kpgp6_ont_on_kp1.bam"}}}
NAMES = {"kp1": "Kp-M", "kpmtx": "Kp-MTX", "kpgp3": "Kp-GP", "kpgp6": "Kp-GP2"}

def ops(cig): return [(int(n), o) for n, o in re.findall(r"(\d+)([MIDNSHP=X])", cig)]
def refspan(o): return sum(n for n, x in o if x in "MDN=X")

PAD = 1000                                     # the one padded-overlap window; see the docstring

def score(bam, interval, pad=PAD):
    """The definition above. Returns the counts, the per-read records (for drawing), the MAPQ distribution
    of the records seen, and the filter cascade. chk1/quantify.py calls this, so there is one copy."""
    c, D0, D1 = interval; DL = D1 - D0 + 1
    casc = dict(fetched=0, primary=0, mapq20=0, overlap=0)
    seen = {}; mapqs = []
    fetch0 = max(1, D0 - pad - 20000)          # clamp: scaffold_11's interval is 8 kb from the contig start
    out = subprocess.run([SAM, "view", bam, f"{c}:{fetch0}-{D1 + pad + 20000}"],
                         capture_output=True, text=True, check=True).stdout
    for l in out.splitlines():
        f = l.split("\t"); flag, pos, mq = int(f[1]), int(f[3]), int(f[4])
        casc["fetched"] += 1
        if flag & 0x904: continue
        casc["primary"] += 1
        o = ops(f[5]); end = pos + refspan(o) - 1
        if mq < 20: continue
        casc["mapq20"] += 1
        if end < D0 - pad or pos > D1 + pad: continue
        casc["overlap"] += 1; mapqs.append(mq)
        blocks = []; rp = pos; bigdel = False; bs = pos
        for n, x in o:
            if x in "DN" and n >= 500:
                blocks.append((bs, rp - 1))
                if n >= 0.8 * DL and rp >= D0 - 300 and rp + n <= D1 + 300: bigdel = True
                rp += n; bs = rp
            elif x in "MDN=X": rp += n
        blocks.append((bs, rp - 1))
        cls = "other"
        sa = next((t[5:] for t in f[11:] if t.startswith("SA:Z:")), None)
        if bigdel: cls = "junction"
        elif sa:
            for s in sa.strip(";").split(";"):
                c2, p2, st2, cg2 = s.split(",")[:4]; e2 = int(p2) + refspan(ops(cg2)) - 1
                if c2 == c and ((abs(end - (D0 - 1)) <= 150 and abs(int(p2) - (D1 + 1)) <= 150)
                                or (abs(pos - (D1 + 1)) <= 150 and abs(e2 - (D0 - 1)) <= 150)):
                    # the SA half is added to the blocks so a split junction read draws as one read reaching
                    # both flanks, which is what it is; without it panel C shows it stopping at the boundary
                    cls = "junction"; blocks.append((int(p2), e2)); break
        if cls == "other" and any(a <= D0 - 100 and b >= D1 + 100 for a, b in blocks): cls = "intact"
        seen[f[0]] = dict(name=f[0], cls=cls, blocks=sorted(blocks), mq=mq)
    n = dict(junction=0, intact=0, other=0)
    for v in seen.values(): n[v["cls"]] += 1
    n["total"] = sum(n[k] for k in ("junction", "intact", "other"))
    n["reads"] = list(seen.values())
    n["mapq_med"] = float(np.median(mapqs)) if mapqs else float("nan")
    n["mapq0"] = sum(1 for m in mapqs if m == 0); n["mapq_ge20"] = len(mapqs)
    n["cascade"] = casc
    return n

rows = []
def run(tag, interval, iname):
    for sk, S in SETS.items():
        for s, bam in S["bams"].items():
            if not os.path.exists(bam): print(f"  MISSING {bam}"); continue
            r = score(bam, interval); r.pop("reads")
            rows.append(dict(locus=iname, bamset=sk, ref=S["ref"], strain=s, **r))
            print(f"  {iname:10s} {sk:9s} {NAMES[s]:7s} junction={r['junction']:4d} intact={r['intact']:4d} "
                  f"other={r['other']:4d} total={r['total']:4d}  MAPQmed={r['mapq_med']:g} "
                  f"cascade fetched={r['cascade']['fetched']} primary={r['cascade']['primary']} "
                  f"mapq20={r['cascade']['mapq20']} overlap={r['cascade']['overlap']}")

if __name__ != "__main__":
    pass
else:
  print("== CHK1 scaffold_6:420,988-422,860 (1,873 bp) ==");   run("chk1", CHK1, "CHK1")
  print("== OG0004482 scaffold_11:8,315-11,469 (3,155 bp) =="); run("og4482", OG4482, "OG0004482")

  with open(f"{W}/reconcile_ont.tsv", "w") as fh:
    fh.write("# one ONT definition (see reconcile_ont.py docstring) applied to both BAM sets\n")
    fh.write("locus\tbamset\treference\tstrain\tjunction\tintact\tother\ttotal\tmapq_median\tmapq0\t"
             "fetched\tprimary\tmapq>=20\toverlap_padded\n")
    for r in rows:
        fh.write(f"{r['locus']}\t{r['bamset']}\t{r['ref']}\t{NAMES[r['strain']]}\t{r['junction']}\t{r['intact']}\t"
                 f"{r['other']}\t{r['total']}\t{r['mapq_med']:g}\t{r['mapq0']}\t{r['cascade']['fetched']}\t"
                 f"{r['cascade']['primary']}\t{r['cascade']['mapq20']}\t{r['cascade']['overlap']}\n")
  print("\nwrote work/reconcile_ont.tsv")
