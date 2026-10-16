#!/usr/bin/env python3
"""Read extents of the ONT deletion reads, for drawing only. Same BAM, same filters and the same classification rule as
ahp1_cnv/ont_evidence.py (MAPQ >= 20, primary + supplementary, anchors 79,000 / 94,500); nothing is re-estimated.
Check: the (type, left breakpoint, right breakpoint) list must equal ahp1_cnv/ont_<strain>.json exactly, else abort.
Output: ont_reads_<strain>.json with, per deletion read, the aligned segments on scaffold_18."""
import sys, os, re, json, subprocess
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__)); CNV = os.path.dirname(HERE)
SAM = os.path.expanduser("~/miniconda3/envs/buscoenv/bin/samtools")
ctg, LA, RA = "scaffold_18", 79000, 94500
for strain in ("kpgp3", "kpgp6"):
    bam = f"{CNV}/bams/{strain}_ont_on_kpgp3.bam"
    out = subprocess.run([SAM, "view", "-F", "0x104", "-q", "20", bam, f"{ctg}:{LA - 20000}-{RA + 20000}"], capture_output=True, text=True).stdout
    aln = defaultdict(list)
    for ln in out.splitlines():
        f = ln.split("\t"); flag = int(f[1]); pos = int(f[3]); rp = pos; dels = []
        for n, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]):
            n = int(n)
            if op in "M=X": rp += n
            elif op in "DN":
                if n >= 1000: dels.append((rp, rp + n - 1))
                rp += n
        aln[f[0]].append(dict(start=pos, end=rp - 1, strand="-" if flag & 16 else "+", dels=dels))
    reads, bps = [], []
    for r, al in aln.items():
        if any(a["start"] <= LA and a["end"] >= RA and not any(d0 < RA and d1 > LA for d0, d1 in a["dels"]) for a in al): continue
        big = [(d0, d1) for a in al for d0, d1 in a["dels"] if d1 - d0 >= 5000 and d0 >= LA - 5000 and d1 <= RA + 5000]
        if big and any(a["start"] < LA for a in al) and any(a["end"] > RA for a in al):
            bps.append(["cigar", big[0][0], big[0][1]])
            a = [x for x in al if any(d == big[0] for d in x["dels"])][0]
            reads.append(dict(read=r, type="cigar", segments=[[a["start"], big[0][0] - 1], [big[0][1] + 1, a["end"]]], gap=list(big[0]))); continue
        left = [a for a in al if a["start"] < LA and LA - 3000 < a["end"] < RA - 3000]
        right = [a for a in al if a["end"] > RA and LA + 3000 < a["start"] < RA + 3000]
        if left and right and any(x["strand"] == y["strand"] for x in left for y in right):
            x = max(left, key=lambda a: a["end"]); y = min(right, key=lambda a: a["start"])
            if y["start"] - x["end"] >= 5000:
                bps.append(["split", x["end"], y["start"]])
                reads.append(dict(read=r, type="split", segments=[[x["start"], x["end"]], [y["start"], y["end"]]], gap=[x["end"], y["start"]]))
    ref = json.load(open(f"{CNV}/ont_{strain}.json"))["breakpoints"]
    ok = sorted(map(tuple, bps)) == sorted(map(tuple, ref))
    print(strain, "reads", len(reads), "matches ont_%s.json:" % strain, ok)
    if not ok: sys.exit(f"breakpoint list differs from ont_{strain}.json - stop")
    json.dump(dict(strain=strain, source_bam=bam, reads=reads), open(f"ont_reads_{strain}.json", "w"), indent=1)
