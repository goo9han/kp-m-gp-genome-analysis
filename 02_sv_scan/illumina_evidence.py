#!/usr/bin/env python3
"""Read-level evidence around a low-depth interval (Illumina, existing BAM).
usage: illumina_evidence.py <label> <bam> <contig> <start> <end> <dip_start> <dip_end>
Per 100-bp bin in [start, end]: reads (primary, MAPQ>=0), reads MAPQ>=20, soft-clip >= 20 bp breakpoints (MAPQ>=20),
reads with SA tag (split), pairs with mate on another contig or |TLEN| > 2000 (MAPQ>=20).
Positions with depth 0 (coverage gaps, MAPQ>=0) in the dip interval. Output: evidence_<label>.json"""
import sys, subprocess, re, json, os
from collections import Counter, defaultdict
label, bam, ctg, S0, E0, D0, D1 = sys.argv[1], sys.argv[2], sys.argv[3], *map(int, sys.argv[4:8])
SAM = os.path.expanduser("~/miniconda3/envs/buscoenv/bin/samtools")
out = subprocess.run([SAM, "view", "-F", "0x904", bam, f"{ctg}:{S0}-{E0}"], capture_output=True, text=True).stdout
b = lambda p: (p - S0) // 100
n_all, n_q = Counter(), Counter(); clip = Counter(); sa = Counter(); disc = Counter()
sa_targets = []; disc_targets = []
for ln in out.splitlines():
    f = ln.split("\t"); flag = int(f[1]); pos = int(f[3]); mq = int(f[4]); cig = f[5]
    if pos < S0 or pos > E0: continue
    n_all[b(pos)] += 1
    if mq >= 20: n_q[b(pos)] += 1
    ops = re.findall(r"(\d+)([MIDNSHP=X])", cig)
    rlen = sum(int(n) for n, o in ops if o in "MDN=X")
    if mq >= 20:
        if ops and ops[0][1] == "S" and int(ops[0][0]) >= 20: clip[pos] += 1
        if ops and ops[-1][1] == "S" and int(ops[-1][0]) >= 20: clip[pos + rlen - 1] += 1
        tags = f[11:]
        sat = [t for t in tags if t.startswith("SA:Z:")]
        if sat:
            sa[b(pos)] += 1
            for s in sat[0][5:].strip(";").split(";"):
                c2, p2 = s.split(",")[:2]; sa_targets.append((pos, c2, int(p2)))
        if flag & 1 and not flag & 8:
            if f[6] != "=" or abs(int(f[8])) > 2000:
                disc[b(pos)] += 1; disc_targets.append((pos, ctg if f[6] == "=" else f[6], int(f[7])))
# coverage gaps in the dip (any MAPQ)
d = subprocess.run([SAM, "depth", "-a", "-r", f"{ctg}:{D0}-{D1}", bam], capture_output=True, text=True).stdout
dep = [int(l.split("\t")[2]) for l in d.splitlines()]
dq = subprocess.run([SAM, "depth", "-a", "-Q", "20", "-r", f"{ctg}:{D0}-{D1}", bam], capture_output=True, text=True).stdout
depq = [int(l.split("\t")[2]) for l in dq.splitlines()]
res = dict(label=label, contig=ctg, window=[S0, E0], dip=[D0, D1],
           dip_min_depth_any=min(dep), dip_positions_depth0=sum(x == 0 for x in dep), dip_positions_depth_lt5_q20=sum(x < 5 for x in depq),
           dip_mean_depth_any=sum(dep) / len(dep), dip_mean_depth_q20=sum(depq) / len(depq),
           dip_frac_reads_mapq_lt20=1 - sum(n_q[k] for k in n_all if D0 <= S0 + 100 * k <= D1) / max(1, sum(n_all[k] for k in n_all if D0 <= S0 + 100 * k <= D1)),
           bins=[dict(start=S0 + 100 * k, reads=n_all[k], reads_q20=n_q[k], sa=sa[k], disc=disc[k]) for k in range((E0 - S0) // 100 + 1)],
           clip_positions=sorted(clip.items()), sa_targets=sa_targets, disc_targets=disc_targets)
json.dump(res, open(f"evidence_{label}.json", "w"))
print({k: v for k, v in res.items() if k not in ("bins", "clip_positions", "sa_targets", "disc_targets")})
top = sorted(clip.items(), key=lambda x: -x[1])[:12]
print("top soft-clip positions (pos, reads):", top)
print("SA reads:", len(sa_targets), "discordant:", len(disc_targets))
