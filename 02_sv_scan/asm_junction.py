#!/usr/bin/env python3
"""For reverse-scan intervals (kp1 sequence missing from GP) with GP assembly coverage < 0.6: blastn kp1 [start-6 kb, end+6 kb]
against the kpgp3 and kpgp6 assemblies (-dust no) and find two collinear HSPs on one GP contig that flank a gap in the kp1 query.
The gap = the kp1 segment missing from that GP assembly (junction coordinates on the GP contig reported). Also reports whether the
two kp1 ends of the missing segment share sequence (repeat at both ends; blastn of +-300 bp around each end).
Output: work/asm_junctions.json"""
import os, json, subprocess, collections
HERE = os.path.dirname(os.path.abspath(__file__)); R = os.path.dirname(HERE); W = f"{HERE}/work"
B = os.path.expanduser("~/miniconda3/envs/buscoenv/bin")
K1 = f"{R}/qpcr_primer_design/validation/genomes/kp1_combined.fa"
NUC = f"{R}/rnaseq_prep/ani/nuclear_clean"
FL = 6000


def sh(*a): return subprocess.run(list(a), capture_output=True, text=True).stdout


out = []
for r in json.load(open(f"{W}/polarized_reverse_gp_absent.json")):
    if min(r["assembly_coverage_ge90"][g] for g in ("kpgp3", "kpgp6")) >= 0.6: continue
    c, s, e = r["contig"], r["start0"], r["end"]
    q0 = max(1, s - FL + 1); q1 = min(r["contig_len"], e + FL)
    open(f"{W}/_j.fa", "w").write(sh(f"{B}/samtools", "faidx", K1, f"{c}:{q0}-{q1}"))
    rec = dict(interval=f"{c}:{s + 1}-{e}", query=f"{c}:{q0}-{q1}", outgroup_cov80=r["outgroup_cov80"], junctions={})
    for a in ("kpgp3", "kpgp6"):
        hs = [l.split("\t") for l in sh(f"{B}/blastn", "-query", f"{W}/_j.fa", "-subject", f"{NUC}/{a}_nuclear_clean.fa", "-dust", "no",
                                        "-evalue", "1e-20", "-outfmt", "6 sseqid pident length qstart qend sstart send").splitlines()]
        hs = [(h[0], float(h[1]), int(h[2]), int(h[3]), int(h[4]), int(h[5]), int(h[6])) for h in hs if int(h[2]) >= 500]
        best = None
        by = collections.defaultdict(list)
        for h in hs: by[h[0]].append(h)
        for sid, H in by.items():
            H.sort(key=lambda h: h[3])
            for x in H:
                for y in H:
                    if y[3] <= x[4]: continue
                    gap = y[3] - x[4] - 1
                    # the gap must cover the scan interval; subject ends must be adjacent (within 50 bp) and collinear
                    if not (x[4] + q0 - 1 <= s + 200 and y[3] + q0 - 1 >= e - 200): continue
                    sdist = min(abs(y[5] - x[6]), abs(y[6] - x[5]))
                    if sdist > 50: continue
                    cand = (gap, sid, x, y)
                    if best is None or cand[0] < best[0]: best = cand
        if best:
            gap, sid, x, y = best
            rec["junctions"][a] = dict(gp_contig=sid, junction_on_gp=x[6], kp1_missing=f"{c}:{x[4] + q0}-{y[3] + q0 - 2}", missing_bp=gap,
                                       flank_identity=[x[1], y[1]])
        else:
            rec["junctions"][a] = None
    j = rec["junctions"].get("kpgp3") or rec["junctions"].get("kpgp6")
    if j:
        a, b = [int(v) for v in j["kp1_missing"].split(":")[1].split("-")]
        open(f"{W}/_e1.fa", "w").write(sh(f"{B}/samtools", "faidx", K1, f"{c}:{max(1, a - 300)}-{a + 300}"))
        open(f"{W}/_e2.fa", "w").write(sh(f"{B}/samtools", "faidx", K1, f"{c}:{b - 300}-{b + 300}"))
        o = sh(f"{B}/blastn", "-query", f"{W}/_e1.fa", "-subject", f"{W}/_e2.fa", "-dust", "no", "-evalue", "1e-5", "-outfmt", "6 pident length")
        rec["repeat_at_both_ends"] = [tuple(x.split("\t")) for x in o.splitlines()][:2]
    out.append(rec); print(json.dumps(rec))
json.dump(out, open(f"{W}/asm_junctions.json", "w"), indent=1)
