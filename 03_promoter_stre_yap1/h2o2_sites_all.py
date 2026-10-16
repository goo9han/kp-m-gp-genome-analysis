#!/usr/bin/env python3
"""Count every YAP1 / STRE (MSN2, MSN4) site per strain in each H2O2-set promoter (not only the best hit).

Merge rule = ahp1_yap1_analysis/ahp1_step0_checks.py:
  within one strain, one family (YAP1; STRE = MSN2 or MSN4) and one strand, FIMO hits whose intervals
  share >=1 bp -> one site (union span); a + and a - site with identical span -> one palindromic site
  ("+/-"); sites in different strains are the same site when their MAFFT columns overlap (same family
  and strand).
Presence per strain: 1 = hit at p<1e-3; 0 = aligned columns present, no hit; NA = columns are N (not
assembled) or all gap. For every 0 the best p of a relaxed scan (--thresh 1, same background) over
the same columns and strand is reported, to separate a sequence change from a near-threshold call.

usage: h2o2_sites_all.py <set>   set = new (true-ATG promoters) | old (meme_fimo/all_promoters.fa)
"""
import os, sys, subprocess
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
H2O2 = os.path.dirname(HERE)
FIMO = os.path.expanduser("~/miniconda3/envs/memeenv/bin/fimo")
MAFFT = os.path.expanduser("~/miniconda3/envs/buscoenv/bin/mafft")
STRAINS = ["kp1", "kpmtx", "C89-2", "ID1LA", "IDNYC", "QC1", "W26G", "kpgp3", "kpgp6"]
FAMILY = {"MSN2": "STRE", "MSN4": "STRE", "YAP1": "YAP1"}
PROM_LEN = 700
rc = lambda s: s.translate(str.maketrans("ACGTN-", "TGCAN-"))[::-1]

SET = sys.argv[1]
if SET == "new":
    FA, OUT = f"{HERE}/H2O2_promoters_trueATG_9strains.fa", f"{HERE}/new"
    from_atg = lambda p: p - (PROM_LEN + 1)
    COORD = "from ATG"
else:
    FA, OUT = f"{H2O2}/meme_fimo/all_promoters.fa", f"{HERE}/old"
    from_atg = lambda p: p                       # old windows are not ATG anchored: keep promoter positions
    COORD = "old promoter pos"
os.makedirs(OUT, exist_ok=True)
split = lambda name: (name.split("__")[0], name.split("__")[1].replace("kp_", ""))


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        elif name is not None:
            seqs[name].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}


prom = defaultdict(dict)
for k, v in read_fasta(FA).items():
    gname, s = split(k)
    prom[gname][s] = v
GENES = list(prom)

# ---- FIMO (strict + relaxed) and per-gene MAFFT
for tag, th in (("fimo_out", "1e-3"), ("fimo_out_relaxed", "1.0")):
    if not os.path.exists(f"{OUT}/{tag}/fimo.tsv"):
        subprocess.run([FIMO, "--thresh", th, "--bfile", f"{H2O2}/meme_fimo/bg.txt", "--oc", f"{OUT}/{tag}",
                        f"{H2O2}/meme_fimo/motifs.meme", FA], check=True, capture_output=True)
os.makedirs(f"{OUT}/aln", exist_ok=True)
aln = {}
for gname in GENES:
    fa, af = f"{OUT}/aln/{gname}.fa", f"{OUT}/aln/{gname}.aln.fa"
    if not os.path.exists(af):
        open(fa, "w").write("".join(f">{gname}__{s}\n{prom[gname][s]}\n" for s in STRAINS))
        open(af, "w").write(subprocess.run([MAFFT, "--auto", "--quiet", fa], capture_output=True, text=True, check=True).stdout)
    aln[gname] = {split(k)[1]: v for k, v in read_fasta(af).items()}
    for s in STRAINS:
        assert aln[gname][s].replace("-", "") == prom[gname][s], (gname, s)


def load(path):
    rows = []
    for line in open(path):
        if line.startswith("#") or line.startswith("motif_id") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        gname, s = split(f[2])
        rows.append(dict(gene=gname, strain=s, pwm=f[1], family=FAMILY[f[1]], start=int(f[3]), stop=int(f[4]),
                         strand=f[5], score=float(f[6]), pvalue=float(f[7]), matched=f[9]))
    return rows


raw, relaxed = load(f"{OUT}/fimo_out/fimo.tsv"), load(f"{OUT}/fimo_out_relaxed/fimo.tsv")


def same_strand(a, b):
    """+/- (palindromic, hit on both strands) is the same site as a one-strand hit at the same columns."""
    return a == b or "+/-" in (a, b)


def variant_type(ref, alt):
    """Classify the aligned site columns of a strain without a hit against the carrier (ref)."""
    keep = [i for i in range(len(ref)) if not (ref[i] == "-" and alt[i] == "-")]   # gap columns made by other strains
    ref, alt = "".join(ref[i] for i in keep), "".join(alt[i] for i in keep)
    if ref == alt:
        return "identical seq (score/threshold)"
    kind = "indel" if "-" in ref + alt else "subst"
    if kind == "indel":
        # homopolymer: the gap sits in a run of >=4 identical bases in either row
        for row in (ref, alt):
            u = row.replace("-", "")
            if any(b * 4 in u for b in "ACGT"):
                return "indel in homopolymer (check reads)"
    return kind


def pos2col(a):
    m, i = {}, 0
    for c, ch in enumerate(a):
        if ch != "-":
            i += 1; m[i] = c
    return m


site_rows, count_rows = [], []
for gname in GENES:
    A = aln[gname]
    p2c = {s: pos2col(A[s]) for s in STRAINS}
    c2kp1 = {c: p for p, c in p2c["kp1"].items()}
    sites = {s: [] for s in STRAINS}
    for s in STRAINS:
        for fam in ("YAP1", "STRE"):
            per_strand = []
            for strand in ("+", "-"):
                hits = sorted([r for r in raw if r["gene"] == gname and r["strain"] == s and r["family"] == fam
                               and r["strand"] == strand], key=lambda r: r["start"])
                cur = None
                for h in hits:
                    if cur and h["start"] <= cur["stop"]:
                        cur["stop"] = max(cur["stop"], h["stop"]); cur["hits"].append(h)
                    else:
                        cur = dict(start=h["start"], stop=h["stop"], strand=strand, family=fam, hits=[h]); per_strand.append(cur)
            merged = []
            for st in per_strand:
                twin = next((m for m in merged if (m["start"], m["stop"]) == (st["start"], st["stop"]) and m["strand"] != st["strand"]), None)
                if twin:
                    twin["strand"] = "+/-"; twin["hits"] += st["hits"]
                else:
                    merged.append(st)
            sites[s] += merged
    intervals = sorted((p2c[s][st["start"]], p2c[s][st["stop"]], st["family"], st["strand"], s, st) for s in STRAINS for st in sites[s])
    canon = []
    for c0, c1, fam, strand, s, st in intervals:
        for cs in canon:
            if cs["family"] == fam and same_strand(cs["strand"], strand) and c0 <= cs["c1"] and c1 >= cs["c0"]:
                cs["c0"] = min(cs["c0"], c0); cs["c1"] = max(cs["c1"], c1); cs["members"][s] = st
                if cs["strand"] != strand:
                    cs["strand"] = "+/-" if "+/-" in (cs["strand"], strand) else cs["strand"]
                break
        else:
            canon.append(dict(c0=c0, c1=c1, family=fam, strand=strand, members={s: st}))
    canon.sort(key=lambda x: (x["c0"], x["family"]))

    def kp1_at(col, step):
        while col not in c2kp1:
            col += step
        return c2kp1[col]

    for i, cs in enumerate(canon, 1):
        lo, hi = kp1_at(cs["c0"], 1), kp1_at(cs["c1"], -1)
        ref = "kp1" if "kp1" in cs["members"] else next(s for s in STRAINS if s in cs["members"])
        st = cs["members"][ref]
        seq = prom[gname][ref][st["start"] - 1:st["stop"]]
        pres, relax = {}, {}
        for s in STRAINS:
            if s in cs["members"]:
                pres[s] = "1"
                continue
            cols = A[s][cs["c0"]:cs["c1"] + 1]
            if "N" in cols:
                pres[s] = "NA"
                continue
            pres[s] = "0"
            c2p = {c: p for p, c in p2c[s].items()}
            ps = [c2p[c] for c in range(cs["c0"], cs["c1"] + 1) if c in c2p]
            fam_strands = ("+", "-") if cs["strand"] == "+/-" else (cs["strand"],)
            ref = "kp1" if "kp1" in cs["members"] else next(x for x in STRAINS if x in cs["members"])
            cand = [r for r in relaxed if r["gene"] == gname and r["strain"] == s and r["family"] == cs["family"]
                    and r["strand"] in fam_strands and ps and r["start"] <= max(ps) and r["stop"] >= min(ps)]
            b = min(cand, key=lambda r: r["pvalue"]) if cand else None
            vt = variant_type(A[ref][cs["c0"]:cs["c1"] + 1], cols)
            relax[s] = (f"{b['pwm']} {b['matched']} p={b['pvalue']:.1e}" if b else "no overlapping hit") + f" [{vt}]"
        best = min((h for m in cs["members"].values() for h in m["hits"]), key=lambda h: h["pvalue"])
        n1 = sum(v == "1" for v in pres.values())
        site_rows.append(dict(gene=gname, site=f"S{i}", family=cs["family"], strand=cs["strand"], kp1_start=lo, kp1_stop=hi,
                              ref_pos=f"{from_atg(st['start'])}..{from_atg(st['stop'])}",
                              pos=f"{from_atg(lo)}..{from_atg(hi)}", seq=seq if cs["strand"] != "-" else rc(seq), ref=ref,
                              pwms=",".join(sorted({h["pwm"] for m in cs["members"].values() for h in m["hits"]})),
                              best_p=best["pvalue"], pres=pres, n=n1, n_na=sum(v == "NA" for v in pres.values()),
                              relax=relax, aln={s: A[s][cs["c0"]:cs["c1"] + 1] for s in STRAINS}))
    for s in STRAINS:
        count_rows.append((gname, s, sum(st["family"] == "YAP1" for st in sites[s]), sum(st["family"] == "STRE" for st in sites[s])))

with open(f"{OUT}/sites_all.tsv", "w") as fh:
    fh.write("\t".join(["gene", "site", "family", "strand", f"kp1_{COORD.replace(' ', '_')}", "kp1_promoter_pos", "matched_seq",
                        "seq_from", f"seq_from_{COORD.replace(' ', '_')}", "pwms", "best_p"] + STRAINS + ["n_present", "n_NA", "absent_strains_relaxed_best"]) + "\n")
    for r in site_rows:
        fh.write("\t".join(map(str, [r["gene"], r["site"], r["family"], r["strand"], r["pos"], f"{r['kp1_start']}-{r['kp1_stop']}",
                                     r["seq"], r["ref"], r["ref_pos"], r["pwms"], f"{r['best_p']:.2e}"] + [r["pres"][s] for s in STRAINS]
                                    + [r["n"], r["n_na"], "; ".join(f"{s}: {v}" for s, v in r["relax"].items())])) + "\n")
with open(f"{OUT}/site_counts_per_strain.tsv", "w") as fh:
    fh.write("gene\t" + "\t".join(f"{s}_YAP1\t{s}_STRE" for s in STRAINS) + "\n")
    for gname in GENES:
        d = {s: (y, t) for g_, s, y, t in count_rows if g_ == gname}
        fh.write(gname + "\t" + "\t".join(f"{d[s][0]}\t{d[s][1]}" for s in STRAINS) + "\n")

# ---- console summary
for gname in GENES:
    rows = [r for r in site_rows if r["gene"] == gname]
    d = {s: (y, t) for g_, s, y, t in count_rows if g_ == gname}
    print(f"## {gname}: {sum(r['family']=='YAP1' for r in rows)} YAP1 + {sum(r['family']=='STRE' for r in rows)} STRE sites; "
          f"per strain Y/S: " + " ".join(f"{s}={d[s][0]}/{d[s][1]}" for s in STRAINS))
    for r in rows:
        if r["n"] + r["n_na"] < 9:
            pat = "".join({"1": "+", "0": ".", "NA": "?"}[r["pres"][s]] for s in STRAINS)
            print(f"   VAR {r['site']} {r['family']} {r['strand']} {COORD} {r['pos']} ({r['ref']} {r['ref_pos']}) {r['seq']} best p {r['best_p']:.1e} "
                  f"[{pat}] absent: " + "; ".join(f"{s}: {v}" for s, v in r["relax"].items()))
            for s in STRAINS:
                print(f"        {s:6s} {r['aln'][s]}")
