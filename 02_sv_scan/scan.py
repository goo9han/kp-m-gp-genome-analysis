#!/usr/bin/env python3
"""Genome-wide scan for GP-specific copy-number changes on the kpgp3 assembly (primary scan) and for sequence missing from
GP on the kp1 assembly (reverse scan). Input: work/bins_<label>.tsv from depth_bins.py (1-kb bins, MAPQ >= 20, normalized
to each strain's own genome median). No read is remapped and no depth is recomputed here.

Windows: sliding windows of 1, 5, 9 and 20 kb (1-kb step) inside each contig; value = width-weighted mean of the 1-kb bins.
Reference distribution per strain and scale: windows on contigs >= 100 kb, > 20 kb from either end (the window_noise.py rule).
  GP-specific deletion   : kpgp3_self AND kpgp6_on3 < own 1st percentile, AND all 7 mouse >= own 5th percentile
  GP-specific gain       : kpgp3_self AND kpgp6_on3 > own 99th percentile, AND all 7 mouse <= own 95th percentile
  mouse-unmappable (separate category, reference bias): all 7 mouse < 0.15 AND both GP >= 0.30 (1-kb scale)
Qualifying windows are merged when they overlap or are separated by <= 1 kb; copy-class labels are never used. A single
qualifying window is kept (its width = the window size). Intervals from all scales are then united.
Reverse scan (kp1 coordinates): kpgp3_on_k1 AND kpgp6_on_k1 < 0.15 while ALL SEVEN mouse strains are >= their own
5th percentile (1-kb scale) - the same mouse condition as the primary deletion rule.
Outputs: work/scan_windows.npz, work/percentiles.json, work/intervals_*.json"""
import os, json, collections
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); W = f"{HERE}/work"
GP = ["kpgp3_self", "kpgp6_on3"]
MOUSE = ["kp1_on3", "kpmtx_on3", "C89-2_on3", "ID1LA_on3", "IDNYC_on3", "QC1_on3", "W26G_on3"]
ONT = ["kpgp3_ont", "kpgp6_ont", "kp1_ont", "kpmtx_ont"]
GP_K1 = ["kpgp3_on_k1", "kpgp6_on_k1"]
MOUSE_K1 = ["kp1_self_k1", "kpmtx_on_k1", "C89-2_on_k1", "ID1LA_on_k1", "IDNYC_on_k1", "QC1_on_k1", "W26G_on_k1"]
SCALES = [1, 5, 9, 20]
QS = [1, 5, 50, 95, 99]


def load(label):
    rows = [l.rstrip("\n").split("\t") for l in open(f"{W}/bins_{label}.tsv")]
    return [(r[0], int(r[1]), int(r[2]), float(r[3]), float(r[4])) for r in rows]


def windows(bins, k):
    """sliding k-bin windows per contig -> list of (contig, start0, end, value)"""
    by = collections.OrderedDict()
    for c, s, e, m, n in bins: by.setdefault(c, []).append((s, e, n))
    out = []
    for c, B in by.items():
        s = np.array([b[0] for b in B]); e = np.array([b[1] for b in B]); n = np.array([b[2] for b in B]); w = e - s
        cs = np.concatenate([[0], np.cumsum(n * w)]); cw = np.concatenate([[0], np.cumsum(w)])
        for i in range(0, max(1, len(B) - k + 1)):
            j = min(len(B), i + k)
            out.append((c, int(s[i]), int(e[j - 1]), float((cs[j] - cs[i]) / (cw[j] - cw[i]))))
    return out


def ref_mask(win, clen):
    return np.array([clen[c] >= 100000 and s >= 20000 and e <= clen[c] - 20000 for c, s, e, _ in win])


def merge(iv, gap=1000):
    iv = sorted(iv); out = []
    for c, s, e in iv:
        if out and out[-1][0] == c and s <= out[-1][2] + gap: out[-1][2] = max(out[-1][2], e)
        else: out.append([c, s, e])
    return [tuple(x) for x in out]


def scan(labels_gp, labels_mouse, extra, tag):
    data = {l: load(l) for l in labels_gp + labels_mouse + extra}
    clen = collections.Counter()
    for c, s, e, *_ in data[labels_gp[0]]: clen[c] = max(clen[c], e)
    res = {"windows": {}, "pct": {}}
    for k in SCALES:
        W_ = {l: windows(data[l], k) for l in data}
        coords = [(c, s, e) for c, s, e, _ in W_[labels_gp[0]]]
        for l in data: assert [(c, s, e) for c, s, e, _ in W_[l]] == coords, l
        V = {l: np.array([v for *_, v in W_[l]]) for l in data}
        m = ref_mask(W_[labels_gp[0]], clen)
        P = {l: dict(zip(map(str, QS), np.percentile(V[l][m], QS).tolist())) for l in data}
        res["windows"][k] = (coords, V); res["pct"][k] = P
    return res, clen


def main():
    out = {}
    # ---------- primary scan: kpgp3 coordinates ----------
    R, clen = scan(GP, MOUSE, ONT, "k3")
    json.dump({str(k): v for k, v in R["pct"].items()}, open(f"{W}/percentiles.json", "w"), indent=1)
    cats = {"deletion": [], "gain": [], "mouse_unmappable": []}
    for k in SCALES:
        coords, V = R["windows"][k]; P = R["pct"][k]
        dl = np.all([V[g] < P[g]["1"] for g in GP], axis=0) & np.all([V[mm] >= P[mm]["5"] for mm in MOUSE], axis=0)
        gn = np.all([V[g] > P[g]["99"] for g in GP], axis=0) & np.all([V[mm] <= P[mm]["95"] for mm in MOUSE], axis=0)
        for cat, mk in (("deletion", dl), ("gain", gn)):
            cats[cat] += [(coords[i][0], coords[i][1], coords[i][2], k) for i in np.where(mk)[0]]
        if k == 1:
            mu = np.all([V[mm] < 0.15 for mm in MOUSE], axis=0) & np.all([V[g] >= 0.30 for g in GP], axis=0)
            cats["mouse_unmappable"] += [(coords[i][0], coords[i][1], coords[i][2], k) for i in np.where(mu)[0]]
    for cat, wl in cats.items():
        iv = merge([(c, s, e) for c, s, e, _ in wl])
        rec = []
        for c, s, e in iv:
            sc = sorted({k for c2, s2, e2, k in wl if c2 == c and s2 < e and e2 > s})
            rec.append(dict(contig=c, start0=s, end=e, width=e - s, scales_kb=sc, contig_len=clen[c],
                            near_contig_end=bool(s < 5000 or e > clen[c] - 5000)))
        out[cat] = rec
        json.dump(rec, open(f"{W}/intervals_{cat}.json", "w"), indent=1)
        print(cat, len(rec), "intervals;", sum(r["width"] for r in rec), "bp")
    np.savez_compressed(f"{W}/scan_k3_1kb.npz", coords=np.array([f"{c}:{s}:{e}" for c, s, e in R["windows"][1][0]]),
                        **{l: v for l, v in R["windows"][1][1].items()})
    for k in (5, 9, 20):
        np.savez_compressed(f"{W}/scan_k3_{k}kb.npz", coords=np.array([f"{c}:{s}:{e}" for c, s, e in R["windows"][k][0]]),
                            **{l: v for l, v in R["windows"][k][1].items()})
    # ---------- reverse scan: kp1 coordinates ----------
    R1, clen1 = scan(GP_K1, MOUSE_K1, [], "k1")
    coords, V = R1["windows"][1]; P = R1["pct"][1]
    mk = np.all([V[g] < 0.15 for g in GP_K1], axis=0) & np.all([V[mm] >= P[mm]["5"] for mm in MOUSE_K1], axis=0)
    iv = merge([coords[i] for i in np.where(mk)[0]])
    rec = [dict(contig=c, start0=s, end=e, width=e - s, contig_len=clen1[c], near_contig_end=bool(s < 5000 or e > clen1[c] - 5000)) for c, s, e in iv]
    json.dump(rec, open(f"{W}/intervals_reverse_gp_absent.json", "w"), indent=1)
    json.dump({str(k): v for k, v in R1["pct"].items()}, open(f"{W}/percentiles_kp1.json", "w"), indent=1)
    np.savez_compressed(f"{W}/scan_k1_1kb.npz", coords=np.array([f"{c}:{s}:{e}" for c, s, e in coords]), **V)
    print("reverse (GP absent on kp1)", len(rec), "intervals;", sum(r["width"] for r in rec), "bp")


if __name__ == "__main__":
    main()
