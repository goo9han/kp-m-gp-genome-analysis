"""Shared code for the STRE null model: PWM scorer (FIMO-equivalent), alignments, events."""
import itertools, json, os, re, subprocess
from collections import defaultdict, Counter
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MAFFT = os.path.expanduser("~/miniconda3/envs/buscoenv/bin/mafft")
STRAINS = ["kp1", "kpmtx", "C89-2", "ID1LA", "IDNYC", "QC1", "W26G", "kpgp3", "kpgp6"]
SEVEN, GP = STRAINS[:7], ["kpgp3", "kpgp6"]
IDX = {b: i for i, b in enumerate("ACGT")}
COMP = str.maketrans("ACGTN-", "TGCAN-")
rc = lambda s: s.translate(COMP)[::-1]
P_THRESH = 1e-3
BG = None


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        elif name is not None:
            seqs[name].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}


def load_bg(path=f"{ROOT}/h2o2_analysis/meme_fimo/bg.txt"):
    bg = {}
    for line in open(path):
        t = line.split()
        if len(t) == 2 and t[0] in "ACGT" and len(t[0]) == 1:
            bg[t[0]] = float(t[1])
    s = sum(bg.values())
    return np.array([bg[b] / s for b in "ACGT"])


def read_meme(path):
    mats, name, nsites, rows = {}, None, None, None
    for line in open(path):
        t = line.split()
        if line.startswith("MOTIF"):
            name = t[1]; rows = []; mats[name] = dict(rows=rows, nsites=20)
        elif "letter-probability" in line and name:
            m = re.search(r"nsites=\s*([\d.]+)", line)
            if m:
                mats[name]["nsites"] = float(m.group(1))
        elif name and len(t) == 4:
            try:
                rows.append([float(x) for x in t])
            except ValueError:
                pass
    return {k: (np.array(v["rows"]), v["nsites"]) for k, v in mats.items()}


class PWM:
    """FIMO scoring: p' = (p*nsites + 0.1*bg)/(nsites + 0.1); score = sum log2(p'/bg).
    p-value of a score: exact, by enumerating the score distribution of all k-mers under the 0-order bg."""

    def __init__(self, name, probs, nsites, bg, pseudo=0.1):
        self.name, self.w = name, probs.shape[0]
        p = (probs * nsites + pseudo * bg) / (nsites + pseudo)
        self.lo = np.log2(p / bg)                        # w x 4
        # exact score distribution
        scores = np.zeros(1); probs_ = np.ones(1)
        for i in range(self.w):
            scores = (scores[:, None] + self.lo[i][None, :]).ravel()
            probs_ = (probs_[:, None] * bg[None, :]).ravel()
        # p(score) = P(S >= score), ties included: group k-mers by (rounded) score
        key = np.round(scores, 6)
        uniq, inv = np.unique(key, return_inverse=True)
        mass = np.bincount(inv, weights=probs_)
        uniq, mass = uniq[::-1], mass[::-1]              # descending score
        cum = np.cumsum(mass)
        ok = cum < P_THRESH
        self.cut = uniq[ok].min() if ok.any() else np.inf
        self._s_sorted, self._cum = uniq, cum

    def pvalue(self, score):
        i = np.searchsorted(-self._s_sorted, -score + 1e-5, side="right") - 1
        return self._cum[i] if i >= 0 else 1.0

    def hits(self, seq):
        """[(start0, end0_inclusive, strand, score)] with p < P_THRESH on both strands."""
        out, w = [], self.w
        for strand, s in (("+", seq), ("-", rc(seq))):
            L = len(s)
            for i in range(L - w + 1):
                k = s[i:i + w]
                if any(c not in IDX for c in k):
                    continue
                sc = sum(self.lo[j, IDX[c]] for j, c in enumerate(k))
                if sc >= self.cut - 1e-5:
                    a = i if strand == "+" else L - i - w
                    out.append((a, a + w - 1, strand, sc))
        return out


def site_list(seq, pwms):
    out = []
    for strand in "+-":
        iv = sorted((a, b) for p in pwms for a, b, st, _ in p.hits(seq) if st == strand)
        cur = None
        for a, b in iv:
            if cur and a <= cur[1]:
                cur[1] = max(cur[1], b)
            else:
                cur = [a, b]; out.append((strand, cur))
    return [(c[0], c[1], s) for s, c in out]


def window_hits(seq, pwms, lo, hi):
    """hits (any pwm) whose window overlaps [lo, hi] -> [(a, b, strand)]"""
    out = []
    for p in pwms:
        w = p.w
        a0, a1 = max(0, lo - w + 1), min(len(seq) - w, hi)
        for i in range(a0, a1 + 1):
            k = seq[i:i + w]
            if any(c not in IDX for c in k):
                continue
            sf = sum(p.lo[j, IDX[c]] for j, c in enumerate(k))
            if sf >= p.cut - 1e-5:
                out.append((i, i + w - 1, "+"))
            kr = rc(k)
            sr = sum(p.lo[j, IDX[c]] for j, c in enumerate(kr))
            if sr >= p.cut - 1e-5:
                out.append((i, i + w - 1, "-"))
    return out


def loss_gain(anc, der, pwms, mut_pos, anc_sites=None):
    """STRE site losses and gains between two equal-length sequences differing at mut_pos."""
    if not mut_pos:
        return 0, 0, [], []
    anc_sites = site_list(anc, pwms) if anc_sites is None else anc_sites
    maxw = max(p.w for p in pwms)
    lo, hi = min(mut_pos) - maxw, max(mut_pos) + maxw
    # only sites touching a mutated base can change
    touched = lambda a, b: any(a <= m <= b for m in mut_pos)
    a_s = [s for s in anc_sites if touched(s[0], s[1])]
    d_h = [h for h in window_hits(der, pwms, lo, hi) if touched(h[0], h[1])]
    a_h = [h for h in window_hits(anc, pwms, lo, hi) if touched(h[0], h[1])]
    losses = [s for s in a_s if not any(h[2] == s[2] and h[0] <= s[1] and h[1] >= s[0] for h in d_h)]
    # gains: derived hits touching a mutation, on a strand/interval with no ancestral hit; merge into sites
    new = [h for h in d_h if not any(a[2] == h[2] and a[0] <= h[1] and a[1] >= h[0] for a in a_h)]
    merged = []
    for strand in "+-":
        cur = None
        for a, b, st in sorted(h for h in new if h[2] == strand):
            if cur and a <= cur[1]:
                cur[1] = max(cur[1], b)
            else:
                cur = [a, b, strand]; merged.append(cur)
    return len(losses), len(merged), losses, merged


def align(unit, seqs, cache_dir):
    os.makedirs(cache_dir, exist_ok=True)
    fn = unit.replace("/", "_")
    af = f"{cache_dir}/{fn}.aln.fa"
    if not os.path.exists(af):
        fa = f"{cache_dir}/{fn}.fa"
        open(fa, "w").write("".join(f">{s}\n{seqs[s]}\n" for s in STRAINS))
        open(af, "w").write(subprocess.run([MAFFT, "--auto", "--quiet", fa], capture_output=True, text=True, check=True).stdout)
    a = read_fasta(af)
    for s in STRAINS:
        assert a[s].replace("-", "") == seqs[s], (unit, s)
    return a


def load_set(fasta, cache_dir):
    d = defaultdict(dict)
    for k, v in read_fasta(fasta).items():
        u, s = k.split("__")
        d[u][s] = v
    return {u: align(u, ss, cache_dir) for u, ss in d.items() if all(s in ss for s in STRAINS)}
