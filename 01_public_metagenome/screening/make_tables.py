#!/usr/bin/env python3
"""Province and two-way prevalence tables from the screen run against the deposited Kp-GP genome.

Replaces cavia_screen/pass1_provinces.tsv and cavia_screen/lanes/prevalence_two_ways.tsv as the source
the figure reads, so the figure's provenance points at the re-run rather than at the earlier analysis.
Column names and ordering are kept so make_panels.py needs only a path change.

(a) every pass-1 positive. (b) conservative, a positive carrying 20 or fewer Kp-confident reads is
removed when it shares a sequencing lane with a run carrying 100 or more. Both thresholds were fixed
before tabulation, and both figures are lower bounds because only part of each run was read.

Writes provinces_v2.tsv and prevalence_two_ways_v2.tsv
"""
import os, csv, json, collections
from scipy.stats import beta

HERE = os.path.dirname(os.path.abspath(__file__))
KP = os.path.dirname(os.path.dirname(HERE))
SCREEN = json.load(open(f"{HERE}/screen_v2.json"))
PREV = json.load(open(f"{HERE}/prevalence_v2.json"))
RATE = 0.661          # Kp-confident reads recovered per spiked read, from the spikes of this same run
DEPTH = "1M pairs per run"


def meta():
    f = open(f"{KP}/public_survey/cavia_screen/pass1_runs.tsv")
    while True:
        p = f.tell()
        if not f.readline().startswith("#"):
            f.seek(p); break
    return {r["run"]: r for r in csv.DictReader(f, delimiter="\t")}


def cp(k, n):
    lo = beta.ppf(0.025, k, n - k + 1) if k else 0.0
    hi = beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return f"{lo:.2f}-{hi:.2f}"


if __name__ == "__main__":
    M = meta()
    reads = PREV["reads_per_1M"]
    runs = list(reads)
    prov = {r: M[r]["province"] for r in runs}
    a = set(PREV["positives"]); b = a - set(PREV["dirty"])
    screened = collections.Counter(prov[r] for r in runs)
    libs = collections.Counter(prov[r] for r in runs if M[r].get("library"))
    da = collections.Counter(prov[r] for r in a)
    db = collections.Counter(prov[r] for r in b)
    mx = collections.defaultdict(int)
    for r in runs:
        mx[prov[r]] = max(mx[prov[r]], reads[r])

    with open(f"{HERE}/provinces_v2.tsv", "w") as fh:
        fh.write("# Province-level result (primary unit). Prevalence = detected / screened runs, a LOWER "
                 "bound (1 M read pairs per run). Reference: the deposited Kp-GP genome.\n")
        fh.write("province\tlibraries\tscreened_runs\tdetected_runs\tprevalence_lower_bound\tCP95"
                 "\tmax_reads_in_a_run\tpooled_upper_bound_if_none\tscreening_depth\n")
        for p in sorted(screened):
            n, d = screened[p], da[p]
            ub = "" if d else f"{3 / RATE / (n * 2e6):.1e}"
            fh.write(f"{p}\t{libs[p]}\t{n}\t{d}\t{d / n:.2f}\t{cp(d, n)}\t{mx[p]}\t{ub}\t{DEPTH}\n")

    with open(f"{HERE}/prevalence_two_ways_v2.tsv", "w") as fh:
        fh.write("# (a) all pass-1 positives; (b) conservative: low positives (<= 20 reads) sharing a lane "
                 "with a >= 100-read run removed. Both are lower bounds.\n")
        fh.write("# Reference: the deposited Kp-GP genome. Every run was screened at 1,000,000 read pairs, "
                 "so these counts are lower bounds set by the screening depth, not by the biology.\n")
        fh.write("unit\tscreened\ta_detected\tb_detected\ta_CP95\tb_CP95\tscreening_depth\n")
        fh.write(f"runs\t{len(runs)}\t{len(a)}\t{len(b)}\t{cp(len(a), len(runs))}\t{cp(len(b), len(runs))}\t{DEPTH}\n")
        pa = sum(1 for p in screened if da[p]); pb = sum(1 for p in screened if db[p])
        fh.write(f"provinces\t{len(screened)}\t{pa}\t{pb}\t\t\t{DEPTH}\n")
        for p in sorted(screened):
            fh.write(f"{p}\t{screened[p]}\t{da[p]}\t{db[p]}\t\t\t{DEPTH}\n")
    print(f"  wrote provinces_v2.tsv and prevalence_two_ways_v2.tsv")
    print(f"  runs {len(runs)}  (a) {len(a)}  (b) {len(b)}   provinces {len(screened)}  (a) {pa}  (b) {pb}")
