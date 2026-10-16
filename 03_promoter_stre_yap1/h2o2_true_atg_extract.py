#!/usr/bin/env python3
"""Re-extract the H2O2-set promoters (700 bp immediately upstream of the verified start codon) in 9 strains.

Generalises ahp1_yap1_analysis/true_atg/ahp1_true_atg_extract.py to all loci of the H2O2 set.

Why the old promoters could not simply be shifted: meme_fimo/all_promoters.fa was cut from tblastn hit
coordinates. For genes on the - strand of the kp1 assembly the window was taken on the + strand at the
low end of the hit, i.e. 3' of the gene (or inside the CDS for YAP1) and in antisense orientation. For
SKN7 the tblastn hit began 120 codons into the CDS, so the window is CDS. Only AHP1, CTA1, TRR1 and
TRX1/TRX2 windows are real upstream sequence (end 9-30 bp before the ATG). The first-ATG offsets in
the old windows are therefore not used; each locus is re-derived from the kp1 funannotate model.

Procedure
  kp1: start codon = funannotate CDS start (rnaseq_prep/reference/kp1_funannotate.gtf), except where
       KP1_ATG_SHIFT says otherwise. Accepted only if the ORF has no stop for >= MIN_CODONS and the
       protein aligns to the S. cerevisiae protein (h2o2_query.fa) with local identity > 40%.
  other strains: the kp1 protein (from the verified ATG) is mapped with miniprot; the rank-1 hit must
       start at kp1 residue 1 (the ATG) with identity >= 0.90, and the same ORF/identity checks apply.
  promoter = genome[-700..-1] of the ATG on the gene strand.
"""
import os, re, subprocess
from Bio.Seq import Seq
from Bio import Align
from Bio.Align import substitution_matrices

HERE = os.path.dirname(os.path.abspath(__file__))
H2O2 = os.path.dirname(HERE)
ROOT = os.path.dirname(H2O2)
GENOME_DIR = f"{ROOT}/rnaseq_prep/ani/nuclear_clean"
GTF = f"{ROOT}/rnaseq_prep/reference/kp1_funannotate.gtf"
MINIPROT = os.path.expanduser("~/miniconda3/envs/buscoenv/bin/miniprot")
STRAINS = ["kp1", "kpmtx", "C89-2", "ID1LA", "IDNYC", "QC1", "W26G", "kpgp3", "kpgp6"]
OLD_STRAIN = {s: f"{'kp_' if s in ('C89-2', 'ID1LA', 'IDNYC', 'QC1', 'W26G') else ''}{s}" for s in STRAINS}
PROM_LEN = 700
rc = lambda s: s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]

# label, Sc query, old promoter set name(s), kp1 funannotate gene, min ORF codons, note
LOCI = [
    ("AHP1", "AHP1_TPx2", ["AHP1_TPx2"], "FUN_003455", 150, "control: must reproduce ahp1_true_atg_extract.py"),
    ("YAP1", "YAP1_TF", ["YAP1_TF"], "FUN_002921", 580, ""),
    ("SKN7", "SKN7_TF", ["SKN7_TF"], "FUN_003997", 720, "GTF ATG = 5'-most in ORF (in-frame stop -21); alt in-frame ATG +243; Sc N-term low complexity"),
    ("CTA1", "CTA1_catalase_A", ["CTA1_catalase_A"], "FUN_000769", 460, "PTS1 ...SKF kept"),
    ("CTT1", "CTT1_catalase_T", ["CTT1_catalase_T"], "FUN_002314", 460, "73% id to Sc CTA1 vs 48% to CTT1: CTA1 paralog, no PTS1"),
    ("SOD1", "SOD1_CuZnSOD", ["SOD1_CuZnSOD"], "FUN_003967", 140, ""),
    ("SOD2", "SOD2_MnSOD", ["SOD2_MnSOD"], "FUN_000528", 210, ""),
    ("TSA1", "TSA1_TPx", ["TSA1_TPx", "TSA2_TPx"], "FUN_000514", 180, "next to YOX1 homolog = Sc TSA1 synteny; old TSA1 and TSA2 windows both at this locus"),
    ("TSA2", "TSA2_TPx", [], "FUN_001514", 180, "scaffold_3 TSA paralog; not in old set; TSA2 orthology not resolved by synteny"),
    ("GPX2", "GPX2_glutathione_peroxidase", ["GPX2_glutathione_peroxidase"], "FUN_003437", 145, ""),
    ("TRX2", "TRX2_thioredoxin", ["TRX1_thioredoxin", "TRX2_thioredoxin"], "FUN_000958", 95,
     "next to SER2 homolog = Sc TRX2 synteny; old TRX1 and TRX2 windows are the same sequence; no TRX1 ortholog"),
    ("TRR1", "TRR1_thioredoxin_reductase", ["TRR1_thioredoxin_reductase"], "FUN_001968", 290, ""),
    ("GLR1", "GLR1_glutathione_reductase", ["GLR1_glutathione_reductase"], "FUN_001390", 430,
     "ATG shifted -36 to 5'-most in-frame ATG (= Sc Glr1 M1, mito presequence); GTF ATG = Sc M17 cytosolic start"),
    ("GSH1", "GSH1_glutamate_cysteine_ligase", ["GSH1_glutamate_cysteine_ligase"], "FUN_000135", 640, ""),
]
KP1_ATG_SHIFT = {"GLR1": -36}


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}


aligner = Align.PairwiseAligner(mode="local", open_gap_score=-10, extend_gap_score=-0.5)
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")


def identity(p, q):
    aln = aligner.align(p, q)[0]
    ident = sum(1 for (a0, a1), (b0, b1) in zip(*aln.aligned) for x, y in zip(p[a0:a1], q[b0:b1]) if x == y)
    alen = sum(a1 - a0 for a0, a1 in aln.aligned[0])
    return ident / alen, alen, aln.aligned[0][0][0] + 1, aln.aligned[1][0][0] + 1


def orf(seq, k):
    s = seq[k:k + 9000]
    return str(Seq(s[:len(s) // 3 * 3]).translate()).split("*")[0]


def old_window_relation(seq, k, oldp):
    """Where the old promoter sits relative to the ATG (index k on the gene strand)."""
    if oldp is None:
        return "not in old set"
    i, j = seq.find(oldp), seq.find(rc(oldp))
    if i >= 0:
        return f"same strand, {i - k}..{i - k + PROM_LEN - 1} from ATG"
    if j >= 0:
        return f"ANTISENSE, gene-strand {j - k}..{j - k + PROM_LEN - 1} from ATG"
    return "other contig"


sc = read_fasta(f"{H2O2}/h2o2_query.fa")
old = read_fasta(f"{H2O2}/meme_fimo/all_promoters.fa")
cds = {}
for line in open(GTF):
    f = line.rstrip("\n").split("\t")
    if len(f) > 8 and f[2] == "CDS":
        gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        cds.setdefault(gid, []).append((f[0], int(f[3]), int(f[4]), f[6]))

genomes = {s: read_fasta(f"{GENOME_DIR}/{s}_nuclear_clean.fa") for s in STRAINS}
TMP = f"{HERE}/miniprot"
os.makedirs(TMP, exist_ok=True)

log = ["locus\tstrain\tcontig\tgene_strand\tATG_A_contig_pos\tORF_codons_no_stop\tprotein_Nterm\tSc_local_identity\t"
       "Sc_aligned_len\tKp_aln_start\tSc_aln_start\tkp1_protein_identity(miniprot)\told_window_vs_ATG\tnew_promoter_contig_range\tnote"]
fa = []


def record(lab, s, contig, strand, seq, k, prot, sc_q, mp_id, oldname, note, minc):
    assert seq[k:k + 3] == "ATG", (lab, s, seq[k:k + 3])
    if len(prot) < minc:
        # accepted only when miniprot aligned the whole kp1 protein from residue 1 across a frameshift
        assert "frameshift" in note and len(prot) >= 150, (lab, s, len(prot), note)
    r, alen, kp0, sc0 = identity(prot, sc[sc_q])
    assert r > 0.4, (lab, s, r)
    L = len(seq)
    up = min(k, PROM_LEN)                       # contig may end before -700: pad the 5' end with N
    if up < PROM_LEN:
        note = (note + "; " if note else "") + f"contig ends {up} bp upstream of ATG: 5' {PROM_LEN - up} bp padded with N"
    if strand == "+":
        atg_pos, rng = k + 1, f"{contig}:{k - up + 1}-{k}(+)"
    else:
        atg_pos, rng = L - k, f"{contig}:{L - k + 1}-{L - k + up}(-)"
    oldp = old.get(f"{oldname}__{OLD_STRAIN[s]}") if oldname else None
    log.append(f"{lab}\t{s}\t{contig}\t{strand}\t{atg_pos}\t{len(prot)}\t{prot[:20]}\t{r:.3f}\t{alen}\t{kp0}\t{sc0}\t"
               f"{mp_id}\t{old_window_relation(seq, k, oldp)}\t{rng}\t{note}")
    fa.append(f">{lab}__{s}\n{'N' * (PROM_LEN - up)}{seq[k - up:k]}")
    print(f"{lab:5s} {s:6s} {contig}:{atg_pos}({strand}) ORF {len(prot)} aa {prot[:10]}.. Sc id {r:.1%}/{alen} "
          f"| old: {old_window_relation(seq, k, oldp)}")


for lab, sc_q, oldnames, gid, minc, note in LOCI:
    oldname = oldnames[0] if oldnames else None
    # ---- kp1 from the annotation
    contig, _, _, strand = cds[gid][0]
    g = genomes["kp1"]
    seq = g[contig] if strand == "+" else rc(g[contig])
    k = (min(x[1] for x in cds[gid]) - 1) if strand == "+" else len(seq) - max(x[2] for x in cds[gid])
    k += KP1_ATG_SHIFT.get(lab, 0)
    kp1_prot = orf(seq, k)
    record(lab, "kp1", contig, strand, seq, k, kp1_prot, sc_q, "1.000", oldname, note, minc)
    qfa = f"{TMP}/{lab}_kp1.faa"
    open(qfa, "w").write(f">{lab}\n{kp1_prot}\n")
    # ---- other strains via miniprot of the kp1 protein
    for s in STRAINS[1:]:
        out = subprocess.run([MINIPROT, "-t", "4", "--gff", f"{GENOME_DIR}/{s}_nuclear_clean.fa", qfa],
                             capture_output=True, text=True, check=True).stdout
        open(f"{TMP}/{lab}_{s}.gff", "w").write(out)
        m = [l.split("\t") for l in out.splitlines() if "\tmRNA\t" in l and "Rank=1;" in l]
        assert len(m) == 1, (lab, s, len(m))
        f = m[0]
        tgt = re.search(r"Target=\S+ (\d+) (\d+)", f[8]); mp_id = float(re.search(r"Identity=([\d.]+)", f[8]).group(1))
        assert tgt.group(1) == "1" and mp_id >= 0.90, (lab, s, tgt.group(0), mp_id)
        contig, strand = f[0], f[6]
        g = genomes[s]
        seq = g[contig] if strand == "+" else rc(g[contig])
        k = int(f[3]) - 1 if strand == "+" else len(seq) - int(f[4])
        fs = re.search(r"Frameshift=(\d+)", f[8])
        full = int(tgt.group(2)) == len(kp1_prot)
        snote = f"miniprot frameshift={fs.group(1)}, kp1 protein aligned 1-{tgt.group(2)} (assembly indel?)" if fs and full else ""
        record(lab, s, contig, strand, seq, k, orf(seq, k), sc_q, f"{mp_id:.3f}", oldname, snote, minc)

open(f"{HERE}/H2O2_promoters_trueATG_9strains.fa", "w").write("\n".join(fa) + "\n")
open(f"{HERE}/true_atg_extraction_log.tsv", "w").write("\n".join(log) + "\n")
