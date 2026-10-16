#!/bin/bash
# PRJDB38973 selection, screened against the deposited Kp-GP genome.
#
# Why this re-run. The earlier screen mapped to a local file that joined the Kp-GP assembly to a 13 kb
# contig of the Kp-M assembly carrying CHK1. That contig served a second question, whether a sample's Kp
# retains CHK1, which this paper does not report. For the question the paper does ask, whether Kp is
# present, the contig was never counted, and 11 kb of it duplicated sequence already in the counting set,
# which can only split mapping quality and suppress reads in that region. The reference is now the
# deposited genome alone. Its contigs of at least 20 kb are identical in name, length and sequence to the
# counting set used before, so the analysis space does not change, only its provenance and its cleanliness.
#
# The deposited nuclear assembly also drops three short contigs that the earlier local assembly carried,
# scaffold_119, scaffold_126 and scaffold_137, which are mitochondrion-derived. All three are far below
# 20 kb and were never in the counting set, so removing them changes nothing there, and it removes three
# spurious nuclear homes for mitochondrial reads.
#
# Streaming and mapping are as before, 1 M read pairs per run, bwa mem -p, mapped reads kept.
# Resumable: a run with done/<run>.ok is skipped. The SRA cache is deleted per run.
set -uo pipefail
KP=/Users/geongoohan/SynologyDrive/Brown/006_Fungal_WGS/kp_wgs_pangenome
D=$KP/public_survey/screen_v2
REF=$KP/public_survey/rigp_ref/rigp.fa
CACHE=$KP/sra_downloads/sra_cache
FQD=~/miniconda3/envs/sra/bin/fastq-dump
BWA=~/miniconda3/envs/pangenome/bin/bwa
SAM=~/miniconda3/envs/buscoenv/bin/samtools
TMP=/tmp/screening_work/screen_v2
mkdir -p $TMP; cd $D
log() { echo "$(date +%H:%M:%S) $*" >> logs/progress.log; }
log "screening started, reference $(basename $REF)"
N=0; T=$(wc -l < runs.txt | tr -d ' ')
while read R; do
  N=$((N+1))
  [ -s done/$R.ok ] && continue
  FREE=$(df -g . | awk 'NR==2{print $4}')
  if [ "$FREE" -lt 15 ]; then log "stopping, ${FREE} GB free"; break; fi
  $FQD --split-spot --stdout -X 1000000 $R 2> logs/fqd_$R.log > $TMP/$R.fq
  if [ ! -s $TMP/$R.fq ]; then log "$R streaming failed"; rm -f $TMP/$R.fq; continue; fi
  $BWA mem -p -t 6 $REF $TMP/$R.fq 2> logs/bwa_$R.log | $SAM view -u -F 4 - \
    | $SAM sort -@2 -m 1G -o bams/${R}_mapped.bam - && $SAM index bams/${R}_mapped.bam && echo ok > done/$R.ok
  rm -f $TMP/$R.fq $CACHE/sra/$R.sra.cache $CACHE/sra/$R.sra
  log "$R ($N/$T) $(grep -o 'Written [0-9]* spots' logs/fqd_$R.log | tail -1)"
done < runs.txt
log "screening mapping finished"
