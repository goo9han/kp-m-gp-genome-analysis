#!/bin/bash
# The 21 positive runs taken through the full-run path against the deposited Kp-GP genome.
# Whole run streamed and piped straight into bwa, so no
# large FASTQ is written. Waits for the screen and the spikes so nothing streams at the same time.
set -uo pipefail
KP=/Users/geongoohan/SynologyDrive/Brown/006_Fungal_WGS/kp_wgs_pangenome
D=$KP/public_survey/typing_v2; REF=$KP/public_survey/rigp_ref/rigp.fa
CACHE=$KP/sra_downloads/sra_cache
FQD=~/miniconda3/envs/sra/bin/fastq-dump; BWA=~/miniconda3/envs/pangenome/bin/bwa
SAM=~/miniconda3/envs/buscoenv/bin/samtools; PY=~/miniconda3/envs/rnaseq_prep/bin/python
cd $D
log() { echo "$(date +%H:%M:%S) $*" | tee -a logs/progress.log; }
while pgrep -f "bash run_screen.sh|bash run_spikes.sh" > /dev/null; do sleep 60; done
log "lineage typing started, reference rigp.fa"
for R in DRR848895 DRR849046 DRR848899 DRR849092 DRR848835 DRR848897 DRR849095 DRR848808 DRR848896 \
         DRR848894 DRR848753 DRR848795 DRR848939 DRR849048 DRR849047 DRR849049 DRR848839 DRR848838 \
         DRR848802 DRR848804 DRR848801; do
  [ -s work/$R.fq ] && { log "$R already done"; continue; }
  FREE=$(df -g . | awk 'NR==2{print $4}')
  if [ "$FREE" -lt 15 ]; then log "stopping, ${FREE} GB free"; break; fi
  $FQD --split-spot --stdout $R 2> logs/fqd_$R.log \
    | $BWA mem -p -t 6 $REF /dev/stdin 2> logs/bwa_$R.log \
    | $SAM view -u -F 4 - | $SAM sort -@2 -m 1G -o bams/${R}_mapped.bam -
  ST=("${PIPESTATUS[@]}")
  rm -rf $CACHE/sra/$R.sra.cache $CACHE/sra/$R.sra
  if [ "${ST[0]}" != "0" ] || [ "${ST[3]}" != "0" ]; then log "$R failed"; rm -f bams/${R}_mapped.bam; continue; fi
  $SAM index bams/${R}_mapped.bam || { log "$R index failed"; continue; }
  log "$R mapped $(grep -o 'Written [0-9]* spots' logs/fqd_$R.log | tail -1)"
  $PY extract_kp_confident.py $R bams/${R}_mapped.bam work/$R.fq 2>&1 | grep -v Runtime | tee -a logs/progress.log
done
log "lineage typing mapping finished"
