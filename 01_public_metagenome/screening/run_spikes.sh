#!/bin/bash
# Pass-1 spike controls, against the deposited Kp-GP genome.
#
# Reproduces cavia_screen/spike_pass1.sh with the new reference. N read pairs (25, 100, 500) of Kp-GP or
# Kp-M Illumina reads, taken after skipping the first 1 M pairs, added to the first 1 M pairs of a run in
# which Kp was not detected, then processed exactly as a pass-1 run. The spikes are what calibrate
# sensitivity, so they have to be measured on the same reference as the runs.
#
# Waits for run_screen.sh so the two do not stream at the same time.
set -uo pipefail
W=/Users/geongoohan/SynologyDrive/Brown/006_Fungal_WGS; KP=$W/kp_wgs_pangenome
D=$KP/public_survey/screen_v2; REF=$KP/public_survey/rigp_ref/rigp.fa
CACHE=$KP/sra_downloads/sra_cache
FQD=~/miniconda3/envs/sra/bin/fastq-dump; BWA=~/miniconda3/envs/pangenome/bin/bwa; SAM=~/miniconda3/envs/buscoenv/bin/samtools
TMP=/tmp/screening_work/screen_v2
mkdir -p $TMP; cd $D
log() { echo "$(date +%H:%M:%S) $*" >> logs/progress.log; }
while pgrep -f "bash run_screen.sh" > /dev/null; do sleep 60; done
log "spikes started"
[ -s $TMP/bg.fq ] || $FQD --split-spot --stdout -X 1000000 DRR848535 2> logs/fqd_spikebg.log > $TMP/bg.fq
rm -f $CACHE/sra/DRR848535.sra.cache
il() { gzip -cd "$1" | tail -n +4000001 | head -n $(( $3 * 4 )) | paste - - - - > $TMP/a.t; gzip -cd "$2" | tail -n +4000001 | head -n $(( $3 * 4 )) | paste - - - - > $TMP/b.t; paste $TMP/a.t $TMP/b.t | tr '\t' '\n'; rm -f $TMP/a.t $TMP/b.t; }
for N in 25 100 500; do for s in gp3 kp1; do
  L=SPIKE1M_${s}n$N; [ -s done/$L.ok ] && continue
  if [ $s = gp3 ]; then R1="$W/kpgp3/Illumina DNA Reads/kpgp3_S203_R1_001.fastq.gz"; R2="$W/kpgp3/Illumina DNA Reads/kpgp3_S203_R2_001.fastq.gz";
  else R1=$W/kp1/Illumina_reads/kp1_S52_R1_001.fastq.gz; R2=$W/kp1/Illumina_reads/kp1_S52_R2_001.fastq.gz; fi
  ( cat $TMP/bg.fq; il "$R1" "$R2" $N | awk -v L=$L 'NR%4==1{print "@" L "." int((NR+7)/8); next}{print}' ) > $TMP/$L.fq
  $BWA mem -p -t 6 $REF $TMP/$L.fq 2> logs/bwa_$L.log | $SAM view -u -F 4 - \
    | $SAM sort -@2 -m 1G -o bams/${L}_mapped.bam - && $SAM index bams/${L}_mapped.bam && echo spike > done/$L.ok
  echo "Written $((1000000 + N)) spots for $L" > logs/fqd_$L.log
  rm -f $TMP/$L.fq; log "$L done"
done; done
rm -f $TMP/bg.fq; log "spikes finished"
