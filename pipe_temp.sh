#!/bin/bash


#### EXTRACT FASTQ FROM SRA FILE
echo -e "FASTQ-dump from SRA data $1 . . .\n"
prefix=$1
fastq_file="$1.fastq.gz"

fastq-dump -Z --gzip $1 > $fastq_file
echo -e "$fastq_file was generated.\n"



### DETECT AND REMOVE ADAPTERS
echo -e "Detecting the adapter sequence . . .\n"
log="$fastq_file.adapter.txt"
adapter_err_rate=0.2
python detect_adapter.py $fastq_file > $log
echo -e "Parsing The adapter sequence . . .\n"

trimming_command="tr -d  \040\011\012\015"
adapter_from="$(tail -1 $log | $trimming_command)"

echo -e "Adapter from :: $adapter_from \n"
adapter_seq=$(awk -v adapt="$adapter_from" '($1==adapt&&$2!=NULL){print $3}' $log)
echo -e "Adapter sequence: $adapter_seq \n"
trimmed_fastq="$prefix.trim.fastq.gz"
cutadapt -m 5 -e 0.2 -a $adapter_seq $fastq_file | gzip -c > $trimmed_fastq

rm $fastq_file



## ALIGNMENT WITH BOWTIE2 (FOR PAIRED END)
multimapping=4
core_num=28
bwt2_idx="hg19"
bam="$prefix.bam"
log_align="$prefix.align.log"

echo -e "Sequence alignment with bowtie2 . . . \n"

bowtie2 -k $multimapping --local -x $bwt2_idx --threads $core_num -U <(zcat -f $trimmed_fastq) 2> $log_align | samtools view -Su /dev/stdin | samtools sort -o $prefix.bam



### REMOVE UNMAPPED, MATE UNMAPPED
qname_sort_bam="$prefix.qnmsrt.bam"
qname_sort_prefix="$prefix.qnmsrt"
qname_sort_bam_filt="$prefix.qnmsrt.filt.bam"

filt_bam_prefix="$prefix.filt.srt"
filt_bam="$filt_bam_prefix.bam"


echo -e "REMOVE UNMAPPED READS . . . \n"
samtools sort -n $bam -T $qname_sort_prefix -o $qname_sort_bam
samtools view -h $qname_sort_bam | python2 assign_multimappers.py -k $multimapping | samtools view -bS > $qname_sort_bam_filt

samtools view -F 1804 -b $qname_sort_bam_filt > $filt_bam

rm -f $qname_sort_bam $qname_sort_bam_filt
rm $bam
rm $trimmed_fastq


### MARK DUMPLICATES

tmp_filt_bam_file="$filt_bam_prefix.dupmakr.bam"
dup_file_qc="$filt_bam_prefix.dup.qc"

final_bam_prefix="$prefix.final"
final_bam="$final_bam_prefix.bam"

bamtools sort -in $filt_bam > "$filt_bam.tmpsrt"

rm $filt_bam

java -jar picard.jar MarkDuplicates I="$filt_bam.tmpsrt" O=$tmp_filt_bam_file M=$dup_file_qc ASSUME_SORTED=true

rm "$filt_bam.tmpsrt"

samtools view -F 1804 -b $tmp_filt_bam_file > $final_bam
rm $tmp_filt_bam_file
