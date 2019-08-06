#!/bin/bash


#### EXTRACT FASTQ FROM SRA FILE
echo -e "FASTQ-dump from SRA data $1 . . .\n"
prefix=$1
fastq_file="$1.fastq.gz"

#fastq-dump -Z --gzip $1 > $fastq_file
echo -e "$fastq_file was generated.\n"



### DETECT AND REMOVE ADAPTERS
echo -e "Detecting the adapter sequence . . .\n"
log="$fastq_file.adapter.txt"
adapter_err_rate=0.2
#python detect_adapter.py $fastq_file > $log
echo -e "Parsing The adapter sequence . . .\n"

trimming_command="tr -d  \040\011\012\015"
adapter_from="$(tail -1 $log | $trimming_command)"

echo -e "Adapter from :: $adapter_from \n"
#adapter_seq=$(awk -v adapt="$adapter_from" '($1==adapt&&$2!=NULL){print $3}' $log)
echo -e "Adapter sequence: $adapter_seq \n"
trimmed_fastq="$prefix.trim.fastq.gz"
#cutadapt -m 5 -e 0.2 -a $adapter_seq $fastq_file | gzip -c > $trimmed_fastq

#rm $fastq_file



## ALIGNMENT WITH BOWTIE2 (FOR PAIRED END)
multimapping=4
core_num=28
bwt2_idx="hg19"
bam="$prefix.bam"
log_align="$prefix.align.log"

echo -e "Sequence alignment with bowtie2 . . . \n"

bowtie2 -k $multimapping --local -x $bwt2_idx --threads $core_num -U <(zcat -f $trimmed_fastq) 2> $log_align | samtools view -Su /dev/stdin | samtools sort -o $prefix.bam
