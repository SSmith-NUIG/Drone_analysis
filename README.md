# Drone_analysis

Set up conda environment
```bash
conda create --name wgs_env
```

install the following packages
```bash
conda install -c bioconda fastqc
conda install -c bioconda samtools
conda install -c bioconda picard
conda install -c bioconda multiqc
conda install -c bioconda bbmap
conda install -c bioconda qualimap
conda install -c bioconda bwa-mem2
conda install -c bioconda gatk
conda install -c bioconda snpeff
conda install -c bioconda vcftools
```
download the fastq files from novogene (.zip format)
```bash
wget -P /path/to/local/foler/ novogene.url.here
```

unzip the folder
```bash
unzip /path/to/local/folder/files.zip
```
You can rename the files to something which is easier to read if you wish with renaming.sh

```bash
#!/bin/sh

#this script will loop through each directory and rename all files with the extension .gz
#it keeps the S[Sample Number]_ and the read number (read 1 or 2)
#NOTE - this command uses the perl rename program, not the built in linux-utils rename
#perl rename can be installed using anaconda. conda install -c bioconda rename


#USE -n AFTER rename TO DO A TEST RUN WHICH WILL SHOW YOU WHAT THE FILE NAME WILL BECOME
#REMOVE THIS -n WHEN YOU ARE SURE THE PATTERN IS CORRECT



dirPath=("/data/ssmith/drone_data/raw_data/clumped_data")

for dir in $dirPath/*/; do
	
	firstName=$(ls $dir | sed -n 1p)
	secondName=$(ls $dir | sed -n 2p)
	rename 's{(S[0-9]_)[^.]*([12])[^.]*}{$1$2}' $dir$firstName
	rename 's{(S[0-9]_)[^.]*([12])[^.]*}{$1$2}' $dir$secondName
done
```
create a tab separated file called adaptors.txt which contains the novogene adaptor sequences
```bash
first	AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
second	GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
```

This analysis was run on a HPC using SLURM job submissions.
All jobs were run using an array format where "$SLURM_ARRAY_TASK_ID" will relate to a novogene sample number

Run FastQC on the files to get a sense of the quality you are dealing with.
```bash
fastqc S"$SLURM_ARRAY_TASK_ID"*.fq.gz --noextract -t 6 -a /path/to/adaptors/adaptors.txt --outdir=Metrics
```
bbduk.sh from BBmap to trim adaptors

```bash
mkdir ./cleaned
```

```bash
bbduk.sh \
in1=A_"$SLURM_ARRAY_TASK_ID"*1.fq.gz \
in2=A_"$SLURM_ARRAY_TASK_ID"*2.fq.gz \
out1=/cleaned/A_"$SLURM_ARRAY_TASK_ID"_1.fq.gz \
out2=/cleaned/A_"$SLURM_ARRAY_TASK_ID"_2.fq.gz \
ref=/adaptors.txt \
mink=11 \
hdist=1 \
tbo=t
```

bbduk.sh to trim low quality reads

```bash
mkdir ./trimmed
```

```bash
bbduk.sh \
in1=/cleaned/A_"$SLURM_ARRAY_TASK_ID"_1.fq.gz \
in2=/cleaned/A_"$SLURM_ARRAY_TASK_ID"_2.fq.gz \
out1=/trimmed/A_"$SLURM_ARRAY_TASK_ID"_1.fq.gz \
out2=//trimmed/A_"$SLURM_ARRAY_TASK_ID"_2.fq.gz \
qtrim=r \
trimq=10
```
Clean as you go, get rid of intermediate files

```bash
rm /cleaned/A_"$SLURM_ARRAY_TASK_ID"_1.fq.gz
rm /cleaned/A_"$SLURM_ARRAY_TASK_ID"_2.fq.gz
```
Run clumpify to get rid of optical duplicates in your files

```bash
mkdir deopt
```

```bash
clumpify.sh \
in1=/trimmed/A_"$SLURM_ARRAY_TASK_ID"_1.fq.gz \
in2=/trimmed/A_"$SLURM_ARRAY_TASK_ID"_2.fq.gz \
out1=/deopt/A_"$SLURM_ARRAY_TASK_ID"_clumped_1.fq.gz \
out2=/deopt/A_"$SLURM_ARRAY_TASK_ID"_clumped_2.fq.gz \
dedupe \
optical;
```

Clean as you go

```bash
rm /trimmed/A_"$SLURM_ARRAY_TASK_ID"_1.fq.gz
rm /trimmed/A_"$SLURM_ARRAY_TASK_ID"_2.fq.gz
```

Run FastQC to check effects of adaptor trimming and low quality read removal
```bash

mkdir After_QC_Metrics
```

```bash
fastqc deopt/A_"$SLURM_ARRAY_TASK_ID"*.fq.gz --noextract -t 6 -a /adaptors.txt --outdir=/After_QC_Metrics
```

Alignment
Download the genome in fastq format

```bash
wget -P /path/to/genome/ url.of.genome
```

Index the genome using samtools and bwa-mem2

```bash
samtools faidx /path/to/genome/genome.fa

bwa-mem2 index /path/to/genome/genome.fa
```

Align the reads to the genome using bwa-mem2
Pipe the output to samtools sort which will sort the reads by chromosome location 
```bash
mkdir sorted
```

```bash
bwa-mem2 mem -t 12 /data/ssmith/c_l_genome/apis_c_l_genome.fa  \
deopt/A_"$SLURM_ARRAY_TASK_ID"_clumped_1.fq.gz \
deopt/A_"$SLURM_ARRAY_TASK_ID"_clumped_2.fq.gz \
| samtools sort -@ 6 -o sorted/A_"$SLURM_ARRAY_TASK_ID".sorted.bam
```
Clean as you go
```bash
rm deopt/A_"$SLURM_ARRAY_TASK_ID"_clumped_1.fq.gz
rm deopt/A_"$SLURM_ARRAY_TASK_ID"_clumped_2.fq.g
```
Use Samtools flagstat to get information on the reads such as duplication levels, alignment, read pairs matching etc.
```bash
samtools flagstat sorted/A_"$SLURM_ARRAY_TASK_ID".sorted.bam \
> /Metrics/A_"$SLURM_ARRAY_TASK_ID"_flagstat.txt
```

Use picard MarkDuplicates to find and mark PCR duplicates in the reads

```bash
picard MarkDuplicates \
I=sorted/A_"$SLURM_ARRAY_TASK_ID".sorted.bam \
M=sorted/After_QC_Metrics/A_"$SLURM_ARRAY_TASK_ID"_dup_metrics.txt \
O=sorted/A_"$SLURM_ARRAY_TASK_ID".sorted.dupM.bam
```

Clean as you go
```bash
rm sorted/A_"$SLURM_ARRAY_TASK_ID".sorted.bam
```

Get new samtools flagstat information to see the effects of Picard on the data
```bash 
samtools flagstat sorted/A_"$SLURM_ARRAY_TASK_ID".sorted.dupM.bam \
> sorted/After_QC_Metrics/A_"$SLURM_ARRAY_TASK_ID"_flagstat.txt
```

Use Samtools index to index the sorted bams
```bash
samtools index "$sortedPath"/A_"$SLURM_ARRAY_TASK_ID".sorted.dupM.bam
```

Samtools to view number of alignments:
```bash
samtools view "$sortedPath"/dupMarked/S"$SLURM_ARRAY_TASK_ID".sorted.dupM.bam \
| wc -l > "$sortedPath"/dupNumAln/S"$SLURM_ARRAY_TASK_ID"_num_alignments.txt
```
SAMTOOLS VIEW HEADER INFO
```bash
samtools view -H "$sortedPath"/dupMarked/S"$SLURM_ARRAY_TASK_ID".sorted.dupM.bam \
> "$sortedPath"/dupHeaders/S"$SLURM_ARRAY_TASK_ID"_header_info.txt
 ```
PICARD ADD OR REPLACE GROUPS

```bash
picard AddOrReplaceReadGroups \
I="$sortedPath"/dupMarked/S"$SLURM_ARRAY_TASK_ID".sorted.dupM.bam \
O="$clumPath"/c_l_dupMarkedGrouped/S"$SLURM_ARRAY_TASK_ID".nobqsr.grpd.bam \
RGID=4 \
RGLB=lib1 \
RGPL=ILLUMINA \
RGPU=unit1 \
RGSM=S"$SLURM_ARRAY_TASK_ID"
```

```bash
rm "$sortedPath"/dupMarked/S"$SLURM_ARRAY_TASK_ID".sorted.dupM.bam
```

BASE RECALIBRATOR
```bash
gatk BaseRecalibrator \
-I "$clumPath"/c_l_dupMarkedGrouped/S"$SLURM_ARRAY_TASK_ID".nobqsr.grpd.bam \
-R /data/ssmith/c_l_genome/apis_c_l_genome.fa \
--known-sites $featureFile \
-O /data/ssmith/c_l_genome/recal_data.table \
--bqsr-baq-gap-open-penalty 30.0
```

BASE RECALIBRATOR APPLY BQSR
```bash
gatk ApplyBQSR \
-R /data/ssmith/c_l_genome/apis_c_l_genome.fa \
-I "$clumPath"/c_l_dupMarkedGrouped/S"$SLURM_ARRAY_TASK_ID".nobqsr.grpd.bam \
--bqsr-recal-file /data/ssmith/c_l_genome/recal_data.table \
-O "$clumPath"/c_l_dupMarkedGrouped/S"$SLURM_ARRAY_TASK_ID".grpd.bam
```

```bash
rm "$clumPath"/c_l_dupMarkedGrouped/S"$SLURM_ARRAY_TASK_ID".nobqsr.grpd.bam
```


QUALIMAP
```bash
qualimap bamqc \
-bam "$clumPath"/c_l_dupMarkedGrouped/S"$SLURM_ARRAY_TASK_ID".grpd.bam \
-sd -os -c -oc "$dirPath"/quali/S"$SLURM_ARRAY_TASK_ID"/S"$SLURM_ARRAY_TASK_ID".txt \
-gff /data/ssmith/c_l_genome/cl.gtf \
-outdir "$dirPath"/quali/S"$SLURM_ARRAY_TASK_ID"/ \
-outformat HTML -outfile S"$SLURM_ARRAY_TASK_ID".html 
```

SAMTOOLS INDEX
```bash
samtools index "$clumPath"/c_l_dupMarkedGrouped/S"$SLURM_ARRAY_TASK_ID".grpd.bam
```

HAPLOTYPECALLER

```bash
gatk HaplotypeCaller \
-I "$clumPath"/c_l_dupMarkedGrouped/S"$SLURM_ARRAY_TASK_ID".grpd.bam \
-R /data/ssmith/gatk_ref/apis_c_l_genome.fa \
-ERC GVCF \
-O /data/ssmith/drone_data/cldata/GATK_c_l_output/S"$SLURM_ARRAY_TASK_ID".g.vcf
```

SNPEFF
```bash
snpEff ann apiscl /data/ssmith/drone_data/cldata/GATK_c_l_output/S"$SLURM_ARRAY_TASK_ID".g.vcf \
-s /data/ssmith/drone_data/cldata/annotation_results/S"$SLURM_ARRAY_TASK_ID"_summary.html \
> /data/ssmith/drone_data/cldata/annotation_results/S"$SLURM_ARRAY_TASK_ID".vcf
```
