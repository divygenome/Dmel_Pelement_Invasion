1 - Prep & mapping
================
Matthew Beaumont
2023-07-17

``` bash
knitr::opts_chunk$set(echo = TRUE)
```

We obtained the following mRNA reads:

``` bash
cd /Volumes/Data/Projects/dmelR2_p-ele/rna/run1/raw
ls *.fastq.gz

cd /Volumes/Data/Projects/dmelR2_p-ele/rna/run2/raw
ls *.fastq.gz
```

    ## dmel_rna_R1_G15_run1_R1.fastq.gz
    ## dmel_rna_R1_G15_run1_R2.fastq.gz
    ## dmel_rna_R1_G21_run1_R1.fastq.gz
    ## dmel_rna_R1_G21_run1_R2.fastq.gz
    ## dmel_rna_R1_G30_run1_R1.fastq.gz
    ## dmel_rna_R1_G30_run1_R2.fastq.gz
    ## dmel_rna_R1_G40_run1_R1.fastq.gz
    ## dmel_rna_R1_G40_run1_R2.fastq.gz
    ## dmel_rna_R1_G6_run1_R1.fastq.gz
    ## dmel_rna_R1_G6_run1_R2.fastq.gz
    ## dmel_rna_R2_G15_run1_R1.fastq.gz
    ## dmel_rna_R2_G15_run1_R2.fastq.gz
    ## dmel_rna_R2_G21_run1_R1.fastq.gz
    ## dmel_rna_R2_G21_run1_R2.fastq.gz
    ## dmel_rna_R2_G30_run1_R1.fastq.gz
    ## dmel_rna_R2_G30_run1_R2.fastq.gz
    ## dmel_rna_R2_G40_run1_R1.fastq.gz
    ## dmel_rna_R2_G40_run1_R2.fastq.gz
    ## dmel_rna_R2_G6_run1_R1.fastq.gz
    ## dmel_rna_R2_G6_run1_R2.fastq.gz
    ## dmel_rna_R3_G15_run1_R1.fastq.gz
    ## dmel_rna_R3_G15_run1_R2.fastq.gz
    ## dmel_rna_R3_G21_run1_R1.fastq.gz
    ## dmel_rna_R3_G21_run1_R2.fastq.gz
    ## dmel_rna_R3_G30_run1_R1.fastq.gz
    ## dmel_rna_R3_G30_run1_R2.fastq.gz
    ## dmel_rna_R3_G40_run1_R1.fastq.gz
    ## dmel_rna_R3_G40_run1_R2.fastq.gz
    ## dmel_rna_R3_G6_run1_R1.fastq.gz
    ## dmel_rna_R3_G6_run1_R2.fastq.gz
    ## dmel_rna_R1_G15_run2_R1.fastq.gz
    ## dmel_rna_R1_G15_run2_R2.fastq.gz
    ## dmel_rna_R1_G21_run2_R1.fastq.gz
    ## dmel_rna_R1_G21_run2_R2.fastq.gz
    ## dmel_rna_R1_G30_run2_R1.fastq.gz
    ## dmel_rna_R1_G30_run2_R2.fastq.gz
    ## dmel_rna_R1_G40_run2_R1.fastq.gz
    ## dmel_rna_R1_G40_run2_R2.fastq.gz
    ## dmel_rna_R1_G6_run2_R1.fastq.gz
    ## dmel_rna_R1_G6_run2_R2.fastq.gz
    ## dmel_rna_R2_G15_run2_R1.fastq.gz
    ## dmel_rna_R2_G15_run2_R2.fastq.gz
    ## dmel_rna_R2_G21_run2_R1.fastq.gz
    ## dmel_rna_R2_G21_run2_R2.fastq.gz
    ## dmel_rna_R2_G30_run2_R1.fastq.gz
    ## dmel_rna_R2_G30_run2_R2.fastq.gz
    ## dmel_rna_R2_G40_run2_R1.fastq.gz
    ## dmel_rna_R2_G40_run2_R2.fastq.gz
    ## dmel_rna_R2_G6_run2_R1.fastq.gz
    ## dmel_rna_R2_G6_run2_R2.fastq.gz
    ## dmel_rna_R3_G15_run2_R1.fastq.gz
    ## dmel_rna_R3_G15_run2_R2.fastq.gz
    ## dmel_rna_R3_G21_run2_R1.fastq.gz
    ## dmel_rna_R3_G21_run2_R2.fastq.gz
    ## dmel_rna_R3_G30_run2_R1.fastq.gz
    ## dmel_rna_R3_G30_run2_R2.fastq.gz
    ## dmel_rna_R3_G40_run2_R1.fastq.gz
    ## dmel_rna_R3_G40_run2_R2.fastq.gz
    ## dmel_rna_R3_G6_run2_R1.fastq.gz
    ## dmel_rna_R3_G6_run2_R2.fastq.gz

We then ran fastQC on all files to assess their quality.

``` bash
fastqc --outdir /Volumes/Data/Projects/dmelR2_p-ele/rna/raw/run1/fastQC /Volumes/Data/Projects/dmelR2_p-ele/rna/run1/raw/*fastq.gz
fastqc --outdir /Volumes/Data/Projects/dmelR2_p-ele/rna/raw/run2/fastQC /Volumes/Data/Projects/dmelR2_p-ele/rna/run2/raw/*fastq.gz
```

``` bash
cd /Volumes/Data/Projects/dmelR2_p-ele/rna/run1/fastQC/
ls

cd /Volumes/Data/Projects/dmelR2_p-ele/rna/run2/fastQC/
ls
```

    ## dmel_rna_R1_G15_run1_R1_fastqc.html
    ## dmel_rna_R1_G15_run1_R1_fastqc.zip
    ## dmel_rna_R1_G15_run1_R2_fastqc.html
    ## dmel_rna_R1_G15_run1_R2_fastqc.zip
    ## dmel_rna_R1_G21_run1_R1_fastqc.html
    ## dmel_rna_R1_G21_run1_R1_fastqc.zip
    ## dmel_rna_R1_G21_run1_R2_fastqc.html
    ## dmel_rna_R1_G21_run1_R2_fastqc.zip
    ## dmel_rna_R1_G30_run1_R1_fastqc.html
    ## dmel_rna_R1_G30_run1_R1_fastqc.zip
    ## dmel_rna_R1_G30_run1_R2_fastqc.html
    ## dmel_rna_R1_G30_run1_R2_fastqc.zip
    ## dmel_rna_R1_G40_run1_R1_fastqc.html
    ## dmel_rna_R1_G40_run1_R1_fastqc.zip
    ## dmel_rna_R1_G40_run1_R2_fastqc.html
    ## dmel_rna_R1_G40_run1_R2_fastqc.zip
    ## dmel_rna_R1_G6_run1_R1_fastqc.html
    ## dmel_rna_R1_G6_run1_R1_fastqc.zip
    ## dmel_rna_R1_G6_run1_R2_fastqc.html
    ## dmel_rna_R1_G6_run1_R2_fastqc.zip
    ## dmel_rna_R2_G15_run1_R1_fastqc.html
    ## dmel_rna_R2_G15_run1_R1_fastqc.zip
    ## dmel_rna_R2_G15_run1_R2_fastqc.html
    ## dmel_rna_R2_G15_run1_R2_fastqc.zip
    ## dmel_rna_R2_G21_run1_R1_fastqc.html
    ## dmel_rna_R2_G21_run1_R1_fastqc.zip
    ## dmel_rna_R2_G21_run1_R2_fastqc.html
    ## dmel_rna_R2_G21_run1_R2_fastqc.zip
    ## dmel_rna_R2_G30_run1_R1_fastqc.html
    ## dmel_rna_R2_G30_run1_R1_fastqc.zip
    ## dmel_rna_R2_G30_run1_R2_fastqc.html
    ## dmel_rna_R2_G30_run1_R2_fastqc.zip
    ## dmel_rna_R2_G40_run1_R1_fastqc.html
    ## dmel_rna_R2_G40_run1_R1_fastqc.zip
    ## dmel_rna_R2_G40_run1_R2_fastqc.html
    ## dmel_rna_R2_G40_run1_R2_fastqc.zip
    ## dmel_rna_R2_G6_run1_R1_fastqc.html
    ## dmel_rna_R2_G6_run1_R1_fastqc.zip
    ## dmel_rna_R2_G6_run1_R2_fastqc.html
    ## dmel_rna_R2_G6_run1_R2_fastqc.zip
    ## dmel_rna_R3_G15_run1_R1_fastqc.html
    ## dmel_rna_R3_G15_run1_R1_fastqc.zip
    ## dmel_rna_R3_G15_run1_R2_fastqc.html
    ## dmel_rna_R3_G15_run1_R2_fastqc.zip
    ## dmel_rna_R3_G21_run1_R1_fastqc.html
    ## dmel_rna_R3_G21_run1_R1_fastqc.zip
    ## dmel_rna_R3_G21_run1_R2_fastqc.html
    ## dmel_rna_R3_G21_run1_R2_fastqc.zip
    ## dmel_rna_R3_G30_run1_R1_fastqc.html
    ## dmel_rna_R3_G30_run1_R1_fastqc.zip
    ## dmel_rna_R3_G30_run1_R2_fastqc.html
    ## dmel_rna_R3_G30_run1_R2_fastqc.zip
    ## dmel_rna_R3_G40_run1_R1_fastqc.html
    ## dmel_rna_R3_G40_run1_R1_fastqc.zip
    ## dmel_rna_R3_G40_run1_R2_fastqc.html
    ## dmel_rna_R3_G40_run1_R2_fastqc.zip
    ## dmel_rna_R3_G6_run1_R1_fastqc.html
    ## dmel_rna_R3_G6_run1_R1_fastqc.zip
    ## dmel_rna_R3_G6_run1_R2_fastqc.html
    ## dmel_rna_R3_G6_run1_R2_fastqc.zip
    ## raw
    ## trimmed

After this, we then realised that there was significant proportion of
the 3’ read ends (~10-30%) which were adapter sequences. To alleviate
this, we decided to trim all of the 150bp reads down to 100bp, using the
following for loop.

``` bash
for file in *gz; do
    gzip -cd "$file" | cut -c-100 | awk 'NR%4==2 && length($0)>=100{print prev; print $0; getline; print $0; getline; print $0} {prev=$0}' | gzip -c > "trimmed/${file}"
done
```

# Reference genome - preparation

First, we obtained the reference D. mel fasta file from flybase, then
run the following command to remove everything but the Flybase ID from
the identifier line, while also attaching “\_mRNA” for downstream
analysis.

``` bash
less dmel-all-transcript-r6.52.fasta | cut -f1 -d";" | gsed 's/ type=/_/' > dmel-transcriptome-r6.52.fasta
```

We then merged it with the list of consensus D. mel TEs and indexed it.

``` bash
cd /Volumes/Data/Tools/RefGenomes/dmel/rna/dmel_TEs
ls
```

    ## dmel-transcriptome-r6.52-TEs.fasta
    ## dmel-transcriptome-r6.52-TEs.fasta.amb
    ## dmel-transcriptome-r6.52-TEs.fasta.ann
    ## dmel-transcriptome-r6.52-TEs.fasta.bwt
    ## dmel-transcriptome-r6.52-TEs.fasta.fai
    ## dmel-transcriptome-r6.52-TEs.fasta.pac
    ## dmel-transcriptome-r6.52-TEs.fasta.sa

Then we used bwa to map the forward and reverse reads to the reference.

``` bash
nohup zsh dmel_RNA_mapping_bwamem.sh > ../logs/dmel_rna_map.log
```

dmel_RNA_mapping_bwamem.sh -

``` bash
ref="/Volumes/Data/Tools/RefGenomes/dmel/rna/dmel_TEs/dmel-transcriptome-r6.52-TEs.fasta"
if="/Volumes/Data/Projects/dmelR2_p-ele/rna/raw/run2/trimmed"
of="/Volumes/Data/Projects/dmelR2_p-ele/rna/raw/run2/map-bwamem"

bwa mem -t 12 $ref $if/dmel_rna_R1_G6_run2_R1.fastq.gz $if/dmel_rna_R1_G6_run2_R2.fastq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_R1G6_run2.sort.bam
bwa mem -t 12 $ref $if/dmel_rna_R2_G6_run2_R1.fastq.gz $if/dmel_rna_R2_G6_run2_R2.fastq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_R2G6_run2.sort.bam
bwa mem -t 12 $ref $if/dmel_rna_R3_G6_run2_R1.fastq.gz $if/dmel_rna_R3_G6_run2_R2.fastq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_R3G6_run2.sort.bam

bwa mem -t 12 $ref $if/dmel_rna_R1_G15_run2_R1.fastq.gz $if/dmel_rna_R1_G15_run2_R2.fastq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_R1G15_run2.sort.bam
bwa mem -t 12 $ref $if/dmel_rna_R2_G15_run2_R1.fastq.gz $if/dmel_rna_R2_G15_run2_R2.fastq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_R2G15_run2.sort.bam
bwa mem -t 12 $ref $if/dmel_rna_R3_G15_run2_R1.fastq.gz $if/dmel_rna_R3_G15_run2_R2.fastq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_R3G15_run2.sort.bam

bwa mem -t 12 $ref $if/dmel_rna_R1_G21_run2_R1.fastq.gz $if/dmel_rna_R1_G21_run2_R2.fastq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_R1G21_run2.sort.bam
bwa mem -t 12 $ref $if/dmel_rna_R2_G21_run2_R1.fastq.gz $if/dmel_rna_R2_G21_run2_R2.fastq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_R2G21_run2.sort.bam
bwa mem -t 12 $ref $if/dmel_rna_R3_G21_run2_R1.fastq.gz $if/dmel_rna_R3_G21_run2_R2.fastq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_R3G21_run2.sort.bam

bwa mem -t 12 $ref $if/dmel_rna_R1_G30_run2_R1.fastq.gz $if/dmel_rna_R1_G30_run2_R2.fastq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_R1G30_run2.sort.bam
bwa mem -t 12 $ref $if/dmel_rna_R2_G30_run2_R1.fastq.gz $if/dmel_rna_R2_G30_run2_R2.fastq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_R2G30_run2.sort.bam
bwa mem -t 12 $ref $if/dmel_rna_R3_G30_run2_R1.fastq.gz $if/dmel_rna_R3_G30_run2_R2.fastq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_R3G30_run2.sort.bam

bwa mem -t 12 $ref $if/dmel_rna_R1_G40_run2_R1.fastq.gz $if/dmel_rna_R1_G40_run2_R2.fastq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_R1G40_run2.sort.bam
bwa mem -t 12 $ref $if/dmel_rna_R2_G40_run2_R1.fastq.gz $if/dmel_rna_R2_G40_run2_R2.fastq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_R2G40_run2.sort.bam
bwa mem -t 12 $ref $if/dmel_rna_R3_G40_run2_R1.fastq.gz $if/dmel_rna_R3_G40_run2_R2.fastq.gz | samtools sort -m 2G --output-fmt BAM --threads 2 -o $of/dmel_R3G40_run2.sort.bam
```

# GSNAP

Run on Vetgrid06

To map via GSNAP, we first need to create the database with gmap_build
using the reference genome.

``` bash
gmap_build -D . -d mel-transcriptome /Volumes/Temp2/Matt/rna/refgenomes/dmel/dmel_TEs/dmel-transcriptome-r6.52-TEs.fasta
```

Then we run the following script for alignment, assessing for novel
spicing events.

``` bash
#!/bin/bash

refdir="/Volumes/Temp2/Matt/rna/refgenomes/GMAP/mel-transcriptome"
refname="mel-transcriptome"

input_dir="/Volumes/Temp2/Matt/rna/trimmed/trimmed"
output_dir="/Volumes/Temp2/Matt/rna/map_GSNAP/output"

gsnap -D "$refdir" -d "$refname" -t 12 -N 1 --format=sam --gunzip "$input_dir/dmel_rna_R1_G6_run2_R1.fastq.gz" "$input_dir/dmel_rna_R1_G6_run2_R2.fastq.gz" | samtools sort -m 2G --output-fmt BAM --threads 2 -o "$output_dir/gt_R1G6.sort.bam"
gsnap -D "$refdir" -d "$refname" -t 12 -N 1 --format=sam --gunzip "$input_dir/dmel_rna_R2_G6_run2_R1.fastq.gz" "$input_dir/dmel_rna_R2_G6_run2_R2.fastq.gz" | samtools sort -m 2G --output-fmt BAM --threads 2 -o "$output_dir/gt_R2G6.sort.bam"
gsnap -D "$refdir" -d "$refname" -t 12 -N 1 --format=sam --gunzip "$input_dir/dmel_rna_R3_G6_run2_R1.fastq.gz" "$input_dir/dmel_rna_R3_G6_run2_R2.fastq.gz" | samtools sort -m 2G --output-fmt BAM --threads 2 -o "$output_dir/gt_R3G6.sort.bam"

gsnap -D "$refdir" -d "$refname" -t 12 -N 1 --format=sam --gunzip "$input_dir/dmel_rna_R1_G15_run2_R1.fastq.gz" "$input_dir/dmel_rna_R1_G15_run2_R2.fastq.gz" | samtools sort -m 2G --output-fmt BAM --threads 2 -o "$output_dir/gt_R1G15.sort.bam"
gsnap -D "$refdir" -d "$refname" -t 12 -N 1 --format=sam --gunzip "$input_dir/dmel_rna_R2_G15_run2_R1.fastq.gz" "$input_dir/dmel_rna_R2_G15_run2_R2.fastq.gz" | samtools sort -m 2G --output-fmt BAM --threads 2 -o "$output_dir/gt_R2G15.sort.bam"
gsnap -D "$refdir" -d "$refname" -t 12 -N 1 --format=sam --gunzip "$input_dir/dmel_rna_R3_G15_run2_R1.fastq.gz" "$input_dir/dmel_rna_R3_G15_run2_R2.fastq.gz" | samtools sort -m 2G --output-fmt BAM --threads 2 -o "$output_dir/gt_R3G15.sort.bam"

gsnap -D "$refdir" -d "$refname" -t 12 -N 1 --format=sam --gunzip "$input_dir/dmel_rna_R1_G21_run2_R1.fastq.gz" "$input_dir/dmel_rna_R1_G21_run2_R2.fastq.gz" | samtools sort -m 2G --output-fmt BAM --threads 2 -o "$output_dir/gt_R1G21.sort.bam"
gsnap -D "$refdir" -d "$refname" -t 12 -N 1 --format=sam --gunzip "$input_dir/dmel_rna_R2_G21_run2_R1.fastq.gz" "$input_dir/dmel_rna_R2_G21_run2_R2.fastq.gz" | samtools sort -m 2G --output-fmt BAM --threads 2 -o "$output_dir/gt_R2G21.sort.bam"
gsnap -D "$refdir" -d "$refname" -t 12 -N 1 --format=sam --gunzip "$input_dir/dmel_rna_R3_G21_run2_R1.fastq.gz" "$input_dir/dmel_rna_R3_G21_run2_R2.fastq.gz" | samtools sort -m 2G --output-fmt BAM --threads 2 -o "$output_dir/gt_R3G21.sort.bam"

gsnap -D "$refdir" -d "$refname" -t 12 -N 1 --format=sam --gunzip "$input_dir/dmel_rna_R1_G30_run2_R1.fastq.gz" "$input_dir/dmel_rna_R1_G30_run2_R2.fastq.gz" | samtools sort -m 2G --output-fmt BAM --threads 2 -o "$output_dir/gt_R1G30.sort.bam"
gsnap -D "$refdir" -d "$refname" -t 12 -N 1 --format=sam --gunzip "$input_dir/dmel_rna_R2_G30_run2_R1.fastq.gz" "$input_dir/dmel_rna_R2_G30_run2_R2.fastq.gz" | samtools sort -m 2G --output-fmt BAM --threads 2 -o "$output_dir/gt_R2G30.sort.bam"
gsnap -D "$refdir" -d "$refname" -t 12 -N 1 --format=sam --gunzip "$input_dir/dmel_rna_R3_G30_run2_R1.fastq.gz" "$input_dir/dmel_rna_R3_G30_run2_R2.fastq.gz" | samtools sort -m 2G --output-fmt BAM --threads 2 -o "$output_dir/gt_R3G30.sort.bam"

gsnap -D "$refdir" -d "$refname" -t 12 -N 1 --format=sam --gunzip "$input_dir/dmel_rna_R1_G40_run2_R1.fastq.gz" "$input_dir/dmel_rna_R1_G40_run2_R2.fastq.gz" | samtools sort -m 2G --output-fmt BAM --threads 2 -o "$output_dir/gt_R1G40.sort.bam"
gsnap -D "$refdir" -d "$refname" -t 12 -N 1 --format=sam --gunzip "$input_dir/dmel_rna_R2_G40_run2_R1.fastq.gz" "$input_dir/dmel_rna_R2_G40_run2_R2.fastq.gz" | samtools sort -m 2G --output-fmt BAM --threads 2 -o "$output_dir/gt_R2G40.sort.bam"
gsnap -D "$refdir" -d "$refname" -t 12 -N 1 --format=sam --gunzip "$input_dir/dmel_rna_R3_G40_run2_R1.fastq.gz" "$input_dir/dmel_rna_R3_G40_run2_R2.fastq.gz" | samtools sort -m 2G --output-fmt BAM --threads 2 -o "$output_dir/gt_R3G40.sort.bam"
```

Which provides the following sorted BAM files.

``` bash
cd /Volumes/Data/Projects/dmelR2_p-ele/rna/run2/map-GSNAP/output_genome
ls
```
