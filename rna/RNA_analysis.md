RNA analysis
================
Matthew Beaumont
2023-06-05

``` bash
knitr::opts_chunk$set(echo = TRUE)
```

We obtained the following mRNA reads:

``` bash
cd /Volumes/Data/Projects/dmelR2_p-ele/rna/run1/raw
ls -lh *.fastq.gz

cd /Volumes/Data/Projects/dmelR2_p-ele/rna/run2/raw
ls -lh *.fastq.gz
```

    ## -rw-r--r--  1 mbeaumont  staff   1.1G May 25 15:21 dmel_rna_R1_G15_run1_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.2G May 25 15:22 dmel_rna_R1_G15_run1_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.1G May 25 15:22 dmel_rna_R1_G21_run1_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.1G May 25 15:22 dmel_rna_R1_G21_run1_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.0G May 25 15:23 dmel_rna_R1_G30_run1_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.1G May 25 15:23 dmel_rna_R1_G30_run1_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.1G May 25 15:23 dmel_rna_R1_G40_run1_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.2G May 25 15:23 dmel_rna_R1_G40_run1_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   957M May 25 15:24 dmel_rna_R1_G6_run1_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.0G May 25 15:24 dmel_rna_R1_G6_run1_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.1G May 25 15:24 dmel_rna_R2_G15_run1_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.2G May 25 15:24 dmel_rna_R2_G15_run1_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   985M May 25 15:25 dmel_rna_R2_G21_run1_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.0G May 25 15:25 dmel_rna_R2_G21_run1_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   987M May 25 15:25 dmel_rna_R2_G30_run1_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.0G May 25 15:25 dmel_rna_R2_G30_run1_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.1G May 25 15:26 dmel_rna_R2_G40_run1_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.2G May 25 15:26 dmel_rna_R2_G40_run1_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.0G May 25 15:26 dmel_rna_R2_G6_run1_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.0G May 25 15:27 dmel_rna_R2_G6_run1_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.0G May 25 15:27 dmel_rna_R3_G15_run1_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.0G May 25 15:27 dmel_rna_R3_G15_run1_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.1G May 25 15:27 dmel_rna_R3_G21_run1_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.1G May 25 15:28 dmel_rna_R3_G21_run1_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.1G May 25 15:28 dmel_rna_R3_G30_run1_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.2G May 25 15:28 dmel_rna_R3_G30_run1_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.1G May 25 15:29 dmel_rna_R3_G40_run1_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.2G May 25 15:29 dmel_rna_R3_G40_run1_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.0G May 25 15:29 dmel_rna_R3_G6_run1_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.1G May 25 15:30 dmel_rna_R3_G6_run1_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.1G May 25 15:31 dmel_rna_R1_G15_run2_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.2G May 25 15:32 dmel_rna_R1_G15_run2_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.1G May 25 15:32 dmel_rna_R1_G21_run2_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.2G May 25 15:33 dmel_rna_R1_G21_run2_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.8G May 25 15:33 dmel_rna_R1_G30_run2_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.0G May 25 15:34 dmel_rna_R1_G30_run2_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.1G May 25 15:34 dmel_rna_R1_G40_run2_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.2G May 25 15:35 dmel_rna_R1_G40_run2_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.8G May 25 15:35 dmel_rna_R1_G6_run2_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.9G May 25 15:36 dmel_rna_R1_G6_run2_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.2G May 25 15:36 dmel_rna_R2_G15_run2_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.3G May 25 15:37 dmel_rna_R2_G15_run2_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.9G May 25 15:38 dmel_rna_R2_G21_run2_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.0G May 25 15:38 dmel_rna_R2_G21_run2_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.8G May 25 15:38 dmel_rna_R2_G30_run2_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.9G May 25 15:39 dmel_rna_R2_G30_run2_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.2G May 25 15:39 dmel_rna_R2_G40_run2_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.3G May 25 15:40 dmel_rna_R2_G40_run2_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.9G May 25 15:40 dmel_rna_R2_G6_run2_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.0G May 25 15:41 dmel_rna_R2_G6_run2_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.9G May 25 15:41 dmel_rna_R3_G15_run2_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.0G May 25 15:42 dmel_rna_R3_G15_run2_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.0G May 25 15:42 dmel_rna_R3_G21_run2_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.2G May 25 15:43 dmel_rna_R3_G21_run2_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.2G May 25 15:43 dmel_rna_R3_G30_run2_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.4G May 25 15:44 dmel_rna_R3_G30_run2_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.2G May 25 15:45 dmel_rna_R3_G40_run2_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.4G May 25 15:45 dmel_rna_R3_G40_run2_R2.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   1.9G May 25 15:46 dmel_rna_R3_G6_run2_R1.fastq.gz
    ## -rw-r--r--  1 mbeaumont  staff   2.0G May 25 15:46 dmel_rna_R3_G6_run2_R2.fastq.gz

We then ran fastQC on all files to assess their quality.

``` bash
fastqc --outdir /Volumes/Data/Projects/dmelR2_p-ele/rna/raw/run1/fastQC /Volumes/Data/Projects/dmelR2_p-ele/rna/run1/raw/*fastq.gz
fastqc --outdir /Volumes/Data/Projects/dmelR2_p-ele/rna/raw/run2/fastQC /Volumes/Data/Projects/dmelR2_p-ele/rna/run2/raw/*fastq.gz
```

``` bash
cd /Volumes/Data/Projects/dmelR2_p-ele/rna/run1/fastQC/
ls -lh

cd /Volumes/Data/Projects/dmelR2_p-ele/rna/run2/fastQC/
ls -lh 
```

    ## total 80048
    ## -rw-r--r--  1 mbeaumont  staff   708K May 31 13:10 dmel_rna_R1_G15_run1_R1_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   624K May 31 13:10 dmel_rna_R1_G15_run1_R1_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   701K May 31 13:12 dmel_rna_R1_G15_run1_R2_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   616K May 31 13:12 dmel_rna_R1_G15_run1_R2_fastqc.zip
    ## -rw-r--r--@ 1 mbeaumont  staff   710K May 31 13:15 dmel_rna_R1_G21_run1_R1_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   623K May 31 13:15 dmel_rna_R1_G21_run1_R1_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   706K May 31 13:17 dmel_rna_R1_G21_run1_R2_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   619K May 31 13:17 dmel_rna_R1_G21_run1_R2_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   705K May 31 13:20 dmel_rna_R1_G30_run1_R1_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   616K May 31 13:20 dmel_rna_R1_G30_run1_R1_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   706K May 31 13:22 dmel_rna_R1_G30_run1_R2_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   622K May 31 13:22 dmel_rna_R1_G30_run1_R2_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   706K May 31 13:24 dmel_rna_R1_G40_run1_R1_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   619K May 31 13:24 dmel_rna_R1_G40_run1_R1_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   702K May 31 13:27 dmel_rna_R1_G40_run1_R2_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   620K May 31 13:27 dmel_rna_R1_G40_run1_R2_fastqc.zip
    ## -rw-r--r--@ 1 mbeaumont  staff   705K May 31 13:29 dmel_rna_R1_G6_run1_R1_fastqc.html
    ## -rw-r--r--@ 1 mbeaumont  staff   618K May 31 13:29 dmel_rna_R1_G6_run1_R1_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   704K May 31 13:31 dmel_rna_R1_G6_run1_R2_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   617K May 31 13:31 dmel_rna_R1_G6_run1_R2_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   710K May 31 13:34 dmel_rna_R2_G15_run1_R1_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   623K May 31 13:34 dmel_rna_R2_G15_run1_R1_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   706K May 31 13:36 dmel_rna_R2_G15_run1_R2_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   620K May 31 13:36 dmel_rna_R2_G15_run1_R2_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   707K May 31 13:39 dmel_rna_R2_G21_run1_R1_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   617K May 31 13:39 dmel_rna_R2_G21_run1_R1_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   701K May 31 13:41 dmel_rna_R2_G21_run1_R2_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   614K May 31 13:41 dmel_rna_R2_G21_run1_R2_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   707K May 31 13:43 dmel_rna_R2_G30_run1_R1_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   619K May 31 13:43 dmel_rna_R2_G30_run1_R1_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   703K May 31 13:45 dmel_rna_R2_G30_run1_R2_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   616K May 31 13:45 dmel_rna_R2_G30_run1_R2_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   710K May 31 13:48 dmel_rna_R2_G40_run1_R1_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   623K May 31 13:48 dmel_rna_R2_G40_run1_R1_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   704K May 31 13:51 dmel_rna_R2_G40_run1_R2_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   621K May 31 13:51 dmel_rna_R2_G40_run1_R2_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   706K May 31 13:53 dmel_rna_R2_G6_run1_R1_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   618K May 31 13:53 dmel_rna_R2_G6_run1_R1_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   703K May 31 13:55 dmel_rna_R2_G6_run1_R2_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   615K May 31 13:55 dmel_rna_R2_G6_run1_R2_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   707K May 31 13:58 dmel_rna_R3_G15_run1_R1_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   621K May 31 13:58 dmel_rna_R3_G15_run1_R1_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   703K May 31 14:00 dmel_rna_R3_G15_run1_R2_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   617K May 31 14:00 dmel_rna_R3_G15_run1_R2_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   713K May 31 14:02 dmel_rna_R3_G21_run1_R1_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   624K May 31 14:02 dmel_rna_R3_G21_run1_R1_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   705K May 31 14:05 dmel_rna_R3_G21_run1_R2_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   621K May 31 14:05 dmel_rna_R3_G21_run1_R2_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   707K May 31 14:08 dmel_rna_R3_G30_run1_R1_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   623K May 31 14:08 dmel_rna_R3_G30_run1_R1_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   705K May 31 14:10 dmel_rna_R3_G30_run1_R2_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   620K May 31 14:10 dmel_rna_R3_G30_run1_R2_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   712K May 31 14:13 dmel_rna_R3_G40_run1_R1_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   629K May 31 14:13 dmel_rna_R3_G40_run1_R1_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   706K May 31 14:15 dmel_rna_R3_G40_run1_R2_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   625K May 31 14:15 dmel_rna_R3_G40_run1_R2_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   708K May 31 14:17 dmel_rna_R3_G6_run1_R1_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   618K May 31 14:17 dmel_rna_R3_G6_run1_R1_fastqc.zip
    ## -rw-r--r--  1 mbeaumont  staff   702K May 31 14:20 dmel_rna_R3_G6_run1_R2_fastqc.html
    ## -rw-r--r--  1 mbeaumont  staff   618K May 31 14:20 dmel_rna_R3_G6_run1_R2_fastqc.zip
    ## total 0
    ## drwxr-xr-x  62 mbeaumont  staff   2.1K Jun  1 18:23 raw
    ## drwxr-xr-x  62 mbeaumont  staff   2.1K Jun  1 19:58 trimmed

After this, we then realised that there was significant proportion of
the 3’ read ends (\~10-30%) which were adapter sequences. To alleviate
this, we decided to trim all of the 150bp reads down to 100bp, using the
following for loop.

``` bash
for file in *gz; do
    gzip -cd "$file" | cut -c-100 | awk 'NR%4==2 && length($0)>=100{print prev; print $0; getline; print $0; getline; print $0} {prev=$0}' | gzip -c > "trimmed/${file}"
done
```

We then obtained the reference D. mel transcriptome and merged it with
the list of consensus D. mel TEs, then indexed it.

``` bash
cd /Volumes/Data/Tools/RefGenomes/dmel/rna/dmel_TEs
ls -lh
```

    ## total 518704
    ## -rw-r--r--@ 1 mbeaumont  staff    93M May 30 17:15 dmel_1215.4_rel_6_ISO1_consensusTEs.fasta
    ## -rw-r--r--  1 mbeaumont  staff   642B Jun  1 11:13 dmel_1215.4_rel_6_ISO1_consensusTEs.fasta.amb
    ## -rw-r--r--  1 mbeaumont  staff   3.6M Jun  1 11:13 dmel_1215.4_rel_6_ISO1_consensusTEs.fasta.ann
    ## -rw-r--r--  1 mbeaumont  staff    89M Jun  1 11:13 dmel_1215.4_rel_6_ISO1_consensusTEs.fasta.bwt
    ## -rw-r--r--  1 mbeaumont  staff   1.1M Jun  2 17:25 dmel_1215.4_rel_6_ISO1_consensusTEs.fasta.fai
    ## -rw-r--r--  1 mbeaumont  staff    22M Jun  1 11:13 dmel_1215.4_rel_6_ISO1_consensusTEs.fasta.pac
    ## -rw-r--r--  1 mbeaumont  staff    44M Jun  1 11:14 dmel_1215.4_rel_6_ISO1_consensusTEs.fasta.sa

Then we used bwa to map the forward and reverse reads to the reference.

``` bash
nohup zsh dmel_RNA_mapping_bwamem.sh > ../logs/dmel_rna_map.log
```

dmel_RNA_mapping_bwamem.sh -

``` bash
ref="/Volumes/Data/Tools/RefGenomes/dmel/rna/dmel_TEs/dmel_1215.4_rel_6_ISO1_consensusTEs.fasta"
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

To map via GMAP, we first need to create the database with the reference
genome.

# Coverage

First we indexed the D. mel + TE fasta file using samtools

``` bash
samtools faidx dmel_1215.4_rel_6_ISO1_consensusTEs.fasta
```

# Coverage

We ran the following script to assess coverage

``` bash
nohup zsh expression-splicing.sh > ../logs/splicing-expression.log
```

``` bash
#!/bin/bash

fai="/Volumes/Data/Tools/RefGenomes/dmel/rna/dmel_TEs/dmel_1215.4_rel_6_ISO1_consensusTEs.fasta.fai"
pyscript="/Volumes/Data/Projects/dmelR2_p-ele/scripts/mRNA-coverage-senseantisense.py"
input_dir="/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/map-bwamem"
output_dir="/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/splicing-expression/raw_expression"
seqs="PPI251,FBtr0141238_mRNA,FBtr0141090_mRNA,FBtr0132023_mRNA,FBtr0137012_mRNA,FBtr0140282_mRNA,FBtr0141095_mRNA,FBtr0142660_mRNA,FBtr0140302_mRNA,FBtr0142342_mRNA,FBtr0143902_mRNA,FBtr0130551_mRNA,FBtr0145198_mRNA,FBtr0130342_mRNA,FBtr0135961_mRNA,FBtr0135197_mRNA,FBtr0130390_mRNA,FBtr0141576_mRNA,FBtr0139614_mRNA,FBtr0130391_mRNA,FBtr0143034_mRNA"

samtools view $input_dir/dmel_R1G6_run2.sort.bam | python $pyscript --sam - --sample-id R1-G6-run2 --seqs $seqs --fai $fai > $output_dir/dmel_R1G6_run2.txt
samtools view $input_dir/dmel_R2G6_run2.sort.bam | python $pyscript --sam - --sample-id R2-G6-run2 --seqs $seqs --fai $fai > $output_dir/dmel_R2G6_run2.txt
samtools view $input_dir/dmel_R3G6_run2.sort.bam | python $pyscript --sam - --sample-id R3-G6-run2 --seqs $seqs --fai $fai > $output_dir/dmel_R3G6_run2.txt

samtools view $input_dir/dmel_R1G15_run2.sort.bam | python $pyscript --sam - --sample-id R1-G15-run2 --seqs $seqs --fai $fai > $output_dir/dmel_R1G15_run2.txt
samtools view $input_dir/dmel_R2G15_run2.sort.bam | python $pyscript --sam - --sample-id R2-G15-run2 --seqs $seqs --fai $fai > $output_dir/dmel_R2G15_run2.txt
samtools view $input_dir/dmel_R3G15_run2.sort.bam | python $pyscript --sam - --sample-id R3-G15-run2 --seqs $seqs --fai $fai > $output_dir/dmel_R3G15_run2.txt

samtools view $input_dir/dmel_R1G21_run2.sort.bam | python $pyscript --sam - --sample-id R1-G21-run2 --seqs $seqs --fai $fai > $output_dir/dmel_R1G21_run2.txt
samtools view $input_dir/dmel_R2G21_run2.sort.bam | python $pyscript --sam - --sample-id R2-G21-run2 --seqs $seqs --fai $fai > $output_dir/dmel_R2G21_run2.txt
samtools view $input_dir/dmel_R3G21_run2.sort.bam | python $pyscript --sam - --sample-id R3-G21-run2 --seqs $seqs --fai $fai > $output_dir/dmel_R3G21_run2.txt

samtools view $input_dir/dmel_R1G30_run2.sort.bam | python $pyscript --sam - --sample-id R1-G30-run2 --seqs $seqs --fai $fai > $output_dir/dmel_R1G30_run2.txt
samtools view $input_dir/dmel_R2G30_run2.sort.bam | python $pyscript --sam - --sample-id R2-G30-run2 --seqs $seqs --fai $fai > $output_dir/dmel_R2G30_run2.txt
samtools view $input_dir/dmel_R3G30_run2.sort.bam | python $pyscript --sam - --sample-id R3-G30-run2 --seqs $seqs --fai $fai > $output_dir/dmel_R3G30_run2.txt

samtools view $input_dir/dmel_R1G40_run2.sort.bam | python $pyscript --sam - --sample-id R1-G40-run2 --seqs $seqs --fai $fai > $output_dir/dmel_R1G40_run2.txt
samtools view $input_dir/dmel_R2G40_run2.sort.bam | python $pyscript --sam - --sample-id R2-G40-run2 --seqs $seqs --fai $fai > $output_dir/dmel_R2G40_run2.txt
samtools view $input_dir/dmel_R3G40_run2.sort.bam | python $pyscript --sam - --sample-id R3-G40-run2 --seqs $seqs --fai $fai > $output_dir/dmel_R3G40_run2.txt

cat *.txt| perl -pe 's/-/\t/'|perl -pe 's/-/\t/' > expr-spli.forr
```

# P-element exrpession

We ran the following script to extract out all expression levels.

``` bash
nohup zsh /Volumes/Data/Projects/dmelR2_p-ele/scripts/dmel_pele_expression.sh > ../../../../../logs/pele-expression.log
```

``` bash
#!/bin/bash

fai="/Volumes/Data/Tools/RefGenomes/dmel/rna/dmel_TEs/dmel_1215.4_rel_6_ISO1_consensusTEs.fasta.fai"
pyscript="/Volumes/Data/Projects/dmelR2_p-ele/scripts/mRNA-expression.py"
input_dir="/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/map-bwamem"
output_dir="/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/splicing-expression/raw_expression/all-expressionlevel"

samtools view $input_dir/dmel_R1G6_run2.sort.bam | python $pyscript --sam - --sample-id R1-G6-run2 --fai $fai > $output_dir/expr_dmel_R1G6_run2.txt
samtools view $input_dir/dmel_R2G6_run2.sort.bam | python $pyscript --sam - --sample-id R2-G6-run2 --fai $fai > $output_dir/expr_dmel_R2G6_run2.txt
samtools view $input_dir/dmel_R3G6_run2.sort.bam | python $pyscript --sam - --sample-id R3-G6-run2 --fai $fai > $output_dir/expr_dmel_R3G6_run2.txt

samtools view $input_dir/dmel_R1G15_run2.sort.bam | python $pyscript --sam - --sample-id R1-G15-run2 --fai $fai > $output_dir/expr_dmel_R1G15_run2.txt
samtools view $input_dir/dmel_R2G15_run2.sort.bam | python $pyscript --sam - --sample-id R2-G15-run2 --fai $fai > $output_dir/expr_dmel_R2G15_run2.txt
samtools view $input_dir/dmel_R3G15_run2.sort.bam | python $pyscript --sam - --sample-id R3-G15-run2 --fai $fai > $output_dir/expr_dmel_R3G15_run2.txt

samtools view $input_dir/dmel_R1G21_run2.sort.bam | python $pyscript --sam - --sample-id R1-G21-run2 --fai $fai > $output_dir/expr_dmel_R1G21_run2.txt
samtools view $input_dir/dmel_R2G21_run2.sort.bam | python $pyscript --sam - --sample-id R2-G21-run2 --fai $fai > $output_dir/expr_dmel_R2G21_run2.txt
samtools view $input_dir/dmel_R3G21_run2.sort.bam | python $pyscript --sam - --sample-id R3-G21-run2 --fai $fai > $output_dir/expr_dmel_R3G21_run2.txt

samtools view $input_dir/dmel_R1G30_run2.sort.bam | python $pyscript --sam - --sample-id R1-G30-run2 --fai $fai > $output_dir/expr_dmel_R1G30_run2.txt
samtools view $input_dir/dmel_R2G30_run2.sort.bam | python $pyscript --sam - --sample-id R2-G30-run2 --fai $fai > $output_dir/expr_dmel_R2G30_run2.txt
samtools view $input_dir/dmel_R3G30_run2.sort.bam | python $pyscript --sam - --sample-id R3-G30-run2 --fai $fai > $output_dir/expr_dmel_R3G30_run2.txt

samtools view $input_dir/dmel_R1G40_run2.sort.bam | python $pyscript --sam - --sample-id R1-G40-run2 --fai $fai > $output_dir/expr_dmel_R1G40_run2.txt
samtools view $input_dir/dmel_R2G40_run2.sort.bam | python $pyscript --sam - --sample-id R2-G40-run2 --fai $fai > $output_dir/expr_dmel_R2G40_run2.txt
samtools view $input_dir/dmel_R3G40_run2.sort.bam | python $pyscript --sam - --sample-id R3-G40-run2 --fai $fai > $output_dir/expr_dmel_R3G40_run2.txt

cat *.txt| perl -pe 's/-/\t/'|perl -pe 's/-/\t/' > expr.forr
```

Then we visualise it using ggplot.

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
    ## ✔ tibble  3.1.7     ✔ dplyr   1.0.9
    ## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
    ## ✔ readr   2.1.2     ✔ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(RColorBrewer)
theme_set(theme_bw())
tresrep<-c("#e41a1c","#377eb8","#4daf4a")

t=read_delim("/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/splicing-expression/raw_expression/all-expressionlevel/expr.forr",delim="\t",col_names=FALSE,comment="#")
```

    ## Rows: 303364 Columns: 10

    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (4): X1, X2, X3, X4
    ## dbl (6): X5, X6, X7, X8, X9, X10
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
names(t)<-c("replicate","generation","run","gene","rawse","rawase","genlen","sense","antisense","total")
t<-subset(t,gene=="PPI251")
t$generation<-as.numeric(substring(t$generation, 2))


s<-ggplot()+geom_line(data=t,aes(x=generation,y=sense,color=replicate),size=1)+
  geom_point(data=t,aes(x=generation,y=sense,color=replicate))+
  theme(strip.text=element_blank(),legend.position=c(0.1,0.9))+
  ylab("expression [rpkm]")+scale_colour_manual(values=tresrep)+xlim(0,48)

plot(s)
```

![](RNA_analysis_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->
