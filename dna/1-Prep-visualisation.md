1 - DNA analysis
================
Matthew Beaumont
2023-10-19

``` bash
knitr::opts_chunk$set(echo = TRUE)
```

# Read trimming

All of the following .fq.gz read files were trimmed to the minimum read
length of the shortest reads from all files (100) for uniformity.

``` bash
cd /Volumes/Data/Projects/DmelR2_P-ele/dna/fastq/raw 
ls
```

    ## dmel_naive_CATTTT_1.fq.gz
    ## dmel_naive_CATTTT_2.fq.gz
    ## mel_R1G01_CGATGT_1.fq.gz
    ## mel_R1G01_CGATGT_2.fq.gz
    ## mel_R1G10_ACAGTG_1.fq.gz
    ## mel_R1G10_ACAGTG_2.fq.gz
    ## mel_R1G20_ACTTGA_1.fq.gz
    ## mel_R1G20_ACTTGA_2.fq.gz
    ## mel_R1G34_GTAGAG_1.fq.gz
    ## mel_R1G34_GTAGAG_2.fq.gz
    ## mel_R1G40_GTAGAG_1.fq.gz
    ## mel_R1G40_GTAGAG_2.fq.gz
    ## mel_R1G48_GTAGAG_1.fq.gz
    ## mel_R1G48_GTAGAG_2.fq.gz
    ## mel_R1G63_ATCACG_1.fq.gz
    ## mel_R1G63_ATCACG_2.fq.gz
    ## mel_R2G01_TTAGGC_1.fq.gz
    ## mel_R2G01_TTAGGC_2.fq.gz
    ## mel_R2G10_GCCAAT_1.fq.gz
    ## mel_R2G10_GCCAAT_2.fq.gz
    ## mel_R2G20_GGCTAC_1.fq.gz
    ## mel_R2G20_GGCTAC_2.fq.gz
    ## mel_R2G34_GGTAGC_1.fq.gz
    ## mel_R2G34_GGTAGC_2.fq.gz
    ## mel_R2G40_GGTAGC_1.fq.gz
    ## mel_R2G40_GGTAGC_2.fq.gz
    ## mel_R2G48_GGTAGC_1.fq.gz
    ## mel_R2G48_GGTAGC_2.fq.gz
    ## mel_R2G63_CGATGT_1.fq.gz
    ## mel_R2G63_CGATGT_2.fq.gz
    ## mel_R3G01_TGACCA_1.fq.gz
    ## mel_R3G01_TGACCA_2.fq.gz
    ## mel_R3G10_CAGATC_1.fq.gz
    ## mel_R3G10_CAGATC_2.fq.gz
    ## mel_R3G20_CTTGTA_1.fq.gz
    ## mel_R3G20_CTTGTA_2.fq.gz
    ## mel_R3G34_ATGAGC_1.fq.gz
    ## mel_R3G34_ATGAGC_2.fq.gz
    ## mel_R3G40_ATGAGC_1.fq.gz
    ## mel_R3G40_ATGAGC_2.fq.gz
    ## mel_R3G48_ATGAGC_1.fq.gz
    ## mel_R3G48_ATGAGC_2.fq.gz
    ## mel_R3G63_TTAGGC_1.fq.gz
    ## mel_R3G63_TTAGGC_2.fq.gz

``` bash
nohup bash scripts/trim.sh > trim.log &
```

``` bash
for i in fastq/*.fq.gz
        do 
                date
                j=`basename $i`
                j=${j%.fq.gz}
                zless $i | awk '{print substr($1,1,100)}' | gzip -9 > fastq/trimmed/${j}_trimmed.fq.gz
                date
        done
```

The following trimmed .fq.gz files were generated.

``` bash
cd /Volumes/Data/Projects/DmelR2_P-ele/dna/fastq/trimmed
ls
```

    ## dmel_naive_CATTTT_1_trimmed.fq.gz
    ## dmel_naive_CATTTT_2_trimmed.fq.gz
    ## mel_R1G01_CGATGT_1_trimmed.fq.gz
    ## mel_R1G01_CGATGT_2_trimmed.fq.gz
    ## mel_R1G10_ACAGTG_1_trimmed.fq.gz
    ## mel_R1G10_ACAGTG_2_trimmed.fq.gz
    ## mel_R1G20_ACTTGA_1_trimmed.fq.gz
    ## mel_R1G20_ACTTGA_2_trimmed.fq.gz
    ## mel_R1G34_GTAGAG_1_trimmed.fq.gz
    ## mel_R1G34_GTAGAG_2_trimmed.fq.gz
    ## mel_R1G40_GTAGAG_1_trimmed.fq.gz
    ## mel_R1G40_GTAGAG_2_trimmed.fq.gz
    ## mel_R1G48_GTAGAG_1_trimmed.fq.gz
    ## mel_R1G48_GTAGAG_2_trimmed.fq.gz
    ## mel_R1G63_ATCACG_1_trimmed.fq.gz
    ## mel_R1G63_ATCACG_2_trimmed.fq.gz
    ## mel_R2G01_TTAGGC_1_trimmed.fq.gz
    ## mel_R2G01_TTAGGC_2_trimmed.fq.gz
    ## mel_R2G10_GCCAAT_1_trimmed.fq.gz
    ## mel_R2G10_GCCAAT_2_trimmed.fq.gz
    ## mel_R2G20_GGCTAC_1_trimmed.fq.gz
    ## mel_R2G20_GGCTAC_2_trimmed.fq.gz
    ## mel_R2G34_GGTAGC_1_trimmed.fq.gz
    ## mel_R2G34_GGTAGC_2_trimmed.fq.gz
    ## mel_R2G40_GGTAGC_1_trimmed.fq.gz
    ## mel_R2G40_GGTAGC_2_trimmed.fq.gz
    ## mel_R2G48_GGTAGC_1_trimmed.fq.gz
    ## mel_R2G48_GGTAGC_2_trimmed.fq.gz
    ## mel_R2G63_CGATGT_1_trimmed.fq.gz
    ## mel_R2G63_CGATGT_2_trimmed.fq.gz
    ## mel_R3G01_TGACCA_1_trimmed.fq.gz
    ## mel_R3G01_TGACCA_2_trimmed.fq.gz
    ## mel_R3G10_CAGATC_1_trimmed.fq.gz
    ## mel_R3G10_CAGATC_2_trimmed.fq.gz
    ## mel_R3G20_CTTGTA_1_trimmed.fq.gz
    ## mel_R3G20_CTTGTA_2_trimmed.fq.gz
    ## mel_R3G34_ATGAGC_1_trimmed.fq.gz
    ## mel_R3G34_ATGAGC_2_trimmed.fq.gz
    ## mel_R3G40_ATGAGC_1_trimmed.fq.gz
    ## mel_R3G40_ATGAGC_2_trimmed.fq.gz
    ## mel_R3G48_ATGAGC_1_trimmed.fq.gz
    ## mel_R3G48_ATGAGC_2_trimmed.fq.gz
    ## mel_R3G63_TTAGGC_1_trimmed.fq.gz
    ## mel_R3G63_TTAGGC_2_trimmed.fq.gz

# TE Mapping

The fast-miner.sh script was altered to accommodate the altered file
name for the trimmed read files.

``` bash
nohup zsh fastq-miner-trimmed.sh mel /Volumes/Data/Projects/DmelR2_P-ele/fastq/trimmed > /Volumes/Data/Projects/DmelR2_P-ele/logs/mel.log &
```

``` bash
if [ $# -lt 2 ]
  # "$#" is number of parameters- here we test
   then
   echo "Usage $0 sampleid inputdir"
exit 2
fi

set -o shwordsplit

source ~/.zshrc

# Files and folders
sampleid=$1
inputdir=$2
outputdir=$2
outabu="../results/$sampleid/rpm"
outbamall="../results/$sampleid/bam"
tmpdir="../tmp"

## the REFERENCE 
refg="../refg/Dmel_tes_3scg/teseqs-3scg.fasta"
hier="../refg/Dmel_tes_3scg/teseqs.hier"

# software
samtools="../bin/samtools"
bwa="../bin/bwa"
popte2="../bin/popte2.jar"
readtorpm="../bin/readstat2rpm_all.py"

mkdir -p $outabu
mkdir -p $outbamall
# paths
for read1 in $inputdir/*_1_trimmed.fq.gz
do
    # variables defined to store the sra ids of the files from the input directory
        
    read2=${read1%_1_trimmed.fq.gz}
    read2=${read2}_2_trimmed.fq.gz
    tmp=`basename $read1`
        readid=${tmp%_1.fq.gz}

    tmpfile=$tmpdir/${readid}.fq.gz
    
    gzip -cd $read1 $read2 | paste - - - - |awk '{print "@" NR,$2,"+" NR,$4}'|tr " " "\n" |gzip -c > $tmpfile
    #gzip -cd $read1 $read2 |gzip -c > $tmpfile
 
        # mapping
    bamfile="$tmpdir/$readid.sort.bam"
    command="$bwa bwasw -t 8 -M $refg $tmpfile | $samtools view -Sb - | $samtools sort -T $sraid.nnnn.bam -O bam -m 4G > $bamfile"
    echo "executing mapping $command"
    eval $command
 
    # TE bam
    allbamfile="$outbamall/$readid.allte.sort.bam"
    samtools view -b -F 0x004 $bamfile > $allbamfile
    
    # PopoolationTE2
    opteabu="$outabu/$readid.rawabu"
    opterpm="$outabu/$readid.rpm"
    echo "Writting raw abundance to $opteabu"
    java -jar $popte2 stat-reads --bam $bamfile --map-qual 10 --hier $hier --output $opteabu
    echo "Writting rpm to $opterpm"
    python $readtorpm --rs $opteabu > $opterpm 

    # Cleanup tmp
        rm $bamfile
    rm $tmpfile
done 
```

This altered fast-miner.sh TEMiner script was then used to generate the
following .bam and .rpm files from the trimmed .fq.gz read files.

``` bash
cd /Volumes/Data/Projects/DmelR2_P-ele/software/teminer-code-r10/results/mel/bam
ls

cd /Volumes/Data/Projects/DmelR2_P-ele/software/teminer-code-r10/results/mel/rpm
ls 
```

    ## dmel_naive_CATTTT_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R1G01_CGATGT_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R1G10_ACAGTG_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R1G20_ACTTGA_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R1G34_GTAGAG_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R1G40_GTAGAG_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R1G48_GTAGAG_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R1G63_ATCACG_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R2G01_TTAGGC_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R2G10_GCCAAT_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R2G20_GGCTAC_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R2G34_GGTAGC_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R2G40_GGTAGC_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R2G48_GGTAGC_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R2G63_CGATGT_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R3G01_TGACCA_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R3G10_CAGATC_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R3G20_CTTGTA_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R3G34_ATGAGC_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R3G40_ATGAGC_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R3G48_ATGAGC_1_trimmed.fq.gz.allte.sort.bam
    ## mel_R3G63_TTAGGC_1_trimmed.fq.gz.allte.sort.bam
    ## dmel_naive_CATTTT_1_trimmed.fq.gz.rawabu
    ## dmel_naive_CATTTT_1_trimmed.fq.gz.rpm
    ## mel_R1G01_CGATGT_1_trimmed.fq.gz.rawabu
    ## mel_R1G01_CGATGT_1_trimmed.fq.gz.rpm
    ## mel_R1G10_ACAGTG_1_trimmed.fq.gz.rawabu
    ## mel_R1G10_ACAGTG_1_trimmed.fq.gz.rpm
    ## mel_R1G20_ACTTGA_1_trimmed.fq.gz.rawabu
    ## mel_R1G20_ACTTGA_1_trimmed.fq.gz.rpm
    ## mel_R1G34_GTAGAG_1_trimmed.fq.gz.rawabu
    ## mel_R1G34_GTAGAG_1_trimmed.fq.gz.rpm
    ## mel_R1G40_GTAGAG_1_trimmed.fq.gz.rawabu
    ## mel_R1G40_GTAGAG_1_trimmed.fq.gz.rpm
    ## mel_R1G48_GTAGAG_1_trimmed.fq.gz.rawabu
    ## mel_R1G48_GTAGAG_1_trimmed.fq.gz.rpm
    ## mel_R1G63_ATCACG_1_trimmed.fq.gz.rawabu
    ## mel_R1G63_ATCACG_1_trimmed.fq.gz.rpm
    ## mel_R2G01_TTAGGC_1_trimmed.fq.gz.rawabu
    ## mel_R2G01_TTAGGC_1_trimmed.fq.gz.rpm
    ## mel_R2G10_GCCAAT_1_trimmed.fq.gz.rawabu
    ## mel_R2G10_GCCAAT_1_trimmed.fq.gz.rpm
    ## mel_R2G20_GGCTAC_1_trimmed.fq.gz.rawabu
    ## mel_R2G20_GGCTAC_1_trimmed.fq.gz.rpm
    ## mel_R2G34_GGTAGC_1_trimmed.fq.gz.rawabu
    ## mel_R2G34_GGTAGC_1_trimmed.fq.gz.rpm
    ## mel_R2G40_GGTAGC_1_trimmed.fq.gz.rawabu
    ## mel_R2G40_GGTAGC_1_trimmed.fq.gz.rpm
    ## mel_R2G48_GGTAGC_1_trimmed.fq.gz.rawabu
    ## mel_R2G48_GGTAGC_1_trimmed.fq.gz.rpm
    ## mel_R2G63_CGATGT_1_trimmed.fq.gz.rawabu
    ## mel_R2G63_CGATGT_1_trimmed.fq.gz.rpm
    ## mel_R3G01_TGACCA_1_trimmed.fq.gz.rawabu
    ## mel_R3G01_TGACCA_1_trimmed.fq.gz.rpm
    ## mel_R3G10_CAGATC_1_trimmed.fq.gz.rawabu
    ## mel_R3G10_CAGATC_1_trimmed.fq.gz.rpm
    ## mel_R3G20_CTTGTA_1_trimmed.fq.gz.rawabu
    ## mel_R3G20_CTTGTA_1_trimmed.fq.gz.rpm
    ## mel_R3G34_ATGAGC_1_trimmed.fq.gz.rawabu
    ## mel_R3G34_ATGAGC_1_trimmed.fq.gz.rpm
    ## mel_R3G40_ATGAGC_1_trimmed.fq.gz.rawabu
    ## mel_R3G40_ATGAGC_1_trimmed.fq.gz.rpm
    ## mel_R3G48_ATGAGC_1_trimmed.fq.gz.rawabu
    ## mel_R3G48_ATGAGC_1_trimmed.fq.gz.rpm
    ## mel_R3G63_TTAGGC_1_trimmed.fq.gz.rawabu
    ## mel_R3G63_TTAGGC_1_trimmed.fq.gz.rpm

# P-element copy number estimates

The deviate-family.sh script was run on the resulting .bam and .rpm
files, to assess and visualise P-element presence, using the family ID
“PPI251”.

``` bash
nohup zsh deviate-family.sh mel PPI251 > /Volumes/Data/Projects/DmelR2_P-ele/logs/deviate_mel &
```

``` bash

if [ $# -lt 2 ]
  # "$#" is number of parameters- here we test
   then
   echo "Usage $0 sampleid family"
exit 2
fi

set -o shwordsplit

# Files and folders
sampleid=$1
tefamily=$2
output="../results/$sampleid/deviate/$tefamily"
inputbam="../results/$sampleid/bam"
# scgold="Dmel_rpl32,Dmel_piwi,Dmel_Act5C"
scg="Dmel_rhi,Dmel_rpl32,Dmel_tj"

## the REFERENCE 
refg="../refg/Dmel_tes_3scg/teseqs-3scg.fasta"
anno="../refg/Dmel_tes_3scg/teseqs.gff"
hier="../refg/Dmel_tes_3scg/teseqs.hier"


# software
samtools="../bin/samtools"
bwa="../bin/bwa"

mkdir -p $output

# paths
for bam in $inputbam/*.allte.sort.bam
do
    
    n=`basename $bam` 
    sampleid=${n%.allte.sort.bam}
    echo $sampleid


    com="deviate --input_bam $bam --library $refg --annotation $anno --single_copy_genes $scg --families $tefamily --minID 1"
    echo $com
    eval $com
    mv ${bam}.${tefamily} $output/$sampleid.$tefamily
    mv ${bam}.${tefamily}.pdf  $output/${sampleid}.${tefamily}.pdf
    mv ${bam}.${tefamily}.raw $output/${sampleid}.${tefamily}.raw
    rm ${bam}.fused.sort.bam
    rm ${bam}.fused.sort.bam.bai

    
done 
```

Generating the following .fq.gz.PPI251 files.

``` bash
cd /Volumes/Data/Projects/DmelR2_P-ele/software/teminer-code-r10/results/mel/deviate/PPI251
ls *fq.gz.*
```

    ## dmel_naive_CATTTT_1_trimmed.fq.gz.PPI251
    ## dmel_naive_CATTTT_1_trimmed.fq.gz.PPI251.pdf
    ## dmel_naive_CATTTT_1_trimmed.fq.gz.PPI251.raw
    ## mel_R1G01_CGATGT_1_trimmed.fq.gz.PPI251
    ## mel_R1G01_CGATGT_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R1G01_CGATGT_1_trimmed.fq.gz.PPI251.raw
    ## mel_R1G10_ACAGTG_1_trimmed.fq.gz.PPI251
    ## mel_R1G10_ACAGTG_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R1G10_ACAGTG_1_trimmed.fq.gz.PPI251.raw
    ## mel_R1G20_ACTTGA_1_trimmed.fq.gz.PPI251
    ## mel_R1G20_ACTTGA_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R1G20_ACTTGA_1_trimmed.fq.gz.PPI251.raw
    ## mel_R1G34_GTAGAG_1_trimmed.fq.gz.PPI251
    ## mel_R1G34_GTAGAG_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R1G34_GTAGAG_1_trimmed.fq.gz.PPI251.raw
    ## mel_R1G40_GTAGAG_1_trimmed.fq.gz.PPI251
    ## mel_R1G40_GTAGAG_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R1G40_GTAGAG_1_trimmed.fq.gz.PPI251.raw
    ## mel_R1G48_GTAGAG_1_trimmed.fq.gz.PPI251
    ## mel_R1G48_GTAGAG_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R1G48_GTAGAG_1_trimmed.fq.gz.PPI251.raw
    ## mel_R1G63_ATCACG_1_trimmed.fq.gz.PPI251
    ## mel_R1G63_ATCACG_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R1G63_ATCACG_1_trimmed.fq.gz.PPI251.raw
    ## mel_R2G01_TTAGGC_1_trimmed.fq.gz.PPI251
    ## mel_R2G01_TTAGGC_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R2G01_TTAGGC_1_trimmed.fq.gz.PPI251.raw
    ## mel_R2G10_GCCAAT_1_trimmed.fq.gz.PPI251
    ## mel_R2G10_GCCAAT_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R2G10_GCCAAT_1_trimmed.fq.gz.PPI251.raw
    ## mel_R2G20_GGCTAC_1_trimmed.fq.gz.PPI251
    ## mel_R2G20_GGCTAC_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R2G20_GGCTAC_1_trimmed.fq.gz.PPI251.raw
    ## mel_R2G34_GGTAGC_1_trimmed.fq.gz.PPI251
    ## mel_R2G34_GGTAGC_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R2G34_GGTAGC_1_trimmed.fq.gz.PPI251.raw
    ## mel_R2G40_GGTAGC_1_trimmed.fq.gz.PPI251
    ## mel_R2G40_GGTAGC_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R2G40_GGTAGC_1_trimmed.fq.gz.PPI251.raw
    ## mel_R2G48_GGTAGC_1_trimmed.fq.gz.PPI251
    ## mel_R2G48_GGTAGC_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R2G48_GGTAGC_1_trimmed.fq.gz.PPI251.raw
    ## mel_R2G63_CGATGT_1_trimmed.fq.gz.PPI251
    ## mel_R2G63_CGATGT_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R2G63_CGATGT_1_trimmed.fq.gz.PPI251.raw
    ## mel_R3G01_TGACCA_1_trimmed.fq.gz.PPI251
    ## mel_R3G01_TGACCA_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R3G01_TGACCA_1_trimmed.fq.gz.PPI251.raw
    ## mel_R3G10_CAGATC_1_trimmed.fq.gz.PPI251
    ## mel_R3G10_CAGATC_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R3G10_CAGATC_1_trimmed.fq.gz.PPI251.raw
    ## mel_R3G20_CTTGTA_1_trimmed.fq.gz.PPI251
    ## mel_R3G20_CTTGTA_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R3G20_CTTGTA_1_trimmed.fq.gz.PPI251.raw
    ## mel_R3G34_ATGAGC_1_trimmed.fq.gz.PPI251
    ## mel_R3G34_ATGAGC_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R3G34_ATGAGC_1_trimmed.fq.gz.PPI251.raw
    ## mel_R3G40_ATGAGC_1_trimmed.fq.gz.PPI251
    ## mel_R3G40_ATGAGC_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R3G40_ATGAGC_1_trimmed.fq.gz.PPI251.raw
    ## mel_R3G48_ATGAGC_1_trimmed.fq.gz.PPI251
    ## mel_R3G48_ATGAGC_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R3G48_ATGAGC_1_trimmed.fq.gz.PPI251.raw
    ## mel_R3G63_TTAGGC_1_trimmed.fq.gz.PPI251
    ## mel_R3G63_TTAGGC_1_trimmed.fq.gz.PPI251.pdf
    ## mel_R3G63_TTAGGC_1_trimmed.fq.gz.PPI251.raw

From the resulting .fq.gz.PPI251 files, we extracted and assembled the
P-element copy number values into a single file, then removed the “G”
from the generation coloumn to aid visualisation.

``` bash

for i in *.PPI251; do echo $i | cut -f1 -d "." | cut -f2 -d "_"  >> mel_pcopies.txt; grep "or " $i | cut -f5 -d " "  | cut -b 1-8 >> mel_pcopies.txt; done

 mel_pcopies.txt |grep  '^R' -A2 | xargs -n3 > mel_pcopies_aligned.txt

done
```

``` bash
cd /Volumes/Data/Projects/DmelR2_P-ele/software/teminer-code-r10/results/mel/deviate/PPI251
ls *pcopies*
```

    ## mel_pcopies.txt
    ## mel_pcopies_aligned.tps
    ## mel_pcopies_aligned.txt

# Visualisation

## P-element copy number

The copy number data was then visualised using ggplot2.

``` r
library (ggplot2)
theme_set(theme_bw())
tresrep<-c("firebrick", "skyblue3", "chartreuse4")
  
d = read.table("mel_pcopies_values.txt")

names(d)<-c ("replicate", "generation", "copy number")

print(d)
```

    ##    replicate generation copy number
    ## 1         R1          1    0.229964
    ## 2         R1         10    1.007938
    ## 3         R1         20    4.177100
    ## 4         R1         34   24.724080
    ## 5         R1         40   26.013780
    ## 6         R1         48   25.789190
    ## 7         R1         63   26.365490
    ## 8         R2          1    0.446720
    ## 9         R2         10    0.999447
    ## 10        R2         20    2.383754
    ## 11        R2         34    6.281913
    ## 12        R2         40    6.159023
    ## 13        R2         48    7.119456
    ## 14        R2         63    6.490166
    ## 15        R3          1    0.581259
    ## 16        R3         10    1.307322
    ## 17        R3         20    3.947311
    ## 18        R3         34   15.565970
    ## 19        R3         40   18.922240
    ## 20        R3         48   20.435890
    ## 21        R3         63   20.984870

``` r
pele_dna <- ggplot(d, aes(x = generation, y = `copy number`, group= replicate)) + 
  geom_line(aes(color=replicate),linewidth=1) +
  geom_point(aes(color=replicate)) +
  ggtitle("P-element invasions in D. melanogaster") +
  scale_colour_manual(values=tresrep) +
  xlim(0, 65)

ggsave("figs/pele_dna.png", pele_dna, width = 8, height = 5, dpi = 600)

knitr::include_graphics("figs/pele_dna.png")
```

<img src="figs/pele_dna.png" width="4800" />

``` r
pele_dna_2 <- ggplot(d, aes(x = generation, y = `copy number`, group = replicate)) + 
  geom_point(aes(colour = replicate), size = 3.5, shape = 18) +
  geom_path(aes(colour = replicate), linewidth = 1.5) +
  scale_colour_manual(values = tresrep) +
  scale_size(guide = "none") +  
  xlim(0, 65) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        plot.margin = margin(10, 10, 10, 10, "pt"),
        legend.position = c(0.1, 0.8),
        legend.box.background = element_rect(size = 0.5))
```

    ## Warning: The `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
pele_dna_sl <- pele_dna_2 + guides(color = guide_legend(title = NULL))

ggsave("figs/pele_dna_sl.png", pele_dna_sl, width = 6, height = 4, dpi = 600)

knitr::include_graphics("figs/pele_dna_sl.png")
```

<img src="figs/pele_dna_sl.png" width="3600" />

## DeviaTE

We also generated DeviaTE plots for each replicate for each measured
generation time-point.

``` r
knitr::include_graphics("figs/pele_coverage.png")
```

<img src="figs/pele_coverage.png" width="13464" />

``` r
knitr::include_graphics("figs/pele_coverage_s.png")
```

<img src="figs/pele_coverage_s.png" width="13493" />

``` r
knitr::include_graphics("figs/pele_coverage_xs.png")
```

<img src="figs/pele_coverage_xs.png" width="13493" />
