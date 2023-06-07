Data overview of Dmel P-element invasion
================
Divya Selvaraju
6/7/2023

``` bash
knitr::opts_chunk$set(echo = TRUE)
```

We are counting the number of reads as well as read length using the
following commands. The number of reads for each sample should be same
in both R1 and R2 of a given paired-end sample.

## Genomic data overview - Processing

``` bash
cd /Volume/Temp/Divya/invasion/fastq/mel/

for i in *_1.fq.gz;do echo $i >> dmel_dna_len_1.txt; zless $i | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c >> dmel_dna_len_1.txt; done &

for i in *_2.fq.gz;do echo $i >> dmel_dna_len_2.txt; zless $i | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c >> dmel_dna_len_2.txt; done &

less dmel_dna_len_1.txt| xargs -n3 | sed 's/_/ /g' | sed 's/1\.fq\.gz//g' | sed 's/[A-Z][A-Z][A-Z][A-Z][A-Z][A-Z]//g' | sed 's/G/ G/'| sed 's/$/ paired/'|  sed 's/  / /g' | sort -k2 > dmel_dna_rl.txt
```

## Output of overview of the genomic data

``` r
dna_data= read.table("dmel_dna_rl.txt",header=TRUE)
print(dna_data)
```

    ##    species replicate generation totalreads readlength readtype
    ## 1      mel        R1        G01   26714134        125   paired
    ## 2      mel        R1        G10   27630152        125   paired
    ## 3      mel        R1        G20   30589764        125   paired
    ## 4      mel        R1        G34   22727883        125   paired
    ## 5      mel        R1        G40   27593549        125   paired
    ## 6      mel        R1        G48   26804510        125   paired
    ## 7      mel        R1        G63   35552063        100   paired
    ## 8      mel        R2        G01   29362505        125   paired
    ## 9      mel        R2        G10   27137385        125   paired
    ## 10     mel        R2        G20   28290767        125   paired
    ## 11     mel        R2        G34   20023157        125   paired
    ## 12     mel        R2        G40   26821470        125   paired
    ## 13     mel        R2        G48   26282343        125   paired
    ## 14     mel        R2        G63   44583633        100   paired
    ## 15     mel        R3        G01   28528059        125   paired
    ## 16     mel        R3        G10   30253806        125   paired
    ## 17     mel        R3        G20   23926839        125   paired
    ## 18     mel        R3        G34   24064611        125   paired
    ## 19     mel        R3        G40   38535020        125   paired
    ## 20     mel        R3        G48   24943600        125   paired
    ## 21     mel        R3        G63   38737252        100   paired
    ## 22     mel     naive          -   14935203        125   paired

## Transcriptomic data overview - Processing

# 

``` bash
cd /Volume/Temp/Divya/invasion/RNA/data/dmel/run1

for i in *_R1.fastq.gz;do echo $i >> dmel_rna_len_run1_1.txt; zless $i | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c >> dmel_rna_len_run1_1.txt; done &

for i in *_R2.fastq.gz;do echo $i >> dmel_rna_len_run1_2.txt; zless $i | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c >>
dmel_rna_len_run1_2.txt; done &

dmel_rna_len_run1_1.txt| xargs -n3 | sed 's/_/ /g' | sed 's/R1\.fastq\.gz//g' |sed 's/$/ paired/' | sort -k3 > dmel_rna_run1_rl.txt

cd ../run2

for i in *_R1.fastq.gz;do echo $i >> dmel_rna_len_run2_1.txt; zless $i | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c >> dmel_rna_len_run2_1.txt; done &

for i in *_R2.fastq.gz;do echo $i >> dmel_rna_len_run2_2.txt; zless $i | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c >> dmel_rna_len_run2_2.txt; done &

less dmel_rna_len_run2_1.txt| xargs -n3 | sed 's/_/ /g' | sed 's/R1\.fastq\.gz//g' |sed 's/$/ paired/'  > dmel_rna_run2_rl.txt
```

## Output of overview of transcriptomic data

``` r
rna_data_run1 = read.table("dmel_rna_run1_rl.txt", header = TRUE)
print(rna_data_run1)
```

    ##    species material replicate generation  run totalreads readlength readtype
    ## 1     dmel      rna        R1        G15 run1   17919773        150   paired
    ## 2     dmel      rna        R1        G21 run1   17415757        150   paired
    ## 3     dmel      rna        R1        G30 run1   16116118        150   paired
    ## 4     dmel      rna        R1        G40 run1   17658726        150   paired
    ## 5     dmel      rna        R1         G6 run1   15182477        150   paired
    ## 6     dmel      rna        R2        G15 run1   18195062        150   paired
    ## 7     dmel      rna        R2        G21 run1   15773656        150   paired
    ## 8     dmel      rna        R2        G30 run1   15749983        150   paired
    ## 9     dmel      rna        R2        G40 run1   18425395        150   paired
    ## 10    dmel      rna        R2         G6 run1   16052901        150   paired
    ## 11    dmel      rna        R3        G15 run1   16040102        150   paired
    ## 12    dmel      rna        R3        G21 run1   17560544        150   paired
    ## 13    dmel      rna        R3        G30 run1   18488825        150   paired
    ## 14    dmel      rna        R3        G40 run1   18837267        150   paired
    ## 15    dmel      rna        R3         G6 run1   16297920        150   paired

``` r
rna_data_run2 = read.table("dmel_rna_run1_rl.txt", header = TRUE)
print(rna_data_run2)
```

    ##    species material replicate generation  run totalreads readlength readtype
    ## 1     dmel      rna        R1        G15 run1   17919773        150   paired
    ## 2     dmel      rna        R1        G21 run1   17415757        150   paired
    ## 3     dmel      rna        R1        G30 run1   16116118        150   paired
    ## 4     dmel      rna        R1        G40 run1   17658726        150   paired
    ## 5     dmel      rna        R1         G6 run1   15182477        150   paired
    ## 6     dmel      rna        R2        G15 run1   18195062        150   paired
    ## 7     dmel      rna        R2        G21 run1   15773656        150   paired
    ## 8     dmel      rna        R2        G30 run1   15749983        150   paired
    ## 9     dmel      rna        R2        G40 run1   18425395        150   paired
    ## 10    dmel      rna        R2         G6 run1   16052901        150   paired
    ## 11    dmel      rna        R3        G15 run1   16040102        150   paired
    ## 12    dmel      rna        R3        G21 run1   17560544        150   paired
    ## 13    dmel      rna        R3        G30 run1   18488825        150   paired
    ## 14    dmel      rna        R3        G40 run1   18837267        150   paired
    ## 15    dmel      rna        R3         G6 run1   16297920        150   paired

## small RNA data overview - Processing

``` bash
cd /Volume/Temp/Divya/invasion/piRNA/data/dmel/batch1

for i in *.fastq.gz;do echo $i >> dmel_smallrna_batch1.txt; zless $i | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c >> dmel_smallrna_batch1.txt; done &

less dmel_smallrna_batch1.txt| xargs -n3 | sed 's/_/ /g' | sed 's/\.fastq\.gz//g' |sed 's/$/ ovaries/' > dmel_smallrna_batch1_rl.txt

cd ../batch2

for i in *.fastq.gz;do echo $i >> dmel_smallrna_len_batch2.txt; zless $i | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c >> dmel_smallrna_len_batch2.txt; done &

less dmel_smallrna_len_batch2.txt| xargs -n3 | sed 's/_/ /g' | sed 's/\.fastq\.gz//g' |gsed 's/$/ wholebody/' > dmel_smallrna_batch2_rl.txt

cat ../batch1/dmel_smallrna_batch1_rl.txt dmel_smallrna_batch2_rl.txt| sort -k2 > dmel_smallrna_rl.txt
```

## Output of overview of the small RNA data

``` r
smallrna_data= read.table("dmel_smallrna_rl.txt",header=TRUE)
print(smallrna_data)
```

    ##    species replicate generation  run totalreads readlength    tissue
    ## 1     dmel        R1         G1 run1   12200215         50   ovaries
    ## 2     dmel        R1        G15 run1   14332444         75 wholebody
    ## 3     dmel        R1        G15 run2   12444828         75 wholebody
    ## 4     dmel        R1        G21 run1   13497060         75 wholebody
    ## 5     dmel        R1        G21 run2   11643714         75 wholebody
    ## 6     dmel        R1        G25 run1   14663976         75 wholebody
    ## 7     dmel        R1        G25 run2   12574018         75 wholebody
    ## 8     dmel        R1        G30 run1   14825830         75 wholebody
    ## 9     dmel        R1        G30 run2   12690099         75 wholebody
    ## 10    dmel        R1        G35 run1   12605340         50   ovaries
    ## 11    dmel        R1        G40 run1   12920000         75 wholebody
    ## 12    dmel        R1        G40 run2   11173097         75 wholebody
    ## 13    dmel        R1        G45 run1   11801273         75 wholebody
    ## 14    dmel        R1        G45 run2   10070017         75 wholebody
    ## 15    dmel        R1         G6 run1   15634888         75 wholebody
    ## 16    dmel        R1         G6 run2   13546426         75 wholebody
    ## 17    dmel        R2         G1 run1   12485241         50   ovaries
    ## 18    dmel        R2        G15 run1   12816885         75 wholebody
    ## 19    dmel        R2        G15 run2   11073322         75 wholebody
    ## 20    dmel        R2        G21 run1   15224282         75 wholebody
    ## 21    dmel        R2        G21 run2   13042586         75 wholebody
    ## 22    dmel        R2        G25 run1   13313137         75 wholebody
    ## 23    dmel        R2        G25 run2   11577780         75 wholebody
    ## 24    dmel        R2        G30 run1   14634092         75 wholebody
    ## 25    dmel        R2        G30 run2   12589509         75 wholebody
    ## 26    dmel        R2        G35 run1   11675753         50   ovaries
    ## 27    dmel        R2        G40 run1   13267775         75 wholebody
    ## 28    dmel        R2        G40 run2   11397755         75 wholebody
    ## 29    dmel        R2        G45 run1   12525298         75 wholebody
    ## 30    dmel        R2        G45 run2   10756270         75 wholebody
    ## 31    dmel        R2         G6 run1   13786454         75 wholebody
    ## 32    dmel        R2         G6 run2   11904680         75 wholebody
    ## 33    dmel        R3         G1 run1   13610358         50   ovaries
    ## 34    dmel        R3        G15 run1   14654721         75 wholebody
    ## 35    dmel        R3        G15 run2   12567287         75 wholebody
    ## 36    dmel        R3        G21 run1   15301439         75 wholebody
    ## 37    dmel        R3        G21 run2   13158327         75 wholebody
    ## 38    dmel        R3        G25 run1   12761202         75 wholebody
    ## 39    dmel        R3        G25 run2   10897447         75 wholebody
    ## 40    dmel        R3        G30 run1   15222316         75 wholebody
    ## 41    dmel        R3        G30 run2   13036716         75 wholebody
    ## 42    dmel        R3        G35 run1    9274309         50   ovaries
    ## 43    dmel        R3        G40 run1   12811052         75 wholebody
    ## 44    dmel        R3        G40 run2   11064679         75 wholebody
    ## 45    dmel        R3        G45 run1   12921264         75 wholebody
    ## 46    dmel        R3        G45 run2   11101567         75 wholebody
    ## 47    dmel        R3         G6 run1   11193408         75 wholebody
    ## 48    dmel        R3         G6 run2    9993762         75 wholebody
