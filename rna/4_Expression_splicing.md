4 - Expression & splicing
================
Matthew Beaumont
2023-07-17

``` bash
knitr::opts_chunk$set(echo = TRUE)
```

# Visualisation

We first need to extract the expression data from our GSNAP aligned
reads.

``` bash
#!/bin/bash

fai="/Volumes/Data/Tools/RefGenomes/dmel/rna/dmel_TEs/dmel-transcriptome-r6.52-TEs.fasta.fai"
pyscript="/Volumes/Data/Projects/dmelR2_p-ele/scripts/mRNA-coverage-senseantisense.py"
if="/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/map-GSNAP/output"
of="/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/splicing-expression/expression"
seqs="PPI251,FBtr0083183_mRNA,FBtr0088034_mRNA,FBtr0086904_mRNA,FBtr0087984_mRNA,FBtr0087189_mRNA,FBtr0080497_mRNA,FBtr0079489_mRNA,FBtr0445185_mRNA,FBtr0080316_mRNA,FBtr0075559_mRNA,FBtr0100641_mRNA,FBtr0080165_mRNA,FBtr0081502_mRNA,FBtr0073637_mRNA,FBtr0080166_mRNA,FBtr0301669_mRNA,FBtr0086897_mRNA,FBtr0085594_mRNA,FBtr0329922_mRNA,FBtr0081328_mRNA"

# All transcripts are 'RA' variants unless stated otherwise.

samtools view $if/gt_R1G6.sort.bam | python $pyscript --sam - --sample-id R1-G6 --seqs $seqs --fai $fai -a  > $of/R1G6.txt
samtools view $if/gt_R2G6.sort.bam | python $pyscript --sam - --sample-id R2-G6 --seqs $seqs --fai $fai -a  > $of/R2G6.txt
samtools view $if/gt_R3G6.sort.bam | python $pyscript --sam - --sample-id R3-G6 --seqs $seqs --fai $fai -a  > $of/R3G6.txt

samtools view $if/gt_R1G15.sort.bam | python $pyscript --sam - --sample-id R1-G15 --seqs $seqs --fai $fai -a  > $of/R1G15.txt
samtools view $if/gt_R2G15.sort.bam | python $pyscript --sam - --sample-id R2-G15 --seqs $seqs --fai $fai -a  > $of/R2G15.txt
samtools view $if/gt_R3G15.sort.bam | python $pyscript --sam - --sample-id R3-G15 --seqs $seqs --fai $fai -a  > $of/R3G15.txt

samtools view $if/gt_R1G21.sort.bam | python $pyscript --sam - --sample-id R1-G21 --seqs $seqs --fai $fai -a  > $of/R1G21.txt
samtools view $if/gt_R2G21.sort.bam | python $pyscript --sam - --sample-id R2-G21 --seqs $seqs --fai $fai -a  > $of/R2G21.txt
samtools view $if/gt_R3G21.sort.bam | python $pyscript --sam - --sample-id R3-G21 --seqs $seqs --fai $fai -a  > $of/R3G21.txt

samtools view $if/gt_R1G30.sort.bam | python $pyscript --sam - --sample-id R1-G30 --seqs $seqs --fai $fai -a  > $of/R1G30.txt
samtools view $if/gt_R2G30.sort.bam | python $pyscript --sam - --sample-id R2-G30 --seqs $seqs --fai $fai -a  > $of/R2G30.txt
samtools view $if/gt_R3G30.sort.bam | python $pyscript --sam - --sample-id R3-G30 --seqs $seqs --fai $fai -a  > $of/R3G30.txt

samtools view $if/gt_R1G40.sort.bam | python $pyscript --sam - --sample-id R1-G40 --seqs $seqs --fai $fai -a  > $of/R1G40.txt
samtools view $if/gt_R2G40.sort.bam | python $pyscript --sam - --sample-id R2-G40 --seqs $seqs --fai $fai -a  > $of/R2G40.txt
samtools view $if/gt_R3G40.sort.bam | python $pyscript --sam - --sample-id R3-G40 --seqs $seqs --fai $fai -a  > $of/R3G40.txt

# Combine outputs into single file, separating ID.
#cat *.txt| perl -pe 's/-/\t/'|perl -pe 's/-/\t/' > expr.forr 

# Tissue column for 'wf' needed adding back.
#awk 'BEGIN{OFS="\t"}{$2 = $2 "\t" "wf"; print}' expr.forr > expr_wf.forr
```

Then again for splicing.

``` bash
#!/bin/bash

pyscript="/Volumes/Data/Projects/dmelR2_p-ele/scripts/mRNA-splicing-senseantisense.py"
if="/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/map-GSNAP/output"
of="/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/splicing-expression/splicing"
seqs="PPI251,FBtr0083183_mRNA,FBtr0088034_mRNA,FBtr0086904_mRNA,FBtr0087984_mRNA,FBtr0087189_mRNA,FBtr0080497_mRNA,FBtr0079489_mRNA,FBtr0445185_mRNA,FBtr0080316_mRNA,FBtr0075559_mRNA,FBtr0100641_mRNA,FBtr0080165_mRNA,FBtr0081502_mRNA,FBtr0073637_mRNA,FBtr0080166_mRNA,FBtr0301669_mRNA,FBtr0086897_mRNA,FBtr0085594_mRNA,FBtr0329922_mRNA,FBtr0081328_mRNA"

# All transcripts are 'RA' variants unless stated otherwise.

samtools view $if/gt_R1G6.sort.bam | python $pyscript --sam - --sample-id R1-G6-wf --seqs $seqs  > $of/R1G6.txt
samtools view $if/gt_R2G6.sort.bam | python $pyscript --sam - --sample-id R2-G6-wf --seqs $seqs  > $of/R2G6.txt
samtools view $if/gt_R3G6.sort.bam | python $pyscript --sam - --sample-id R3-G6-wf --seqs $seqs  > $of/R3G6.txt

samtools view $if/gt_R1G15.sort.bam | python $pyscript --sam - --sample-id R1-G15-wf --seqs $seqs  > $of/R1G15.txt
samtools view $if/gt_R2G15.sort.bam | python $pyscript --sam - --sample-id R2-G15-wf --seqs $seqs  > $of/R2G15.txt
samtools view $if/gt_R3G15.sort.bam | python $pyscript --sam - --sample-id R3-G15-wf --seqs $seqs  > $of/R3G15.txt

samtools view $if/gt_R1G21.sort.bam | python $pyscript --sam - --sample-id R1-G21-wf --seqs $seqs  > $of/R1G21.txt
samtools view $if/gt_R2G21.sort.bam | python $pyscript --sam - --sample-id R2-G21-wf --seqs $seqs  > $of/R2G21.txt
samtools view $if/gt_R3G21.sort.bam | python $pyscript --sam - --sample-id R3-G21-wf --seqs $seqs  > $of/R3G21.txt

samtools view $if/gt_R1G30.sort.bam | python $pyscript --sam - --sample-id R1-G30-wf --seqs $seqs  > $of/R1G30.txt
samtools view $if/gt_R2G30.sort.bam | python $pyscript --sam - --sample-id R2-G30-wf --seqs $seqs  > $of/R2G30.txt
samtools view $if/gt_R3G30.sort.bam | python $pyscript --sam - --sample-id R3-G30-wf --seqs $seqs  > $of/R3G30.txt

samtools view $if/gt_R1G40.sort.bam | python $pyscript --sam - --sample-id R1-G40-wf --seqs $seqs  > $of/R1G40.txt
samtools view $if/gt_R2G40.sort.bam | python $pyscript --sam - --sample-id R2-G40-wf --seqs $seqs  > $of/R2G40.txt
samtools view $if/gt_R3G40.sort.bam | python $pyscript --sam - --sample-id R3-G40-wf --seqs $seqs  > $of/R3G40.txt

# Combine outputs into single file, separating ID.
# cat *.txt| perl -pe 's/-/\t/'|perl -pe 's/-/\t/' > spli.forr
```

And then we visualised both the P-element expression and the splicing of
its’ three introns together in ggplot2.

``` r
library(ggplot2)
theme_set(theme_bw())

target<-"PPI251"
ttissue<-"wf" # target tissue
sminfreq<-0.1

h<-read.table("/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/splicing-expression/expression/expr_wf.forr")
names(h)<-c("rep","time","tissue","strand","gene","pos","cov")

a<-subset(h,gene==target & tissue==ttissue)
a$key<-paste0(a$rep,"_",a$time,"_",a$pos,"_",a$strand)


a$time <- factor(a$time, levels=c("G6", "G15", "G21", "G30", "G40"))
s<-subset(a,strand=="se")
as<-subset(a,strand=="ase")

spli<-read.table("/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/splicing-expression/splicing/spli.forr")

names(spli)<-c("rep","time","tissue","strand","gene","skey","start","end","rawcount","freq")
aspli<-subset(spli,gene==target & freq>sminfreq & tissue==ttissue)
aspli$time <- factor(aspli$time, levels=c("G6", "G15", "G21", "G30", "G40"))
aspli$rep<- as.factor(aspli$rep)
aspli$start<-aspli$start-1 # position inaccuracy, graph is more appealing
aspli$keystart<-paste0(aspli$rep,"_",aspli$time,"_",aspli$start,"_",aspli$strand)
aspli$keyend<-paste0(aspli$rep,"_",aspli$time,"_",aspli$end,"_",aspli$strand)

aspli<-merge(x=aspli,y=a[,c("key","cov")],by.x="keystart",by.y="key")
aspli<-merge(x=aspli,y=a[,c("key","cov")],by.x="keyend",by.y="key")
aspli$size<-log(aspli$freq+1)

a_s<-subset(aspli,strand=="se")
a_as<-subset(aspli,strand=="ase")

expr_spli_plot <- ggplot() +
    geom_polygon(data=s,mapping=aes(x=pos, y=cov), fill='grey', color='grey') +
    geom_polygon(data=as, aes(x=pos, y=-cov), fill='lightgrey', color='lightgrey')+
    geom_curve(data=a_s, mapping=aes(x=start, y=cov.x, xend=end, yend=cov.y,linewidth=size), curvature=-0.2, ncp=10, show.legend=FALSE)+
    geom_curve(data=a_as, mapping=aes(x=start, y=-cov.x, xend=end, yend=-cov.y,linewidth=size), curvature=0.2, ncp=10, show.legend=FALSE)+
    facet_grid(time~rep)+scale_size(range=c(0.2,2))+xlab("position")+ylab("expression level [rpm]")

ggsave("figs/expr_spli.png", expr_spli_plot, width = 14, height = 14)

knitr::include_graphics("figs/expr_spli.png")
```

<img src="figs/expr_spli.png" width="4200" />

# P-element expression

We ran the following script to extract out expression levels.

``` bash
#!/bin/bash

fai="/Volumes/Data/Tools/RefGenomes/dmel/rna/dmel_TEs/dmel-transcriptome-r6.52-TEs.fasta.fai"
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

# Combine outputs into single file, separating ID.
# cat *.txt| perl -pe 's/-/\t/'|perl -pe 's/-/\t/' > expr.forr
```

Then we visualise it using ggplot.

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.2     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ lubridate 1.9.2     ✔ tibble    3.2.1
    ## ✔ purrr     1.0.1     ✔ tidyr     1.3.0
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(RColorBrewer)
theme_set(theme_bw())
tresrep <- c("firebrick", "skyblue3", "chartreuse4")

t <- read_delim("/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/splicing-expression/raw_expression/all-expressionlevel/expr.forr", delim = "\t", col_names = FALSE, comment = "#", show_col_types = FALSE)
names(t) <- c("replicate", "generation", "run", "gene", "rawse", "rawase", "genlen", "sense", "antisense", "total")
t <- subset(t, gene == "PPI251")
t$generation <- as.numeric(substring(t$generation, 2))

# Add data point for gen 0 with expression level 0
zero_data <- data.frame(replicate = rep(c("R1", "R2", "R3"), each = 1),
                        generation = rep(0, 3),
                        sense = rep(0, 3),
                        stringsAsFactors = FALSE)
t <- bind_rows(t, zero_data)

width <- 16
height <- 12
resolution <- 600

s <- ggplot(data = t, aes(linetype = replicate)) +
  geom_line(data = t, aes(x = generation, y = sense, color = replicate), linewidth = 1) +
  geom_point(data = t, aes(x = generation, y = sense, color = replicate)) +
  theme(strip.text = element_blank(), 
        legend.position = c(0.1, 0.85),
        legend.box.background = element_rect(color = "black", fill = "transparent"),
        legend.title = element_blank()) + 
  ylab("expression [rpkm]") +
  scale_colour_manual(values = tresrep) +
  xlim(0, 48) +
  xlab("generation")

ggsave("figs/P-ele_expression.png", plot = s, width = 10, height = 6, dpi = 600)

knitr::include_graphics("figs/P-ele_expression.png")
```

<img src="figs/P-ele_expression.png" width="6000" />

## IVS3

``` r
library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())
colrep1 <- c("firebrick", "skyblue3", "chartreuse4")

# Annotation
# PPI251    ensembl exon    153 442 .   +   .   gene_id "pele"; transcript_id "pele1";
# PPI251    ensembl exon    501 1168    .   +   .   gene_id "pele"; transcript_id "pele1";
# PPI251    ensembl exon    1222    1947    .   +   .   gene_id "pele"; transcript_id "pele1";
# PPI251    ensembl exon    2138    2709    .   +   .   gene_id "pele"; transcript_id "pele1";

t <- read.table("/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/splicing-expression/splicing/spli.forr")

# head
# R1    G15 wf  se  PPI251  443-502 443 502 51  0.719474282680966
# R1    G15 wf  se  PPI251  1169-1223   1169    1223    51  0.719474282680966

names <- c("rep", "generation", "tissue", "strand", "gene", "skey", "start", "end", "rawcount", "freq")
names(t) <- names
t <- subset(t, gene == "PPI251" & tissue == "wf" & strand == "se")
t$generation <- as.numeric(substring(t$generation, 2))

t <- subset(t, skey == "1948-2139")  # delta23
missingvalues <- data.frame(rep = c("R1", "R2", "R3"), generation = c(0, 0, 0), tissue = rep("wf", 3),
                            strand = rep("se", 3), gene = rep("PPI251", 3), skey = rep("1948-2139", 3),
                            start = rep(0, 3), end = rep(0, 3), rawcount = rep(0, 3), freq = rep(0, 3))

t <- rbind(t, missingvalues)

# Add data point at gen 6 for R1 with splicing level 0
zero_data <- data.frame(rep = "R1", generation = 6, tissue = "wf", strand = "se", gene = "PPI251",
                        skey = "1948-2139", start = 0, end = 0, rawcount = 0, freq = 0)
t <- rbind(t, zero_data)

s <- ggplot(data = t, aes(linetype = rep)) +
  geom_line(data = t, aes(x = generation, y = freq, color = rep), linewidth = 1) +
  geom_point(data = t, aes(x = generation, y = freq, color = rep)) +
  theme(strip.text = element_blank(),
        legend.position = c(0.1, 0.85),
        legend.box.background = element_rect(color = "black", fill = "transparent"),
        legend.title = element_blank()) +
  ylab("splicing level IVS3 [srpm]") + 
  scale_colour_manual(values = colrep1) +
  xlim(0, 48) + xlab("generation")

ggsave("figs/IVS3.png", plot = s, width = 10, height = 6, dpi = 600)

knitr::include_graphics("figs/IVS3.png")
```

<img src="figs/IVS3.png" width="6000" />

## IVS2

``` r
library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())
colrep2<-c("firebrick", "skyblue3", "chartreuse4")

# Annotation
# PPI251    ensembl exon    153 442 .   +   .   gene_id "pele"; transcript_id "pele1";
# PPI251    ensembl exon    501 1168    .   +   .   gene_id "pele"; transcript_id "pele1";
# PPI251    ensembl exon    1222    1947    .   +   .   gene_id "pele"; transcript_id "pele1";
# PPI251    ensembl exon    2138    2709    .   +   .   gene_id "pele"; transcript_id "pele1";

t<-read.table("/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/splicing-expression/splicing/spli.forr")

# head
# R1    G15 wf  se  PPI251  443-502 443 502 51  0.719474282680966
# R1    G15 wf  se  PPI251  1169-1223   1169    1223    51  0.719474282680966

names<-c("rep","generation","tissue","strand","gene","skey","start","end","rawcount","freq")
names(t)<-names
t<-subset(t,gene=="PPI251" & tissue=="wf" & strand=="se")
t$generation<-as.numeric(substring(t$generation, 2))

t<-subset(t,skey=="1169-1223")
missingvalues=data.frame(rep=c("R1","R2","R3"),generation=c(0,0,0),tissue=rep("wf",3),
                         strand=rep("se",3),gene=rep("PPI251",3),skey=rep("1169-1223",3),
                         start=rep(0,3),end=rep(0,3),rawcount=rep(0,3),freq=rep(0,3))

t<-rbind(t,missingvalues)

s <- ggplot(data = t, aes(linetype = rep)) +
  geom_line(data = t, aes(x = generation, y = freq, color = rep), linewidth = 1) +
  geom_point(data = t, aes(x = generation, y = freq, color = rep)) +
  theme(strip.text = element_blank(),
        legend.position = c(0.1, 0.85),
        legend.box.background = element_rect(color = "black", fill = "transparent"),
        legend.title = element_blank()) +
  ylab("splicing level IVS2 [srpm]") +
  scale_colour_manual(values=colrep2)+xlim(0,48) +
  xlab("generation")

ggsave("figs/IVS2.png", plot = s, width = 10, height = 6, dpi = 600)

knitr::include_graphics("figs/IVS2.png")
```

<img src="figs/IVS2.png" width="6000" />

## IVS1

``` r
library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())
colrep3<-c("firebrick", "skyblue3", "chartreuse4")

# Annotation
# PPI251    ensembl exon    153 442 .   +   .   gene_id "pele"; transcript_id "pele1";
# PPI251    ensembl exon    501 1168    .   +   .   gene_id "pele"; transcript_id "pele1";
# PPI251    ensembl exon    1222    1947    .   +   .   gene_id "pele"; transcript_id "pele1";
# PPI251    ensembl exon    2138    2709    .   +   .   gene_id "pele"; transcript_id "pele1";

t<-read.table("/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/splicing-expression/splicing/spli.forr")

# head
# R1    G15 wf  se  PPI251  443-502 443 502 51  0.719474282680966
# R1    G15 wf  se  PPI251  1169-1223   1169    1223    51  0.719474282680966

names<-c("rep","generation","tissue","strand","gene","skey","start","end","rawcount","freq")
names(t)<-names
t<-subset(t,gene=="PPI251" & tissue=="wf" & strand=="se")
t$generation<-as.numeric(substring(t$generation, 2))

t<-subset(t,skey=="443-502")
missingvalues=data.frame(rep=c("R1","R2","R3"),generation=c(0,0,0),tissue=rep("wf",3),
                        strand=rep("se",3),gene=rep("PPI251",3),skey=rep("443-502",3),
                        start=rep(0,3),end=rep(0,3),rawcount=rep(0,3),freq=rep(0,3))

t<-rbind(t,missingvalues)

s <- ggplot(data = t, aes(linetype = rep)) +
  geom_line(data = t, aes(x = generation, y = freq, color = rep), linewidth = 1) +
  geom_point(data = t, aes(x = generation, y = freq, color = rep)) +
  theme(strip.text = element_blank(),
        legend.position = c(0.1, 0.85),
        legend.box.background = element_rect(color = "black", fill = "transparent"),
        legend.title = element_blank()) +
  ylab("splicing level IVS1 [srpm]") +
  scale_colour_manual(values=colrep3)+xlim(0,48) +
  xlab("generation")

ggsave("figs/IVS1.png", plot = s, width = 10, height = 6, dpi = 600)

knitr::include_graphics("figs/IVS1.png")
```

<img src="figs/IVS1.png" width="6000" />
