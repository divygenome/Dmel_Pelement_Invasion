4 - Expression & splicing
================
Matthew Beaumont
2024-10-30

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
#cat *.txt| perl -pe 's/-/\mrna_data/'|perl -pe 's/-/\mrna_data/' > expr.forr 

# Tissue column for 'wf' needed adding back.
#awk 'BEGIN{OFS="\mrna_data"}{$2 = $2 "\mrna_data" "wf"; print}' expr.forr > expr_wf.forr
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
# cat *.txt| perl -pe 's/-/\mrna_data/'|perl -pe 's/-/\mrna_data/' > spli.forr
```

And then we visualised both the P-element expression and the splicing of
its’ three introns together in ggplot2.

``` r
library(ggplot2)
library(svglite)
theme_set(theme_bw())

target<-"PPI251"
ttissue<-"wf" # target tissue
sminfreq<-0.1

h<-read.table("/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/splicing-expression/expression/expr_wf.forr")
names(h)<-c("rep","time","tissue","strand","gene","pos","cov")

a<-subset(h,gene==target & tissue==ttissue)
a$key<-paste0(a$rep,"_",a$time,"_",a$pos,"_",a$strand)


a$time <- factor(a$time, levels=c("G0", "G6", "G15", "G21", "G30", "G40", "G50", "G100"))
s<-subset(a,strand=="se")
as<-subset(a,strand=="ase")

spli<-read.table("/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/splicing-expression/splicing/spli.forr")

names(spli)<-c("rep","time","tissue","strand","gene","skey","start","end","rawcount","freq")
aspli<-subset(spli,gene==target & freq>sminfreq & tissue==ttissue)
aspli$time <- factor(aspli$time, levels=c("G0", "G6", "G15", "G21", "G30", "G40", "G50", "G100"))
aspli$rep<- as.factor(aspli$rep)
aspli$start<-aspli$start-1 # position inaccuracy, graph is more appealing
aspli$keystart<-paste0(aspli$rep,"_",aspli$time,"_",aspli$start,"_",aspli$strand)
aspli$keyend<-paste0(aspli$rep,"_",aspli$time,"_",aspli$end,"_",aspli$strand)

aspli<-merge(x=aspli,y=a[,c("key","cov")],by.x="keystart",by.y="key")
aspli<-merge(x=aspli,y=a[,c("key","cov")],by.x="keyend",by.y="key")
aspli$size<-log(aspli$freq+1)

a_s<-subset(aspli,strand=="se")
a_as<-subset(aspli,strand=="ase")

a_s$size_scaled <- a_s$size - 1
a_as$size_scaled <- a_as$size -1 

expr_spli_plot <- ggplot() +
  geom_polygon(data = s, mapping = aes(x = pos, y = cov), fill = 'grey', color = 'grey') +
  geom_polygon(data = as, aes(x = pos, y = -cov), fill = 'lightgrey', color = 'lightgrey') +
  geom_curve(data = a_s, mapping = aes(x = start, y = cov.x, xend = end, yend = cov.y, linewidth = size_scaled), curvature = -0.15, ncp = 10, show.legend = FALSE) +
  geom_curve(data = a_as, mapping = aes(x = start, y = -cov.x, xend = end, yend = -cov.y, linewidth = size_scaled), curvature = 0.15, ncp = 10, show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "black") + # Add black line at Y = 0
  facet_grid(time ~ rep) +
  scale_size(range = c(0.1, 1)) + # Adjust the range for line thickness
  coord_cartesian(ylim = c(min(a$cov) -1, max(a$cov) + 1)) + # Adjust y-axis limits
  xlab("position") +
  ylab("expression [rpm]") +
  theme(
    panel.grid.major = element_blank(), # Remove major gridlines
    panel.grid.minor = element_blank(), # Remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 1), # Add border around each plot
    strip.background = element_rect(fill = "grey80", color = "black"), # Shaded grey box for facet titles
    strip.text = element_text(size = 10),
    panel.spacing = unit(0.4, "lines"), # Reduce space between plots
    legend.position = "bottom", # Move legend to the bottom
    legend.direction = "horizontal" # Make the legend horizontal
  )
```

    ## Warning: The `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
# Save plots
ggsave("figs/expr_spli.svg", expr_spli_plot, width = 14, height = 8)
ggsave("figs/expr_spli.png", expr_spli_plot, width = 14, height = 8, dpi = 600)

# Display the plot
knitr::include_graphics("figs/expr_spli.png")
```

<img src="figs/expr_spli.png" width="8400" />

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
# cat *.txt| perl -pe 's/-/\mrna_data/'|perl -pe 's/-/\mrna_data/' > expr.forr
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

# Read the dataset
mrna_data <- read_delim("/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/splicing-expression/raw_expression/all-expressionlevel/expr_fixed.forr", delim = "\t", col_names = FALSE, comment = "#", show_col_types = FALSE)

# Set the column names as per your data
names(mrna_data) <- c("replicate", "generation", "run", "gene", "rawse", "rawase", "genlen", "sense", "antisense", "total")

# Filter the data for the specific gene
mrna_data <- subset(mrna_data, gene == "PPI251")

# Convert generation to numeric if necessary
mrna_data$generation <- as.numeric(substring(mrna_data$generation, 2))

# Add data point for gen 0 with expression level 0
zero_data <- data.frame(replicate = rep(c("R1", "R2", "R3"), each = 1),
                        generation = rep(0, 3),
                        sense = rep(0, 3),
                        stringsAsFactors = FALSE)
mrna_data <- bind_rows(mrna_data, zero_data)

# Define plot dimensions and resolution
width <- 16
height <- 12
resolution <- 600

# First plot
mrna_1 <- ggplot(data = mrna_data, aes(linetype = replicate)) +
  geom_line(data = mrna_data, aes(x = generation, y = total, color = replicate), linewidth = 1) +
  geom_point(data = mrna_data, aes(x = generation, y = total, color = replicate)) +
  theme(strip.text = element_blank(), 
        legend.position = c(0.1, 0.85),
        legend.box.background = element_rect(color = "black", fill = "transparent"),
        legend.title = element_blank()) + 
  ylab("expression [rpkm]") +
  scale_colour_manual(values = tresrep) +
  xlim(0, 105) +
  xlab("generation")

ggsave("figs/P-ele_expression.png", plot = mrna_1, width = 10, height = 6, dpi = 600)
```

    ## Warning: Removed 3 rows containing missing values (`geom_line()`).

    ## Warning: Removed 3 rows containing missing values (`geom_point()`).

``` r
knitr::include_graphics("figs/P-ele_expression.png")
```

<img src="figs/P-ele_expression.png" width="6000" />

``` r
# Second plot
mrna_2 <- ggplot(data = mrna_data, aes(x = generation, y = total, color = replicate, shape = replicate)) + 
  geom_point(size = 3.5) +
  geom_line(size = 1.5) +
  scale_colour_manual(values = tresrep) +
  scale_size(guide = "none") +  
  ylab("expression [rpkm]") +
  xlim(0, 105) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        plot.margin = margin(10, 10, 10, 10, "pt"),
        legend.position = c(0.08, 0.85),
        legend.box.background = element_rect(fill = NULL, colour = NULL, size = 0, linetype = NULL, color = NULL)) +
  guides(color = guide_legend(title = NULL), shape = guide_legend(title = NULL))
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
mrna_sl <- mrna_2 + guides(color = guide_legend(title = NULL), shape = guide_legend(title = NULL))

ggsave("figs/P-ele_expression_2.png", plot = mrna_sl, width = 6, height = 4, dpi = 600)
```

    ## Warning: Removed 3 rows containing missing values (`geom_point()`).

    ## Warning: Removed 3 rows containing missing values (`geom_line()`).

``` r
ggsave("figs/P-ele_expression_2.svg", plot = mrna_sl,  width = 6, height = 4)
```

    ## Warning: Removed 3 rows containing missing values (`geom_point()`).
    ## Removed 3 rows containing missing values (`geom_line()`).

``` r
knitr::include_graphics("figs/P-ele_expression_2.png")
```

<img src="figs/P-ele_expression_2.png" width="3600" />

``` r
library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())
custom_palette <- c("firebrick", "skyblue3", "chartreuse4")

mrna_data <- read_delim("/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/splicing-expression/raw_expression/all-expressionlevel/expr_fixed.forr", delim = "\t", col_names = FALSE, comment = "#", show_col_types = FALSE)
names(mrna_data) <- c("replicate", "generation", "run", "gene", "rawse", "rawase", "genlen", "sense", "antisense", "total")
mrna_data <- subset(mrna_data, gene == "PPI251")
mrna_data$generation <- as.numeric(substring(mrna_data$generation, 2))

# Add data point for gen 0 with expression level 0
zero_data <- data.frame(replicate = rep(c("R1", "R2", "R3"), each = 1),
                        generation = rep(0, 3),
                        sense = rep(0, 3),
                        stringsAsFactors = FALSE)
mrna_data <- bind_rows(mrna_data, zero_data)

width <- 20
height <- 12
resolution <- 600

mrna_data <- mrna_data %>%
  gather(key = "expression_type", value = "value", sense, antisense) %>%
  mutate(linetype = ifelse(expression_type == "sense", "Sense", "Antisense"))

mrna_sa <- ggplot(data = mrna_data, aes(x = generation, y = value, color = replicate, linetype = linetype, shape = replicate)) + 
  geom_line(size = 1.2) +
  geom_point(size = 3) +  
  scale_colour_manual(values = custom_palette) +
  scale_linetype_manual(values = c(Sense = 1, Antisense = 6), breaks = c("Sense", "Antisense")) + 
  ylab("expression [rpkm]") +
  xlim(0, 105) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        plot.margin = margin(10, 10, 10, 10, "pt"),
        legend.position = c(0.12, 0.728),
        legend.box.background = element_rect(fill = NULL, colour = NULL, size = 0, linetype = NULL, color = NULL)) +
  guides(color = guide_legend(title = NULL), linetype = guide_legend(title = NULL), shape = guide_legend(title = NULL))

ggsave("figs/P-ele_expression_sa.png", plot = mrna_sa, width = 6, height = 4, dpi = 600)
```

    ## Warning: Removed 3 rows containing missing values (`geom_line()`).

    ## Warning: Removed 3 rows containing missing values (`geom_point()`).

``` r
ggsave("figs/P-ele_expression_sa.svg", plot = mrna_sa, width = 6, height = 4)
```

    ## Warning: Removed 3 rows containing missing values (`geom_line()`).
    ## Removed 3 rows containing missing values (`geom_point()`).

``` r
knitr::include_graphics("figs/P-ele_expression_sa.png")
```

<img src="figs/P-ele_expression_sa.png" width="3600" />

Then we visualised the splicing level of all three introns of the
P-element over time.

## IVS3

``` r
library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())
colrep3 <- c("firebrick", "skyblue3", "chartreuse4")

# Read the data
sp3 <- read.table("/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/splicing-expression/splicing/spli.forr")

# Define the column names
names <- c("rep", "generation", "tissue", "strand", "gene", "skey", "start", "end", "rawcount", "freq")
names(sp3) <- names

# Filter the data
sp3 <- subset(sp3, gene == "PPI251" & tissue == "wf" & strand == "se")
sp3$generation <- as.numeric(substring(sp3$generation, 2))
sp3 <- subset(sp3, skey == "1948-2139")  # delta23

# Add zero data points for missing values in dataset
zero_data <- data.frame(rep = "R1", generation = 6, tissue = "wf", strand = "se", gene = "PPI251",
                        skey = "1948-2139", start = 0, end = 0, rawcount = 0, freq = 0)
sp3 <- rbind(sp3, zero_data)

zero_data_2 <- data.frame(rep = c("R1", "R2", "R3"), generation = 0, tissue = "wf", strand = "se", gene = "PPI251",
                          skey = "1948-2139", start = 0, end = 0, rawcount = 0, freq = 0)
sp3 <- rbind(sp3, zero_data_2)

zero_data_3 <- data.frame(rep = c("R1", "R3"), generation = 100, tissue = "wf", strand = "se", gene = "PPI251",
                          skey = "1948-2139", start = 0, end = 0, rawcount = 0, freq = 0)
sp3 <- rbind(sp3, zero_data_3)

# Ensure 'rep' is treated as a factor
sp3$rep <- as.factor(sp3$rep)

# Plot 1
spli_3 <- ggplot(data = sp3, aes(linetype = rep)) +
  geom_line(data = sp3, aes(x = generation, y = freq, color = rep), linewidth = 1) +
  geom_point(data = sp3, aes(x = generation, y = freq, color = rep)) +
  theme(strip.text = element_blank(),
        legend.position = "none") +
  ylab("splicing level [srpm]") + 
  scale_colour_manual(values = colrep3) +
  xlab("generation")

ggsave("figs/IVS3.png", plot = spli_3, width = 5, height = 3, dpi = 600)
knitr::include_graphics("figs/IVS3.png")
```

<img src="figs/IVS3.png" width="3000" />

``` r
# Plot 2
spli_3_x <- ggplot(data = sp3, aes(color = rep, shape = rep)) + 
  geom_point(aes(x = generation, y = freq), size = 3.5) +
  geom_line(aes(x = generation, y = freq), size = 1.5) +
  scale_colour_manual(values = colrep3) +
  scale_size(guide = "none") +  
  ylab("splicing level [srpm]") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        plot.margin = margin(10, 10, 10, 10, "pt"),
        legend.position = c(0.08, 0.85),
        legend.box.background = element_rect(fill = NULL, colour = NULL, size = 0, linetype = NULL, color = NULL)) +
  guides(color = guide_legend(title = NULL), shape = guide_legend(title = NULL))

spli_3_sl <- spli_3_x + guides(color = guide_legend(title = NULL), shape = guide_legend(title = NULL))

ggsave("figs/IVS3_sl.png", plot = spli_3_sl, width = 6, height = 4, dpi = 600)
knitr::include_graphics("figs/IVS3_sl.png")
```

<img src="figs/IVS3_sl.png" width="3600" />

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

sp2<-read.table("/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/splicing-expression/splicing/spli.forr")

# head
# R1    G15 wf  se  PPI251  443-502 443 502 51  0.719474282680966
# R1    G15 wf  se  PPI251  1169-1223   1169    1223    51  0.719474282680966

names<-c("rep","generation","tissue","strand","gene","skey","start","end","rawcount","freq")
names(sp2)<-names
sp2<-subset(sp2,gene=="PPI251" & tissue=="wf" & strand=="se")
sp2$generation<-as.numeric(substring(sp2$generation, 2))

sp2<-subset(sp2,skey=="1169-1223")

zero_data <- data.frame(rep = c("R1", "R2", "R3"), generation = 0, tissue = "wf", strand = "se", gene = "PPI251",
                          skey = "1169-1223", start = 0, end = 0, rawcount = 0, freq = 0)
sp2 <- rbind(sp2, zero_data)

spli_2 <- ggplot(data = sp2, aes(linetype = rep)) +
  geom_line(data = sp2, aes(x = generation, y = freq, color = rep), linewidth = 1) +
  geom_point(data = sp2, aes(x = generation, y = freq, color = rep)) +
  theme(strip.text = element_blank(),
        legend.position = "none") +
        #legend.position = c(0.1, 0.85),
        #legend.box.background = element_rect(color = "black", fill = "transparent"),
        #legend.title = element_blank()) +
  ylab("splicing level [srpm]") +
  scale_colour_manual(values=colrep2) + 
  #xlim(0,48) +
  xlab("generation")

ggsave("figs/IVS2.png", plot = spli_2, width = 5, height = 3, dpi = 600)

knitr::include_graphics("figs/IVS2.png")
```

<img src="figs/IVS2.png" width="3000" />

``` r
spli_2_x <- ggplot(data = sp2, aes(color = rep, shape = rep)) + 
  geom_point(aes(x = generation, y = freq), size = 3.5) +
  geom_line(aes(x = generation, y = freq), size = 1.5) +
  scale_colour_manual(values = colrep3) +
  scale_size(guide = "none") +  
  ylab("splicing level [srpm]") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        plot.margin = margin(10, 10, 10, 10, "pt"),
        legend.position = c(0.08, 0.85),
        legend.box.background = element_rect(fill = NULL, colour = NULL, size = 0, linetype = NULL, color = NULL)) +
  guides(color = guide_legend(title = NULL), shape = guide_legend(title = NULL))

spli_2_sl <- spli_2_x + guides(color = guide_legend(title = NULL))

ggsave("figs/IVS2_sl.png", plot = spli_2_sl, width = 6, height = 4, dpi = 600)

knitr::include_graphics("figs/IVS2_sl.png")
```

<img src="figs/IVS2_sl.png" width="3600" />

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

sp1<-read.table("/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/splicing-expression/splicing/spli.forr")

# head
# R1    G15 wf  se  PPI251  443-502 443 502 51  0.719474282680966
# R1    G15 wf  se  PPI251  1169-1223   1169    1223    51  0.719474282680966

names<-c("rep","generation","tissue","strand","gene","skey","start","end","rawcount","freq")
names(sp1)<-names
sp1<-subset(sp1,gene=="PPI251" & tissue=="wf" & strand=="se")
sp1$generation<-as.numeric(substring(sp1$generation, 2))

sp1<-subset(sp1,skey=="443-502")

zero_data <- data.frame(rep = c("R1", "R2", "R3"), generation = 0, tissue = "wf", strand = "se", gene = "PPI251",
                          skey = "443-502", start = 0, end = 0, rawcount = 0, freq = 0)
sp1 <- rbind(sp1, zero_data)

spli_1 <- ggplot(data = sp1, aes(linetype = rep)) +
  geom_line(data = sp1, aes(x = generation, y = freq, color = rep), linewidth = 1) +
  geom_point(data = sp1, aes(x = generation, y = freq, color = rep)) +
  theme(strip.text = element_blank(),
        legend.position = c(0.1, 0.75),
        legend.box.background = element_rect(color = "black", fill = "transparent"),
        legend.title = element_blank()) +
  ylab("splicing level [srpm]") +
  scale_colour_manual(values=colrep3) + 
  #xlim(0,48) +
  xlab("generation")

ggsave("figs/IVS1.png", plot = spli_1, width = 5, height = 3, dpi = 600)

knitr::include_graphics("figs/IVS1.png")
```

<img src="figs/IVS1.png" width="3000" />

``` r
spli_1_x <- ggplot(data = sp1, aes(color = rep, shape = rep)) + 
  geom_point(aes(x = generation, y = freq), size = 3.5) +
  geom_line(aes(x = generation, y = freq), size = 1.5) +
  scale_colour_manual(values = colrep3) +
  scale_size(guide = "none") +  
  ylab("splicing level [srpm]") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        plot.margin = margin(10, 10, 10, 10, "pt"),
        legend.position = c(0.08, 0.85),
        legend.box.background = element_rect(fill = NULL, colour = NULL, size = 0, linetype = NULL, color = NULL)) +
  guides(color = guide_legend(title = NULL), shape = guide_legend(title = NULL))

spli_1_sl <- spli_1_x + guides(color = guide_legend(title = NULL))

ggsave("figs/IVS1_sl.png", plot = spli_1_sl, width = 6, height = 4, dpi = 600)

knitr::include_graphics("figs/IVS1_sl.png")
```

<img src="figs/IVS1_sl.png" width="3600" />

``` r
spli_1_xs <- spli_1_sl + ggtitle("IVS1")
spli_2_xs <- spli_2_sl + ggtitle("IVS2") +
  ylab(NULL) +
  theme(legend.position = "none")
spli_3_xs <- spli_3_sl + ggtitle("IVS3") +
  ylab(NULL) +
  theme(legend.position = "none")

spli_comb <- cowplot::plot_grid(spli_1_xs, spli_2_xs, spli_3_xs, ncol = 3)

ggsave("figs/spli_comb.png", plot = spli_comb, width = 21, height = 2, dpi = 600)
ggsave("figs/spli_comb.svg", plot = spli_comb, width = 21, height = 2)

knitr::include_graphics("figs/spli_comb.png")
```

<img src="figs/spli_comb.png" width="12600" />
