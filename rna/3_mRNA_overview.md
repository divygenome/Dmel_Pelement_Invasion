3 - mRNA overview
================
Matthew Beaumont
2023-07-17

``` bash
knitr::opts_chunk$set(echo = TRUE)
```

# mRNA overview

We generated a list of summary statistics for the mRNA reads for each
mapping software output.

``` bash
nohup zsh dmel_mRNA_overview.sh > ../logs/dmel_mRNA_overview.log
```

``` bash
#!/bin/bash

pyscript="/Volumes/Data/Projects/dmelR2_p-ele/scripts/mRNA-overview.py"
input_dir="/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/map-bwamem"
output_dir="/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/mRNA/overview_bwa"

# samtools view $if/gtdmelf_R1G5_run2.sort.bam | python $pyscript --sam - --sample-id R1-G5-wf > $of/wf_R1G5.txt

samtools view $input_dir/dmel_R1G6_run2.sort.bam | python $pyscript --sam - --sample-id R1-G6-wf > $output_dir/wf_R1G6.txt
samtools view $input_dir/dmel_R2G6_run2.sort.bam | python $pyscript --sam - --sample-id R2-G6-wf > $output_dir/wf_R2G6.txt
samtools view $input_dir/dmel_R3G6_run2.sort.bam | python $pyscript --sam - --sample-id R3-G6-wf > $output_dir/wf_R3G6.txt

samtools view $input_dir/dmel_R1G15_run2.sort.bam | python $pyscript --sam - --sample-id R1-G15-wf > $output_dir/wf_R1G15.txt
samtools view $input_dir/dmel_R2G15_run2.sort.bam | python $pyscript --sam - --sample-id R2-G15-wf > $output_dir/wf_R2G15.txt
samtools view $input_dir/dmel_R3G15_run2.sort.bam | python $pyscript --sam - --sample-id R3-G15-wf > $output_dir/wf_R3G15.txt

samtools view $input_dir/dmel_R1G21_run2.sort.bam | python $pyscript --sam - --sample-id R1-G21-wf > $output_dir/wf_R1G21.txt
samtools view $input_dir/dmel_R2G21_run2.sort.bam | python $pyscript --sam - --sample-id R2-G21-wf > $output_dir/wf_R2G21.txt
samtools view $input_dir/dmel_R3G21_run2.sort.bam | python $pyscript --sam - --sample-id R3-G21-wf > $output_dir/wf_R3G21.txt

samtools view $input_dir/dmel_R1G30_run2.sort.bam | python $pyscript --sam - --sample-id R1-G30-wf > $output_dir/wf_R1G30.txt
samtools view $input_dir/dmel_R2G30_run2.sort.bam | python $pyscript --sam - --sample-id R2-G30-wf > $output_dir/wf_R2G30.txt
samtools view $input_dir/dmel_R3G30_run2.sort.bam | python $pyscript --sam - --sample-id R3-G30-wf > $output_dir/wf_R3G30.txt

samtools view $input_dir/dmel_R1G40_run2.sort.bam | python $pyscript --sam - --sample-id R1-G40-wf > $output_dir/wf_R1G40.txt
samtools view $input_dir/dmel_R2G40_run2.sort.bam | python $pyscript --sam - --sample-id R2-G40-wf > $output_dir/wf_R2G40.txt
samtools view $input_dir/dmel_R3G40_run2.sort.bam | python $pyscript --sam - --sample-id R3-G40-wf > $output_dir/wf_R3G40.txt

# Combine outputs into single file, separating ID.
# cat *.txt|perl -pe 's/-/\t/'|perl -pe 's/-/\t/' > dmel_all.forr 
```

``` bash
#!/bin/bash

pyscript="/Volumes/Data/Projects/dmelR2_p-ele/scripts/mRNA-overview.py"
input_dir="/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/map-GSNAP/output"
output_dir="/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/mRNA/overview"

samtools view $input_dir/gt_R1G6.sort.bam | python $pyscript --sam - --sample-id R1-G6-wf   > $output_dir/wf_R1G6.txt
samtools view $input_dir/gt_R2G6.sort.bam | python $pyscript --sam - --sample-id R2-G6-wf   > $output_dir/wf_R2G6.txt
samtools view $input_dir/gt_R3G6.sort.bam | python $pyscript --sam - --sample-id R3-G6-wf   > $output_dir/wf_R3G6.txt

samtools view $input_dir/gt_R1G15.sort.bam | python $pyscript --sam - --sample-id R1-G15-wf   > $output_dir/wf_R1G15.txt
samtools view $input_dir/gt_R2G15.sort.bam | python $pyscript --sam - --sample-id R2-G15-wf   > $output_dir/wf_R2G15.txt
samtools view $input_dir/gt_R3G15.sort.bam | python $pyscript --sam - --sample-id R3-G15-wf   > $output_dir/wf_R3G15.txt

samtools view $input_dir/gt_R1G21.sort.bam | python $pyscript --sam - --sample-id R1-G21-wf   > $output_dir/wf_R1G21.txt
samtools view $input_dir/gt_R2G21.sort.bam | python $pyscript --sam - --sample-id R2-G21-wf   > $output_dir/wf_R2G21.txt
samtools view $input_dir/gt_R3G21.sort.bam | python $pyscript --sam - --sample-id R3-G21-wf   > $output_dir/wf_R3G21.txt

samtools view $input_dir/gt_R1G30.sort.bam | python $pyscript --sam - --sample-id R1-G30-wf   > $output_dir/wf_R1G30.txt
samtools view $input_dir/gt_R2G30.sort.bam | python $pyscript --sam - --sample-id R2-G30-wf   > $output_dir/wf_R2G30.txt
samtools view $input_dir/gt_R3G30.sort.bam | python $pyscript --sam - --sample-id R3-G30-wf   > $output_dir/wf_R3G30.txt

samtools view $input_dir/gt_R1G40.sort.bam | python $pyscript --sam - --sample-id R1-G40-wf   > $output_dir/wf_R1G40.txt
samtools view $input_dir/gt_R2G40.sort.bam | python $pyscript --sam - --sample-id R2-G40-wf   > $output_dir/wf_R2G40.txt
samtools view $input_dir/gt_R3G40.sort.bam | python $pyscript --sam - --sample-id R3-G40-wf   > $output_dir/wf_R3G40.txt

# Combine outputs into single file, separating ID.
# cat *.txt|perl -pe 's/-/\t/g'
# cat *.txt|perl -pe 's/-/\t/'|perl -pe 's/-/\t/' > dmel_all.forr 
```

Then visualised the most relevant ones in ggplot2.

## bwa-transcriptome

``` r
library(ggplot2)
library(scales)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
overview <- read.table("/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/mRNA/overview_bwa/dmel_all.forr", 
                       header = FALSE, sep = "\t")
colnames(overview) <- c("Replicate", "Generation", "Tissue", "Reads", "MappedReads", "MappedReadsWithMinMQ", 
                        "SenseGene", "AntisenseGene", "SenseTranscriptExon", "AntisenseTranscriptExon", 
                        "SensePelement", "AntisensePelement")

generation_order <- c("G6", "G15", "G21", "G30", "G40")

overview$Generation <- factor(overview$Generation, levels = generation_order)

overview$AntisenseGene <- -overview$AntisenseGene
overview$AntisensePelement <- -overview$AntisensePelement

mapped_reads_plot <- ggplot(overview, aes(x = Generation)) +
  geom_bar(aes(y = Reads, fill = "Reads"), stat = "identity", position = "dodge") +
  geom_bar(aes(y = MappedReads, fill = "Mapped reads"), stat = "identity", position = "dodge") +
  geom_bar(aes(y = MappedReadsWithMinMQ, fill = "Mapped reads w/ min mq"), stat = "identity", position = "dodge") +
  facet_wrap(. ~ Replicate, ncol = 1) +
  labs(x = "Generation", y = "Read counts") +
  ggtitle("Mapped reads & mapped reads w/ min mq") +
  scale_fill_manual(values = c("Reads" = "azure2", "Mapped reads" = "azure3", "Mapped reads w/ min mq" = "azure4"),
                    name = NULL,
                    breaks = c("Reads", "Mapped reads", "Mapped reads w/ min mq"),
                    labels = c("Reads", "Mapped reads", "Mapped reads w/ min mq"),
                    guide = guide_legend(reverse = TRUE)) +
                    theme_bw() +
                    theme(legend.position = "bottom") +
                    scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M"))

sense_antisense_gene_plot <- ggplot(overview, aes(x = Generation, y = SenseGene, fill = "Sense")) +
  geom_bar(stat = "identity") +
  geom_bar(aes(y = AntisenseGene, fill = "Antisense"), stat = "identity") +
  facet_wrap(. ~ Replicate, ncol = 1) +
  labs(x = "Generation", y = "Gene Totals") +
  ggtitle("Sense/antisense gene read counts") +
  scale_fill_manual(values = c("Sense" = "lightblue", "Antisense" = "darksalmon"),
                    name = NULL,
                    breaks = c("Sense", "Antisense"),
                    labels = c("Sense", "Antisense"),
                    guide = guide_legend(reverse = TRUE)) +
                    theme_bw() +
                    theme(legend.position = "bottom") +
                    scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M"))

sense_antisense_pelement_plot <- ggplot(overview, aes(x = Generation, y = SensePelement, fill = "Sense")) +
  geom_bar(stat = "identity") +
  geom_bar(aes(y = AntisensePelement, fill = "Antisense"), stat = "identity") +
  facet_wrap(. ~ Replicate, ncol = 1) +
  labs(x = "Generation", y = NULL) +
  ggtitle("Sense/antisense P-element read counts") +
  scale_fill_manual(values = c("Sense" = "lightblue", "Antisense" = "darksalmon"),
                    name = NULL,
                    breaks = c("Sense", "Antisense"),
                    labels = c("Sense", "Antisense"),
                    guide = guide_legend(reverse = TRUE)) +
                    theme_bw() +
                    theme(legend.position = "bottom") +
                    scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k"))

combined_plot <- cowplot::plot_grid(mapped_reads_plot, sense_antisense_gene_plot, sense_antisense_pelement_plot, ncol = 3)

ggsave("figs/mRNA_SS_bwa.png", combined_plot, width = 20, height = 10, dpi = 600)

knitr::include_graphics("figs/mRNA_SS_bwa.png")
```

<img src="figs/mRNA_SS_bwa.png" width="12000" />

## bwa-genome

``` r
library(ggplot2)
library(scales)
library(dplyr)

overview <- read.table("/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/mRNA/overview_bwa_genome/dmel_all.forr", 
                       header = FALSE, sep = "\t")
colnames(overview) <- c("Replicate", "Generation", "Tissue", "Reads", "MappedReads", "MappedReadsWithMinMQ", 
                        "SenseGene", "AntisenseGene", "SenseTranscriptExon", "AntisenseTranscriptExon", 
                        "SensePelement", "AntisensePelement")

generation_order <- c("G6", "G15", "G21", "G30", "G40")

overview$Generation <- factor(overview$Generation, levels = generation_order)

overview$AntisenseGene <- -overview$AntisenseGene
overview$AntisensePelement <- -overview$AntisensePelement

mapped_reads_plot <- ggplot(overview, aes(x = Generation)) +
  geom_bar(aes(y = Reads, fill = "Reads"), stat = "identity", position = "dodge") +
  geom_bar(aes(y = MappedReads, fill = "Mapped reads"), stat = "identity", position = "dodge") +
  geom_bar(aes(y = MappedReadsWithMinMQ, fill = "Mapped reads w/ min mq"), stat = "identity", position = "dodge") +
  facet_wrap(. ~ Replicate, ncol = 1) +
  labs(x = "Generation", y = "Read counts") +
  ggtitle("Mapped reads & mapped reads w/ min mq") +
  scale_fill_manual(values = c("Reads" = "azure2", "Mapped reads" = "azure3", "Mapped reads w/ min mq" = "azure4"),
                    name = NULL,
                    breaks = c("Reads", "Mapped reads", "Mapped reads w/ min mq"),
                    labels = c("Reads", "Mapped reads", "Mapped reads w/ min mq"),
                    guide = guide_legend(reverse = TRUE)) +
                    theme_bw() +
                    theme(legend.position = "bottom") +
                    scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M"))

sense_antisense_pelement_plot <- ggplot(overview, aes(x = Generation, y = SensePelement, fill = "Sense")) +
  geom_bar(stat = "identity") +
  geom_bar(aes(y = AntisensePelement, fill = "Antisense"), stat = "identity") +
  facet_wrap(. ~ Replicate, ncol = 1) +
  labs(x = "Generation", y = NULL) +
  ggtitle("Sense/antisense P-element read counts") +
  scale_fill_manual(values = c("Sense" = "lightblue", "Antisense" = "darksalmon"),
                    name = NULL,
                    breaks = c("Sense", "Antisense"),
                    labels = c("Sense", "Antisense"),
                    guide = guide_legend(reverse = TRUE)) +
                    theme_bw() +
                    theme(legend.position = "bottom") +
                    scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k"))

combined_plot <- cowplot::plot_grid(mapped_reads_plot, sense_antisense_pelement_plot, ncol = 2)

ggsave("figs/mRNA_SS_bwa_g.png", combined_plot, width = 15, height = 10, dpi = 600)

knitr::include_graphics("figs/mRNA_SS_bwa_g.png")
```

<img src="figs/mRNA_SS_bwa_g.png" width="9000" />

## GSNAP-transcriptome

``` r
library(ggplot2)
library(scales)
library(dplyr)

overview <- read.table("/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/mRNA/overview/dmel_all.forr", 
                       header = FALSE, sep = "\t")
colnames(overview) <- c("Replicate", "Generation", "Tissue", "Reads", "MappedReads", "MappedReadsWithMinMQ", 
                        "SenseGene", "AntisenseGene", "SenseTranscriptExon", "AntisenseTranscriptExon", 
                        "SensePelement", "AntisensePelement")

generation_order <- c("G6", "G15", "G21", "G30", "G40")

overview$Generation <- factor(overview$Generation, levels = generation_order)

overview$AntisenseGene <- -overview$AntisenseGene
overview$AntisensePelement <- -overview$AntisensePelement

mapped_reads_plot <- ggplot(overview, aes(x = Generation)) +
  geom_bar(aes(y = Reads, fill = "Reads"), stat = "identity", position = "dodge") +
  geom_bar(aes(y = MappedReads, fill = "Mapped reads"), stat = "identity", position = "dodge") +
  geom_bar(aes(y = MappedReadsWithMinMQ, fill = "Mapped reads w/ min mq"), stat = "identity", position = "dodge") +
  facet_wrap(. ~ Replicate, ncol = 1) +
  labs(x = "Generation", y = "Read counts") +
  ggtitle("Mapped reads & mapped reads w/ min mq") +
  scale_fill_manual(values = c("Reads" = "azure2", "Mapped reads" = "azure3", "Mapped reads w/ min mq" = "azure4"),
                    name = NULL,
                    breaks = c("Reads", "Mapped reads", "Mapped reads w/ min mq"),
                    labels = c("Reads", "Mapped reads", "Mapped reads w/ min mq"),
                    guide = guide_legend(reverse = TRUE)) +
                    theme_bw() +
                    theme(legend.position = "bottom") +
                    scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M"))

sense_antisense_gene_plot <- ggplot(overview, aes(x = Generation, y = SenseGene, fill = "Sense")) +
  geom_bar(stat = "identity") +
  geom_bar(aes(y = AntisenseGene, fill = "Antisense"), stat = "identity") +
  facet_wrap(. ~ Replicate, ncol = 1) +
  labs(x = "Generation", y = "Gene Totals") +
  ggtitle("Sense/antisense gene read counts") +
  scale_fill_manual(values = c("Sense" = "lightblue", "Antisense" = "darksalmon"),
                    name = NULL,
                    breaks = c("Sense", "Antisense"),
                    labels = c("Sense", "Antisense"),
                    guide = guide_legend(reverse = TRUE)) +
                    theme_bw() +
                    theme(legend.position = "bottom") +
                    scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M"))

sense_antisense_pelement_plot <- ggplot(overview, aes(x = Generation, y = SensePelement, fill = "Sense")) +
  geom_bar(stat = "identity") +
  geom_bar(aes(y = AntisensePelement, fill = "Antisense"), stat = "identity") +
  facet_wrap(. ~ Replicate, ncol = 1) +
  labs(x = "Generation", y = NULL) +
  ggtitle("Sense/antisense P-element read counts") +
  scale_fill_manual(values = c("Sense" = "lightblue", "Antisense" = "darksalmon"),
                    name = NULL,
                    breaks = c("Sense", "Antisense"),
                    labels = c("Sense", "Antisense"),
                    guide = guide_legend(reverse = TRUE)) +
                    theme_bw() +
                    theme(legend.position = "bottom") +
                    scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k"))

combined_plot <- cowplot::plot_grid(mapped_reads_plot, sense_antisense_gene_plot, sense_antisense_pelement_plot, ncol = 3)

ggsave("figs/mRNA_SS_GSNAP_tm.png", combined_plot, width = 20, height = 10, dpi = 600)

knitr::include_graphics("figs/mRNA_SS_GSNAP_tm.png")
```

<img src="figs/mRNA_SS_GSNAP_tm.png" width="12000" />

## GSNAP-genome

``` r
library(ggplot2)
library(scales)
library(dplyr)

overview <- read.table("/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/mRNA/overview_genome/dmel_all.forr", 
                       header = FALSE, sep = "\t")
colnames(overview) <- c("Replicate", "Generation", "Tissue", "Reads", "MappedReads", "MappedReadsWithMinMQ", 
                        "SenseGene", "AntisenseGene", "SenseTranscriptExon", "AntisenseTranscriptExon", 
                        "SensePelement", "AntisensePelement")

generation_order <- c("G6", "G15", "G21", "G30", "G40")

overview$Generation <- factor(overview$Generation, levels = generation_order)

overview$AntisenseGene <- -overview$AntisenseGene
overview$AntisensePelement <- -overview$AntisensePelement

mapped_reads_plot <- ggplot(overview, aes(x = Generation)) +
  geom_bar(aes(y = Reads, fill = "Reads"), stat = "identity", position = "dodge") +
  geom_bar(aes(y = MappedReads, fill = "Mapped reads"), stat = "identity", position = "dodge") +
  geom_bar(aes(y = MappedReadsWithMinMQ, fill = "Mapped reads w/ min mq"), stat = "identity", position = "dodge") +
  facet_wrap(. ~ Replicate, ncol = 1) +
  labs(x = "generation", y = "read counts") +
  ggtitle("Mapped reads & mapped reads w/ min mq") +
  scale_fill_manual(values = c("Reads" = "azure2", "Mapped reads" = "azure3", "Mapped reads w/ min mq" = "azure4"),
                    name = NULL,
                    breaks = c("Reads", "Mapped reads", "Mapped reads w/ min mq"),
                    labels = c("Reads", "Mapped reads", "Mapped reads w/ min mq"),
                    guide = guide_legend(reverse = TRUE)) +
                    theme_bw() +
                    theme(legend.position = "bottom") +
                    scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M"))

sense_antisense_pelement_plot <- ggplot(overview, aes(x = Generation, y = SensePelement, fill = "Sense")) +
  geom_bar(stat = "identity") +
  geom_bar(aes(y = AntisensePelement, fill = "Antisense"), stat = "identity") +
  facet_wrap(. ~ Replicate, ncol = 1) +
  labs(x = "generation", y = NULL) +
  ggtitle("Sense/antisense P-element read counts") +
  scale_fill_manual(values = c("Sense" = "lightblue", "Antisense" = "darksalmon"),
                    name = NULL,
                    breaks = c("Sense", "Antisense"),
                    labels = c("Sense", "Antisense"),
                    guide = guide_legend(reverse = TRUE)) +
                    theme_bw() +
                    theme(legend.position = "bottom") +
                    scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k"))

combined_plot <- cowplot::plot_grid(mapped_reads_plot, sense_antisense_pelement_plot, ncol = 2)

ggsave("figs/mRNA_SS_GSNAP_g.png", combined_plot, width = 15, height = 10, dpi = 600)

knitr::include_graphics("figs/mRNA_SS_GSNAP_g.png")
```

<img src="figs/mRNA_SS_GSNAP_g.png" width="9000" />
