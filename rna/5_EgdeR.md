5 - Edge-R - differential expression
================
Matthew Beaumont
2023-07-17

``` bash
knitr::opts_chunk$set(echo = TRUE)
```

Prepped for addition of naive datasets.

``` r
library(tidyverse)
library(RColorBrewer)
library(dplyr)
library(edgeR)
library(data.table)
library(gridExtra)

theme_set(theme_bw())
pallette <- c("darkgrey","orange","red")
colour<-c("grey50","#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628","#f781bf")
cbbPalette<-colour

d<-read_delim("/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/splicing-expression/raw_expression/all-expressionlevel/expr.forr",delim="\t",col_names=FALSE,comment="#")
names(d)<-c("replicate","generation","run","gene","rawse","rawase","genlen","se","ase","total")


dgeme <- function(cou,gro) {
  y <- DGEList(counts=cou,group=gro)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~gro)
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y,design)
  lrt <- glmQLFTest(fit,coef=2)
  return(lrt)
  
}

edgeme <- function(n1,n2,n3,e1,e2) {
  merge<- n1 %>%  inner_join(n2, by = "gene")  %>%  inner_join(n3, by = "gene") %>%  inner_join(e1, by = "gene") %>%  inner_join(e2, by = "gene")
  s<-setDT(merge)
  dm<-data.matrix(s[,2:6])
  rownames(dm)<-merge$gene
  group <- factor(c(1,1,1,2,2))
  lrt<-dgeme(dm,group)
  return(lrt)
}

edgemeext <- function(n1,n2,n3,e5,e1,e2) {
  merge<- n1 %>%  inner_join(n2, by = "gene")  %>%  inner_join(n3, by = "gene") %>% inner_join(e5, by = "gene") %>%  inner_join(e1, by = "gene") %>%  inner_join(e2, by = "gene")
  s<-setDT(merge)
  dm<-data.matrix(s[,2:7])
  rownames(dm)<-merge$gene
  group <- factor(c(1,1,1,1,2,2))
  lrt<-dgeme(dm,group)
  return(lrt)
}

# Need naive data for diff expression.

naive<-subset(d,generation=="naive")
nR1<-subset(naive,replicate=="R1") %>% select(gene,rawse) %>% rename(naiveR1=rawse)
nR2<-subset(naive,replicate=="R2") %>% select(gene,rawse) %>% rename(naiveR2=rawse)
nR3<-subset(naive,replicate=="R3") %>% select(gene,rawse) %>% rename(naiveR3=rawse)


rep1<-subset(d,replicate=="R1")
eg6R1<-subset(rep1,generation=="G6") %>% select(gene,rawse) %>% rename(eg6R1=rawse)
eg15R1<-subset(rep1,generation=="G15") %>% select(gene,rawse) %>% rename(eg15R1=rawse)
eg21R1<-subset(rep1,generation=="G21") %>% select(gene,rawse) %>% rename(eg21R1=rawse)
eg30R1<-subset(rep1,generation=="G30") %>% select(gene,rawse) %>% rename(eg30R1=rawse)
eg40R1<-subset(rep1,generation=="G40") %>% select(gene,rawse) %>% rename(eg40R1=rawse)


rep2<-subset(d,replicate=="R2")
eg6R2<-subset(rep2,generation=="G6") %>% select(gene,rawse) %>% rename(eg6R2=rawse)
eg15R2<-subset(rep2,generation=="G15") %>% select(gene,rawse) %>% rename(eg15R2=rawse)
eg21R2<-subset(rep2,generation=="G21") %>% select(gene,rawse) %>% rename(eg21R2=rawse)
eg30R2<-subset(rep2,generation=="G30") %>% select(gene,rawse) %>% rename(eg30R2=rawse)
eg40R2<-subset(rep2,generation=="G40") %>% select(gene,rawse) %>% rename(eg40R2=rawse)

rep3<-subset(d,replicate=="R3")
eg6R3<-subset(rep3,generation=="G6") %>% select(gene,rawse) %>% rename(eg6R3=rawse)
eg15R3<-subset(rep3,generation=="G15") %>% select(gene,rawse) %>% rename(eg15R3=rawse)
eg21R3<-subset(rep3,generation=="G21") %>% select(gene,rawse) %>% rename(eg21R3=rawse)
eg30R3<-subset(rep3,generation=="G30") %>% select(gene,rawse) %>% rename(eg30R3=rawse)
eg40R3<-subset(rep3,generation=="G40") %>% select(gene,rawse) %>% rename(eg40R3=rawse)

#lm1<-edgeme(nR1,nR2,nR3,eg30R1,eg40R1)
#lm2<-edgeme(nR1,nR2,nR3,eg30R2,eg40R2)
#lm3<-edgeme(nR1,nR2,nR3,eg30R3,eg40R3)
#dflm1<-data.frame(gene=rownames(lm1$table),logfc=lm1$table$logFC, pval=lm1$table$PValue, fdr=p.adjust(lm1$table$PValue,method="BH"))
#dflm2<-data.frame(gene=rownames(lm2$table),logfc=lm2$table$logFC, pval=lm2$table$PValue, fdr=p.adjust(lm2$table$PValue,method="BH"))
#dflm3<-data.frame(gene=rownames(lm3$table),logfc=lm3$table$logFC, pval=lm3$table$PValue, fdr=p.adjust(lm3$table$PValue,method="BH"))

lmvs2<-edgemeext(eg30R1,eg40R1,eg30R2,eg40R2,eg30R3,eg40R3)
dflmvs2<-data.frame(gene=rownames(lmvs2$table),logfc=lmvs2$table$logFC, pval=lmvs2$table$PValue,fdr=p.adjust(as.numeric(lmvs2$table$PValue),method="BH"))

#### Colouring FDR 
fdr<-0.01
fdr2<-0.001

#dflm1$col<-0
#dflm1[dflm1$fdr<fdr,]$col<-1
#dflm1[dflm1$fdr<fdr2,]$col<-2
#dflm1$col<-as.factor(dflm1$col)

#dflm2$col<-0
#dflm2[dflm2$fdr<fdr,]$col<-1
#dflm2[dflm2$fdr<fdr2,]$col<-2
#dflm2$col<-as.factor(dflm2$col)

#dflm3$col<-0
#dflm3[dflm3$fdr<fdr,]$col<-1
#dflm3[dflm3$fdr<fdr2,]$col<-2
#dflm3$col<-as.factor(dflm3$col)

dflmvs2$col <- 0
dflmvs2[dflmvs2$fdr < fdr, "col"] <- 1
dflmvs2[dflmvs2$fdr < fdr2, "col"] <- 2
dflmvs2$col <- as.factor(dflmvs2$col)

pvs2<-ggplot(data = dflmvs2, aes(y = -log10(pval), x = logfc,color=col )) +
  geom_point(alpha = 1, size = 2)+scale_colour_manual(values=pallette)+theme(legend.position = "none")+
  ylab("-log10(p)")+xlab("log2(fold-change)")

#p1<-ggplot(data = dflm1, aes(y = -log10(pval), x = logfc,color=col )) +
  #geom_point(alpha = 1, size = 2)+scale_colour_manual(values=palete)+theme(legend.position = "none")+
  #ylab("-log10(p)")+xlab("log2(fold-change)")

#p2<-ggplot(data = dflm2, aes(y = -log10(pval), x = logfc,color=col )) +
 #geom_point(alpha = 1, size = 2)+scale_colour_manual(values=palete)+theme(legend.position = "none")+
  #ylab("-log10(p)")+xlab("log2(fold-change)")

#p3<-ggplot(data = dflm3, aes(y = -log10(pval), x = logfc,color=col )) +
  #geom_point(alpha = 1, size = 2)+scale_colour_manual(values=palete)+theme(legend.position = "none")+
  #ylab("-log10(p)")+xlab("log2(fold-change)")


#grid.arrange(p1, p2, p3, ncol=3)

plot(pvs2)
```

``` r
library(edgeR)
```

    ## Warning: package 'edgeR' was built under R version 4.2.1

    ## Loading required package: limma

    ## Warning: package 'limma' was built under R version 4.2.1

``` r
library(ggplot2)

counts <- read.table("/Volumes/Data/Projects/dmelR2_p-ele/rna/run2/map-bwamem/genome/dmel_fc_counts_bwa.txt", header = TRUE, skip = 1, stringsAsFactors = FALSE)
counts <- counts[, c("Geneid", "dmel_R1G30_run2.sort.bam", "dmel_R2G30_run2.sort.bam", "dmel_R3G30_run2.sort.bam", "dmel_R1G40_run2.sort.bam", "dmel_R2G40_run2.sort.bam", "dmel_R3G40_run2.sort.bam")]

# Prepare counts matrix
rownames(counts) <- counts$Geneid
counts <- counts[, -1]  # Remove the Geneid column

# Create sample information matrix
samples <- colnames(counts)
group <- c(rep("R1", 2), rep("R2", 2), rep("R3", 2))
subgroup <- c("G30", "G40")
group <- factor(paste(group, subgroup, sep = "_"), levels = unique(paste(group, subgroup, sep = "_")))

# Create separate DGEList objects
dge_R1_R2 <- DGEList(counts = counts[, c("dmel_R1G30_run2.sort.bam", "dmel_R1G40_run2.sort.bam", "dmel_R2G30_run2.sort.bam", "dmel_R2G40_run2.sort.bam")], group = group[c(1, 4, 2, 5)])
dge_R2_R3 <- DGEList(counts = counts[, c("dmel_R2G30_run2.sort.bam", "dmel_R2G40_run2.sort.bam", "dmel_R3G30_run2.sort.bam", "dmel_R3G40_run2.sort.bam")], group = group[c(2, 5, 3, 6)])
dge_R1_R3 <- DGEList(counts = counts[, c("dmel_R1G30_run2.sort.bam", "dmel_R1G40_run2.sort.bam", "dmel_R3G30_run2.sort.bam", "dmel_R3G40_run2.sort.bam")], group = group[c(1, 4, 3, 6)])

dge_R1_R2 <- calcNormFactors(dge_R1_R2)
dge_R2_R3 <- calcNormFactors(dge_R2_R3)
dge_R1_R3 <- calcNormFactors(dge_R1_R3)

dge_R1_R2 <- estimateGLMCommonDisp(dge_R1_R2)
dge_R2_R3 <- estimateGLMCommonDisp(dge_R2_R3)
dge_R1_R3 <- estimateGLMCommonDisp(dge_R1_R3)

dge_R1_R2 <- estimateGLMTrendedDisp(dge_R1_R2)
dge_R2_R3 <- estimateGLMTrendedDisp(dge_R2_R3)
dge_R1_R3 <- estimateGLMTrendedDisp(dge_R1_R3)

dge_R1_R2 <- estimateGLMTagwiseDisp(dge_R1_R2)
dge_R2_R3 <- estimateGLMTagwiseDisp(dge_R2_R3)
dge_R1_R3 <- estimateGLMTagwiseDisp(dge_R1_R3)

fit_R1_R2 <- glmFit(dge_R1_R2)
fit_R2_R3 <- glmFit(dge_R2_R3)
fit_R1_R3 <- glmFit(dge_R1_R3)

lrt_R1_R2 <- glmLRT(fit_R1_R2)
lrt_R2_R3 <- glmLRT(fit_R2_R3)
lrt_R1_R3 <- glmLRT(fit_R1_R3)

# Extract the differential expression results and create volcano plots
de_results_R1_R2 <- topTags(lrt_R1_R2, n = Inf)
de_table_R1_R2 <- as.data.frame(de_results_R1_R2$table)
de_table_R1_R2$logFC <- -de_table_R1_R2$logFC  # Flip the logFC column for better visualization
de_table_R1_R2$FDR <- p.adjust(de_table_R1_R2$PValue, method = "BH")
volcano_plot_R1_R2 <- ggplot(de_table_R1_R2, aes(x = logFC, y = -log10(PValue))) +
  geom_point(aes(color = FDR < 0.2), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red"), labels = c("FDR >= 0.2", "FDR < 0.2")) +
  theme_minimal() +
  labs(x = "log2 Fold Change", y = "-log10 p-value", title = "R1 vs R2") +
  theme(legend.position = "bottom")

de_results_R2_R3 <- topTags(lrt_R2_R3, n = Inf)
de_table_R2_R3 <- as.data.frame(de_results_R2_R3$table)
de_table_R2_R3$logFC <- -de_table_R2_R3$logFC
de_table_R2_R3$FDR <- p.adjust(de_table_R2_R3$PValue, method = "BH")
volcano_plot_R2_R3 <- ggplot(de_table_R2_R3, aes(x = logFC, y = -log10(PValue))) +
  geom_point(aes(color = FDR < 0.2), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red"), labels = c("FDR >= 0.2", "FDR < 0.2")) +
  theme_minimal() +
  labs(x = "log2 Fold Change", y = "-log10 p-value", title = "R2 vs R3") +
  theme(legend.position = "bottom")

de_results_R1_R3 <- topTags(lrt_R1_R3, n = Inf)
de_table_R1_R3 <- as.data.frame(de_results_R1_R3$table)
de_table_R1_R3$logFC <- -de_table_R1_R3$logFC
de_table_R1_R3$FDR <- p.adjust(de_table_R1_R3$PValue, method = "BH")
volcano_plot_R1_R3 <- ggplot(de_table_R1_R3, aes(x = logFC, y = -log10(PValue))) +
  geom_point(aes(color = FDR < 0.2), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red"), labels = c("FDR >= 0.8", "FDR < 0.8")) +
  theme_minimal() +
  labs(x = "log2 Fold Change", y = "-log10 p-value", title = "R1 vs R3") +
  theme(legend.position = "bottom")

print(volcano_plot_R1_R2)
```

![](5_EgdeR_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
print(volcano_plot_R2_R3)
```

![](5_EgdeR_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
print(volcano_plot_R1_R3)
```

![](5_EgdeR_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->
