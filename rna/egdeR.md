Edge-R
================
Matthew Beaumont
2023-07-11

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
#library(gridExtra)
theme_set(theme_bw())
palette <- c("darkgrey","orange","red")
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
  #fit <- glmFit(y,design)
  #lrt <- glmLRT(fit,coef=2)
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

lm1<-edgeme(nR1,nR2,nR3,eg30R1,eg40R1)
lm2<-edgeme(nR1,nR2,nR3,eg30R2,eg40R2)
lm3<-edgeme(nR1,nR2,nR3,eg30R3,eg40R3)
dflm1<-data.frame(gene=rownames(lm1$table),logfc=lm1$table$logFC, pval=lm1$table$PValue, fdr=p.adjust(lm1$table$PValue,method="BH"))
dflm2<-data.frame(gene=rownames(lm2$table),logfc=lm2$table$logFC, pval=lm2$table$PValue, fdr=p.adjust(lm2$table$PValue,method="BH"))
dflm3<-data.frame(gene=rownames(lm3$table),logfc=lm3$table$logFC, pval=lm3$table$PValue, fdr=p.adjust(lm3$table$PValue,method="BH"))

lmvs2<-edgemeext(eg30R1,eg40R1,eg30R2,eg40R2,eg30R3,eg40R3)
dflmvs2<-data.frame(gene=rownames(lmvs2$table),logfc=lmvs2$table$logFC, pval=lmvs2$table$PValue,fdr=p.adjust(as.numeric(lmvs2$table$PValue),method="BH"))

#### Colouring FDR 
fdr<-0.01
fdr2<-0.001

dflm1$col<-0
dflm1[dflm1$fdr<fdr,]$col<-1
dflm1[dflm1$fdr<fdr2,]$col<-2
dflm1$col<-as.factor(dflm1$col)

dflm2$col<-0
dflm2[dflm2$fdr<fdr,]$col<-1
dflm2[dflm2$fdr<fdr2,]$col<-2
dflm2$col<-as.factor(dflm2$col)

dflm3$col<-0
dflm3[dflm3$fdr<fdr,]$col<-1
dflm3[dflm3$fdr<fdr2,]$col<-2
dflm3$col<-as.factor(dflm3$col)

dflmvs2$col<-0
dflmvs2[dflmvs2$fdr<fdr,]$col<-1
#dflmvs2[dflmvs2$fdr<fdr2,]$col<-2
dflmvs2$col<-as.factor(dflmvs2$col)

pvs2<-ggplot(data = dflmvs2, aes(y = -log10(pval), x = logfc,color=col )) +
  geom_point(alpha = 1, size = 2)+scale_colour_manual(values=palete)+theme(legend.position = "none")+
  ylab("-log10(p)")+xlab("log2(fold-change)")

p1<-ggplot(data = dflm1, aes(y = -log10(pval), x = logfc,color=col )) +
  geom_point(alpha = 1, size = 2)+scale_colour_manual(values=palete)+theme(legend.position = "none")+
  ylab("-log10(p)")+xlab("log2(fold-change)")

p2<-ggplot(data = dflm2, aes(y = -log10(pval), x = logfc,color=col )) +
  geom_point(alpha = 1, size = 2)+scale_colour_manual(values=palete)+theme(legend.position = "none")+
  ylab("-log10(p)")+xlab("log2(fold-change)")

p3<-ggplot(data = dflm3, aes(y = -log10(pval), x = logfc,color=col )) +
  geom_point(alpha = 1, size = 2)+scale_colour_manual(values=palete)+theme(legend.position = "none")+
  ylab("-log10(p)")+xlab("log2(fold-change)")

 
#pdf(file="",width=8,height=3)

#grid.arrange(p1, p2, p4, ncol=3)
#dev.off()

#pdf(file="",width=3,height=3)
plot(pvs2)
#dev.off()
```
