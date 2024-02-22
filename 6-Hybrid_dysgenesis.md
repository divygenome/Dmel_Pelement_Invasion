6-Hybrid_dysgenesis
================
Matthew Beaumont
2024-02-22

``` bash
knitr::opts_chunk$set(echo = TRUE)
```

We wanted to assess the dysgenic inducibility of the different invaded
replicates. To do so, we set up reciprocal crosses of each replicate
with the original P-element naive strain, DM68. We also used a
reciprocal cross between DM68 and the Harwich strain, known to carry the
P-element, as a control.

``` bash
less hybrid_dysgenesis/dmel_HDassay.tsv
```

    ## line replicate   dysgenic    intermediate    non-dysgenic    total_ovaries   percentage_dysgenic
    ## R1x68    1   13  5   32  100 31
    ## R1x68    2   12  2   36  100 26
    ## R1x68    3   10  1   39  100 21
    ## R2x68    1   26  4   20  100 56
    ## R2x68    2   30  2   18  100 62
    ## R2x68    3   18  11  21  100 47
    ## R3x68    1   10  8   32  100 28
    ## R3x68    2   14  2   34  100 30
    ## R3x68    3   16  4   30  100 36
    ## 68xR1    1   50  0   0   100 100
    ## 68xR1    2   48  0   2   100 96
    ## 68xR1    3   47  3   0   100 97
    ## 68xR2    1   18  8   24  100 44
    ## 68xR2    2   33  4   13  100 70
    ## 68xR2    3   31  2   17  100 64
    ## 68xR3    1   44  1   5   100 89
    ## 68xR3    2   46  2   2   100 94
    ## 68xR3    3   41  4   5   100 86
    ## Harx68   1   50  0   0   100 100
    ## Harx68   2   50  0   0   100 100
    ## Harx68   3   50  0   0   100 100
    ## 68xHar   1   0   0   50  100 0
    ## 68xHar   2   0   0   50  100 0
    ## 68xHar   3   0   0   50  100 0

``` r
library(ggplot2)

data <- read.table("hybrid_dysgenesis/dmel_HDassay.tsv", header = TRUE)

data$group <- NA  
data$replicate <- NA  

data$group[data$line %in% c("68xR1", "68xR2", "68xR3", "68xHar")] <- "Maternal"
data$group[data$line %in% c("R1x68", "R2x68", "R3x68", "Harx68")] <- "Paternal"

data$replicate[data$line %in% c("68xR1", "R1x68")] <- "R1"
data$replicate[data$line %in% c("68xR2", "R2x68")] <- "R2"
data$replicate[data$line %in% c("68xR3", "R3x68")] <- "R3"
data$replicate[data$line %in% c("Harx68", "68xHar")] <- "control"

# Order of replicates
replicate_order <- c("R1", "R2", "R3", "control")

# Order of groups
#group_order <- c("Maternal", "Paternal")

# Create a grouped boxplot with manually nudged points
HD <- ggplot(data, aes(x = replicate, y = percentage_dysgenic, fill = group)) +
        geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8, color = "black") +
        #geom_jitter(position = position_nudge(x = 0.2), size = 2, alpha = 0.7) +
        labs(y = "ovarian dysgenesis (%)", x = NULL) +
        scale_fill_manual(values = c("Maternal" = "steelblue", "Paternal" = "coral")) +
        theme_minimal() +
        theme(legend.position = "right",
              legend.title = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, size = 0.7)) +
        ylim(0, 100) +
        scale_x_discrete(limits = replicate_order)
```

    ## Warning: The `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
ggsave("figs/HD.png", HD, width = 6, height = 4, dpi = 600)

knitr::include_graphics("figs/HD.png")
```

<img src="figs/HD.png" width="3600" />
