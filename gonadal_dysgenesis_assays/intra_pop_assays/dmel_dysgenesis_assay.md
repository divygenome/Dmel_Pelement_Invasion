D. melanogaster intrapopulation hybrid dysgenesis assay
================
Divya Selvaraju
2024-12-12

## Generate hybrid dysgenesis plot for D. melanogaster

``` r
library(ggplot2)
theme_set(theme_bw())

dys_dm=read.table("/home/divya/Desktop/invasion/mel/Dmel_3R_3TP.tsv",header=TRUE)
print(dys_dm)
```

    ##    Species Replicate Generation Doer Dysgenic Intermediate Normal Total
    ## 1     Dmel        R1         G1   Di        4            3     63    70
    ## 2     Dmel        R2         G1   Di        0            6     64    70
    ## 3     Dmel        R3         G1   Ro        2            0     65    67
    ## 4     Dmel        R1         G5   Di        0            0     82    82
    ## 5     Dmel        R2         G5   Di        3            0     69    72
    ## 6     Dmel        R3         G5   Di        0            0     70    70
    ## 7     Dmel        R1        G10   Di       13            6     72    91
    ## 8     Dmel        R2        G10   Di        7            3     79    89
    ## 9     Dmel        R3        G10   Ro        3            0     53    56
    ## 10    Dmel        R1        G15   Di        5            4     61    70
    ## 11    Dmel        R2        G15   Di        2            1     67    70
    ## 12    Dmel        R3        G15   Di       10            3     58    71
    ## 13    Dmel        R1        G20   Di       41            4     32    77
    ## 14    Dmel        R2        G20   Di        3            4     63    70
    ## 15    Dmel        R3        G20   Di       20            3     60    83
    ## 16    Dmel        R1        G25   Di       75            0      3    78
    ## 17    Dmel        R2        G25   Di       10            3     80    93
    ## 18    Dmel        R3        G25   Di       69            2     18    89
    ## 19    Dmel        R1        G30   Di       38            1     42    81
    ## 20    Dmel        R2        G30   Di        6            1     64    71
    ## 21    Dmel        R3        G30   Di       41            0     32    73
    ## 22    Dmel        R1        G34   Di        8            4     63    75
    ## 23    Dmel        R2        G34   Di        8            0     72    80
    ## 24    Dmel        R3        G34   Di       29            0     45    74
    ## 25    Dmel        R1        G41   Di        0            3     78    81
    ## 26    Dmel        R2        G41   Di        3            7     87    97
    ## 27    Dmel        R3        G41   Di       19            8     95   122
    ## 28    Dmel        R1        G45   Di        0            2     78    80
    ## 29    Dmel        R2        G45   Di        3            1     94    98
    ## 30    Dmel        R3        G45   Di        5            4    102   111
    ## 31    Dmel        R1        G50   Di        4            8     88   100
    ## 32    Dmel        R2        G50   Di        4            2     97   103
    ## 33    Dmel        R3        G50   Di        0            3     98   101
    ##    Dys_percentile
    ## 1           7.86%
    ## 2           4.29%
    ## 3           2.99%
    ## 4           0.00%
    ## 5           4.17%
    ## 6           0.00%
    ## 7          17.58%
    ## 8           9.55%
    ## 9           5.36%
    ## 10         10.00%
    ## 11          3.57%
    ## 12         16.20%
    ## 13         55.84%
    ## 14          7.14%
    ## 15         25.90%
    ## 16         96.15%
    ## 17         12.37%
    ## 18         78.65%
    ## 19         47.53%
    ## 20          9.15%
    ## 21         56.16%
    ## 22         13.33%
    ## 23         10.00%
    ## 24         39.19%
    ## 25          1.85%
    ## 26          6.70%
    ## 27         18.85%
    ## 28          1.25%
    ## 29          3.57%
    ## 30          6.31%
    ## 31          8.00%
    ## 32          4.85%
    ## 33          1.49%

``` r
dys_dmel=data.frame(dys_dm$Replicate, dys_dm$Generation,dys_dm$Dys_percentile)
names(dys_dmel)[1]<-"replicate"
names(dys_dmel)[2]<-"generation"
names(dys_dmel)[3]<-"dys_percent"
dys_dmel$dys_percent=as.numeric(gsub("%", "", dys_dmel$dys_percent))
dys_dmel$generation=as.numeric(gsub("G", "", dys_dmel$generation))

print(dys_dmel)
```

    ##    replicate generation dys_percent
    ## 1         R1          1        7.86
    ## 2         R2          1        4.29
    ## 3         R3          1        2.99
    ## 4         R1          5        0.00
    ## 5         R2          5        4.17
    ## 6         R3          5        0.00
    ## 7         R1         10       17.58
    ## 8         R2         10        9.55
    ## 9         R3         10        5.36
    ## 10        R1         15       10.00
    ## 11        R2         15        3.57
    ## 12        R3         15       16.20
    ## 13        R1         20       55.84
    ## 14        R2         20        7.14
    ## 15        R3         20       25.90
    ## 16        R1         25       96.15
    ## 17        R2         25       12.37
    ## 18        R3         25       78.65
    ## 19        R1         30       47.53
    ## 20        R2         30        9.15
    ## 21        R3         30       56.16
    ## 22        R1         34       13.33
    ## 23        R2         34       10.00
    ## 24        R3         34       39.19
    ## 25        R1         41        1.85
    ## 26        R2         41        6.70
    ## 27        R3         41       18.85
    ## 28        R1         45        1.25
    ## 29        R2         45        3.57
    ## 30        R3         45        6.31
    ## 31        R1         50        8.00
    ## 32        R2         50        4.85
    ## 33        R3         50        1.49

``` r
p1=ggplot(dys_dmel,aes(x=generation,y=dys_percent))+geom_line(aes(color=replicate))+geom_point(aes(color=replicate))+ ggtitle("dmel_dysgenesis")+ xlab("generations")+ylab("fraction of dysgenesis") + scale_color_manual(values=c ("firebrick", "skyblue3","chartreuse4"))

plot(p1)
```

![](dmel_dysgenesis_assay_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
## Generate Plots for manuscript

png("dmel_dys.png", width = 30, height = 20, units = "cm", res = 100)
print(p1)
dev.off()
```

    ## png 
    ##   2

``` r
postscript("dmel_dys.eps", family = "ArialMT", width = 30, height = 20)
print(p1)
dev.off()
```

    ## png 
    ##   2
