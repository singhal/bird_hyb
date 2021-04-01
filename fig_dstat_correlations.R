library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

setwd("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/results/")
d1 = read.csv("dstat_af_subsample100.csv", stringsAsFactors = F)
d2 = read.csv("dstat_af_subsample250.csv", stringsAsFactors = F)
d3 = read.csv("dstat_af_subsample500.csv", stringsAsFactors = F)

d4 = left_join(d1, d2, by = c("sp1" = "sp1",
                              "sp2" = "sp2",
                              "sp3" = "sp3",
                              "out" = "out"))
d5 = left_join(d4, d3, by = c("sp1" = "sp1",
                              "sp2" = "sp2",
                              "sp3" = "sp3",
                              "out" = "out"))

a = ggplot(d5) + geom_point(aes(D.y, D.x), alpha = 0.6) +
  ylab('D-statistic, 100 SNPs') +
  xlab('D-statistic, 250 SNPs') +
  geom_abline(slope = 1, intercept = 0, col = "navyblue")
b = ggplot(d5 %>% filter(complete.cases(D))) + 
  geom_point(aes(D.y, D), alpha = 0.6) +
  ylab('D-statistic, 500 SNPs') +
  xlab('D-statistic, 250 SNPs') +
  geom_abline(slope = 1, intercept = 0, col = "navyblue")
ab = plot_grid(a, b, labels = c("A", "B"))
save_plot("../figures/dstat_subsample.png", ab,
          ncol = 2, base_height = 3, base_width = 4)

cor.test(d5$D.x, d5$D.y)

cor.test(d5$D.y, d5$D)
