library(ggplot2)
library(ape)
library(cowplot)
library(nlme)
library(geiger)
library(phytools)
theme_set(theme_cowplot())

d = read.csv("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/results/dstat_af_subsample250.expanded.csv",
             stringsAsFactors = F)
xx = ggplot(d) + 
  geom_point(aes(abs(latitude),log(abs(hyb_stability)))) +
  xlab("absolute Latitude") +
  ylab("Log absolute climate instability")
save_plot("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/figures/lat_instability.png", xx)

setwd("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/")
t = read.tree("data/T400F_AOS_HowardMoore.tre")
d = d[complete.cases(d$hyb_stability), ]
t1 = keep.tip(t, d$sp1)
d1 = d[match(t1$tip.label, d$sp1), ]

bm = corBrownian(1, t1)
xx= gls(log(abs(hyb_stability)) ~ abs(latitude), data = d1, correlation = bm)
summary(xx)

d1$stab = log(abs(d1$all_stability))
d1$lat = abs(d1$latitude)
obj = phyl.vcv(as.matrix(d1[, c("stab", "lat")]), vcv(t1), 1)
corvar = cov2cor(obj$R)[1, 2]

