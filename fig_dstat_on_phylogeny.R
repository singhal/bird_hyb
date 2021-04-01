library(ape)
library(phytools)

setwd("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/")

t = read.tree("data/T400F_AOS_HowardMoore.tre")
d = read.csv("results/dstat_af_subsample250.csv",
                stringsAsFactors = F)

t2 = keep.tip(t, d$sp1)
d1 = d[match(t2$tip.label, d$sp1), ]
d1$sig = ifelse(d1$p < 0.05, "#fc8d62", "gray40")

edgecols = rep("gray40", length(t2$edge.length))
red = match(d1[d1$sig == "#fc8d62", "sp1"], t2$tip.label)
edgecols[ which(t2$edge[,2] %in% red) ] = "#fc8d62"

pdf("figures/Dstat_on_phylogeny.pdf", height = 5.5, width = 4)
layout(matrix(c(1,2), nrow = 1), widths = c(0.5, 0.5))
par(mar = c(2, 0, 0, 0))
plot.phylo(t2, show.tip.label = F, edge.color = edgecols)
plot.new()
plot.window(xlim = c(-1, 1), ylim = c(1, Ntip(t2)))
points(d1$D, 1:nrow(d1), pch = 21, bg = d1$sig)
ticks = c(-1, -0.5, 0, 0.5, 1)
laxis = -0.2
axis(1, ticks, tck=-0.02, labels=NA, line = laxis)
mtext(side = 1, l = laxis + 0.1, text = ticks, 
      at = ticks, las = 1, cex = 0.7)
mtext("D-statistic", 
      side = 1, line = laxis + 1, cex = 0.8)
dev.off()

