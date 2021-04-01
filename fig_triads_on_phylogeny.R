library(ape)
library(phytools)

setwd("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/")
t = read.tree("data/T400F_AOS_HowardMoore.tre")

# not all triads identified
# could be included
all = read.csv("metadata/triads-2020-07-07.csv",
               stringsAsFactors = F)
used = read.csv("results/dstat_af_subsample250.csv",
                stringsAsFactors = F)
all$used = ifelse(all$sp1 %in% used$sp1, "#fc8d62", "#8da0cb")

edge.cols = rep("gray60", length(t$edge.length))
edge.widths = rep(0.1, length(t$edge.length))
for (i in 1:nrow(all)) {
  sps = as.character(all[i, 1:3])
  node = findMRCA(t, sps)
  desc = getDescendants(t, node)
  edge.cols[ which(t$edge[,2] %in% desc) ] = all[i, "used"]
  edge.widths[ which(t$edge[,2] %in% desc) ] = 1
}  
 
pdf("figures/triads_on_phylogeny.pdf", height = 6, width = 6)
par(mar = c(0, 0, 0, 0))
plot.phylo(t, show.tip.label = F, edge.color = edge.cols,
           edge.width = edge.widths, label.offset = 1,
           type = "fan")
dev.off()
