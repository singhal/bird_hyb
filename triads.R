library(ape)
library(phytools)
library(phangorn)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
theme_set(theme_cowplot())

mdir = "/Users/Sonal/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/"
setwd(mdir)
loccts = read.csv("metadata/sample_contig_counts.csv", stringsAsFactors = F)

# load tree
t = read.tree("data/T400F_AOS_HowardMoore.tre")

# climb through all nodes
nstart = Ntip(t) + 1
nodes = seq(nstart, nstart + Nnode(t))
desc = lapply(nodes, function(x) {getDescendants(t, x)})
desctips = lapply(desc, function(x) {x[x <= Ntip(t)]})
cts = lapply(desctips, length)

# what are our triads?
triads = desctips[which(cts == 3)]

# get outgroups
ntriads = nodes[which(cts == 3)]
parent = lapply(ntriads, function(x) {Ancestors(t, x, type = "parent")})
parentdesc = lapply(parent, function(x) {getDescendants(t, x)})
outs = rep(NA, length(ntriads))
for (i in 1:length(parentdesc)) {
  trio = triads[[i]]
  alldesc = parentdesc[[i]]
  possout = alldesc[!(alldesc %in% trio)]
  possout = possout[possout <= Ntip(t)]
  possout = t$tip.label[possout]
  possdf = loccts[which(loccts$sample %in% possout), ] %>% arrange(-contig)
  outs[i] = possdf$sample[1]
}

nametriads = lapply(ntriads, function(x) {read.tree(text = 
                                                     write.tree(ladderize(extract.clade(t,
                                                                                  x), right = T)))$tip.label})
df = as.data.frame(do.call(rbind, nametriads))
names(df) = c("sp1", "sp2", "sp3")
df$outgroup = outs 
write.csv(df, paste0("triads-", Sys.Date(), ".csv"), row.names = F)

# get MRCA
nodenums = nodes[which(cts == 3)]
mrcas = max(nodeHeights(t)) - 
  unlist(lapply(nodenums, function(x) {nodeheight(t, x)}))
  
# get internode length
## get the internal node
internodes1 = lapply(nodenums, function(x) {getDescendants(t, x)})
internodes2 = unlist(lapply(internodes1, function(x) {x[x > Ntip(t)]}))
intermrcas = max(nodeHeights(t)) - 
  unlist(lapply(internodes2, function(x) {nodeheight(t, x)}))
interlengths = mrcas - intermrcas

# get distance to outgroup
outnodes = unlist(lapply(nodenums, 
                  function(x) {Ancestors(t, x, type = "parent")}))
outmrcas = max(nodeHeights(t)) - 
  unlist(lapply(outnodes, function(x) {nodeheight(t, x)}))

# summarize all
d = data.frame(df, nodenums, mrcas, interlengths, outmrcas)
used = read.csv("results/dstat_af_subsample250.csv",
                stringsAsFactors = F)
d$used = ifelse(d$sp1 %in% used$sp1, TRUE, FALSE)
                
              #   "#fc8d62", "#8da0cb")

a = ggplot(d, aes(mrcas, fill = used)) + 
  geom_histogram(bins = 30) + 
  xlab("MRCA of triad (myr)") +
  scale_fill_manual(values = c("#8da0cb", "#fc8d62"))
b = ggplot(d, aes(interlengths, fill = used)) + 
  geom_histogram(bins = 30) + 
  xlab("internode length of triad (myr)") +
  scale_fill_manual(values = c("#8da0cb", "#fc8d62"))
c = ggplot(d, aes(outmrcas, fill = used)) + 
  geom_histogram(bins = 30) + 
  xlab("MRCA of triad + outgroup (myr)") +
  scale_fill_manual(values = c("#8da0cb", "#fc8d62"))

prow <- plot_grid(
  b + theme(legend.position="none"),
  a + theme(legend.position="none"),
  c + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C"),
  hjust = -1,
  nrow = 1
)
save_plot("figures/triad_descriptions.pdf", 
          prow, ncol = 3, base_width = 4, 
          base_height = 3)

g = ggplot(d, aes(mrcas, interlengths)) +
  geom_point() + xlab("MRCA of triad") +
  ylab("internode length of triad")
e = ggplot(d, aes(mrcas, outmrcas)) +
  geom_point() + xlab("MRCA of triad") +
  ylab("MRCA of triad + outgroup")
f = ggplot(d, aes(interlengths, outmrcas)) +
  geom_point() + xlab("internode length of triad") +
  ylab("MRCA of triad + outgroup")
abc = (g | e | f)
save_plot("~/Desktop/triad_description_correlations.pdf", 
          abc, ncol = 3)