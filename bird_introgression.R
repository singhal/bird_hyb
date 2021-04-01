library(ape)
library(corrplot)
library(cowplot)
library(dplyr)
library(ggplot2)
library(phangorn)
library(phytools)
library(rgdal)
library(rgeos)
library(sf)
library(tidyr)
theme_set(theme_cowplot())

##############################
# define the data set
##############################

# load data
setwd("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/")
d = read.csv("results/dstat_af_subsample250.csv", stringsAsFactors = F)

# define the hybridizing species
d$hyb_sp1 = d$sp3
# positive D - introgression btn 2 & 3
# negative D - introgression btn 1 & 3
d$hyb_sp2 = ifelse(d$D > 0, d$sp2, d$sp1)

##############################
# read in the div time data
##############################

t = read.tree("data/T400F_AOS_HowardMoore.tre")
times = branching.times(t)

# divtime btn 2
get_time2 <- function(row, tree) {
  return(getMRCA(tree, as.character(row[c("sp1", "sp2")])))
}
nodes2 = apply(d, 1, get_time2, t)
d$divtime2 = times[nodes2 - Ntip(t)]

# divtime btn 3
# this will also be the divergence 
# time between hybridizing pairs
# so perhaps the most relevant?
get_time3 <- function(row, tree) {
  return(getMRCA(tree, as.character(row[c("sp1", "sp2", "sp3")])))
}
nodes3 = apply(d, 1, get_time3, t)
d$divtime3 = times[nodes3 - Ntip(t)]

# internode
d$internode = d$divtime3 - d$divtime2

# len to outgroup
get_time4 <- function(row, tree) {
  return(getMRCA(tree, as.character(row[c("sp1", "sp2", "sp3", "out")])))
}
nodes4 = apply(d, 1, get_time4, t)
d$length_to_out = times[nodes4 - Ntip(t)] - d$divtime3

##############################
# read in the ILS data
##############################

ils = read.tree("data/phyparts/out_filter/out80.concon.tre")
# conc
con = ils[[1]]
con$node.label = as.numeric(con$node.label)
cf = ils[[2]]
cf$node.label = as.numeric(cf$node.label)

# a tree with levels of conflict
cf$node.label = (cf$node.label) / (con$node.label + cf$node.label)
mean(cf$node.label, na.rm = T)
cfnodes = cf$node.label

# ILS at node 1
nodes2 = apply(d, 1, get_time2, cf)
d$ils2 = cfnodes[nodes2 - Ntip(t)]

# ILS at node 2
nodes3 = apply(d, 1, get_time3, cf)
d$ils3 = cfnodes[nodes3 - Ntip(t)]

# ILS at node 3
nodes4 = apply(d, 1, get_time4, cf)
d$ils4 = cfnodes[nodes4 - Ntip(t)]

##############################
# geography
##############################

# old way of defining latitude
# lat = readRDS("data/latitude_reconstruction.Rds")
# lattree = readRDS("data/latitude_reconstruction_tree.Rds")

geo = vector("list", nrow(d))

# skip row 62 bc
# it is has Machae_reglus
skip = grep("Machae_reglus", d$sp1)
rows = seq(1:nrow(d))
rows = rows[rows != skip]

for (i in rows) {
  r1 = paste0("data/range_maps_shape_AOSHM/",
              d[i, "sp1"], "/",
              d[i, "sp1"], ".shp")
  r2 = paste0("data/range_maps_shape_AOSHM/",
              d[i, "sp2"], "/",
              d[i, "sp2"], ".shp")
  r3 = paste0("data/range_maps_shape_AOSHM/",
              d[i, "sp3"], "/",
              d[i, "sp3"], ".shp")
  
  r1 = readOGR(r1)
  # r1 = spTransform(r1, projstring)
  r2 = readOGR(r2)
  # r2 = spTransform(r2, projstring)
  r3 = readOGR(r3)
  # r3 = spTransform(r3, projstring)
  ranges = list(r1, r2, r3)
  
  # range size (as proxy for pop size)
  rangesize = mean(unlist(lapply(ranges, gArea))) / 1e6
  
  # geographic overlap
  # btn hyb species
  if (d[i, "D"] < 0) {
    hybranges = list(r1, r3)
  } else {
    hybranges = list(r2, r3)
  }
  inter = gIntersection(hybranges[[1]], hybranges[[2]])
  minarea = min(unlist(lapply(hybranges, gArea)))
  if (is.null(inter)) {
    overlap = 0
  } else {
    overlap = gArea(inter) / minarea
  }
  
  # geographic distance
  # btn hyb species
  distance = gDistance(hybranges[[1]], hybranges[[2]])
  
  # average latitude
  hyb1 = spTransform(hybranges[[1]], 
                     CRS('+proj=longlat +datum=WGS84 +no_defs'))
  hyb2 = spTransform(hybranges[[2]],
                     CRS('+proj=longlat +datum=WGS84 +no_defs'))
  splat = mean(c(gCentroid(hyb1)@coords[2],
               gCentroid(hyb2)@coords[2]))
  
  # latitude from ancestral reconstruction
  # sps = as.character(d[i, 1:3])
  # sps = sps[sps %in% lattree$tip.label]
  # if (length(sps) > 1) {
  #  latnode = findMRCA(lattree, sps)
  #  splat = lat[as.character(latnode)]
  # } else {
  #   splat = NA
  # }
  
  geo[[i]] = c(rangesize, overlap, distance, splat)
}
geo[[skip]] = c(NA, NA, NA, NA)
geo2 = as.data.frame(do.call(rbind, geo))
names(geo2) = c("rangesize", "overlap", "distance", "latitude")

# geographic stability
s = read.csv("data/bird.climate_change.csv", stringsAsFactors = F)
s1 = spread(s, type, value)
velocity = s1$mean_vel
names(velocity) = s1$species

all_stability = rep(NA, nrow(d))
hyb_stability = rep(NA, nrow(d))
for (i in rows) {
  sps = as.character(d[i, 1:3])
  sps = sps[sps %in% names(velocity)]
  all_stability[i] = mean(velocity[sps], na.rm = T)
  
  sps = as.character(d[i, c("hyb_sp1", "hyb_sp2")])
  sps = sps[sps %in% names(velocity)]
  hyb_stability[i] = mean(velocity[sps], na.rm = T)
}

  
d2 = cbind(d, geo2, all_stability, hyb_stability)
write.csv(d2, "results/dstat_af_subsample250.expanded.csv",
          row.names = F, quote = F)

# divtime2, divtime3, internode, length_to_out, 
# ils2, ils3, ils4, D / Z
# rangesize, overlap, distance, latitude
# all_stability, hyb_stability

vars1 = c("divtime2", "divtime3", "internode",
         "length_to_out")
vars2 = c("distance", "overlap")
vars3 = c("all_stability", "hyb_stability", "latitude")
vars4 = c("ils2", "ils3", "ils4")
vars5 = c("D", "Z")

d3 = d2 %>% dplyr::select(vars1) %>% drop_na(any_of(vars1))
m1 = cor(d3)
corrplot.mixed(m1, upper='ellipse', lower="number",
               tl.col="black", tl.cex=0.5, tl.srt=45)
# divtime3 is very correlated with divtime2 and internode

d4 = d2 %>% dplyr::select(vars2) %>% drop_na(any_of(vars2))
m2 = cor(d4)
corrplot.mixed(m2, upper='ellipse', lower="number",
               tl.col="black", tl.cex=0.5, tl.srt=45)
# shouldn't drop anything

d5 = d2 %>% dplyr::select(vars3) %>% drop_na(any_of(vars3))
m3 = cor(d5)
corrplot.mixed(m3, upper='ellipse', lower="number",
               tl.col="black", tl.cex=0.5, tl.srt=45)
# should drop all_stability
cor.test(abs(d2$hyb_stability), abs(d2$latitude))
# as expected, higher latitudes show less stability

d6 = d2 %>% dplyr::select(vars4) %>% drop_na(any_of(vars4))
m4 = cor(d6)
corrplot.mixed(m4, upper='ellipse', lower="number",
               tl.col="black", tl.cex=0.5, tl.srt=45)
# shouldn't drop any of them 