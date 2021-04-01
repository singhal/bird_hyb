library(ape)
library(phytools)
library(rgdal)
library(rgeos)
library(sf)
library(raster)
library(pbapply)
library(hadsstr)

setwd("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/")
t = read.tree("data/T400F_AOS_HowardMoore.tre")

shps = list.files("data/range_maps_shape_AOSHM/",
                   recursive = T, full.names = T, 
                   pattern = ".shp")
sps = gsub(".*\\/", "", shps)
sps = gsub(".shp", "", sps)
ranges = lapply(shps, readOGR)
 
names(ranges) = sps
saveRDS(ranges, file = "data/range_maps_shape_AOSHM.Rds")

# ranges = readRDS("data/range_maps_shape_AOSHM.Rds")

# load raster file to get the right projection
r = raster("data/spatial_data/CHELSA_bio10_01.tif")
ranges2 = pblapply(ranges, spTransform, crs(r))
names(ranges2) = sps
saveRDS(ranges2, file = "data/range_maps_shape_AOSHM.projected.Rds")

ranges2 = readRDS("data/range_maps_shape_AOSHM.projected.Rds")

get_lat <- function(range) {
  lonlat = gCentroid(range)
  return(lonlat@coords[2])
}
lats1 = pblapply(ranges, get_lat)
lats2 = pblapply(ranges2, get_lat)

# now do ancestral reconstruction
# the two latitudes are very highly correlated
# so not that concerned about projection issues
cor.test(unlist(lats1), unlist(lats2))
# we will do actual latitude and take absolute later
lats = unlist(lats2)

# downsample tree 
# to match lats
to_drop = t$tip.label[! (t$tip.label %in% names(lats)) ]
t1 = drop.tip(t, to_drop)
lat_recon = fastAnc(t1, lats[t1$tip.label])
saveRDS(lat_recon, "data/latitude_reconstruction.Rds")
saveRDS(t1, "data/latitude_reconstruction_tree.Rds")
