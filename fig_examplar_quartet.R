setwd("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/")
getwd()

library(ape)
library(rgdal)
library(maptools)
library(raster)

# Get subtree
tree <- read.tree('./data/T400F_AOS_HowardMoore.tre')
targets <- c("Pipreo_rief_IAvHBTAMC910", "Pipreo_inter_L7982", "Pipreo_arcu_L31988", "Pipreo_pulhra_L8055")
subtree <- drop.tip(tree, setdiff(tree$tip.label, targets))
pdf("./figures/examplar_quartet_subtree.pdf", width=5, height=5)
plot(ladderize(subtree), direction="downwards", label.offset=0.4, edge.width=2)
dev.off()

# Pull discordance fraction
phyparts.trees <- read.tree("./data/phyparts/out_filter/out80.concon.tre")
con.tree <- phyparts.trees[1][[1]]
disc.tree <- phyparts.trees[2][[1]]
con.subtree <- drop.tip(con.tree, setdiff(con.tree$tip.label, targets))
disc.subtree <- drop.tip(disc.tree, setdiff(disc.tree$tip.label, targets))
plot(con.subtree, show.node.label=TRUE)
plot(disc.subtree, show.node.label=TRUE)

# Loop to get range maps
ranges <- vector()
for(i in 1:length(targets)) {
	print(targets[i])
	range.file <- paste0('./data/range_maps_shape_AOSHM/', targets[i], '/', targets[i], '.shp')
	range <- readOGR(dsn=range.file, layer=targets[i], verbose=FALSE)	
	range.dissolved <- try(unionSpatialPolygons(range, rep(as.character(i), length(range@polygons))))
	if(class(range.dissolved) == "try-error") {
		print(paste("Range error for: ", names(targets[i])))
	} else {
		ranges <- c(ranges, range)
	}
}

# Get base map
alt <- raster('/Users/eebuser/Documents/Harvey/research/seabirds/sampling_labwork/map/Vert Lunch_ R Maps Tutorial/alt.NewWorld.30arcsec.grd') # Too large to put in Dropbox
projection(alt)
options(scipen=999)
ex <- c(-83, -59, -24, 13)
alt.crop <- crop(alt, extent(ex))  
mapcol = gray.colors(50, start = 0.9, end = 0, alpha = 1)
rivers <- readShapeSpatial('./data/other_map_data/ne_50m_rivers_lake_centerlines/ne_50m_rivers_lake_centerlines.shp')
coastline <- readShapeSpatial('./data/other_map_data/ne_50m_coastline/ne_50m_coastline.shp')

pdf("./figures/examplar_quartet_maps.pdf", width=8, height=5)
par(mfrow=c(1,4))
for(i in 1:length(ranges)) {
	plot(alt.crop, col = mapcol, axes = FALSE, box = FALSE, legend = FALSE)
	plot(rivers, add=T, lwd=0.5)
	plot(coastline, add=T, lwd=0.5)
	range <- spTransform(ranges[i][[1]], CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
	plot(range, add=T, col=alpha("red", 0.75), border=NA)
}
dev.off()
