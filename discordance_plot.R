setwd("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/")
getwd()

library(ape)
library(maptools)
library(raster)
library(rgdal)
library(rworldmap)
library(RColorBrewer)
library(colourschemes)

# Get data
phyparts.trees <- read.tree("./data/phyparts/out_filter/out80.concon.tre")
con.tree <- phyparts.trees[1][[1]]
disc.tree <- phyparts.trees[2][[1]]
time.tree.complete <- read.tree("./data/T400F_AOS_HowardMoore.tre")
time.tree <- drop.tip(time.tree.complete, setdiff(time.tree.complete$tip.label, con.tree$tip.label))


# Discordance tree plot

# Get node values
ntip <- Ntip(time.tree)
comp <- all.equal.phylo(con.tree, time.tree, use.edge.length=FALSE, index.return=TRUE)
node.map <- comp[,2] # con.tree nodes
names(node.map) <- as.character(comp[,1]) # time.tree nodes
internal.node.map <- as.numeric(node.map[as.character((ntip+1):((ntip*2)-1))])-ntip
prop.disc <- as.numeric(disc.tree$node.label[internal.node.map])/(as.numeric(disc.tree$node.label[internal.node.map])+as.numeric(con.tree$node.label[internal.node.map]))
min(prop.disc, na.rm=TRUE)
max(prop.disc, na.rm=TRUE)
mean(prop.disc, na.rm=TRUE)

prop.disc[1] <- 0 # Assign value of zero to root node
rs <- rampInterpolate(limits=c(0,1), ramp=c("cornflowerblue", "cornflowerblue", "cornflowerblue", "bisque", "tomato1"))
node.cols <- rs(prop.disc)
#colfunc <- colorRampPalette(c("purple4", "mediumorchid", "bisque", "darkorange"))
#node.cols <- colfunc(length(prop.disc))[cut(prop.disc, breaks=length(prop.disc))]

# Calculate tip values
rootnode <- length(time.tree$tip.label) + 1
tip.disc <- vector()
for (i in 1:length(time.tree$tip.label)){
	vals <- vector()
	weights <- vector()
	node <- i
	index <- 1
	while (node != rootnode){
		node <- time.tree$edge[time.tree$edge[,2] == node][1]
		vals <- c(vals, prop.disc[node-length(time.tree$tip.label)])
		weights <- c(weights, 1/index)
		index <- index*2
		
	}
	tip.disc <- c(tip.disc, weighted.mean(vals, weights))
}		
names(tip.disc) <- time.tree$tip.label
#write.table(tip.disc, "./data/phyparts/discordance_tipvals.txt")
min(tip.disc, na.rm=TRUE)
max(tip.disc, na.rm=TRUE)
mean(tip.disc, na.rm=TRUE)

# Add to tree plot
tip.cols <- rs(tip.disc)
# Radial version
pdf(file="./Figures/Discordance_tree.pdf", width=6, height=6)
plot(ladderize(time.tree), type="fan", show.tip.label=FALSE)
nodelabels(col=node.cols, pch=16, cex=0.45)
theta <- seq(from=0, to=2*pi, by=2*pi/ntip)
segments(46*cos(theta), 46*sin(theta), 48*cos(theta), 48*sin(theta), col=tip.cols)
dev.off()

# Make a color bar
pdf(file="./Figures/Discordance_tree_colorbar.pdf", width=2, height=10)
color.bar <- function(ramp, min, max, nticks=11, ticks=seq(min, max, len=nticks)) { 
	scale = (length(ramp)-1)/(max-min) 
	plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='') 
	axis(2, ticks, las=1) 
	for (i in 1:(length(ramp)-1)) { 
		y = (i-1)/scale + min 
		rect(0,y,10,y+1/scale, col=ramp[i], border=NA) 
	} 
}
color.bar(colorRampPalette(c("cornflowerblue", "cornflowerblue", "cornflowerblue", "bisque", "tomato1"))(1000), 0, 1)
dev.off()


# Discordance map

# Loop to get range maps
file.NAs <- vector()
polygons <- vector()
my.disc <- vector()
for(i in 1:length(tip.disc)) {
	print(names(tip.disc[i]))
	my.file <- paste0('./data/range_maps_shape_AOSHM/', names(tip.disc[i]), '/', names(tip.disc[i]), '.shp')
	if(file.exists(my.file)) {
		range <- readShapeSpatial(my.file)	
		range.dissolved <- try(unionSpatialPolygons(range, rep(as.character(i), length(range@polygons))))
		if(class(range.dissolved) == "try-error") {
			print(paste("Range error for: ", names(tip.disc[i])))
		} else {
			polygons <- c(polygons, range.dissolved)			
			my.disc <- c(my.disc, tip.disc[i])
		}
	} else {
		print(paste("No range data for: ", names(tip.disc[i])))
		file.NAs <- c(file.NAs, names(tip.disc[i]))
	}
}
print(file.NAs)
joined = SpatialPolygons(lapply(polygons, function(x){x@polygons[[1]]}))

# Make an empty raster with 200 km resolution
r <- raster(xmn = -12000000, xmx = 16000000, ymn = -7000000, ymx = 8000000, vals=0, crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs") # whole world

# For each species, make raster for discordance and another for weights
disc.rasters <- stack()
count.rasters <- stack()
for(i in 1:length(joined)) {
	print(i)
	cell.count <- lengths(extract(r, joined[i]))
	disc.raster <- r
	count.raster <- r
	disc.raster[joined[i],] <- my.disc[i]
	count.raster[joined[i],] <- 1/cell.count
	disc.rasters <- stack(disc.rasters, disc.raster)
	count.rasters <- stack(count.rasters, count.raster)
}
wm <- weighted.mean(x=disc.rasters, w=count.rasters, na.rm=TRUE)
Pal <- colorRampPalette(c("cornflowerblue", "cornflowerblue", "cornflowerblue", "bisque", "tomato1"))(11)
min(values(wm), na.rm=TRUE)
max(values(wm), na.rm=TRUE)
mean(values(wm), na.rm=TRUE)

# Plot
pdf(file="./figures/Discordance_map.pdf", height=5, width=9)
par(mar=c(1,1,1,5))
proj4string(wm) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
coastline <- readOGR('./data/other_map_data/ne_50m_coastline/ne_50m_coastline.shp')
coastline2 <- spTransform(coastline, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
sPDF <- getMap()
sPDF <- spTransform(sPDF, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
mapCountryData(sPDF, nameColumnToPlot="continent", colourPalette=c("white", "white"), addLegend=FALSE, borderCol="white", mapTitle="")
polygon((bbox(coastline2)[1,2]*1.024)*sin(seq(0,2*pi,length.out=100)),(bbox(coastline2)[2,2]*1.024)*cos(seq(0,2*pi,length.out=100)), lwd=1, border="black", col="gray")
mapCountryData(sPDF, nameColumnToPlot="continent", colourPalette=c("white", "white"), addLegend=FALSE, borderCol="white", mapTitle="", add=TRUE)
plot(wm, col=Pal, add=TRUE)
plot(coastline2, lwd=0.5, col="black", add=TRUE)
polygon((bbox(coastline2)[1,2]*1.024)*sin(seq(0,2*pi,length.out=100)),(bbox(coastline2)[2,2]*1.024)*cos(seq(0,2*pi,length.out=100)), lwd=1, border="black")
dev.off()

