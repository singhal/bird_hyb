setwd("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/")
getwd()

library(rgeos)

# Get Dstat data
data = read.csv("results/dstat_af_subsample250.csv", stringsAsFactors = F)

# Get triads involved in significant comparisons
sig.rows <- data[data$p < 0.05,]
nonsig.rows <- data[data$p >= 0.05,] 
nrow(sig.rows)
nrow(nonsig.rows)

pdf(file="./figures/sig_nonsig_dstat_ranges.pdf", height=8, width=7)
par(mfrow=c(2,1))

# Map significant triads
polygons <- vector()
for(i in 1:nrow(sig.rows)) {
	
	print(sig.rows[i,]$sp1)
	sp1.file <- paste0('./data/range_maps_shape_AOSHM/', sig.rows[i,]$sp1, '/', sig.rows[i,]$sp1, '.shp')
	sp2.file <- paste0('./data/range_maps_shape_AOSHM/', sig.rows[i,]$sp2, '/', sig.rows[i,]$sp2, '.shp')
	sp3.file <- paste0('./data/range_maps_shape_AOSHM/', sig.rows[i,]$sp3, '/', sig.rows[i,]$sp3, '.shp')
	sp1.range <- readShapeSpatial(sp1.file)	
	sp2.range <- readShapeSpatial(sp2.file)	
	sp3.range <- readShapeSpatial(sp3.file)		
	
	sp1.range.dissolved <- try(unionSpatialPolygons(sp1.range, rep(as.character(i), length(sp1.range@polygons))))
	sp2.range.dissolved <- try(unionSpatialPolygons(sp2.range, rep(as.character(i), length(sp2.range@polygons))))
	sp3.range.dissolved <- try(unionSpatialPolygons(sp3.range, rep(as.character(i), length(sp3.range@polygons))))
	
	range.dissolved <- gUnion(sp1.range.dissolved, sp2.range.dissolved, byid=TRUE)
	range.dissolved <- gUnion(range.dissolved, sp3.range.dissolved, byid=TRUE)

	#plot(range.dissolved, col=alpha("red", 0.5))
	#plot(sp1.range.dissolved, col=alpha("blue", 0.5), add=TRUE)
	#plot(sp2.range.dissolved, col=alpha("blue", 0.5), add=TRUE)
	#plot(sp3.range.dissolved, col=alpha("blue", 0.5), add=TRUE)

	polygons <- c(polygons, range.dissolved)			
}
joined = SpatialPolygons(lapply(polygons, function(x){x@polygons[[1]]}))

# Make an empty raster with 200 km resolution
r <- raster(xmn = -12000000, xmx = 16000000, ymn = -7000000, ymx = 8000000, vals=0, crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs") # whole world

# For each species, make raster for discordance and another for weights
count.rasters <- stack()
for(i in 1:length(joined)) {
	print(i)
	cell.count <- lengths(extract(r, joined[i]))
	count.raster <- r
	count.raster[joined[i],] <- 1
	count.rasters <- stack(count.rasters, count.raster)
}
wm <- sum(count.rasters, na.rm=TRUE)
values(wm)[values(wm)==0] <- NA
Pal <- colorRampPalette(c("cornflowerblue", "bisque", "tomato1"))(11)
min(values(wm), na.rm=TRUE)
max(values(wm), na.rm=TRUE)
mean(values(wm), na.rm=TRUE)

# Plot
par(mar=c(1,1,1,5))
proj4string(wm) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
coastline <- readOGR('./data/other_map_data/ne_50m_coastline/ne_50m_coastline.shp')
coastline2 <- spTransform(coastline, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
sPDF <- getMap()
sPDF <- spTransform(sPDF, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
mapCountryData(sPDF, nameColumnToPlot="continent", colourPalette=c("white", "white"), addLegend=FALSE, borderCol="white", mapTitle="")
polygon((bbox(coastline2)[1,2]*1.024)*sin(seq(0,2*pi,length.out=100)),(bbox(coastline2)[2,2]*1.024)*cos(seq(0,2*pi,length.out=100)), lwd=1, border="black", col="gray")
mapCountryData(sPDF, nameColumnToPlot="continent", colourPalette=c("white", "white"), addLegend=FALSE, borderCol="white", mapTitle="", add=TRUE)
plot(wm, col=Pal, add=TRUE, main="Significant D-statistics")
plot(coastline2, lwd=0.5, col="black", add=TRUE)
polygon((bbox(coastline2)[1,2]*1.024)*sin(seq(0,2*pi,length.out=100)),(bbox(coastline2)[2,2]*1.024)*cos(seq(0,2*pi,length.out=100)), lwd=1, border="black")


# Map significant triads
polygons <- vector()
for(i in 1:nrow(nonsig.rows)) {
	
	print(nonsig.rows[i,]$sp1)
	sp1.file <- paste0('./data/range_maps_shape_AOSHM/', nonsig.rows[i,]$sp1, '/', nonsig.rows[i,]$sp1, '.shp')
	sp2.file <- paste0('./data/range_maps_shape_AOSHM/', nonsig.rows[i,]$sp2, '/', nonsig.rows[i,]$sp2, '.shp')
	sp3.file <- paste0('./data/range_maps_shape_AOSHM/', nonsig.rows[i,]$sp3, '/', nonsig.rows[i,]$sp3, '.shp')
	sp1.range <- readShapeSpatial(sp1.file)	
	sp2.range <- readShapeSpatial(sp2.file)	
	sp3.range <- readShapeSpatial(sp3.file)		
	
	sp1.range.dissolved <- try(unionSpatialPolygons(sp1.range, rep(as.character(i), length(sp1.range@polygons))))
	sp2.range.dissolved <- try(unionSpatialPolygons(sp2.range, rep(as.character(i), length(sp2.range@polygons))))
	sp3.range.dissolved <- try(unionSpatialPolygons(sp3.range, rep(as.character(i), length(sp3.range@polygons))))
	
	range.dissolved <- gUnion(sp1.range.dissolved, sp2.range.dissolved, byid=TRUE)
	range.dissolved <- gUnion(range.dissolved, sp3.range.dissolved, byid=TRUE)

	#plot(range.dissolved, col=alpha("red", 0.5))
	#plot(sp1.range.dissolved, col=alpha("blue", 0.5), add=TRUE)
	#plot(sp2.range.dissolved, col=alpha("blue", 0.5), add=TRUE)
	#plot(sp3.range.dissolved, col=alpha("blue", 0.5), add=TRUE)

	polygons <- c(polygons, range.dissolved)			
}
joined = SpatialPolygons(lapply(polygons, function(x){x@polygons[[1]]}))

# Make an empty raster with 200 km resolution
r <- raster(xmn = -12000000, xmx = 16000000, ymn = -7000000, ymx = 8000000, vals=0, crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs") # whole world

# For each species, make raster for discordance and another for weights
count.rasters <- stack()
for(i in 1:length(joined)) {
	print(i)
	cell.count <- lengths(extract(r, joined[i]))
	count.raster <- r
	count.raster[joined[i],] <- 1
	count.rasters <- stack(count.rasters, count.raster)
}
wm <- sum(count.rasters, na.rm=TRUE)
values(wm)[values(wm)==0] <- NA
Pal <- colorRampPalette(c("cornflowerblue", "bisque", "tomato1"))(11)
min(values(wm), na.rm=TRUE)
max(values(wm), na.rm=TRUE)
mean(values(wm), na.rm=TRUE)

# Plot
par(mar=c(1,1,1,5))
proj4string(wm) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
coastline <- readOGR('./data/other_map_data/ne_50m_coastline/ne_50m_coastline.shp')
coastline2 <- spTransform(coastline, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
sPDF <- getMap()
sPDF <- spTransform(sPDF, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
mapCountryData(sPDF, nameColumnToPlot="continent", colourPalette=c("white", "white"), addLegend=FALSE, borderCol="white", mapTitle="")
polygon((bbox(coastline2)[1,2]*1.024)*sin(seq(0,2*pi,length.out=100)),(bbox(coastline2)[2,2]*1.024)*cos(seq(0,2*pi,length.out=100)), lwd=1, border="black", col="gray")
mapCountryData(sPDF, nameColumnToPlot="continent", colourPalette=c("white", "white"), addLegend=FALSE, borderCol="white", mapTitle="", add=TRUE)
plot(wm, col=Pal, add=TRUE, main="Non-significant D-statistics")
plot(coastline2, lwd=0.5, col="black", add=TRUE)
polygon((bbox(coastline2)[1,2]*1.024)*sin(seq(0,2*pi,length.out=100)),(bbox(coastline2)[2,2]*1.024)*cos(seq(0,2*pi,length.out=100)), lwd=1, border="black")
dev.off()
?plot.raster()