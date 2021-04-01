library(hadsstr)

ranges = readRDS("~/Desktop/range_maps_shape_AOSHM.projected.Rds")
get_extent <- function(r) {
  as.vector(extent(r))
}
extents = pblapply(ranges, get_extent)
ex = do.call("rbind", extents)
ex = data.frame(ex)

# for buffer
to_crop = extent(range(ex$X1, ex$X2) * 1.02,
                 range(ex$X3, ex$X4) * 1.02)
# 30 arc seconds, ~1km data
cur_t = raster("data/spatial_data/CHELSA_bio10_01.tif")
cur_t = crop(cur_t, to_crop)
past_t = raster("data/spatial_data/CHELSA_PMIP_CCSM4_BIO_01.tif")
past_t = crop(past_t, to_crop)
past_t1 = (past_t / 10)

# https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12184
# https://science.sciencemag.org/content/334/6056/660.full
# https://science.sciencemag.org/content/sci/suppl/2011/11/03/334.6056.652.DC1/Burrows.SOM.pdf
# https://www.nature.com/articles/nature08649
# https://www.nature.com/articles/nclimate3127
# Climate‐change velocity is thus a 
# measure of the local temporal rate of
# geographical displacement of climatic conditions, 
# calculated by dividing the temperature change over time
# by the local temperature change across space, 
# which is lower where elevation gradients are present. 
for (i in 1:length(ranges)) {
  range1 = ranges[[i]]
  
  cur1 = crop(cur_t, range1)
  past1 = crop(past_t1, range1)
  
  ns = get_NS_diffs(cur1)
  we = get_WE_diffs(cur1)
  
  tg = cur1 - past1
  sg = get_spatial_gradient(ns, we)
  # all values 0.01 °C/km were replaced with 0.01
  sg2 = sg[[1]]
  sg2[sg2 < 0.01] = 0.01
  
  vel = tg / sg2
  vel2 = raster::extract(vel, range1)[[1]]
  tdiff = raster::extract(tg, range1)[[1]]
  
  df = data.frame(velocity = vel2, tempdiff = tdiff)
  write.csv(df, paste0("~/Desktop/", names(ranges)[i], ".csv"))
}