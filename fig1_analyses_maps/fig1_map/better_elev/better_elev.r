# Elevation data source:
# https://www.sciencebase.gov/catalog/item/5f7783a482ce1d74e7d6c123
# https://www.sciencebase.gov/catalog/item/5f7783a582ce1d74e7d6c125

# Rainfall data source:
# http://rainfall.geography.hawaii.edu/downloads.html

library(raster)


# read in elev (masl) raster
b <- raster("USGS_13_n22w158.tif")
a <- raster("USGS_13_n22w159.tif")
elev_r <- merge(a,b)
rm(a,b)

# read in mean annual rainfall (mm) raster
rain_r <- raster("OahuASCIIGrids_mm/rfgrid_mm_oahu_ann.txt")
crs(rain_r) <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"

# subset for waimea area
samppoints <- read.delim("../sampling_site_points.csv", sep=",", stringsAsFactors=FALSE)

plot(elev_r)
points(x=samppoints$longitude_deg, y=samppoints$latitude_deg)

# extent: w,e,s,n
e <- as(extent(-158.11, -157.9, 21.55, 21.68), 'SpatialPolygons')
crs(e) <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"
plot(e, add=T)
# r <- crop(worldpopcount, e)

elev_r_wai <- crop(elev_r, e)
plot(elev_r_wai)
points(x=samppoints$longitude_deg, y=samppoints$latitude_deg)

rain_r_wai <- crop(rain_r, e)
plot(rain_r_wai)
points(x=samppoints$longitude_deg, y=samppoints$latitude_deg)

# lazily save rasters as rdata
save(list=c("elev_r_wai", "rain_r_wai"), file="elev_rain.rdata")