
library(raster)
library(sf) # for reading shape files; read_sf
library(scales) # for viridis
library(stars) # for st_rasterize

# read in data
# elev <- raster("elevation_raster.asc")
# rain <- raster("rain_raster.asc")
load("better_elev/elev_rain.rdata")
elev <- elev_r_wai
rain <- rain_r_wai
sitepoints <- read.delim("sampling_site_points.csv", sep=",", stringsAsFactors=FALSE)
river <- read_sf("waimea_river.kml")
coast <- read_sf("oahu_cost.dbf")

# subset points for marine vs not marine
sitepoints_mar <- sitepoints[sitepoints$habitat == "Marine", ]
sitepoints <- sitepoints[sitepoints$habitat != "Marine", ]

# calculate ocean
ocean <- st_rasterize(coast, st_as_stars(elev))
ocean[ocean > 0] <- NA
ocean[!is.na(ocean)] <- 1

# define colors
elevcols <- viridis_pal(begin=0.3, end=1)(10)
watercol <- "cyan"

# make draft plot
plot(elev, col=elevcols)
plot(ocean, col=watercol, add=T, legend=FALSE)
plot(coast[1], add=T, col=rgb(0,0,0,0), border=watercol, lwd=2)
plot(river[1], lwd=2, col=watercol, add=TRUE)
points(sitepoints$latitude_deg ~ sitepoints$longitude_deg, col="white", pch=20, cex=2)
points(sitepoints$latitude_deg ~ sitepoints$longitude_deg, col="red", pch=17, cex=1)
points(sitepoints_mar$latitude_deg ~ sitepoints_mar$longitude_deg, col="white", pch=20, cex=2)
points(sitepoints_mar$latitude_deg ~ sitepoints_mar$longitude_deg, col="blue", pch=17, cex=1)

# rectangle to figure out cropping
xlims <- c(-158.077, -157.95)
ylims <- c(21.586, 21.6509)
cropbox <- as(extent(xlims[1], xlims[2], ylims[1], ylims[2]), 'SpatialPolygons')
crs(cropbox) <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"
plot(cropbox, add=TRUE)

# crop stuff
elev_crop <- crop(elev, cropbox)
ocean_crop <- st_crop(ocean, st_as_sf(cropbox))
st_crs(coast) <- crs(cropbox) # will give warning, doesn't matter since projections are the same really
coast_crop <- st_crop(coast, st_as_sf(cropbox))

# make map plot
pdf("map.pdf", useDingbats = FALSE)
plot(elev_crop, col=elevcols)
plot(ocean_crop, col=watercol, add=T, legend=FALSE)
plot(coast_crop[1], add=T, col=rgb(0,0,0,0), border=watercol, lwd=2)
plot(river[1], lwd=2, col=watercol, add=TRUE)
points(sitepoints$latitude_deg ~ sitepoints$longitude_deg, col="white", pch=20, cex=3)
points(sitepoints$latitude_deg ~ sitepoints$longitude_deg, col="red", pch=17, cex=1.4)
points(sitepoints_mar$latitude_deg ~ sitepoints_mar$longitude_deg, col="white", pch=20, cex=3)
points(sitepoints_mar$latitude_deg ~ sitepoints_mar$longitude_deg, col="blue", pch=17, cex=1.4)
dev.off()

# plot elev and rainfall for terrestrial sites
sitepoints <- sitepoints[order(sitepoints$longitude_deg),]
sitepoints$rain <- extract(rain, data.frame(sitepoints$longitude_deg, sitepoints$latitude_deg))
pdf("rain_elev.pdf", useDingbats = FALSE)
par(mfrow=c(2,1))
plot(sitepoints$elevation ~ sitepoints$longitude_deg, type="o", col="red", pch=20, 
	cex=2, lwd=2, xlim=xlims, xaxs="i")
plot(sitepoints$rain ~ sitepoints$longitude_deg, type="o", col="blue", pch=20, 
	cex=2, lwd=2, xlim=xlims, xaxs="i")
dev.off()

# make inset maps

# NOTE: because of github space limitations, the folder "better_elev"
# must be downloaded from FigShare and unzipped here.
# https://figshare.com/articles/dataset/Hawaii_elevation_data_for_waimea_microbe_paper/14538249
# https://doi.org/10.6084/m9.figshare.14538249.v1 

elev_oahu <- merge( raster("better_elev/USGS_13_n22w158.tif"), raster("better_elev/USGS_13_n22w159.tif"))
oahx <- c(-158.33, -157.6)
oahy <- c(21.22, 21.73)
oahu_crop <- as(extent(oahx[1], oahx[2], oahy[1], oahy[2]), 'SpatialPolygons')
crs(oahu_crop) <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs "
elev_oahu <- crop(elev_oahu, oahu_crop)

# make vector for ocean
oahu_crop_st <- as.matrix(data.frame(x=c(oahx, rev(oahx)),y=rep(oahy,each=2)))
oahu_crop_st <- rbind(oahu_crop_st, oahu_crop_st[1,])
oahu_crop_st <- st_polygon(list(oahu_crop_st))
oahu_crop_st <- st_geometry(oahu_crop_st)
st_crs(oahu_crop_st) <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"
oceanpoly <- st_difference( oahu_crop_st, coast[1])[1]



pdf("map_inset.pdf", useDingbats = FALSE)
plot(elev_oahu, col=elevcols)
plot(oceanpoly, col=watercol, border=watercol, add=T)
plot(cropbox, border="black", add=T, lwd=4)
plot(cropbox, border="yellow", add=T, lwd=2)
dev.off()
