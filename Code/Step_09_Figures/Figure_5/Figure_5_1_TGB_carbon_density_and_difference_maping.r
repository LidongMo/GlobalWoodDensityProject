# load the library
library(raster)
library(RColorBrewer)
library(rgeos)
library(sp)
library(rgdal)
# load the forestcover map from hansen of the year 2010

tgbMerged = raster("MergedMaps/WD_derived_TGB_Carbon_Density_Merged.tif")

# change the values larger than 500 to 500
tgbMergedCarbon = tgbMerged
tgbMergedCarbon[tgbMergedCarbon>250] = 250
tgbMergedCarbon[tgbMergedCarbon<=0] = NA


# load the difference percentage map
differenceMap = raster("MergedMaps/WD_derived_Carbon_Density_Maps_Difference_Merged.tif")
differenceMap[differenceMap>0.2] = 0.2
differenceMap[differenceMap< -0.2] = -0.2
differenceMap = differenceMap*100 #transfer to percentage
# > quantile(tt)
#          0%         25%         50%         75%        100% 
# -1.00000000 -0.12675101 -0.01429734  0.15061748  0.29027703 
# pal= colorRampPalette(brewer.pal(9, "YlGnBu"))

# define the equal earth projection
equalEarthProj = "+proj=eqearth +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# project raster to equal earth
extentInfo = projectExtent(tgbMerged, equalEarthProj)
# reproject the rasters
tgbCarbonProjected = projectRaster(tgbMergedCarbon,to=extentInfo,over=T)
differenceMap_Projected = projectRaster(differenceMap, to=extentInfo,over=T)


# since the differenc was already scaled into t/ha, we just directly plot them

worldBorderInit = shapefile("WORLD_BORDERS/TM_WORLD_BORDERS_SIMPL-0.3.shp")
worldBorderRaw = crop(worldBorderInit,c(-180, 180, -60, 84))
# set the resolution of how the polygon plot process
worldBorder = gSimplify(worldBorderRaw, tol = 0.01, topologyPreserve = TRUE)
# reproject the world border
worldEqualEarth = spTransform(worldBorder,to=extentInfo, CRS(equalEarthProj))

pdf("Plots/Fig_5_1_TGB_carbon_density_and_Difference.pdf",width =12, height=12)
par(mfrow=c(2,1))
pal= colorRampPalette((c("gray80","#A1D99B" ,"#74C476", "#41AB5D" ,"#238B45", "#006D2C", "#00441B","#004521")))
plot(worldEqualEarth,border="gray70",col="gray70",axes=T,main="Tree carbon density",lwd=0.2)
plot(tgbCarbonProjected,col=pal(100),add=T,breaks=c(seq(0,250,2.5)),legend=T,maxpixels = 40000000)

pal= colorRampPalette(c('#F46D43', '#FDAE61', '#FEE090', '#FFFFFF','#E0F3F8', '#ABD9E9','#74ADD1'))

plot(worldEqualEarth,border="gray70",col="gray70",axes=T,main="Difference of biomass in percentage",lwd=0.2)
plot(differenceMap_Projected,col=pal(100),add=T,breaks=c(seq(-20,20,0.4)),legend=T,maxpixels = 40000000)
dev.off()

