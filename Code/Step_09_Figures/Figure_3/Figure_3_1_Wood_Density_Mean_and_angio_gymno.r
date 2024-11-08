# this is the code for Figure 3 making

# load the libraries
#  the wood density predication map
library(raster)
library(RColorBrewer)
library(rgeos)
library(viridis)
library(sp)
# library(rgdal)
# library(sf)

# maps coould be downloade from zenodo composite

# load the world border
# worldBorderRaw = shapefile("WORLD_BORDERS/TM_WORLD_BORDERS_SIMPL-0.3.shp")
# worldBorderRaw = crop(worldBorderRaw,c(-180, 180, -60, 84)) 
# # set the resolution of how the polygon plot process
# worldBorder = gSimplify(worldBorderRaw, tol = 0.01, topologyPreserve = TRUE)
# load the wood density map
woodDensityMerged = raster("MergedMaps/Community_Wood_Density_Map_Merged.tif")
woodDensityMerged[woodDensityMerged<=0] = NA
woodDensityMerged[woodDensityMerged<0.4] =0.4
woodDensityMerged[woodDensityMerged>0.75] =0.75

angioMerged = raster("MergedMaps/Angio_Wood_Density_Map_Merged.tif")
angioMerged[angioMerged<=0] = NA
angioMerged[angioMerged<0.4] =0.4
angioMerged[angioMerged>0.75] =0.75

gymnoMerged = raster("MergedMaps/Gymno_Wood_density_Map_Merged.tif")
gymnoMerged[gymnoMerged<=0] = NA
gymnoMerged[gymnoMerged<0.4] =0.4
gymnoMerged[gymnoMerged>0.6] =0.6

# # load the uncertainty maps for wood density and biomass estimations
# modelUncertainMapWD = raster("MergedMaps/Wood_density_variation_coefficients_Map_Merged.tif")
# modelUncertainMapWD[modelUncertainMapWD>=0.05] =0.05
# modelUncertainMapWD[modelUncertainMapWD <= 0] =NA
 # define the equal earth projection
equalEarthProj = "+proj=eqearth +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# project raster to equal earth
extentInfo = projectExtent(woodDensityMerged, equalEarthProj)
# reproject the rasters
woodDensityProjected = projectRaster(woodDensityMerged, to=extentInfo,over=T)
angioProjected = projectRaster(angioMerged, to=extentInfo,over=T)
gymnoProjected = projectRaster(gymnoMerged, to=extentInfo,over=T)
# modelUncertainWD_Projected = projectRaster(modelUncertainMapWD, to=extentInfo,over=T)


worldBorderInit = shapefile("WORLD_BORDERS/TM_WORLD_BORDERS_SIMPL-0.3.shp")
worldBorderRaw = crop(worldBorderInit,c(-180, 180, -60, 84))
# set the resolution of how the polygon plot process
worldBorder = gSimplify(worldBorderRaw, tol = 0.01, topologyPreserve = TRUE)
# reproject the world border
worldEqualEarth = spTransform(worldBorder,to=extentInfo, CRS(equalEarthProj),over=T)

# make the plot and write into local PDF
pdf("Plots/Fig_3_Community_wood_density_Mean_Angio_and_Gymono.pdf",width =12, height=18)
par(mfrow=c(3,1))

rbPal1= colorRampPalette(c("#FDE725FF","#3CBB75FF","#2D708EFF","#481567FF"))

plot(worldEqualEarth,border="gray70",col="gray70",axes=T,main="Community wood density mean",lwd=0.2)
plot(woodDensityProjected,col=rbPal1(100),add=T,breaks=c(seq(0.4,0.75,0.0035)),legend=T,maxpixels = 40000000)

plot(worldEqualEarth,border="gray70",col="gray70",axes=T,main="Angio wood density mean",lwd=0.2)
plot(angioProjected,col=rbPal1(100),add=T,breaks=c(seq(0.4,0.75,0.0035)),legend=T,maxpixels = 40000000)

plot(worldEqualEarth,border="gray70",col="gray70",axes=T,main="Gymno wood density mean",lwd=0.2)
plot(gymnoProjected,col=rbPal1(100),add=T,breaks=c(seq(0.4,0.55,0.0015)),legend=T,maxpixels = 40000000)

dev.off() 

# rbPal2= colorRampPalette(c("#F7F7F7","#E6F598","#ABDDA4","#66C2A5","#3288BD","#5E4FA2"))

# plot(worldEqualEarth,border="gray70",col="gray25",axes=T,main="Uncertainty",lwd=0.2)
# plot(modelUncertainWD_Projected,col=rbPal2(100),add=T,breaks=c(seq(0,0.05,0.0005)),legend=T,maxpixels = 40000000)


