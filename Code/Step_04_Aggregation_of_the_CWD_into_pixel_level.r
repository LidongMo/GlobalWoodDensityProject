###*####*####*####*####*####*####*####*####*####*####*####*####*####*####*##
# STEP 4_1. aggregate the plots in the same pixel(1km) into one
###*####*####*####*####*####*####*####*####*####*####*####*####*####*####*##
# one raster for wood density mean and one for variation
library(raster)
library(data.table)
library(parallel)
library(vegan)
library(dplyr)
library(ape)
library(geiger)
library(phytools)
setwd("~/WoodDensityProject/")
# set the temperary folder to save the temperary files
# rasterOptions(tmpdir = "/Volumes/Scratch/Lidong_Mo/Temp/")

# RASTERIZE the points to a raster of biodiveristy
plotWDtable = fread("Data/Plot_Level_Wood_Density_Data_Frame_GFBI_Metadata_Conservative_with_Diversity.csv")[,-1]  %>% mutate(Richness = SpeciesNumber/PlotArea)
# rasterize into a raster
StartRaster = raster(nrows=17400, ncols=43200, xmn=-180.0001, xmx=179.9999, ymn=-90.00014,ymx=83.99986,crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",resolution=c(0.008333333, 0.008333333))
# backgroud=0, is an option to set the na values in the data absence cells
RecordsRaster = rasterize(plotWDtable[,c("LON","LAT")], StartRaster, subset(plotWDtable,select=c("WeightedWoodDensity", "WoodDensitySD","SpeciesNumber","Richness","ShannonIndex","SimpsonIndex","InvSimpsonIndex","MinWD","MaxWD","FiveQuantile","NintyfiveQuantile","treeDensity","diameterMean","angiospermRatio","GymnoWD","AngioWD")), fun=mean,na.rm=T)
# write raster into the local folder
writeRaster(RecordsRaster,"Data/Plot_Level_Biodiveristy_Rasterized_raster_GFBI_all_Vars.tif",overwrite=T)


###*####*####*####*####*####*####*####*####*####*####*####*####*####*####*##
# STEP 4_2. transfer the pixel level wood density maps into csv format
###*####*####*####*####*####*####*####*####*####*####*####*####*####*####*##
library(dplyr)
library(raster)
library(tidyr)
# stack all the raster layers
dataStack = stack("Data/Plot_Level_Biodiveristy_Rasterized_raster_GFBI_all_Vars.tif")
# transform of the dataStack 
variRasterToTable = as.data.frame(dataStack,xy=TRUE) %>%tidyr::drop_na(c("layer.1")) # drop the rows with the wood density is NA
# 582435 observations left
# rename the data frame
names(variRasterToTable) = c("x","y","WoodDensity", "WoodDensitySD","SpeciesNumber","Richness","ShannonIndex","SimpsonIndex","InvSimpsonIndex","MinWD","MaxWD","FiveQuantile","NintyfiveQuantile","TreeDensity","DiameterMean","angiospermRatio","GymnoWD","AngioWD")

# extract the biome information
biomeLayer = raster("Data/WWFRaster/WWF_Raster_Uniformed.tif")
# extract the wwf biome information
extractedTable = raster::extract(biomeLayer,as.data.frame(variRasterToTable[,c("x","y")]))
variRasterToTable$Biome = extractedTable

write.csv(variRasterToTable,"Data/Pixel_Level_Wood_Density_Meta_Table_with_All_Vars_Biome.csv") #591277   row

###*####*####*####*####*####*####*####*####*####*####*####*####*####*####*##
# STEP 4_3. transfer the aggregated pixel level observations into shapefile
# #***********************************************
# # tranfer this table into shapefiles

# all the three shapefile can be accessed in Goolge earth engine
library(dplyr)
trainDataTable = fread("Data/Pixel_Level_Wood_Density_Meta_Table_with_All_Vars_Biome.csv")[,-1] 
# kick out the outliers based on MAD approach 
MADValue = mad(trainDataTable$WoodDensity)
# define the large value margin
largeRange = 3*MADValue + median(trainDataTable$WoodDensity)
# defein the small value margin 
smallRange = median(trainDataTable$WoodDensity) -3*MADValue 
# subset the data frame by the MAD small and large margin range
cleanedTrainTable = trainDataTable[trainDataTable$WoodDensity<=largeRange&trainDataTable$WoodDensity>=smallRange,]

# metaTable = na.omit(metaTable)
dupilcateTable = cleanedTrainTable
# allote the coordinates
coordinates(cleanedTrainTable) = ~x+y
# allocate projection system
proj4string(cleanedTrainTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
cleanedTrainTable@data = dupilcateTable  
# write the shapefile into local folder
# writeOGR(cleanedTrainTable,dsn="ShapefilesFolder",layer="WoodDensity_Diversity_Pixel_2023", driver = "ESRI Shapefile",overwrite=T) #564197 this is the data for modeling, which was in the GEE already
# the file below is for the later updated table with speciesNumber standardized by plot size 2023/12/08

writeOGR(cleanedTrainTable,dsn="Data/ShapefilesFolder",layer="WoodDensity_Diversity_SpiciesNumber_updated_2023", driver = "ESRI Shapefile",overwrite=T) #564197  


trainDataTable = fread("Data/Pixel_Level_Wood_Density_Meta_Table_with_All_Vars_Biome.csv")[,-1] %>% dplyr::select("x","y","AngioWD") %>% na.omit()
# kick out the outliers based on MAD approach 
MADValue = mad(trainDataTable$AngioWD)
# define the large value margin
largeRange = 3*MADValue + median(trainDataTable$AngioWD)
# defein the small value margin 
smallRange = median(trainDataTable$AngioWD) -3*MADValue 
# subset the data frame by the MAD small and large margin range
cleanedTrainTable = trainDataTable[trainDataTable$AngioWD<=largeRange&trainDataTable$AngioWD>=smallRange,]

# metaTable = na.omit(metaTable)
dupilcateTable = cleanedTrainTable
# allote the coordinates
coordinates(cleanedTrainTable) = ~x+y
# allocate projection system
proj4string(cleanedTrainTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
cleanedTrainTable@data = dupilcateTable  
# write the shapefile into local folder
writeOGR(cleanedTrainTable,dsn="Data/ShapefilesFolder",layer="WoodDensity_Diversity_Pixel_Angio", driver = "ESRI Shapefile",overwrite=T) #450486



library(dplyr)
trainDataTable = fread("Data/Pixel_Level_Wood_Density_Meta_Table_with_All_Vars_Biome.csv")[,-1] %>% dplyr::select("x","y","GymnoWD") %>% na.omit()
# kick out the outliers based on MAD approach 
MADValue = mad(trainDataTable$GymnoWD)
# define the large value margin
largeRange = 3*MADValue + median(trainDataTable$GymnoWD)
# defein the small value margin 
smallRange = median(trainDataTable$GymnoWD) -3*MADValue 
# subset the data frame by the MAD small and large margin range
cleanedTrainTable = trainDataTable[trainDataTable$GymnoWD<=largeRange&trainDataTable$GymnoWD>=smallRange,]

# metaTable = na.omit(metaTable)
dupilcateTable = cleanedTrainTable
# allote the coordinates
coordinates(cleanedTrainTable) = ~x+y
# allocate projection system
proj4string(cleanedTrainTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
cleanedTrainTable@data = dupilcateTable  
# write the shapefile into local folder
writeOGR(cleanedTrainTable,dsn="Data/ShapefilesFolder",layer="WoodDensity_Diversity_Pixel_Gymno", driver = "ESRI Shapefile",overwrite=T) #407136  