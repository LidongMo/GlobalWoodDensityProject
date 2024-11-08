
library(data.table)
library(raster)
library(parallel)
library(rgdal)
library(stringr)

# define an empty data frame to rbind the extractions
randomSampledTable = data.frame()
rSquaredTable = data.frame()
# load the biomass maps
woodDensityDiff = raster("MergedMaps/WD_derived_Carbon_Density_Maps_Difference_Merged.tif")

wwfBiome = raster("WWFRaster/WWF_Biome_Merged.tif")

mapStack = stack(woodDensityDiff,wwfBiome)
# do random sampling
set.seed(1000)
randomSampled_1_km = as.data.frame(sampleRandom(mapStack, 50000000, na.rm=TRUE,xy=T))
# # replace all the 0s to NA
# randomSampled_1_km[randomSampled_1_km==0] = NA
# randomSampled_1_km = na.omit(randomSampled_1_km)
names(randomSampled_1_km) = c("x","y","Diff","Biome")

randomSampledTable = rbind(randomSampledTable,randomSampled_1_km)

write.csv(randomSampled_1_km,"ModelCompare/Random_sample_difference_of_derived_biomass.csv")



