###*####*####*####*####*####*####*####*####*####*####*####*####*####*####*##
# STEP 2. allocate the woode density data to each individual of the GFBI data
#         allocate the angisperm or gymnosperm information
###*####*####*####*####*####*####*####*####*####*####*####*####*####*####*##
library(data.table)
library(parallel)
library(Taxonstand)
setwd("~/WoodDensityProject/")

# due to the data using policy of GFBi database, the raw data will not be shared here, but you will the dervied CWD information
biomassBiomeTPH = fread("****/GFBI_Biomass_Tree_Level_Biomass_New_PLT_Data_Frame_20181024.csv")[,-1]
# get the plot names
plotNames = as.vector(unique(biomassBiomeTPH$NewPLT))
woodDensityReference = fread("Data/Resources/GFBI_Wood_Density_Reference_Table_with_A_G.csv")[,-1]
# in order to get the plot with time series, we just to check the plot with multiple years records
woodDensityAllocatingFunc = function(pn)
{
    # read the data frame for the each plot folder
    perplotDF = fread(paste("****/PerPLT/PLT_",pn,".csv",sep=""))[,-1]
    print(dim(perplotDF))
    # apply a function each row
    rowFunc = function(rn) 
    {
        referedRow = woodDensityReference[woodDensityReference$RawName == gsub("[^[:alnum:] ]", "", perplotDF[rn]$SPCD),]
        return(referedRow[1,c("WoodDensity","WD_Level","A_G","AG_Level")])#$WoodDensity,referedRow$level
    }
    # allocate the wood density to the perplot data frame
    lapplyList  = lapply(1:nrow(perplotDF),rowFunc)
    lapplyResultTable = rbindlist(lapplyList)
    # this allocation approach can avoid multiple running times resulting a duplicated wood density and level columns
    perplotDF$WoodDensity = lapplyResultTable$WoodDensity
    perplotDF$WD_Level = lapplyResultTable$WD_Level
    # this allocation approach can avoid multiple running times resulting a duplicated wood density and level columns
    perplotDF$A_G = lapplyResultTable$A_G
    perplotDF$AG_Level = lapplyResultTable$AG_Level
    print(dim(perplotDF))
    # write the data into the local folder
    write.csv(perplotDF,paste("****/PerPlotWithWoodDensity/",pn,"_GFBI_Data_With_Individual_Wood_Density.csv",sep=""))
    print(pn)
    # return(perplotDF)
}

# paralle running of it
system.time(mclapply(plotNames,woodDensityAllocatingFunc,mc.preschedule=F,mc.cores=32))
