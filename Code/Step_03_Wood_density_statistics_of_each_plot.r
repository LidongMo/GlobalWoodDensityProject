
###*####*####*####*####*####*####*####*####*####*####*####*####*####*####*##
# STEP 3_1. get the plot level wood density metadata which is weighted by basal area
###*####*####*####*####*####*####*####*####*####*####*####*####*####*####*##

library(data.table)
library(parallel)
library(vegan)
library(dplyr)
library(ape)
library(geiger)
library(phytools)
setwd("~/WoodDensityProject/")
# this calcualtion is baseon the data folder PerPlotWithWoodDensity 
# biomassBiomeTPH = fread("/Volumes/CrowtherLabRAID/Lidong_Mo/BiomassEstimation/GFBI_Biomass_Tree_Level_Biomass_New_PLT_Data_Frame_20181024.csv")[,-1]
# # get the plot names
# plotNames = as.vector(unique(biomassBiomeTPH$NewPLT))
# load the full phylogeny tree
fullPhylogenyTree = read.tree("Data/PhylogenyAnalysis/Phylogeny_tree_of_woody_plants.tre")
# load the subset of the wood density rederence table by the species with phylogeny information
woodDensityReferencePhylogeny = fread("Data/PhylogenyAnalysis/Filtered_Wood_Density_Table_for_Species_with_Phylogeny_Info.csv")[,-1]%>% as.data.frame() 
# get the plot names
plotNames = list.files(path="****/PerPlotWithWoodDensity",pattern= "*_GFBI_Data_With_Individual_Wood_Density.csv") %>% substr(1,12) %>% as.vector()

plotLevelStatFunc = function(pn=PlotName,inputPlot = InputPlot)
{
    # get the subset data for each plot by plot name
    perPlotDataFrameRaw = inputPlot
    # perPlotDataFrame2 = inputPlot
    # delete the individuals in the plots have dataset level wood density information
    perPlotDataFrame = perPlotDataFrameRaw[perPlotDataFrameRaw$WD_Level != 'dataset',]
    # if all the level in the data frame is 'dataset', then we will have a data frame with 0 row, then we just allocat an empty data frame with NAs
    if(nrow(perPlotDataFrame) == 0)
    {
        outputRow = data.frame( PLT=NA, LAT = NA, LON=NA, YEAR=NA, PlotArea = NA,IndividualNumber = NA,treeDensity=NA, diameterMean=NA,SpeciesNumber=NA, ShannonIndex=NA, SimpsonIndex = NA, InvSimpsonIndex = NA, BasalArea =NA, WeightedWoodDensity =NA, WoodDensitySD =NA,MinWD =NA,  MaxWD=NA,GymnoWD = NA,AngioWD = NA, FiveQuantile = NA, NintyfiveQuantile =NA, angiospermRatio =NA,gymnospermRatio =NA,nonDefinedAG =NA,dataset=NA,genus=NA,species= NA)
    }else
    {
        # calculate the basal area for each indiviudal
        perPlotDataFrame$basalArea = pi*((perPlotDataFrame$DBH/2)^2)
        # get the diameter mean
        diameterMean = mean(perPlotDataFrame$DBH)
        tphMean = mean(1/perPlotDataFrame$TPH)
        # get the total tree number in per plot
        individualNumber = nrow(perPlotDataFrame)
        # get the species number
        speciesNumber = length(unique(perPlotDataFrame$SPCD))
        # get the individual number in the plot
        individualNumber = nrow(perPlotDataFrame)
        # get the basal area weighted wood density
        WeightedWoodDensity = sum(perPlotDataFrame$basalArea * perPlotDataFrame$WoodDensity)/sum(perPlotDataFrame$basalArea)
        woodDensitySD= sd(perPlotDataFrame$WoodDensity)
        # extract the clumn as a vector
        columnVector = perPlotDataFrameRaw$WD_Level
        # get the ratio of the three levels in the level vector
        levelRatios = data.frame(dataset = length(which(columnVector== "dataset"))/length(columnVector),
                                genus =length(which(columnVector== "genus"))/length(columnVector),
                                species = length(which(columnVector== "species"))/length(columnVector))
        minimumWD = min(perPlotDataFrame$WoodDensity)
        maximumWD = max(perPlotDataFrame$WoodDensity)
        # get the quantile of the wood density at 5% and 95% based on individual information
        quantileResult = as.vector(quantile(perPlotDataFrame$WoodDensity,probs = c(0.05,0.95)))
        # extract the anigosperm and gymnosperm tree based wood density
        angioTable = perPlotDataFrame %>% filter(A_G == "a") %>% mutate(basalArea = pi*((DBH/2)^2))
        gymnoTable = perPlotDataFrame %>% filter(A_G == "g") %>% mutate(basalArea = pi*((DBH/2)^2))
        # calculate the wood density from the pure angiosperm / gymnosperm part in each plot
        if(nrow(angioTable)>0)
        {
            angioWoodDensity = sum(angioTable$basalArea * angioTable$WoodDensity)/sum(angioTable$basalArea)
        }else 
        {
            angioWoodDensity = NA
        }

        if(nrow(gymnoTable)>0)
        {
            gymnoWoodDensity = sum(gymnoTable$basalArea * gymnoTable$WoodDensity)/sum(gymnoTable$basalArea)
        }else 
        {
            gymnoWoodDensity = NA
        }

        # get the angiosperm/gymnosperm individual ratio in the plot
        a_and_g_Length = length(grep('a', perPlotDataFrame$A_G))+ length(grep('g', perPlotDataFrame$A_G))
        a_Frequency = length(grep('a', perPlotDataFrame$A_G))/a_and_g_Length
        g_Frequency = length(grep('g', perPlotDataFrame$A_G))/a_and_g_Length
        # return the information row out
        outputRow = data.frame(PLT=pn,
                                LAT=unique(perPlotDataFrame$LAT),
                                LON=unique(perPlotDataFrame$LON),
                                YEAR=unique(perPlotDataFrame$Year),
                                PlotArea = tphMean, 
                                IndividualNumber = individualNumber,
                                treeDensity = individualNumber*mean(perPlotDataFrame$TPH),
                                diameterMean = diameterMean,
                                SpeciesNumber = speciesNumber,
                                ShannonIndex = diversity(table(perPlotDataFrame$SPCD),'shannon'),
                                SimpsonIndex = diversity(table(perPlotDataFrame$SPCD),'simpson'),
                                InvSimpsonIndex = diversity(table(perPlotDataFrame$SPCD),'invsimpson'),
                                BasalArea = sum(perPlotDataFrame$basalArea),
                                WeightedWoodDensity = WeightedWoodDensity,
                                WoodDensitySD = woodDensitySD,
                                MinWD = minimumWD,
                                MaxWD = maximumWD,
                                GymnoWD = gymnoWoodDensity,
                                AngioWD = angioWoodDensity,
                                FiveQuantile = quantileResult[1],
                                NintyfiveQuantile = quantileResult[2],
                                angiospermRatio = a_Frequency,
                                gymnospermRatio = g_Frequency,
                                nonDefinedAG = 1-(a_Frequency+g_Frequency),
                                levelRatios)
        }
    
    print(paste("--- the metadata for plot ",pn," has been calculated ---"))
    return(outputRow)   
}

# in order to get the plot with time series, we just to check the plot with multiple years records
# for (pn in plotNames)
multipleYearTPHStatistic = function(pn)
{
    # get the per plot data 
    perPlotDF = fread(paste("****/PerPlotWithWoodDensity/",pn,"_GFBI_Data_With_Individual_Wood_Density.csv",sep=""))[,-1]
    # perPlotDF = na.omit(perPlotDF)
    # check is the plot has more than one years observations
    startDatFrame = data.frame()
    if (length(unique(perPlotDF$Year)) > 1)
    {
        # get the years in that plot
        yearsVector = unique(perPlotDF$Year)
        # lopp by year
        for (yr in yearsVector)
        {
            # subset the yearly data frame 
            yearlyDataFrame = perPlotDF[perPlotDF$Year == yr,]
            startDatFrame = rbind(startDatFrame,plotLevelStatFunc(pn,yearlyDataFrame))
        }
    }else
    {
        startDatFrame = rbind(startDatFrame,plotLevelStatFunc(pn,perPlotDF))
    } 
    return(startDatFrame) 
    
}

# do the paralle running
system.time(outputListRaw <- mclapply(plotNames,multipleYearTPHStatistic,mc.cores =28,mc.preschedule =F))
# rbind the output list into a data frame
plotLevelWoodDensityFull = rbindlist(outputListRaw) 

# we only used the latest observation of each plot
# write the latest years wood density information to local folder
outputList = outputListRaw
# clean the data by get the latest year observation
for (i in 1:length(outputList))
{
    perPlot = outputList[[i]]
    outputList[[i]] = perPlot[perPlot$YEAR == max(perPlot$YEAR),]
    print(i)
}

# rbind the output list into a data frame
plotLevelWoodDensity = rbindlist(outputList) 
# write it into local folder
write.csv(plotLevelWoodDensity,"Plot_Level_Wood_Density_Data_Frame_Metadata_With_Diversity_Latest_year.csv")


###*####*####*####*####*####*####*####*####*####*####*####*####*####*####*##
# STEP 3_2. do statistics of the sources of wood density values (species, genus,dataset)
###*####*####*####*####*####*####*####*####*####*####*####*####*####*####*##

library(data.table)
library(parallel)
library(vegan)
setwd("/Volumes/Scratch/Lidong_Mo/WoodDensityProject/")

# clean the data with the highest conservation on genus and species wood density 
latestYearObsTable = fread("Data/Plot_Level_Wood_Density_Data_Frame_Metadata_With_Diversity_Latest_year.csv")[,-1] %>%dplyr::filter(!is.na(LON))%>%dplyr::filter(!is.na(LAT))

# adaptation of the far east russian points
farEastRussiaPlots = latestYearObsTable %>% dplyr::filter(LON>136&LAT>58) %>% dplyr::mutate(LON= LON*(-1))#2014

otherPlots = latestYearObsTable %>% dplyr::filter(LON<=136|LAT<=58) #1186086
# rbind the modified tables together
latestYearObsTableNew = rbind(farEastRussiaPlots,otherPlots)
# add a new column in the data frame  with name ConservativeRatio by add the column genus and species
latestYearObsTableNew$ConservativeRatio = latestYearObsTableNew$genus+latestYearObsTableNew$species 
# make a hist plot
# 1183070/1188100 = 99.5%   0.4%
# here we use the ratio 100%
plotLevelWoodDensityConservative = latestYearObsTableNew[latestYearObsTableNew$ConservativeRatio >=0.75,]
# ACHTUNG! you may found more rows than the plot numbers as many plots has TPH duplicates
write.csv(plotLevelWoodDensityConservative,"Data/Plot_Level_Wood_Density_Data_Frame_GFBI_Metadata_Conservative_with_Diversity.csv")
