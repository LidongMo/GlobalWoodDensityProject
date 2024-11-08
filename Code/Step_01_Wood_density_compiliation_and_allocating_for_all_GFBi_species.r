library(data.table)
library(raster)
library(stringr)
library(parallel)
# library(Taxonstand) #TPL
library(TNRS) #TNRS
library(dplyr)

# set the working directory 
setwd("～/WoodDensityProject/")

###*####*####*####*####*####*####*####*####*####*####*####*####*####*####*##
# STEP 1_1. Compile all the data sets we collected
###*####*####*####*####*####*####*####*####*####*####*####*####*####*####*##

# list the wood density reference data frames
fileList = list.files(path="Data/ReferenceWoodDensityCSVs",pattern= ".csv")
# write a smmall function which only get the columns for Family, SpeciesName,Wood density
columnSubsetFunc = function(filename)
{
    # load the single data frame
    singleDF = fread(paste("Data/ReferenceWoodDensityCSVs/",filename,sep=""))
    # subseting the data frame
    returnDataFrame = singleDF[,c("Family","SpeciesName","WoodDensity")]
    returnDataFrame$DataSource = gsub(".csv", "", filename)
    print(nrow(returnDataFrame))
    # retrun it
    return(returnDataFrame)
}
# lapply to do the subseting process
subsetedDataFrameList = lapply (fileList, columnSubsetFunc)
# rbind them into on data frame
rbindedDataFrame = rbindlist(subsetedDataFrameList)
# remove the "\xa0" in the species names 
rbindedDataFrame$SpeciesName = gsub("[^[:alnum:] ]", "", rbindedDataFrame$SpeciesName)
# valueDataFrame = rbindedDataFrame[!is.na(rbindedDataFrame$Family),]
# NADataFrame = rbindedDataFrame[is.na(rbindedDataFrame$Family),]
# kick out the binomial names longer than two words
# write a function which works by row 


#  This part is using the TNRS (The TNRS database and R package) to do the binomial correction and family name allocation.
# get the unique name in the rbindedDataFrame
uniqueSpeciesTable = data.frame (SpeciesID = c(1:length(unique(rbindedDataFrame$SpeciesName))),Family=NA,SpeciesName = unique(rbindedDataFrame$SpeciesName)) %>% mutate_if(is.factor, as.character)
rowNameFormatingTNRS_Func = function(ro)
# for (ro in 1:nrow(uniqueSpeciesTable))
{
    # get the data for per row
    perRow = uniqueSpeciesTable[ro,]
    print(ro)
    if (str_count(perRow$SpeciesName,"\xa0")==1)
    {
        perRow$SpeciesName = gsub("\xa0"," ",perRow$SpeciesName)
    }
    # check the the space numbers in the speciesName
    if (str_count(perRow$SpeciesName," ") >1)
    {
        # find the second space position
        secondSpacePostion = as.data.frame(str_locate_all(perRow$SpeciesName," "))$start[2]
        # subset the speciename form the start to the second space
        perRow$SpeciesName = substr(perRow$SpeciesName,start=1,stop=(secondSpacePostion-1))
    }
    # check the family name by TNRS
    allocatedInfo = TNRS(perRow$SpeciesName)
    if (nchar(allocatedInfo$Name_matched_accepted_family)>4&!is.na(allocatedInfo$Name_matched_accepted_family)&allocatedInfo$Name_matched_accepted_family!=''&!is.null(allocatedInfo$Name_matched_accepted_family))
    {
        perRow$Family = allocatedInfo$Name_matched_accepted_family
    }else
    {
        perRow$Family = NA
    }
    # allocate the information to the perRow
    perRow$Genus = allocatedInfo$Genus_matched
    perRow$Species = allocatedInfo$Specific_epithet_matched
    perRow$Binomial = paste(perRow$Genus,perRow$Species,sep=" ")
    # show the output
    # print(perRow)
    return(perRow)
    # rbindedDataFrame[ro,] = perRow
}
# do parallel running
system.time(outputList <- mclapply(c(1:nrow(uniqueSpeciesTable)),rowNameFormatingTNRS_Func,mc.cores=32,mc.preschedule = F))
# rbindlist 
nameAllocatedDataFrame = rbindlist(outputList) 
# merge this table with the rbindedDataFrame
mergedWoodDensityTable = merge(rbindedDataFrame %>% select(SpeciesName,WoodDensity), nameAllocatedDataFrame %>% select(Family,SpeciesName,Genus,Species,Binomial), by="SpeciesName") %>% filter(WoodDensity <2)

write.csv(mergedWoodDensityTable,"Data/Resources/Wood_density_Table_with_Family_allocated_TNRS_2024.csv")


###*####*####*####*####*####*####*####*####*####*####*####*####*####*####*##
# STEP 1_2. allocate family name to GFBi species list 
###*####*####*####*####*####*####*####*####*####*####*####*####*####*####*##
# load the GFBI table first 
GFBISpeciesTableRaw = fread("Data/GFBI_raw_and_tnrs_corrected_binomial_df.csv",head=T)[,-1] %>% rename(SpeciesName = RawName,Genus = CorrectedGenus,Species = CorrectedSpecies) %>% mutate(SpeciesName = gsub("[^[:alnum:] ]", "", GFBISpeciesTableRaw$SpeciesName))
# correct the name and allocate the family name use the TNRS database and the function we used above 
rowNameFormatingTNRS_Func = function(ro)
# for (ro in 1:nrow(GFBISpeciesTableRaw))
{
    # get the data for per row
    perRow = GFBISpeciesTableRaw[ro,]
    print(ro)
    if (str_count(perRow$SpeciesName,"\xa0")==1)
    {
        perRow$SpeciesName = gsub("\xa0"," ",perRow$SpeciesName)
    }
    # check the the space numbers in the speciesName
    if (str_count(perRow$SpeciesName," ") >1)
    {
        # find the second space position
        secondSpacePostion = as.data.frame(str_locate_all(perRow$SpeciesName," "))$start[2]
        # subset the speciename form the start to the second space
        perRow$SpeciesName = substr(perRow$SpeciesName,start=1,stop=(secondSpacePostion-1))
    }
    # check the family name by TNRS
    allocatedInfo = TNRS(perRow$SpeciesName)
    if (nchar(allocatedInfo$Name_matched_accepted_family)>4&!is.na(allocatedInfo$Name_matched_accepted_family)&allocatedInfo$Name_matched_accepted_family!=''&!is.null(allocatedInfo$Name_matched_accepted_family))
    {
        perRow$Family = allocatedInfo$Name_matched_accepted_family
    }else
    {
        perRow$Family = NA
    }
    # allocate the information to the perRow
    perRow$Genus = allocatedInfo$Genus_matched
    perRow$Species = allocatedInfo$Specific_epithet_matched
    perRow$Binomial = paste(perRow$Genus,perRow$Species,sep=" ")
    # show the output
    print(perRow)
    return(perRow)
    # rbindedDataFrame[ro,] = perRow
}

system.time(outputList <- mclapply(1:nrow(GFBISpeciesTableRaw),rowNameFormatingTNRS_Func,mc.preschedule=F,mc.cores=12))
GFBISpeciesTable = rbindlist(outputList)
write.csv(GFBISpeciesTable,"Data/Resources/GFBI_raw_and_tnrs_corrected_binomial_df_with_family_name.csv") # this will be mannually corrected for some missing recordings


###*####*####*####*####*####*####*####*####*####*####*####*####*####*####*##
# STEP 1_3. allocate wood density to GFBI species list referes to wood denisty data
###*####*####*####*####*####*####*####*####*####*####*####*####*####*####*##
woodDensityTable = fread("Data/Resources/Wood_density_Table_with_Family_allocated_TNRS_Mannul_corr.csv")[,-1] %>% filter(!is.na(Family))
# > dim(woodDensityTable )
# [1] 80115     6
# > length(unique( woodDensityTable$Binomial))
# [1] 10703
# > length(unique( woodDensityTable$Genus))
# [1] 2026
# > length(unique( woodDensityTable$Family))
# [1] 212

GFBISpeciesTable = fread("Data/Resources/GFBI_raw_and_tnrs_corrected_binomial_df_with_family_name_Mannul_corr.csv",head=T)[,-1]
# define the wood density mathing function
woodDensityMatchingFunc = function(ro)
{
    perRow = GFBISpeciesTable[ro,]
    # allocate the gs and sp
    gs = perRow$Genus
    sp = perRow$Species
    fm = perRow$Family
    # get the binomial name
    binomialSp = perRow$Binomial
    print(ro)
    # check fisrt the binomial in the reference table
    if (binomialSp %in% woodDensityTable$Binomial)
    {
        # subet the wood density table
        subWoodDensityTable = woodDensityTable[woodDensityTable$Binomial == binomialSp,]
        # get the mean of the wd
        speciesWoodDensity = mean(subWoodDensityTable$WoodDensity,na.rm=T)
        return(data.frame(RawName = perRow$SpeciesName,Family=perRow$Family,Genus=perRow$Genus,Species=perRow$Species,Binomial = perRow$Binomial,WoodDensity=speciesWoodDensity,WD_Level ="species"))
    } else
    if (gs %in% woodDensityTable$Genus)
    {
        genusWoodDensityTable = woodDensityTable[woodDensityTable$Genus == gs,]
        speciesWoodDensity = mean(genusWoodDensityTable$WoodDensity,na.rm=T)
        return(data.frame(RawName = perRow$SpeciesName,Family=perRow$Family,Genus=perRow$Genus,Species=perRow$Species,Binomial = perRow$Binomial,WoodDensity=speciesWoodDensity,WD_Level ="genus"))
    }else
    # this part is useless
    # if(fm %in% woodDensityTable$family &is.na(fm))
    # {
    #     familyWoodDensityTable = woodDensityTable[woodDensityTable$family == fm,]
    #     speciesWoodDensity = mean(familyWoodDensityTable$wd)
    #     return(data.frame(WoodDensity=speciesWoodDensity,level ="family"))
    # }else
    {
        return(data.frame(RawName = perRow$SpeciesName,Family=perRow$Family,Genus=perRow$Genus,Species=perRow$Species,Binomial = perRow$Binomial,WoodDensity=mean(woodDensityTable$WoodDensity,na.rm=T),WD_Level ="dataset"))
    }
    
}
# do parallel running 
system.time(outputList <- mclapply(1:nrow(GFBISpeciesTable),woodDensityMatchingFunc,mc.preschedule=F,mc.cores=24))
GFBISpeciesWoodDensityReference = rbindlist(outputList)
write.csv(GFBISpeciesWoodDensityReference,"Data/Resources/GFBI_Wood_Density_Reference_Table.csv")


###*####*####*####*####*####*####*####*####*####*####*####*####*####*####*##
# STEP 1_4. add the angiosperm and gymnosperm information into the wood density table
###*####*####*####*####*####*####*####*####*####*####*####*####*####*####*##
GFBISpeciesTable = fread("Data/Resources/GFBI_Wood_Density_Reference_Table.csv",head=T)[,-1]
taxonomyTable = fread("Data/Resources/Species_taxonomy_reference_list_TNRS_Updated.csv")[,-1]

# define the wood density mathing function
AG_MatchingFunc = function(ro)
{
    perRow = GFBISpeciesTable[ro,]
    # allocate the gs and sp
    gs = perRow$Genus
    sp = perRow$Species
    fm = perRow$Family
    # generate the binomial name
    binomialSp = perRow$Binomial
    print(ro)
    # check fisrt the binomial in the reference table
    if (binomialSp %in% taxonomyTable$Binomial&nchar(binomialSp) >0)
    {
        # subet the wood density table
        subTaxonTable = taxonomyTable[taxonomyTable$Binomial == binomialSp,]
        # get the taxonomic information at species level
        speciesTaxon = subTaxonTable$angio_gymno[1]
        return(data.frame(perRow,A_G=speciesTaxon,AG_Level ="species"))
    } else
    if (gs %in% taxonomyTable$Genus &nchar(gs) >0)
    {
        # get the taxonomic information at genus level
        genusTaxonTable = taxonomyTable[taxonomyTable$Genus == gs,]
        counts = table(genusTaxonTable$angio_gymno)
        speciesTaxon = names(counts)[which.max(counts)]
        return(data.frame(perRow,A_G=speciesTaxon,AG_Level ="genus"))
    }else
    if(fm %in% taxonomyTable$Family &!is.na(fm))
    {
        # get the taxonomic information at family level
        familyTaxonTable = taxonomyTable[taxonomyTable$Family == fm,]
        counts = table(familyTaxonTable$angio_gymno)
        speciesTaxon = names(counts)[which.max(counts)]

        return(data.frame(perRow,A_G=speciesTaxon,AG_Level ="family"))
    }else
    {
        return(data.frame(perRow,A_G=NA,AG_Level ="dataset"))
    }
    
}

# do parallel running 
system.time(outputList <- mclapply(1:nrow(GFBISpeciesTable),AG_MatchingFunc,mc.preschedule=F,mc.cores=4))
GFBISpeciesTaxonReference = rbindlist(outputList)
write.csv(GFBISpeciesTaxonReference,"Data/Resources/GFBI_Wood_Density_Reference_Table_with_A_G.csv")