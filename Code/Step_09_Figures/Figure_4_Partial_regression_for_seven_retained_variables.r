
library(scales)
library(car)
library(data.table)
library(ggpubr)
library(corrplot)
library(gridExtra)
library(ggplot2)
library(raster)
library(purrr)
library(dplyr)
library(h2o)


# kset working directory 
setwd("~/WoodDensityMapingProject/")
source("Functions/sample.grid_from_GSIF_package.r")
##############################################################################################################################################
# Step 1 Partial regression for Angiosperm
##############################################################################################################################################
rawData = fread("Data/20231201_Wood_Density_Project_Merged_sampled_dataset_Diversity_Pixel_2023.csv")[,-1] %>% filter(angsprR>=1) #the csv could be downloaded from step 5 
processTable = data.frame()
# read the subsampled tables
for (ss in seq(100,10000,100))
{
    cleanedTrainTable = rawData %>% select(c("x","y","WdDnsty","CHELSA_Annual_Mean_Temperature","Richnss","cnRatio","Lai","SoilMoisture","Human_Disturbance","Fire_Frequency","ForestAge","angsprR","DimtrMn")) %>% na.omit()

    duplicateTable = cleanedTrainTable
    # tranform the data frame format lat lon into spatial lat long as spatial points
    coordinates(cleanedTrainTable) = ~ x + y
    # allocate the projection
    proj4string(cleanedTrainTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    cleanedTrainTable@data = duplicateTable
    # # set the 
    set.seed(ss)
    gridSubsampledPoints = sample.grid(cleanedTrainTable, cell.size = c(0.25,0.25), n = 1)
    # print(dim(gridSubsampledPoints$subset@data))
    gridSubsampledPoints = gridSubsampledPoints[[1]]
    # get the grid subsample data frame
    gridSubsampledTable = gridSubsampledPoints@data
    # gridSubsampledTable = cleanedTrainTable[sample(nrow(cleanedTrainTable), round(nrow(cleanedTrainTable)*0.05,0)), ]
    # define the vector of independent variables
    independentVars = c("CHELSA_Annual_Mean_Temperature","SoilMoisture","Human_Disturbance","Richnss","cnRatio","Fire_Frequency","Lai","angsprR","ForestAge","DimtrMn")

    # build a subset table
    subsetTable = subset(gridSubsampledTable,select=c(independentVars,"WdDnsty")) %>% na.omit()
    myScale = function(x) (x - mean(x, na.rm=T)) / (2*sd(x, na.rm=T))
    # scale the data
    # subsetTable$WdDnsty = myScale(sqrt(subsetTable$WdDnsty))
    boxTran = caret::BoxCoxTrans(subsetTable$WdDnsty)
    subsetTable$WdDnsty = myScale(predict(boxTran,subsetTable$WdDnsty))

    subsetTable$CHELSA_Annual_Mean_Temperature = myScale(subsetTable$CHELSA_Annual_Mean_Temperature)
    subsetTable$SoilMoisture= myScale(sqrt(subsetTable$SoilMoisture))
    subsetTable$Richnss = myScale(sqrt(subsetTable$Richnss))
    subsetTable$Lai = myScale(sqrt(subsetTable$Lai))
    subsetTable$cnRatio = myScale(subsetTable$cnRatio)
    subsetTable$Human_Disturbance = myScale(sqrt(subsetTable$Human_Disturbance))
    subsetTable$ForestAge= myScale(sqrt(subsetTable$ForestAge))
    subsetTable$DimtrMn= myScale(sqrt(subsetTable$DimtrMn))
    # dont need 
    # build the model
    # > library(ggfortify)
    # > autoplot( lmMulti)

    lmMulti = glm(WdDnsty~CHELSA_Annual_Mean_Temperature+Richnss+cnRatio+Lai+SoilMoisture+Human_Disturbance+angsprR+Fire_Frequency+ForestAge+DimtrMn,data = subsetTable) 
    # summary(lmMulti)
    # # autoplot( lmMulti)
    # summary(lmMulti)$coefficients[1,]

    
    processRow = as.data.frame(t(summary(lmMulti)$coefficients[-1,1]))
    print(processRow)
    processTable = rbind(processTable,processRow)
}

write.csv( processTable,"Data/PartialResults/Partial_Regression_Results_for_Mean_models_Angiosperm.csv")




##############################################################################################################################################
# Step 2 Partial regression for Gymnosperm 
##############################################################################################################################################

rawData = fread("Data/20231201_Wood_Density_Project_Merged_sampled_dataset_Diversity_Pixel_2023.csv")[,-1]  %>% filter(angsprR<=0)
processTable = data.frame()
# read the subsampled tables
for (ss in seq(100,10000,100))
{
    cleanedTrainTable = rawData %>% select(c("x","y","WdDnsty","CHELSA_Annual_Mean_Temperature","Richnss","cnRatio","Lai","SoilMoisture","Human_Disturbance","Fire_Frequency","ForestAge","DimtrMn")) %>% na.omit()
    duplicateTable = cleanedTrainTable
    # tranform the data frame format lat lon into spatial lat long as spatial points
    coordinates(cleanedTrainTable) = ~ x + y
    # allocate the projection
    proj4string(cleanedTrainTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    cleanedTrainTable@data = duplicateTable
    # # set the 
    set.seed(ss)
    gridSubsampledPoints = sample.grid(cleanedTrainTable, cell.size = c(0.25,0.25), n = 1)
    # print(dim(gridSubsampledPoints$subset@data))
    gridSubsampledPoints = gridSubsampledPoints[[1]]
    # get the grid subsample data frame
    gridSubsampledTable = gridSubsampledPoints@data
    # gridSubsampledTable = cleanedTrainTable[sample(nrow(cleanedTrainTable), round(nrow(cleanedTrainTable)*0.05,0)), ]
    # define the vector of independent variables
    independentVars = c("CHELSA_Annual_Mean_Temperature","SoilMoisture","Human_Disturbance","Richnss","cnRatio","Fire_Frequency","Lai","ForestAge","DimtrMn")

    # build a subset table
    subsetTable = subset(gridSubsampledTable,select=c(independentVars,"WdDnsty")) %>% na.omit()
    myScale = function(x) (x - mean(x, na.rm=T)) / (2*sd(x, na.rm=T))
    # scale the data
    # subsetTable$WdDnsty = myScale(sqrt(subsetTable$WdDnsty))
    boxTran = caret::BoxCoxTrans(subsetTable$WdDnsty)
    subsetTable$WdDnsty = myScale(predict(boxTran,subsetTable$WdDnsty))

    subsetTable$CHELSA_Annual_Mean_Temperature = myScale(subsetTable$CHELSA_Annual_Mean_Temperature)
    subsetTable$SoilMoisture= myScale(sqrt(subsetTable$SoilMoisture))
    subsetTable$Richnss = myScale(sqrt(subsetTable$Richnss))
    subsetTable$Lai = myScale(sqrt(subsetTable$Lai))
    subsetTable$cnRatio = myScale(subsetTable$cnRatio)
    subsetTable$Human_Disturbance = myScale(sqrt(subsetTable$Human_Disturbance))
    subsetTable$ForestAge= myScale(sqrt(subsetTable$ForestAge))
    subsetTable$DimtrMn= myScale(sqrt(subsetTable$DimtrMn))
    # dont need 
    # build the model
    # > library(ggfortify)
    # > autoplot( lmMulti)

    lmMulti = glm(WdDnsty~CHELSA_Annual_Mean_Temperature+Richnss+cnRatio+Lai+SoilMoisture+Human_Disturbance+Fire_Frequency+ForestAge+DimtrMn,data = subsetTable) 
    # # autoplot( lmMulti)
    # summary(lmMulti)$coefficients[1,]

    processRow = as.data.frame(t(summary(lmMulti)$coefficients[-1,1]))
    print(processRow)
    processTable = rbind(processTable,processRow)
}

write.csv( processTable,"Data/PartialResults/Partial_Regression_Results_for_Mean_models_Gymnosperm.csv")

library(scales)
library(ggeffects)
library(effects)
library(car)
library(data.table)
library(ggpubr)
library(corrplot)
library(gridExtra)
library(ggplot2)
library(raster)
library(purrr)
library(dplyr)
library(ggfortify)
library(caret)


rawData = fread("Data/20231201_Wood_Density_Project_Merged_sampled_dataset_Diversity_Pixel_2023.csv")[,-1] %>% mutate(Fire_Frequency = ifelse(Fire_Frequency <= 0, 0, 1)) %>% mutate(angsprR = ifelse(angsprR <= 0.5, 0, 1)) %>% filter(CanopyHeight >=2)

processTable = data.frame()
# read the subsampled tables
for (ss in seq(100,10000,100))
{
    cleanedTrainTable = rawData %>% select(c("x","y","WdDnsty","CHELSA_Annual_Mean_Temperature","Richnss","cnRatio","Lai","SoilMoisture","Human_Disturbance","Fire_Frequency","ForestAge","angsprR","DimtrMn")) %>% na.omit()

    duplicateTable = cleanedTrainTable
    # tranform the data frame format lat lon into spatial lat long as spatial points
    coordinates(cleanedTrainTable) = ~ x + y
    # allocate the projection
    proj4string(cleanedTrainTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    cleanedTrainTable@data = duplicateTable
    # # set the 
    set.seed(ss)
    gridSubsampledPoints = sample.grid(cleanedTrainTable, cell.size = c(0.25,0.25), n = 1)
    # print(dim(gridSubsampledPoints$subset@data))
    gridSubsampledPoints = gridSubsampledPoints[[1]]
    # get the grid subsample data frame
    gridSubsampledTable = gridSubsampledPoints@data
    # gridSubsampledTable = cleanedTrainTable[sample(nrow(cleanedTrainTable), round(nrow(cleanedTrainTable)*0.05,0)), ]
    # define the vector of independent variables
    independentVars = c("CHELSA_Annual_Mean_Temperature","SoilMoisture","Human_Disturbance","Richnss","cnRatio","Fire_Frequency","Lai","angsprR","ForestAge","DimtrMn")

    # build a subset table
    subsetTable = subset(gridSubsampledTable,select=c(independentVars,"WdDnsty")) %>% na.omit()
    myScale = function(x) (x - mean(x, na.rm=T)) / (2*sd(x, na.rm=T))
    # scale the data
    # subsetTable$WdDnsty = myScale(sqrt(subsetTable$WdDnsty))
    boxTran = caret::BoxCoxTrans(subsetTable$WdDnsty)
    subsetTable$WdDnsty = myScale(predict(boxTran,subsetTable$WdDnsty))

    subsetTable$CHELSA_Annual_Mean_Temperature = myScale(subsetTable$CHELSA_Annual_Mean_Temperature)
    subsetTable$SoilMoisture= myScale(sqrt(subsetTable$SoilMoisture))
    subsetTable$Richnss = myScale(sqrt(subsetTable$Richnss))
    subsetTable$Lai = myScale(sqrt(subsetTable$Lai))
    subsetTable$cnRatio = myScale(subsetTable$cnRatio)
    subsetTable$Human_Disturbance = myScale(sqrt(subsetTable$Human_Disturbance))
    subsetTable$ForestAge= myScale(sqrt(subsetTable$ForestAge))
    subsetTable$DimtrMn= myScale(sqrt(subsetTable$DimtrMn))
    # define the model 
    lmMulti = glm(WdDnsty~CHELSA_Annual_Mean_Temperature+Richnss+cnRatio+Lai+SoilMoisture+Human_Disturbance+angsprR+Fire_Frequency+ForestAge+DimtrMn,data = subsetTable) 
    # summary(lmMulti)
    # # autoplot( lmMulti)
    # summary(lmMulti)$coefficients[1,]

    
    processRow = as.data.frame(t(summary(lmMulti)$coefficients[-1,1]))
    print(processRow)
    processTable = rbind(processTable,processRow)
}

write.csv( processTable,"Data/PartialResults/Partial_Regression_Results_for_Mean_models.csv")


##############################################################################################################################################
# Step 2 making plots 
##############################################################################################################################################


partialCoeffMeanTable = fread("PartialResults/Partial_Regression_Results_for_Mean_models.csv")[,-1]

partialCoeffMeanFit = apply(partialCoeffMeanTable, 2, mean)
quantile2_5 = function(x) {quantile(x,0.025)}
quantile97_5 = function(x) {quantile(x,0.975)}

partialCoeffMeanLower = apply(partialCoeffMeanTable, 2, quantile2_5)
partialCoeffMeanUpper = apply(partialCoeffMeanTable , 2, quantile97_5)


# plot data frame 
meanMetricTable = data.frame(Variables = c("CHELSA_Annual_Mean_Temperature","Richnss","cnRatio","Lai","SoilMoisture","Human_Disturbance","angsprR","FireFrequency","ForestAge","DimtrMn"),
                             Fit = partialCoeffMeanFit,
                             Lower = partialCoeffMeanLower,
                             Upper = partialCoeffMeanUpper,
                             Model = "Wood density")


finalTable = meanMetricTable


p1 = ggplot(finalTable, aes(x = Variables, y = Fit)) + 
            geom_pointrange(aes(ymax = Upper, ymin = Lower), color = "darkblue") +
            #   geom_text(aes(label = Fit), nudge_x = 0.15) + 
            scale_x_discrete(limits = c("FireFrequency","Richnss","DimtrMn","ForestAge","cnRatio","Lai","SoilMoisture","Human_Disturbance","CHELSA_Annual_Mean_Temperature","angsprR"),
                              labels = c("CHELSA_Annual_Mean_Temperature" = "Mean annual temperature","SoilMoisture" = "Soil moisture","Richnss" = "Biodiversity","cnRatio" = "Soil C:N ratio","Lai" = "Leaf area index","Human_Disturbance" = "Human modification", "angsprR" = "Angiosperm ratio","FireFrequency" = "Fire frequency","ForestAge" = "Forest age","DimtrMn" = "DBH"))+
            geom_hline(yintercept = 0, color = "red") + 
            theme_bw() +
            ylab("Partial regression coefficients") + 
            ylim(-0.3,0.8)+
            xlab(NULL)+
            coord_flip()+
            theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.text=element_text(size=10))
         #+
            #facet_grid(~Model)


##############################################################################################################################################
# Step 3 random forest based  
##############################################################################################################################################


library(data.table)
library(parallel)
library(ggplot2)
library(raster)
library(purrr)
library(dplyr)
library(h2o)



rawData = fread("Data/20231201_Wood_Density_Project_Merged_sampled_dataset_Diversity_Pixel_2023.csv")[,-1] %>% filter(CanopyHeight >=2)

portNumber=54321
# read the subsampled tables
for (ss in 1:100)
{
    cleanedTrainTable = rawData %>% dplyr::select(c("x","y","WdDnsty","CHELSA_Annual_Mean_Temperature","Richnss","cnRatio","Lai","SoilMoisture","Human_Disturbance","Fire_Frequency","ForestAge","angsprR","DimtrMn")) %>% na.omit()

    duplicateTable = cleanedTrainTable
    # tranform the data frame format lat lon into spatial lat long as spatial points
    coordinates(cleanedTrainTable) = ~ x + y
    # allocate the projection
    proj4string(cleanedTrainTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    cleanedTrainTable@data = duplicateTable
    # # set the 
    set.seed(ss)
    gridSubsampledPoints = sample.grid(cleanedTrainTable, cell.size = c(0.25,0.25), n = 1)
    # print(dim(gridSubsampledPoints$subset@data))
    gridSubsampledPoints = gridSubsampledPoints[[1]]
    # get the grid subsample data frame
    gridSubsampledTable = gridSubsampledPoints@data
    # define the vector of independent variables
    trainVaribles = c("CHELSA_Annual_Mean_Temperature","SoilMoisture","Human_Disturbance","Richnss","cnRatio","Fire_Frequency","Lai","angsprR","ForestAge","DimtrMn")
    localH2O = h2o.init(ip = 'localhost', port = portNumber, nthreads= 5,max_mem_size = '8g') #portNumber[portOrder]
    # define the train data into the H@O format
    trainData.hex = as.h2o(gridSubsampledTable, destination_frame = "trainData.hex")
    # finalTable = subsetTails

    # load the parameter table

    fullRandomForestModel = h2o.randomForest(x=trainVaribles, y="WdDnsty",
                                             training_frame=trainData.hex,
                                             ntrees=250,
                                             nfolds = 10,
                                             seed=ss,
                                             max_depth=20,
                                             mtries=3,
                                             min_rows = 20,
                                             sample_rate= 0.632,
                                             keep_cross_validation_predictions=T)
    variableImportance = h2o.varimp(fullRandomForestModel)
    write.csv(variableImportance,paste("Data/RandomForestImportance/Variable_Importance_for_seed_",ss,".csv",sep=""))
    h2o.shutdown(prompt = F)
    portNumber = portNumber+2
}

# merge all the 20 importance table
tableNameList = paste("Data/RandomForestImportance/Variable_Importance_for_seed_",1:100,".csv",sep="")
#
tableBinded = lapply(tableNameList,fread) %>% rbindlist()
tableAggred = aggregate(x = tableBinded$percentage, by = list(tableBinded$variable), FUN = mean) %>% dplyr::rename(variable = Group.1,percentage=x) %>% mutate(percentage = percentage/sum(percentage))

p2 = ggplot(tableAggred,aes(x=variable,y=percentage)) +
        scale_x_discrete(limits = c("Fire_Frequency","Richnss","DimtrMn","ForestAge","cnRatio","Lai","SoilMoisture","Human_Disturbance","CHELSA_Annual_Mean_Temperature","angsprR"),
                              labels = c("CHELSA_Annual_Mean_Temperature" = "Mean annual temperature","Richnss" = "Biodiversity","cnRatio" = "Soil C:N ratio","Lai" = "Leaf area index","SoilMoisture" = "Soil moisture","Human_Disturbance" = "Human modification","DimtrMn" = "DBH",  "angsprR" = "Angiosperm ratio","Fire_Frequency" = "Fire frequency","ForestAge" = "Forest age"))+
        geom_col(fill="cornflowerblue") + 
        coord_flip()+
        ylab("Importance (Percentage)")+
        ylim(0,0.55)+
        xlab(NULL)+
        theme_bw()+
        theme(axis.text=element_text(size=10))
        




##############################################################################################################################################
# Step 4 random forest based for angiosperm/gmynosperm alone
##############################################################################################################################################



library(scales)
library(ggeffects)
library(effects)
library(car)
library(data.table)
library(ggpubr)
library(corrplot)
library(gridExtra)
library(ggplot2)
library(raster)
library(purrr)
library(dplyr)
library(ggfortify)
library(caret)
library(h2o)

source("Functions/sample.grid_from_GSIF_package.r")


rawData = fread("Data/20231201_Wood_Density_Project_Merged_sampled_dataset_Diversity_Pixel_2023.csv")[,-1] %>% filter(angsprR<=0) %>% filter(CanopyHeight >=2)

portNumber=54321
# read the subsampled tables
for (ss in 1:100)
{
    cleanedTrainTable = rawData %>%   dplyr::select(c("x","y","WdDnsty","CHELSA_Annual_Mean_Temperature","Richnss","cnRatio","Lai","SoilMoisture","Human_Disturbance","Fire_Frequency","ForestAge","DimtrMn")) %>% na.omit()

    duplicateTable = cleanedTrainTable
    # tranform the data frame format lat lon into spatial lat long as spatial points
    coordinates(cleanedTrainTable) = ~ x + y
    # allocate the projection
    proj4string(cleanedTrainTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    cleanedTrainTable@data = duplicateTable
    # # set the 
    set.seed(ss)
    gridSubsampledPoints = sample.grid(cleanedTrainTable, cell.size = c(0.25,0.25), n = 1)
    # print(dim(gridSubsampledPoints$subset@data))
    gridSubsampledPoints = gridSubsampledPoints[[1]]
    # get the grid subsample data frame
    gridSubsampledTable = gridSubsampledPoints@data
    # gridSubsampledTable = cleanedTrainTable[sample(nrow(cleanedTrainTable), round(nrow(cleanedTrainTable)*0.05,0)), ]
    # define the vector of independent variables
    trainVaribles = c("CHELSA_Annual_Mean_Temperature","Richnss","cnRatio","Lai","SoilMoisture","Human_Disturbance","Fire_Frequency","ForestAge","DimtrMn")
    localH2O = h2o.init(ip = 'localhost', port = portNumber, nthreads= 12,max_mem_size = '60g') #portNumber[portOrder]
    # define the train data into the H@O format
    trainData.hex = as.h2o(gridSubsampledTable, destination_frame = "trainData.hex")
    # finalTable = subsetTails

    # load the parameter table

    fullRandomForestModel = h2o.randomForest(x=trainVaribles, y="WdDnsty",
                                             training_frame=trainData.hex,
                                             ntrees=250,
                                             nfolds = 10,
                                             seed=ss,
                                             max_depth=20,
                                             mtries=3,
                                             min_rows = 20,
                                             sample_rate= 0.632,
                                             keep_cross_validation_predictions=T)
    variableImportance = h2o.varimp(fullRandomForestModel)
    write.csv(variableImportance,paste("Data/RandomForestImportance/Variable_Importance_for_gymnosperm_seed_",ss,".csv",sep=""))
    h2o.shutdown(prompt = F)
    portNumber = portNumber+2
}

# merge all the 20 importance table
tableNameList3 = paste("Data/RandomForestImportance/Variable_Importance_for_gymnosperm_seed_",1:100,".csv",sep="")
#
tableBinded3 = lapply(tableNameList3,fread) %>% rbindlist()
tableAggred3 = aggregate(x = tableBinded3$percentage, by = list(tableBinded3$variable), FUN = mean) %>% dplyr::rename(variable = Group.1,percentage=x) %>% mutate(Class="Gymnosperm") %>% mutate(percentage = percentage/sum(percentage))

p3 = ggplot(tableAggred3,aes(x=variable,y=percentage)) +
        scale_x_discrete(limits = c("Fire_Frequency","Richnss","DimtrMn","ForestAge","cnRatio","Lai","SoilMoisture","Human_Disturbance","CHELSA_Annual_Mean_Temperature"),labels = c("CHELSA_Annual_Mean_Temperature" = "Mean annual temperature","Richnss" = "Biodiversity","cnRatio" = "Soil C:N ratio","Lai" = "Leaf area index","SoilMoisture" = "Soil moisture","Human_Disturbance" = "Human modification","DimtrMn" = "DBH", "Fire_Frequency" = "Fire frequency","ForestAge" = "Forest age"))+
        geom_col(fill="cornflowerblue") + 
        coord_flip()+
        ylab("Importance (Percentage)")+
        ylim(0,0.55)+
        xlab(NULL)+
        theme_bw()+
        theme(axis.text=element_text(size=10))
        # facet_wrap( ~Class, scales = "free_y", nrow = 1) 


rawData = fread("Data/20231201_Wood_Density_Project_Merged_sampled_dataset_Diversity_Pixel_2023.csv")[,-1]  %>% filter(angsprR>=1) %>% filter(CanopyHeight >=2)

portNumber=54321
# read the subsampled tables
for (ss in 1:100)
{
    cleanedTrainTable = rawData %>% select(c("x","y","WdDnsty","CHELSA_Annual_Mean_Temperature","Richnss","cnRatio","Lai","SoilMoisture","Human_Disturbance","Fire_Frequency","ForestAge","DimtrMn")) %>% na.omit()

    duplicateTable = cleanedTrainTable
    # tranform the data frame format lat lon into spatial lat long as spatial points
    coordinates(cleanedTrainTable) = ~ x + y
    # allocate the projection
    proj4string(cleanedTrainTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    cleanedTrainTable@data = duplicateTable
    # # set the 
    set.seed(ss)
    gridSubsampledPoints = sample.grid(cleanedTrainTable, cell.size = c(0.25,0.25), n = 1)
    # print(dim(gridSubsampledPoints$subset@data))
    gridSubsampledPoints = gridSubsampledPoints[[1]]
    # get the grid subsample data frame
    gridSubsampledTable = gridSubsampledPoints@data
    # gridSubsampledTable = cleanedTrainTable[sample(nrow(cleanedTrainTable), round(nrow(cleanedTrainTable)*0.05,0)), ]
    # define the vector of independent variables
    trainVaribles = c("CHELSA_Annual_Mean_Temperature","Richnss","cnRatio","Lai","SoilMoisture","Human_Disturbance","Fire_Frequency","ForestAge","DimtrMn")
    localH2O = h2o.init(ip = 'localhost', port = portNumber, nthreads= 12,max_mem_size = '60g') #portNumber[portOrder]
    # define the train data into the H@O format
    trainData.hex = as.h2o(gridSubsampledTable, destination_frame = "trainData.hex")
    # finalTable = subsetTails

    # load the parameter table

    fullRandomForestModel = h2o.randomForest(x=trainVaribles, y="WdDnsty",
                                             training_frame=trainData.hex,
                                             ntrees=250,
                                             nfolds = 10,
                                             seed=ss,
                                             max_depth=20,
                                             mtries=3,
                                             min_rows = 20,
                                             sample_rate= 0.632,
                                             keep_cross_validation_predictions=T)
    variableImportance = h2o.varimp(fullRandomForestModel)
    write.csv(variableImportance,paste("Data/RandomForestImportance/Variable_Importance_for_angiosperm_seed_",ss,".csv",sep=""))
    h2o.shutdown(prompt = F)
    portNumber = portNumber+2
}

# merge all the 20 importance table
tableNameList4 = paste("Data/RandomForestImportance/Variable_Importance_for_angiosperm_seed_",1:100,".csv",sep="")
#
tableBinded4 = lapply(tableNameList4,fread) %>% rbindlist()
tableAggred4 = aggregate(x = tableBinded4$percentage, by = list(tableBinded4$variable), FUN = mean) %>% dplyr::rename(variable = Group.1,percentage=x) %>% mutate(Class="Angiosperm") %>% mutate(percentage = percentage/sum(percentage))




p4 = ggplot(tableAggred4,aes(x=variable,y=percentage)) +
        scale_x_discrete(limits = c("Fire_Frequency","Richnss","DimtrMn","ForestAge","cnRatio","Lai","SoilMoisture","Human_Disturbance","CHELSA_Annual_Mean_Temperature"),
                              labels = c("CHELSA_Annual_Mean_Temperature" = "Mean annual temperature","Richnss" = "Biodiversity","cnRatio" = "Soil C:N ratio","Lai" = "Leaf area index","SoilMoisture" = "Soil moisture","Human_Disturbance" = "Human modification","DimtrMn" = "DBH","Fire_Frequency" = "Fire frequency","ForestAge" = "Forest age"))+
        geom_col(fill="cornflowerblue") + 
        coord_flip()+
        ylab("Importance (Percentage)")+
        ylim(0,0.55)+
        xlab(NULL)+
        theme_bw()+
        theme(axis.text=element_text(size=10)) #+
        # facet_wrap( ~Class, scales = "free_y", nrow = 1) 


##############################################################################################################################################
# Step 6 partial regression based for angiosperm/gmynosperm alone
##############################################################################################################################################

library(data.table)
library(parallel)
library(ggplot2)
library(purrr)
library(dplyr)
library(rgdal)



partialCoeffMeanTable = fread("Data/PartialResults/Partial_Regression_Results_for_Mean_models_Angiosperm.csv")[,-1]
partialCoeffMeanFit = apply(partialCoeffMeanTable, 2, mean)
quantile2_5 = function(x) {quantile(x,0.025)}
quantile97_5 = function(x) {quantile(x,0.975)}

partialCoeffMeanLower = apply(partialCoeffMeanTable, 2, quantile2_5)
partialCoeffMeanUpper = apply(partialCoeffMeanTable , 2, quantile97_5)
# plot data frame 
meanMetricTablePartialAng = data.frame(Variables = c("CHELSA_Annual_Mean_Temperature","Richnss","cnRatio","Lai","SoilMoisture","Human_Disturbance","FireFrequency","ForestAge","DimtrMn"),
                             Fit = partialCoeffMeanFit,
                             Lower = partialCoeffMeanLower,
                             Upper = partialCoeffMeanUpper,
                             Model = "Wood density")

finalTablePartialAng = meanMetricTablePartialAng


p5 = ggplot(finalTablePartialAng, aes(x = Variables, y = Fit)) + 
            geom_pointrange(aes(ymax = Upper, ymin = Lower), color = "darkblue") +
            #   geom_text(aes(label = Fit), nudge_x = 0.15) + 
            scale_x_discrete(limits =   c("FireFrequency","Richnss","DimtrMn","ForestAge","cnRatio","Lai","SoilMoisture","Human_Disturbance","CHELSA_Annual_Mean_Temperature"),
                              labels = c("CHELSA_Annual_Mean_Temperature" = "Mean annual temperature","SoilMoisture" = "Soil moisture","Richnss" = "Biodiversity","cnRatio" = "Soil C:N ratio","Lai" = "Leaf area index","Human_Disturbance" = "Human modification","DimtrMn" = "DBH","FireFrequency" = "Fire frequency","ForestAge" = "Forest age"))+
            geom_hline(yintercept = 0, color = "red") + 
            theme_bw() +
            ylab("Partial regression coefficients") + 
            ylim(-0.3,0.8)+
            xlab(NULL)+
            coord_flip()+
            theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.text.x =element_text(size=10))

library(data.table)
library(parallel)
library(ggplot2)
library(purrr)
library(dplyr)
library(rgdal)


partialCoeffMeanTable = fread("Data/PartialResults/Partial_Regression_Results_for_Mean_models_Gymnosperm.csv")[,-1]
partialCoeffMeanFit = apply(partialCoeffMeanTable, 2, mean)

partialCoeffMeanFit = apply(partialCoeffMeanTable, 2, mean)
quantile2_5 = function(x) {quantile(x,0.025)}
quantile97_5 = function(x) {quantile(x,0.975)}

partialCoeffMeanLower = apply(partialCoeffMeanTable, 2, quantile2_5)
partialCoeffMeanUpper = apply(partialCoeffMeanTable , 2, quantile97_5)

# plot data frame 
meanMetricTablePartialGym = data.frame(Variables = c("CHELSA_Annual_Mean_Temperature","Richnss","cnRatio","Lai","SoilMoisture","Human_Disturbance","FireFrequency","ForestAge","DimtrMn"),
                             Fit = partialCoeffMeanFit,
                             Lower = partialCoeffMeanLower,
                             Upper = partialCoeffMeanUpper,
                             Model = "Wood density")

# variMetricTable = data.frame(Variables = c("CHELSA_Annual_Mean_Temperature","CHELSA_Annual_Precipitation","Human_Disturbance","ShnnnIn","cnRatio","Fire_Frequency","Lai"),
#                              Fit = partialCoeffVariFit,
#                              Lower = partialCoeffVariLower,
#                              Upper = partialCoeffVariUpper,
#                              Model = "M2")

finalTablePartialGym = meanMetricTablePartialGym


p6 = ggplot(finalTablePartialGym, aes(x = Variables, y = Fit)) + 
            geom_pointrange(aes(ymax = Upper, ymin = Lower), color = "darkblue") +
            #   geom_text(aes(label = Fit), nudge_x = 0.15) + 
            scale_x_discrete(limits =  c("FireFrequency","Richnss","DimtrMn","ForestAge","cnRatio","Lai","SoilMoisture","Human_Disturbance","CHELSA_Annual_Mean_Temperature"),
                              labels = c("CHELSA_Annual_Mean_Temperature" = "Mean annual temperature","SoilMoisture" = "Soil moisture","Richnss" = "Biodiversity","cnRatio" = "Soil C:N ratio","Lai" = "Leaf area index","Human_Disturbance" = "Human modification","DimtrMn" = "DBH","FireFrequency" = "Fire frequency","ForestAge" = "Forest age"))+
            geom_hline(yintercept = 0, color = "red") + 
            theme_bw() +
            ylab("Partial regression coefficients") + 
            ylim(-0.3,0.8)+
            xlab(NULL)+
            coord_flip()+
            theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.text.x =element_text(size=10))

pdf("Plots/Fig_4_Partial_Regression_Top_8_Variables_With_Human_Modification_and_taxonamy_class_2024_Revision02.pdf",width =8, height=9)
library("ggpubr")
ggarrange(p2,p1,p4,p5,p3,p6,ncol = 2, nrow = 3,align="hv",widths = c(1,1),heights = c(1.125,1,1),labels = c("a","b","c","d","e","f") )
dev.off()   

