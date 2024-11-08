
library(scales)
library(ggeffects)
library(knitr)
library(xtable)
library(printr)
library(effects)
library(car)
library(AER)
library(broom)
library(raster)
library(data.table)
library(dplyr)
library(ggpubr)
library(corrplot)
library(gridExtra)
library(RColorBrewer)
library(sf)
# kset working directory 
setwd("~/Desktop/ETH_Study/WoodDensityMapingProject/")
# load the data
# timeSeriesDataTable = fread("CovariatesExtractionFolder/20201012_Wood_Density_Project_Merged_sampled_dataset_Diversity_TimeSeries.csv")[,-1] %>%
#                              filter(nnDfnAG<0.5) %>% 
#                             filter(CHELSA_Annual_Mean_Temperature<100) %>%
#                             filter(Human_Disturbance>0.2)
# # check how many plots that have the full angisperm and gymnosperm defined
# #  simple linear model of community composition and wood density
# lm(WghtdWD~Human_Disturbance,timeSeriesDataTable)
# plot(WghtdWD~Human_Disturbance,data=timeSeriesDataTable,pch=17)
# abline(lm(timeSeriesDataTable$WghtdWD~timeSeriesDataTable$Human_Disturbance))

fullPlotLevelData = fread("Plot_Level_Wood_Density_Data_Frame_GFBI_Metadata_Conservative_with_Diversity.csv") %>% filter(nonDefinedAG<=0.25) #filter these plots with more than 75% not recognized as angiosperm or gymnosperm
                                                                        #  filter(angiospermRatio>0.02) %>%
                                                                        #  filter(angiospermRatio<0.98)
# add a new column which with rounded angiosperm ratio
fullPlotLevelData$TrueAgiosperm = fullPlotLevelData$angiospermRatio/(fullPlotLevelData$angiospermRatio+fullPlotLevelData$gymnospermRatio)
fullPlotLevelData$RoundedAgiosperm = round(fullPlotLevelData$TrueAgiosperm,1)%>%as.factor()
# define the colour pallete
pal= colorRampPalette(brewer.pal(5, "Greens"))
# make boxplot based on the rounded angiosperm ratio
boxplotFig = ggplot(fullPlotLevelData, aes(x=RoundedAgiosperm,y=WeightedWoodDensity,fill=RoundedAgiosperm))+
             geom_boxplot(outlier.size=-1)+
             ylim(0.15,1.3)+
             scale_fill_manual(values =pal(11))+
             labs(x="Angiosperm ratio",tag = "b")+ #y = expression ("Community wood density ("~g/cm^3~")"),
             theme_classic()+
             scale_y_continuous(sec.axis = dup_axis(),limits=c(0.15,1.3))+
             theme(legend.position='none',#panel.border = element_rect(colour = "black", fill=NA, size=1),
                   plot.tag = element_text(size=32),
                   plot.tag.position = c(0.08, 0.96),
                   axis.title.x = element_text(size=20),
                #    axis.title.y = element_text(size=18),
                   axis.ticks.length=unit(0.2, "cm"),
                   axis.text.x = element_text(colour = "black", size = 18),
                   axis.text.y = element_blank(),
                   axis.title.y = element_blank(),
                   plot.margin = unit(c(0.05,0,0,0), "cm"))
            
# read the wood density reference table with the angiosperm and gymnosperm information
woodDensityReferenceTable = fread("ReferenceTables/GFBI_Wood_Density_Reference_Table_with_A_G.csv") 
angiospermWoodDensityTable = woodDensityReferenceTable %>% filter(A_G=='a') %>% filter(WD_Level != "dataset")
gymnospermWoodDensityTable = woodDensityReferenceTable %>% filter(A_G=='g') %>% filter(WD_Level != "dataset")
# define the species level wood density plot
gymnospermPlot = ggplot(gymnospermWoodDensityTable, aes(x=WoodDensity))+
                  geom_density(color="gray30",fill ="#EDF8E9")+
                  xlim(0.15,1.3)+
                  geom_vline(xintercept = median(gymnospermWoodDensityTable$WoodDensity),color = "darkblue", size=0.8)+
                  annotate("text", x=median(gymnospermWoodDensityTable$WoodDensity),y=3, label=paste("Mean = ",round(mean(gymnospermWoodDensityTable$WoodDensity),2),sep=""), colour="black",vjust = -1,size=6)+
                  labs(x = expression ("Wood density ("~g/cm^3~")"),y="Density",tag = "a")+
                  coord_flip()+
                  scale_y_reverse()+
                  theme_classic()+
                  expand_limits(y = 0)+
                  theme(plot.tag = element_text(size=32), #panel.border = element_rect(colour = "black", fill=NA, size=1),
                        plot.tag.position = c(0.4, 0.96),
                        axis.title.x = element_text(size=20),
                        axis.title.y = element_text(size=20),
                        axis.ticks.length=unit(0.2, "cm"),
                        axis.text.x = element_text(colour = "black", size = 18),
                        axis.text.y = element_text(colour = "black", size = 18),
                        plot.margin= unit(c(0.05,0.2,0,0), "cm"))
                  
angiospermPlot = ggplot(angiospermWoodDensityTable, aes(x=WoodDensity))+
                  geom_density(color="gray30",fill ="#006D2C")+
                  xlim(0.12,1.3)+
                  geom_vline(xintercept = median(angiospermWoodDensityTable$WoodDensity),color = "darkblue", size=0.8)+
                  annotate("text", x=median(angiospermWoodDensityTable$WoodDensity),y=1.5, label=paste("Mean = ",round(mean(angiospermWoodDensityTable$WoodDensity),2),sep=""), colour="black",vjust = -1,size=6)+
                  labs(x = expression ("Wood density ("~g/cm^3~")"),y="Density",tag = "c")+
                  coord_flip()+ 
                  theme_classic()+
                  expand_limits(y = 0)+
                  scale_x_continuous(position = "top",limits=c(0.15,1.3))+ #breaks = c(0.2,0.4,0.6,0.8,1.0,1.2)
                  theme(plot.tag = element_text(size=32), #panel.border = element_rect(colour = "black", fill=NA, size=1),
                        plot.tag.position = c(0.1, 0.96),
                        axis.title.x = element_text(size=20),
                        axis.title.y = element_text(size=20),
                        axis.ticks.length=unit(0.2, "cm"),
                        axis.text.x = element_text(colour = "black", size = 18),
                        axis.text.y = element_text(colour = "black", size = 18),
                        plot.margin= unit(c(0.05,0,0,0.2), "cm"))
                  
# > sd(angiospermWoodDensityTable$WoodDensity)
# [1] 0.1438034
# > sd(gymnospermWoodDensityTable$WoodDensity)
# [1] 0.06877453
# pdf("WritingFolder/NewPlots/Fig_1_ab_Species_and_Community_Level_Wood_density_comparison.pdf",width =14, height=5)
# ggarrange(gymnospermPlot,boxplotFig,angiospermPlot,ncol = 3, nrow = 1,align="h",widths= c(0.28,0.8,0.28))#
# dev.off()


library(raster)
library(RColorBrewer)
library(data.table)
library(ggpubr)
library(maps)
library(gridExtra)
library(dplyr)
library(cowplot)
library(ggrastr)

setwd("~/Desktop/ETH_Study/WoodDensityMapingProject/")

trainDataTableRaw = fread("CovariatesExtractionFolder/Plot_Level_Wood_Density_Data_Frame_Metadata_With_Diversity_Latest_year.csv")[,-1] %>% dplyr::select(LAT,LON,WeightedWoodDensity) %>% na.omit() %>% mutate(x=LON,y=LAT,LON=NULL,LAT=NULL)

# # adaptation of the far east russian points
farEastRussiaPlots = trainDataTableRaw %>% dplyr::filter(x>136&y>58) %>% dplyr::mutate(x= x*(-1))#2014

otherPlots = trainDataTableRaw %>% dplyr::filter(x<=136|y<=58) #1186086
# # rbind the modified tables together
trainDataTable = rbind(farEastRussiaPlots,otherPlots)


# make a copy of the table 
# rbPal= colorRampPalette(c("#ffffcc","#006837"))
rbPal1= colorRampPalette(c("#FDE725FF","#3CBB75FF","#2D708EFF","#481567FF"))
# pal= colorRampPalette(c("#FDE725FF","#3CBB75FF","#2D708EFF","#481567FF")) #
# define the colour for the points
# inputTable$Col = pal(100)[as.numeric(cut(inputTable$mod,breaks = 100))]
duplicateTable = trainDataTable


# load the world border polygon
worldBorderRaw = shapefile("WORLD_BORDERS/TM_WORLD_BORDERS_SIMPL-0.3.shp")
worldBorder = crop(worldBorderRaw,c(-180, 180, -60, 84)) 
# set the resolution of how the polygon plot process
# worldBorder = gSimplify(worldBorderRaw, tol = 0.01, topologyPreserve = TRUE)
# define the equal earth projection
equalEarthProj =  "+proj=eqearth +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# project raster
extentInfo = projectExtent(worldBorder,equalEarthProj)
# reproject the world border
worldEqualEarth = spTransform(worldBorder ,to=extentInfo, CRS(equalEarthProj))


# transform the data into spatial points
coordinates(duplicateTable) = ~x+y
proj4string(duplicateTable) = CRS("+proj=longlat")
transformedPoints = spTransform(duplicateTable, to=extentInfo,CRS(equalEarthProj))
newTable = data.frame(transformedPoints@coords,CWD=trainDataTable$WeightedWoodDensity) %>% mutate(CWD=ifelse(CWD<0.4, 0.4,CWD)) %>%  mutate(CWD=ifelse(CWD>0.75, 0.75,CWD)) %>% mutate(x=coords.x1,y=coords.x2,coords.x2=NULL,coords.x1=NULL)



panelA = ggplot() + geom_polygon(data = worldEqualEarth,aes(x = long, y = lat, group = group),fill = "grey50",color = "grey50",linewidth = 0.1) + 
                    coord_fixed(1)+
                    ggrastr::rasterise(geom_point(data = newTable,aes(x = x, y = y,colour = CWD),size=2),dpi=1000,dev = "cairo") +
                    scale_colour_gradientn(colours = rbPal1(100),breaks= c(0.4,0.55,0.75),labels=c(expression(""<=0.4),"0.5",expression("">=0.75)))+
                    theme_classic()+
                    labs(tag = "d")+
                    guides(fill=guide_colorbar(frame.colour = "black", ticks.colour = "white"))+
                    theme(legend.position = c(0.14, 0.3),
                        legend.title = element_text(size=20),
                        legend.text = element_text(size=18),
                        legend.key.size = unit(1.6,"line"),
                        plot.margin= unit(c(-1,0.2,-1,2.1), "cm"), 
                        # panel.border = element_rect(colour = "black", fill=NA, size=2),
                        axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
                        axis.ticks.x = element_blank(),axis.text.x = element_blank(),
                        axis.title.x = element_blank(),axis.title.y = element_blank(),
                        axis.line = element_blank(),
                        plot.tag = element_text(size=32),
                        plot.tag.position = c(0.03, 0.95))+
                    labs(colour = expression(paste("Wood density (g/",cm^3,')',sep="")))  
                    
# load one of the subsample 
# singleSubSample = fread("ModelOptimization/WoodDensity_one_subsample_Seed_0.csv") %>% select("x","y","WdDnsty")
# singleSubSampleDup = singleSubSample


# # transform the data into spatial points
# coordinates(singleSubSampleDup) = ~x+y
# proj4string(singleSubSampleDup) = CRS("+proj=longlat")
# transformedPointsSub = spTransform(singleSubSampleDup, to=extentInfo,CRS(equalEarthProj))
# newTableSubsample = data.frame(transformedPointsSub@coords,CWD=singleSubSample$WdDnsty) %>% mutate(CWD=ifelse(CWD<0.4, 0.4,CWD)) %>%  mutate(CWD=ifelse(CWD>0.75, 0.75,CWD))%>% mutate(x=coords.x1,y=coords.x2,coords.x2=NULL,coords.x1=NULL)


# PanelB = ggplot() + geom_polygon(data = worldEqualEarth,aes(x = long, y = lat, group = group),fill = "grey50",color = "grey50",size = 0.1) + 
#                     coord_fixed(1)+
#                     geom_point(data = newTableSubsample,aes(x = x, y = y,colour = CWD),size=1.5) +
#                     scale_colour_gradientn(colours = rbPal1(100),breaks= c(0.4,0.55,0.75),labels=c(expression(""<=0.4),"0.5",expression("">=0.75)))+
#                     theme_classic()+
#                     labs(tag = "a")+
#                     guides(fill=guide_colorbar(frame.colour = "black", ticks.colour = "white"))+
#                     theme(legend.position = c(0.14, 0.3),
#                         legend.title = element_text(size=18),
#                         legend.text = element_text(size=14),
#                         legend.key.size = unit(1.6,"line"),
#                         plot.margin= unit(c(-1,0.2,-1,2.1), "cm"), 
#                         # panel.border = element_rect(colour = "black", fill=NA, size=2),
#                         axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
#                         axis.ticks.x = element_blank(),axis.text.x = element_blank(),
#                         axis.title.x = element_blank(),axis.title.y = element_blank(),
#                         axis.line = element_blank(),
#                         plot.tag = element_text(size=28),
#                         plot.tag.position = c(0.03, 0.95))+
#                     labs(colour = expression(paste("Wood density (g/",cm^3,')',sep="")))  


# pdf("WritingFolder/NewPlots/Fig_1_test.pdf",width = 16,height=12)
# ggarrange(panelA, ncol = 1, nrow = 1)
# dev.off()


#########################################################################################################################
################PANEL B and D
#########################################################################################################################
# exclude mangroves since it only has 3 records
trainDataTable = fread("CovariatesExtractionFolder/20231201_Wood_Density_Project_Merged_sampled_dataset_Diversity_Pixel_2023.csv")[,-1] %>% filter(WWF_Biome <14)


# make a copy of the table 
aggregatedTableNew = trainDataTable
# generate a new column to save the forest types not the biome
aggregatedTableNew$ForestType = aggregatedTableNew$WWF_Biome
# change the forest type representative number
aggregatedTableNew$ForestType[aggregatedTableNew$ForestType %in% c(1,2,3,7,9)] <- 30 #Tropical
aggregatedTableNew$ForestType[aggregatedTableNew$ForestType %in% c(4,5,8,10)] <- 40 #Temperate
aggregatedTableNew$ForestType[aggregatedTableNew$ForestType %in% c(6,11)] <- 50 #Boreal
aggregatedTableNew$ForestType[aggregatedTableNew$ForestType %in% c(12,13)] <- 60 # others

# transfer th
aggregatedTableNew$WWF_Biome = as.factor(aggregatedTableNew$WWF_Biome)
aggregatedTableNew$ForestType = as.factor(aggregatedTableNew$ForestType)

# biome level statistics of wood density 
biomeStat = aggregate(x=aggregatedTableNew[,c("WdDnsty")],by=aggregatedTableNew[,c("WWF_Biome")],FUN=mean)  %>% `colnames<-`(c("Biome","BiomeDensityMesn"))

biomeSD = aggregate(x=aggregatedTableNew[,c("WdDnsty")],by=aggregatedTableNew[,c("WWF_Biome")],FUN=sd)  %>% `colnames<-`(c("Biome","WD_SD"))
biomeStat$BiomeDensitySD = biomeSD$WD_SD
#    Biome BiomeDensityMesn BiomeDensitySD
# 1      1        0.5777197     0.08488713
# 2      2        0.5907360     0.09566713
# 3      3        0.5994275     0.14469650
# 4      4        0.5312310     0.09110975
# 5      5        0.4883861     0.06705488
# 6      6        0.4577088     0.04687929
# 7      7        0.5776676     0.11148264
# 8      8        0.5701719     0.09968431
# 9      9        0.4633262     0.07962942
# 10    10        0.5678115     0.10355938
# 11    11        0.5249678     0.06486145
# 12    12        0.6019746     0.08912799
# 13    13        0.5453389     0.08234713

# forest type level wood density 
forestTypeStat = aggregate(x=aggregatedTableNew[,c("WdDnsty")],by=aggregatedTableNew[,c("ForestType")],FUN=mean) %>% `colnames<-`(c("ForestType","TypeDensityMesn")) %>% mutate(ForestType = c("Tropical","Temperate","Boreal","Dryland"))

forestTypeSD = aggregate(x=aggregatedTableNew[,c("WdDnsty")],by=aggregatedTableNew[,c("ForestType")],FUN=sd) %>% `colnames<-`(c("ForestType","WD_SD"))
forestTypeStat$TypeDensitySD = forestTypeSD$WD_SD

#   ForestType TypeDensityMesn TypeDensitySD
# 1   Tropical       0.5717622    0.09801427
# 2  Temperate       0.5235225    0.08921965
# 3     Boreal       0.4623159    0.05122399
# 4    Dryland       0.5891659    0.09078602

write.csv(biomeStat,"BiomeLevelStatistics/Biome_level_wood_density_statistics.csv")
write.csv(forestTypeStat,"BiomeLevelStatistics/Forest_type_level_wood_density_statistics.csv")

rbPal= colorRampPalette(brewer.pal(9,"RdYlBu"))

boxPlotMean = ggplot(data = aggregatedTableNew[aggregatedTableNew$WWF_Biome %in% c(1:13),], aes(x = WWF_Biome, y = WdDnsty,fill=WWF_Biome)) +
            #    geom_jitter(aes(colour = WWF_Biome),alpha = 1/40) +
                geom_boxplot(outlier.size=-1)+ 
            #    scale_fill_brewer(palette="Spectral")+
            scale_fill_manual(values=rbPal(13))+
                # coord_flip()+
                theme_classic()+
                ylim(0.2,0.9)+
            #    scale_color_brewer(palette = "Spectral") +
                # ylab(expression(paste("Wood density (g/",cm^3,')',sep="")))+
                xlab("Biome")+
                labs(tag = "f")+
                scale_x_discrete(labels=c("1" = "Tropical moist",
                                          "2" = "Tropical dry",
                                          "3" = "Tropical coniferous",
                                          "4" = "Temperate broadleaf",
                                          "5" = "Temperate coniferous",
                                          "6" = "Boreal",
                                          "7" = "Tropical savanna",
                                          "8" = "Temperate savanna",
                                          "9" = "Flooded savanna",
                                          "10" = "Montane grassland",
                                          "11" = "Tundra",
                                          "12" = "Mediterranean forest",
                                          "13" = "Desert"),guide = guide_axis(angle = 45))+
                theme(axis.title.x = element_blank(),
                        axis.title.y = element_blank(),
                        axis.ticks.length=unit(0.15, "cm"),
                        axis.text.x = element_text(colour = "black", size = 18),
                        axis.text.y = element_text(colour = "black", size = 18),
                        plot.tag = element_text(size=32),
                        plot.tag.position = c(0.11, 0.95))+
                theme(legend.position="none")#,panel.border = element_rect(colour = "black", fill=NA, size=2)
            #    stat_summary(fun.y = mean, color = "darkred", position = position_dodge(0.75),geom = "point", shape = 18, size = 3,show.legend = FALSE)

boxPlotForestTypeMean = ggplot(data = aggregatedTableNew[aggregatedTableNew$ForestType %in% c(30,40,50,60),], aes(x = ForestType, y = WdDnsty,fill=ForestType)) +
                     #    geom_jitter(aes(colour = WWF_Biome),alpha = 1/40) +
                     geom_boxplot(outlier.size=-1)+ 
                      #    scale_fill_brewer(palette="Spectral")+
            scale_fill_manual(values=rbPal(4))+
                # coord_flip()+
                theme_classic()+
                ylim(0.2,0.9)+
            #    scale_color_brewer(palette = "Spectral") +
                ylab(expression(paste("Wood density (g/",cm^3,')',sep="")))+
                # xlab("Forest type")+
                labs(tag = "e")+
                scale_x_discrete(labels=c("30" = "Tropical forest", "40" = "Temperate forest","50" = "Boreal forest","60" = "Dryland forest"),guide = guide_axis(angle = 45))+
                theme(axis.title.x = element_blank(),
                        axis.title.y = element_text(size=20),
                        axis.ticks.length=unit(0.15, "cm"),
                        axis.text.x = element_text(colour = "black", size = 18),
                        axis.text.y = element_text(colour = "black", size = 18),
                        plot.tag = element_text(size=32),
                        plot.tag.position = c(0.3, 0.95))+
                theme(legend.position="none")#,panel.border = element_rect(colour = "black", fill=NA, size=2)
            #    stat_summary(fun.y = mean, color = "darkred", position = position_dodge(0.75),geom = "point", shape = 18, size = 3,show.legend = FALSE)



pdf("WritingFolder/NewPlots/Fig_1_Wood_Density_Plots_Box+Plot_Panel_a_b_c_d_e_f.pdf",width = 14,height=15)
top_row = plot_grid(ggarrange(gymnospermPlot,boxplotFig,angiospermPlot,ncol = 3, nrow = 1,align="h",widths= c(0.28,0.8,0.28)))
mid_row = plot_grid(panelA, align = "hv")
bottom_row = plot_grid(boxPlotForestTypeMean,boxPlotMean, align = "hv", rel_widths = c(5,13))
plot_grid(top_row,mid_row,bottom_row,ncol=1,rel_heights = c(7,9,6),align = "hv")
dev.off()





