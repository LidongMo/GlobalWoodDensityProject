
library(ggplot2)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(data.table)
# set working directory 
setwd('~/WoodDensityMapingProject')

# exclude mangroves since it only has 3 records
#!!!!
# Due to the large size of the data, pleae run the code "Figure_5_2_Different_wood_density_derived_biomass_map_random_sample.r"
sampledDataTable = fread("Data/ModelCompare/Random_sample_difference_of_derived_biomass.csv")[,-1] %>% filter(Biome %in% c(1:13)) %>% mutate(Diff= Diff*100) %>%na.omit()

# transfer th
sampledDataTable1 = sampledDataTable
sampledDataTable1$Biome = as.factor(sampledDataTable1$Biome)

# biome level statistics of wood density 
biomeStat = aggregate(x=sampledDataTable1[,c("Diff")],by=sampledDataTable1[,c("Biome")],FUN=mean)  %>% `colnames<-`(c("Biome","DiffMean")) 

write.csv(biomeStat,"Data/BiomeLevelStatistics/Biome_level_wood_density_statistics.csv")

rbPal= colorRampPalette(brewer.pal(9,"RdYlBu"))

boxPlotMean = ggplot(data = sampledDataTable, aes(x = as.factor(Biome), y = Diff,fill=as.factor(Biome))) +
                geom_boxplot(outlier.size=-1)+ 
                scale_fill_manual(values=rbPal(13))+
                theme_classic()+
                ylim(-30,30)+
                ylab("Difference (in percentage)")+
                xlab("Biome")+
                labs(tag = "c")+
                geom_hline(yintercept=0,linetype="solid", color="gray30", linewidth=1)+
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
                      axis.title.y = element_text(size=18),
                      axis.ticks.length=unit(0.15, "cm"),
                      axis.text.x = element_text(colour = "black", size = 15),
                      axis.text.y = element_text(colour = "black", size = 15),
                      plot.tag = element_text(size=28),
                      plot.tag.position = c(0.11, 0.95),
                      legend.position="none")

pdf("WPlots/Fig_5_2_Difference_in_WD_derived_biomass_Panel_c.pdf",width = 12,height=5)
plot(boxPlotMean)
dev.off()


# > biomeStat
#    Biome   DiffMean
# 1      1 -11.500628
# 2      2 -17.297711
# 3      3 -12.119646
# 4      4  -3.222271
# 5      5   9.906479
# 6      6  12.605651
# 7      7 -16.966413
# 8      8  -1.421306
# 9      9 -12.471866
# 10    10  -5.708690
# 11    11  -8.318943
# 12    12 -20.954873
# 13    13 -19.408399