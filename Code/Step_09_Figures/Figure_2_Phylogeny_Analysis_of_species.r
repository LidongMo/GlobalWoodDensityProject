# load the package which could provide the phylogeny information for all the species
# here is the link : https://github.com/jinyizju/V.PhyloMaker

library(V.PhyloMaker)
library(data.table)
library(ape)
library(ggtree)
library(tidytree)
library(phytools)
library(geiger) #for matching tree and traits
library(diversitree)
library(RColorBrewer)
library(phylocomr)
library(phylosignal)
library(phylobase)
library(plantlist)
library(stringr)
library(castor)
library(dplyr)
library(purrr)

# set the working directory
setwd("~/Desktop/ETH_Study/WoodDensityMapingProject/")

###################################################################################################################
# STEP1 Phylogeny analysis and plotting for all the species
###################################################################################################################

# load the table of wood density database with species name uniformed in TNRS
woodDensityDataRaw = fread("ReferenceTables/Family_name_Allocated_Wood_Density_Data_Frame_TPL.csv")[,-1] %>% dplyr::select(Family,Genus,Species,WoodDensity) #4289 species retained by using the TPL coreccted data

# woodDensityDataRaw = fread("ReferenceTables/Wood_density_Table_with_Family_allocated_TNRS_Mannul_corr.csv")[,-1] %>% dplyr::select(Family,Genus,Species,WoodDensity) #4276 species retained by using the TNRS coreccted data
woodDensityData = aggregate(x=woodDensityDataRaw[,c("WoodDensity")],by=woodDensityDataRaw [,c("Species","Genus","Family")],FUN=mean,narm=T) %>% rename(family = Family,genus =Genus, species=Species) %>% filter(!is.na(WoodDensity)&WoodDensity<2)


# get the phylogeny information
phyloResult = phylo.maker(woodDensityData,scenarios = "S3",output.tree=T)

# phyloResult$scenario.3


phyloTree = phyloResult$tree.scenario.3
# phyloTree = chronos(phyloTree, lambda = 0, model = "relaxed") 
# write the tree to the local folder
write.tree(phyloTree, "PhylogenyAnalysis/Phylogeny_tree_of_vascular_plants.tre")

# # tt= groupOTU(phyloTree, as.factor(phyloTree$tip.label))
# # phyloTree = chronos(phyloTree, lambda = 0, model = "relaxed") 
# ggtree(phyloTree,layout='circular') + geom_tiplab(size=1, aes(angle=angle))

# woodDensityData$species.name=paste(woodDensityData$genus, woodDensityData$species, sep="_")
rownames(woodDensityData)=paste(woodDensityData$genus, woodDensityData$species, sep="_")

dat1 = geiger::treedata(phyloTree,as.matrix(woodDensityData))
cleanedWoodDensityData = as.data.frame(dat1$data)
newTree = dat1$phy
# write.nexus(newTree,file="nexus.trees")
# newTree = read.nexus(file="nexus.trees")

write.tree(newTree, "PhylogenyAnalysis/Phylogeny_tree_of_woody_plants.tre")
# load the oders and group data from plantlist package
data("orders_dat")
# orders_dat
updatedWoodDensityData = left_join(cleanedWoodDensityData,orders_dat %>% dplyr::select("FAMILY","ORDER","GROUP"),by = c("family" = "FAMILY")) %>% mutate(GROUP=ifelse(GROUP == "", "Angiosperms", GROUP))

# newTree = chronos(chronoTree, lambda = 0, model = "relaxed") 
rownames(updatedWoodDensityData)=paste(updatedWoodDensityData$genus, updatedWoodDensityData$species, sep="_")
newWoodDensityData = updatedWoodDensityData[newTree$tip.label, ] %>% mutate(WoodDensity = as.numeric(as.character(WoodDensity)))

traitVector = newWoodDensityData$WoodDensity
names(traitVector)<-newTree$tip.label

traitTable = data.frame(name= as.vector(newTree$tip.label), WoodDensity = newWoodDensityData$WoodDensity,WoodDensity1 = newWoodDensityData$WoodDensity) %>% purrr::modify_if(is.factor, as.character) 
# run the aot analysis
aotResult = ph_aot(traits=traitTable,phylo=newTree,trait_contrasts = 1, randomizations = 999,ebl_unstconst = F)
# get the conservatism result
conservTable = aotResult$trait_conservatism %>% filter(trait.name == "WoodDensity") %>% filter(ntaxa<574&ntaxa>50)
calculationTable = fread("PhylogenyAnalysis/Order_Level_Lambda_BlombergsK_Table.csv")[,-1]

orderNames = calculationTable$Order
# oder level tree subset and phylogeny analysis
orderLevelNodeFindingFunc = function(ord)
{
    # subet the data frame
    orderLevelTable = calculationTable %>% filter(Order == ord)
    # subset the conservTable bu tip.mn and ntaxa
    orderLevelConservTable = conservTable %>% filter(tip.mn == orderLevelTable$WoodDensityMean& ntaxa==orderLevelTable$SpeciesNumber) %>% dplyr::select(ntaxa,tip.mn,tmn.ranklow,tmn.rankhi,node.mn)
    # get the node position 
    nodePosition = find_root_of_monophyletic_tips(newTree, newWoodDensityData %>% filter(ORDER ==ord) %>% rownames(), as_MRCA=TRUE, is_rooted=FALSE)
    # 
    if (nrow(orderLevelConservTable) ==0)
    {
        ValCompare = NA
    }else
    {
        if (orderLevelConservTable$tmn.ranklow<25|orderLevelConservTable$tmn.rankhi>975)
        {
            ValCompare = "Lower"
        }else 
        if (orderLevelConservTable$tmn.ranklow>975|orderLevelConservTable$tmn.rankhi<25)
        {
            ValCompare = "Higher"
        }else 
        {
            ValCompare = NA
        }
    }
    
    return(data.frame(cbind(orderLevelTable,orderLevelConservTable),NodePosition = nodePosition,ValueCompare = ValCompare))
}

# use the lapply to sun the calculation
calculationList = lapply(orderNames,orderLevelNodeFindingFunc)
# rbind the result
orderPositionTable = rbindlist(calculationList)
write.csv(orderPositionTable,"PhylogenyAnalysis/Order_Level_Phy_AOT_analysis_result_Table.csv")
# lets get the information for each 
# get the higher or lower wood density clades
# lowerTable = conservTable %>% filter(nmn.ranklow<25 &nmn.rankhi>975) %>% filter(ntaxa<600&ntaxa>50)
# higherTable = conservTable %>% filter(nmn.ranklow<975 &nmn.rankhi>25) %>% filter(ntaxa<600&ntaxa>50)

# Lambda
phylosig(newTree,traitVector,method="lambda", test=T)

# Phylogenetic signal lambda : 0.917751 
# logL(lambda) : 2689.15 
# LR(lambda=0) : 2390.33 
# P-value (based on LR test) : 0 

# Blombergs K
phylosig(newTree,traitVector,method="K", test=T)

# Phylogenetic signal K : 0.00510281 
# P-value (based on 1000 randomizations) : 0.005


# define the vector which contains the order names 
habitat = setNames(newWoodDensityData$ORDER,rownames(newWoodDensityData))
n = Ntip(newTree)
colourOrders = calculationTable %>% filter(SpeciesNumber>=50) %>% dplyr::select(Order) %>% unlist()
grayOrders = calculationTable %>% filter(SpeciesNumber<50) %>% dplyr::select(Order)%>% unlist()
# ranked vector of order names 
fullOrderVector = unique(newWoodDensityData$ORDER)
colourOrdersUp = fullOrderVector[fullOrderVector %in% colourOrders]
# allocate the colour by order names
col.hab = setNames(c(colorRampPalette(brewer.pal(8, "Set1"))(length(colourOrdersUp)),rep("gray15",55-length(colourOrdersUp))),c(colourOrdersUp,grayOrders))

# construct the plot object
obj = contMap(newTree, traitVector, plot=F, fsize=0.2, tip.labels=T,lwd=1)#type="fan", 
# define a colour pallete
# rbPal = colorRampPalette(c("black","darkred",brewer.pal(9, "RdYlBu")))(1001) %>% rev()
rbPal = colorRampPalette(c("black",'#67001F', '#B2182B', '#D6604D', '#F4A582','#D1E5F0', '#92C5DE','#4393C3', '#2166AC'))(1001) %>% rev()

#This adds a column of color values
# based on the wood density values
obj$cols = setNames(c(rbPal),0:1000)
plottingOrder = orderPositionTable %>% filter(!is.na(ValueCompare)) 
compareCol = setNames(c(ifelse(plottingOrder$ValueCompare=="Higher", "#cb1b16", "#1368aa")),plottingOrder$ValueCompare)
# ttt= setNames(rbPal(10)[as.numeric(cut(x,breaks = 10))],newTree$tip.label)
pdf("WritingFolder/NewPlots/Figure_3_Phylogeny_and_traits_plot_latest.pdf",width = 40,height=40)
par(bg="gray15",xpd=TRUE)
plotTree.wBars(obj$tree, traitVector, fsize=0.2, scale=40,width=0.5,lwd=2, tip.labels=F,
              method="plotSimmap", colors=obj$cols,type="fan",outline=F,border="transparent",mar=c(2,2,2,25))#type="fan",
nodelabels(text=str_pad(round(plottingOrder$tip.mn,2), 4, pad = "0",side = c("right")),node=plottingOrder$NodePosition,frame="circle",col="white",cex=2.8,bg=compareCol) #+Ntip(obj$tree)
objTT<-get("last_plot.phylo",envir=.PlotPhyloEnv)
for(i in 1:n)
{
    cc = if(objTT$xx[i]>0) 65 else -65
    th = atan(objTT$yy[i]/objTT$xx[i])
    segments(#objTT$xx[i],objTT$yy[i],objTT$xx[i]+cc*cos(th),objTT$yy[i]+cc*sin(th),
    	objTT$xx[i]+cc*cos(th),objTT$yy[i]+cc*sin(th),objTT$xx[i]+1.15*cc*cos(th),objTT$yy[i]+1.15*cc*sin(th),
        lwd=4,lend=2,col=col.hab[habitat[newTree$tip.label[i]]])
}
add.color.bar(160,cols=obj$cols,title=expression("Wood density (g/cm3)"),obj$lims,digits=3,lwd=50,fsize=3,prompt=FALSE,x=340,y=400,subtitle="",outline=T)

legend(x=405,y=370,colourOrdersUp,pch=15,col=c(col.hab[1:length(colourOrdersUp)]),pt.cex=6,cex=4,bty="o",ncol = 1,text.col="white")

dev.off()

###################################################################################################################
# STEP 2 Phylogeny analysis and plotting for each order
###################################################################################################################

library(V.PhyloMaker)
library(data.table)
library(ape)
library(ggtree)
library(tidytree)
library(phytools)
library(geiger) #for matching tree and traits
library(diversitree)
library(RColorBrewer)
library(phylocomr)
library(phylosignal)
library(phylobase)
library(plantlist)

# set the working directory
setwd("~/Desktop/ETH_Study/WoodDensityMapingProject/")

# load the table of wood density database with species name uniformed in TNRS
woodDensityDataRaw = fread("ReferenceTables/Family_name_Allocated_Wood_Density_Data_Frame_TPL.csv")[,-1] %>% dplyr::select(SpeciesName,Family,Genus,Species,WoodDensity) #4289 species retained by using the TPL coreccted data

# woodDensityDataRaw = fread("ReferenceTables/Wood_density_Table_with_Family_allocated_TNRS_Mannul_corr.csv")[,-1] %>% dplyr::select(Family,Genus,Species,WoodDensity) #4276 species retained by using the TNRS coreccted data
woodDensityData = aggregate(x=woodDensityDataRaw[,c("WoodDensity")],by=woodDensityDataRaw [,c("SpeciesName","Species","Genus","Family")],FUN=mean,narm=T) %>% rename(family = Family,genus =Genus, species=Species) %>% filter(!is.na(WoodDensity)&WoodDensity<2)


# get the phylogeny information
phyloResult = phylo.maker(woodDensityData,scenarios = "S3",output.tree=T)

# phyloResult$scenario.3


phyloTree = phyloResult$tree.scenario.3
# phyloTree = chronos(phyloTree, lambda = 0, model = "relaxed") 
# write the tree to the local folder
# write.tree(phyloTree, "PhylogenyAnalysis/Woody_species_phylogeny.tre")

# # tt= groupOTU(phyloTree, as.factor(phyloTree$tip.label))
# # phyloTree = chronos(phyloTree, lambda = 0, model = "relaxed") 
# ggtree(phyloTree,layout='circular') + geom_tiplab(size=1, aes(angle=angle))

# woodDensityData$species.name=paste(woodDensityData$genus, woodDensityData$species, sep="_")
rownames(woodDensityData)=paste(woodDensityData$genus, woodDensityData$species, sep="_")

dat1 = geiger::treedata(phyloTree,as.matrix(woodDensityData))
cleanedWoodDensityData = as.data.frame(dat1$data)
newTree = dat1$phy
# load the oders and group data from plantlist package
data("orders_dat")
# orders_dat
updatedWoodDensityData = left_join(cleanedWoodDensityData,orders_dat %>% dplyr::select("FAMILY","ORDER","GROUP"),by = c("family" = "FAMILY")) %>% mutate(GROUP=ifelse(GROUP == "", "Angiosperms", GROUP))
# write the updated full wood density data with family, order and group 
write.csv(updatedWoodDensityData,"PhylogenyAnalysis/Filtered_Wood_Density_Table_with_Orders_and_Family_Information.csv")
updatedWoodDensityData = fread("PhylogenyAnalysis/Filtered_Wood_Density_Table_with_Orders_and_Family_Information.csv")[,-1] %>% mutate(TaxonName = paste(genus,species,sep="_"))
# generate a table which is the combine of the family, genus and species name
taxonMatchedTable = aggregate(x=woodDensityDataRaw[,c("WoodDensity")],by=woodDensityDataRaw[,c("SpeciesName","Species","Genus","Family")],FUN=mean,narm=T) %>% mutate(TaxonName = paste(Family,Genus,Species,sep="_")) %>% filter(TaxonName %in% unlist(updatedWoodDensityData$TaxonName))
write.csv(taxonMatchedTable,"PhylogenyAnalysis/Filtered_Wood_Density_Table_for_Species_with_Phylogeny_Info.csv")
# use the order names as the identifier to subset the data frame and do anaylsis
orderNames = unique(updatedWoodDensityData$ORDER)
# oder level tree subset and phylogeny analysis
orderLevelFunc = function(ord)
{
    # subet the data frame
    orderLevelTable = updatedWoodDensityData %>% filter(ORDER == ord)
    if(nrow(orderLevelTable)>=5)
    {
        # add the row names 
        rownames(orderLevelTable)=paste(orderLevelTable$genus, orderLevelTable$species, sep="_")
        # use the data frame to get the phylogeny information by the treedata function
        orderLevelTreeData = geiger::treedata(phyloTree,orderLevelTable)
        # get the tree data paired with traits information
        orderWoodDensityData = as.data.frame(orderLevelTreeData$data)
        # get the paired phylogeny tree
        orderTree = orderLevelTreeData$phy
        # write the tree and data frame to the local folder
        write.csv(orderWoodDensityData,paste("PhylogenyAnalysis/OrderLevelWDTables/Order_",ord,"_Subset_wood_Density_data.csv",sep=""))
        # and write the tree
        write.tree(orderTree, paste("PhylogenyAnalysis/OrderLevelTrees/Order_",ord,"_Subset_wood_Density_data.tre",sep=""))
        # do the test for phylogeny
        traitVector = as.numeric(as.character(orderWoodDensityData$WoodDensity))
        names(traitVector) =  orderTree$tip.label
        # apply the test
        lambdaResult = phylosig(orderTree,traitVector,method="lambda", test=T)
        blombergsResult = phylosig(orderTree,traitVector,method="K", test=T)
        meanWD = round(mean(orderLevelTable$WoodDensity),3)


        # define the output table 
        orderLevelOutput = data.frame(Order = ord,
                                      SpeciesNumber = nrow(orderLevelTable),
                                      Lambda = lambdaResult$lambda,
                                      LambdaP_Val = lambdaResult$P,
                                      BlombergsK = blombergsResult$K,
                                      BlombergsP_Val = blombergsResult$P,
                                      WoodDensityMean = meanWD)
        

    }else
    {
        orderLevelOutput = data.frame(Order = ord,
                                      SpeciesNumber = nrow(orderLevelTable),
                                      Lambda = NA,
                                      LambdaP_Val = NA,
                                      BlombergsK = NA,
                                      BlombergsP_Val = NA,
                                      WoodDensityMean = round(mean(orderLevelTable$WoodDensity),3))

    }
    return(orderLevelOutput)
    
}
# use the lapply to sun the calculation
calculationList = lapply(orderNames,orderLevelFunc)
# rbind the result
calculationTable = rbindlist(calculationList)
# write to local folder
write.csv(calculationTable,"PhylogenyAnalysis/Order_Level_Lambda_BlombergsK_Table.csv")




