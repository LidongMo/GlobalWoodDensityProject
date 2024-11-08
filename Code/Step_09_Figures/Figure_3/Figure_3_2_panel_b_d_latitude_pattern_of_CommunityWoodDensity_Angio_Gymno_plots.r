library(ggplot2)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(data.table)
# set working directory 
setwd('~/WoodDensityMapingProject')
##################################################################
# STEP 1 SD 
##################################################################
# load the data table
latitudeTable_WD_VC = fread("Data/Aggregated_CWD_along_latitude_based_on_random_sample.csv")[,-1] %>%mutate_at(vars(c('lowCI')),~ifelse(lowCI <=0, 0, .))


# display.brewer.pal(n, name)
p1 = latitudeTable_WD_VC %>% filter(Type == 'Model_CWD')%>%
    ggplot(aes(x = LatRound, y= mean, color=Type, group=Type)) + 
    geom_ribbon(aes(ymin=lowCI, ymax=highCI), fill="darkcyan",alpha=0.3,color=NA)+
    # scale_fill_manual(values = c("GS_Max1" = "#009e73","GS_Max2" = "#009e73","GS_Mean1" = "#cc79a7","GS_Mean2" = "#cc79a7","HM1"="#56b4e9","HM2"= "#56b4e9","SD1"="#e69f00","SD2"="#e69f00","WK1"="#d55e00","WK2"= "#d55e00"))+
    # geom_hline(yintercept=0.6,linetype="dashed", color="gray30", size=0.5)+
    geom_line(aes(linetype=Type),,size = 1)+
    scale_color_manual(values = c("darkcyan"))+
    scale_linetype_manual(values=c(rep(c("solid","dashed"), 5)))+
    scale_x_continuous(expand = c(0,0),limits = c(-60, 90), breaks = c(-60,-40,-20,0,20,40,60,80,90))+
    scale_y_continuous(expand = c(0,0),limits = c(0.3, 0.9), breaks = c(0.4,0.6,0.8))+
    # scale_alpha_discrete(range = c(0.2, 0.7))+
    ylab("") +
    theme_bw()+
    coord_flip() + #ylim = c(-0.4, 0.4),xlim=c(27,75)
    # ylim(0.3,0.9)+
    xlab("Latitude")+
    # guides( col = guide_legend(ncol=2),linetype = guide_legend())+
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),legend.position="none")

# display.brewer.pal(n, name)
p2 = latitudeTable_WD_VC %>% filter(Type == 'Model_AngioWD')%>%
    ggplot(aes(x = LatRound, y= mean, color=Type, group=Type)) + 
    geom_ribbon(aes(ymin=lowCI, ymax=highCI), fill="darkcyan",alpha=0.3,color=NA)+
    # scale_fill_manual(values = c("GS_Max1" = "#009e73","GS_Max2" = "#009e73","GS_Mean1" = "#cc79a7","GS_Mean2" = "#cc79a7","HM1"="#56b4e9","HM2"= "#56b4e9","SD1"="#e69f00","SD2"="#e69f00","WK1"="#d55e00","WK2"= "#d55e00"))+
    # geom_hline(yintercept=0.6,linetype="dashed", color="gray30", size=0.5)+
    geom_line(aes(linetype=Type),size = 1)+
    scale_color_manual(values = c("darkcyan"))+
    scale_linetype_manual(values=c(rep(c("solid","dashed"), 5)))+
    scale_x_continuous(expand = c(0,0),limits = c(-60, 90), breaks = c(-60,-40,-20,0,20,40,60,80,90))+
    scale_y_continuous(expand = c(0,0),limits = c(0.3, 0.9), breaks = c(0.4,0.6,0.8))+
    # scale_alpha_discrete(range = c(0.2, 0.7))+
    ylab("") +
    theme_bw()+
    coord_flip() + #ylim = c(-0.4, 0.4),xlim=c(27,75)
    # ylim(0.3,0.9)+
    xlab("Latitude")+
    # guides( col = guide_legend(ncol=2),linetype = guide_legend())+
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),legend.position="none")

p3 = latitudeTable_WD_VC %>% filter(Type == 'Model_GymnoWD')%>%
    ggplot(aes(x = LatRound, y= mean, color=Type, group=Type)) + 
    geom_ribbon(aes(ymin=lowCI, ymax=highCI), fill="darkcyan",alpha=0.3,color=NA)+
    # scale_fill_manual(values = c("GS_Max1" = "#009e73","GS_Max2" = "#009e73","GS_Mean1" = "#cc79a7","GS_Mean2" = "#cc79a7","HM1"="#56b4e9","HM2"= "#56b4e9","SD1"="#e69f00","SD2"="#e69f00","WK1"="#d55e00","WK2"= "#d55e00"))+
    # geom_hline(yintercept=0.6,linetype="dashed", color="gray30", size=0.5)+
    geom_line(aes(linetype=Type),size = 1)+
    scale_color_manual(values = c("darkcyan"))+
    scale_linetype_manual(values=c(rep(c("solid","dashed"), 5)))+
    scale_x_continuous(expand = c(0,0),limits = c(-60, 90), breaks = c(-60,-40,-20,0,20,40,60,80,90))+
    scale_y_continuous(expand = c(0,0),limits = c(0.3, 0.9), breaks = c(0.4,0.6,0.8))+
    # scale_alpha_discrete(range = c(0.2, 0.7))+
    ylab("") +
    theme_bw()+
    coord_flip() + #ylim = c(-0.4, 0.4),xlim=c(27,75)
    # ylim(0.3,0.9)+
    xlab("Latitude")+
    # guides( col = guide_legend(ncol=2),linetype = guide_legend())+
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),legend.position="none")

    # ,
    #       panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',colour = "gray65"), 
    #       panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',colour = "gray85")



# px = latitudeTable_WD_VC %>% filter(Type == 'Model_WD_Vari') %>%
#     ggplot(aes(x = LatRound, y= mean, color=Type, group=Type)) + 
#     geom_ribbon(aes(ymin=lowCI, ymax=highCI, fill=Type),alpha=0.1,color=NA)+
#     # scale_fill_manual(values = c("GS_Max1" = "#009e73","GS_Max2" = "#009e73","GS_Mean1" = "#cc79a7","GS_Mean2" = "#cc79a7","HM1"="#56b4e9","HM2"= "#56b4e9","SD1"="#e69f00","SD2"="#e69f00","WK1"="#d55e00","WK2"= "#d55e00"))+
#     # geom_hline(yintercept=0)+
#     geom_line(aes(linetype=Type),lwd = 0.35,)+
#     scale_color_manual(values = c("#009e73"))+
#     scale_linetype_manual(values=c(rep(c("solid","dashed"), 5)))+
#     scale_x_continuous(expand = c(0,0),limits = c(-60, 90), breaks = c(-60,-40,-20,0,20,40,60,80,90))+
#     # scale_alpha_discrete(range = c(0.2, 0.7))+
#     ylab("") +
#     theme_classic()+
#     coord_flip() + #ylim = c(-0.4, 0.4),xlim=c(27,75)
#     ylim(0,0.05)+
#     xlab("Latitude")+
#     guides( col = guide_legend(ncol=2),linetype = guide_legend())+
#     theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),legend.position=c(0.8,0.8))




pdf("Plots/Fig_3_Panel_b_d_f_Wood_density_pattern_along_latitude.pdf",width =3.5, height=17)
library("gridExtra")
ggarrange(p1,p2,p3,ncol = 1, nrow = 3,align="hv",widths = c(1),heights = c(3,3,3) )

dev.off() 

