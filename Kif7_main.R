options(stringsAsFactors = F)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(monocle)
library(Seurat)
library(cowplot)
library(RColorBrewer)
library(pheatmap)

source("https://raw.githubusercontent.com/leezx/Toolsets/master/R/Toolsets.R")
source("https://raw.githubusercontent.com/leezx/Toolsets/master/R/Plot.R")

# lineage
# myColors <- brewer.pal(8,"Set2")
myColors5 <- c('#FDB462','#B3DE69','#FCCDE5',"#ff309a",'#D9D9D9')
names(myColors5) <- c("GP","BP","NPearly","NPlate","ENMFB")


# integrated seurat obj from E9.5 to P21
print(load("B.N.G_P.annotated.E95_P21.seurat.Rdata"))

# only use E13.5 and E16.5
print(load("E135_E165_ens.integrated.Rdata"))

# markers & heatmap
final.genes <- c('Sox10',"Plp1", 'Fabp7','Sparc','Rgcc',
                 'Ret','Gal','Nefm',
                 'Tubb3','Phox2b','Meg3','Stmn2','Stmn3',
                 'Pcsk1n','Rtn1','Cartpt','Prph','Sncg',
                 'Acta2','Mylk','Tagln','Col1a1','Dlk1','Myl9','Tpm2','Tpm1')

options(repr.plot.width=8, repr.plot.height=7)
hp <- DoHeatmap(ens.integrated.E13.16, features = final.genes, group.by = "lineage.sub", hjust=0.5, 
        # cells = rownames(subset(all_umap, stage %in% c("E13.5","Kif7 cKO","E16.5"))),
        group.bar = T, group.colors = myColors5, 
        disp.min = -1.5, disp.max = 1.5, 
        slot = "scale.data", label = T,
        size = 3, angle = 0, raster = F,
        draw.lines = T, lines.width = 100, group.bar.height = 0.03, combine = T) + 
        # NoLegend() +
        scale_fill_gradientn(colors = PurpleAndYellow(), na.value = "white") + # PurpleAndYellow(), c("blue", "white", "red")
        theme(axis.text.y = element_text(size = 15, face="italic"))
hp


# barplot
all_umap$lineage.sub <- all_umap$oneByone
all_umap[all_umap$cluster=="c4" & all_umap$oneByone=="NP",]$lineage.sub <- "late NP"
library(dplyr)
countTable <- all_umap %>% 
				dplyr::group_by(stage) %>% 
				dplyr::count(lineage.sub) %>% 
				dplyr::mutate(prop = n/sum(n))
plot.data <- subset(countTable, stage %in% c("E13.5","Kif7 cKO"))
plot.data$lineage.sub <- factor(plot.data$lineage.sub, levels = c("GP", "BP", "NP", "late NP", "MF"))
plot.data$stage <- factor(plot.data$stage, levels = c("E13.5","Kif7 cKO"))
plot.data <- plot.data[c(3,5,4,1,2, 8,10,9,6,7),]
plot.data$text.y <- c(0.01, 0.01225394+0.13909529/2, 0.01225394+0.13909529+0.18276626/2, 0.01225394+0.32186156+0.39303872/2, 0.01225394+0.32186156+0.39303872+0.27284578/2,
                      0.01, 0.02180000+0.11240000/2, 0.02180000+0.11240000+0.15260000/2, 0.02180000+0.26500000+0.41620000/2, 0.02180000+0.26500000+0.41620000+0.29700000/2)
plot.data$lineage.final <- plyr::mapvalues(plot.data$lineage.sub, 
                                           from=c("GP","BP","NP","late NP","MF"),
                                           to=c("GP","BP","NPearly","NPlate","ENMFB"))
plot.data$lineage.final <- factor(plot.data$lineage.final, levels = c("GP","BP","NPearly","NPlate","ENMFB"))
names(myColors5) <- c("GP","BP","NPearly","NPlate","ENMFB")
options(repr.plot.width=3.5, repr.plot.height=6)
library(scales)
g <- ggplot(data=plot.data, aes(x=stage, y=prop, fill=lineage.final)) +
  geom_bar(stat="identity", position="fill", alpha=1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "",y = "Percentage of Population (%)", title = " ") + 
  theme(axis.text.x  = element_text(face="bold", angle=30, size = 16, color = "black", vjust=0.5),
        axis.text.y  = element_text(face="plain", size = 10, color = "black"),
        axis.title =element_text(size = 14)) +
  scale_y_continuous(labels = percent_format(suffix="")) +
  # coord_flip() +
  # change legend position and size
  theme(legend.title=element_blank(), axis.ticks.x=element_blank(), legend.text=element_text(size=14)) +
  # guides(colour = guide_legend(override.aes = list(size=5))) +
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1), 
                  y = text.y, group = stage), position = position_dodge(width = 0.9), vjust = 0.5, size=5) +
  scale_fill_manual(values=myColors5[c("GP","BP","NPearly","NPlate","ENMFB")], guide = guide_legend(reverse=F)) +
  scale_x_discrete(labels=c("E13.5" = "Control", "Kif7 cKO" = "Kif7 cKO")) 
  
g

# tSNE
centers <- all_tsne %>% 
            dplyr::group_by(cluster2) %>% 
            summarize(X = median(X), Y = median(Y))
options(repr.plot.width=4, repr.plot.height=4)
g1 <- ggplot(all_tsne, aes(x=X, y=Y, color=cluster2)) +
    # facet_wrap( ~group , ncol=3) +
    # facet_grid(cols = vars(group)) +
    geom_point(size=0.3, alpha=1) +
    geom_density_2d(color='black', size=0.1, alpha=0.2) +
    geom_text(data = centers, mapping = aes(label = cluster2), size = 4, color="black") +
    labs(x = "tSNE-1",y = "tSNE-2", title = "") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    # change legend
    theme(legend.position = "none", legend.title=element_blank()) +
    theme(strip.background = element_rect(fill = "gray97", color = NA)) + # strip background color
    # change strip, 17 and bold
    theme(strip.placement = "outside", strip.text.x = element_text(face="bold", size = 17), #italic
          strip.text.y = element_text(face="plain", size = 11)) +
    theme(panel.spacing=unit(.3, "lines"),
          panel.border = element_rect(color = "black", fill = NA, size = 0.2,colour = "black")) + 
    # change xy axis and title text size
    # axis test: 10
    # axis title: 14, with vjust or hjust 
    theme(axis.text  = element_text(face="plain", angle=0, size = 10, color = "black"),
            axis.title.y =element_text(size = 14), axis.title.x =element_text(size = 14, vjust=-1)) +
    scale_color_manual(values=myColors5)
g1

options(repr.plot.width=5.2, repr.plot.height=3.5)
g1 <- ggplot(all_tsne, aes(x=X, y=Y, color=cluster2)) +
    # facet_wrap( ~group , ncol=3) +
    facet_grid(cols = vars(group)) +
    geom_point(size=0.3, alpha=1) +
    geom_density_2d(color='black', size=0.1, alpha=0.2) +
    geom_text(data = facet_centers, mapping = aes(label = cluster2), size = 4, color="black") +
    labs(x = "tSNE-1",y = "tSNE-2", title = "") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    # change legend
    theme(legend.position = "none", legend.title=element_blank()) +
    theme(strip.background = element_rect(fill = "gray97", color = NA)) + # strip background color
    # change strip, 17 and bold
    theme(strip.placement = "outside", strip.text.x = element_text(face="bold", size = 17), #italic
          strip.text.y = element_text(face="plain", size = 11)) +
    theme(panel.spacing=unit(.3, "lines"),
          panel.border = element_rect(color = "black", fill = NA, size = 0.2,colour = "black")) + 
    # change xy axis and title text size
    # axis test: 10
    # axis title: 14, with vjust or hjust 
    theme(axis.text  = element_text(face="plain", angle=0, size = 10, color = "black"),
            axis.title.y =element_text(size = 14), axis.title.x =element_text(size = 14, vjust=-1)) +
    scale_color_manual(values=myColors5)
g1

# trajectory
# brewer.pal(9,"Set3")[6:9]
options(repr.plot.width=4.3, repr.plot.height=3)
p <- plot_cell_trajectory(cds.sub, color_by = "lineage2", cell_size = 0.5, show_state_number=F) + # CellType, lineage
        # facet_wrap(~stage, nrow = 2) +
        # change legend position and size
        theme(legend.position = "right", legend.title=element_blank()) + #"none"
        guides(colour = guide_legend(override.aes = list(size=3))) +
        # scale_color_discrete(name = "", labels = c("Control", "Kif7 cKO")) +
        scale_color_manual(values=myColors5[c("BP", "GP", "NPearly","NPlate")]) +
        geom_vline(xintercept = -3.3, colour=myColors5["NPlate"], linetype="dashed", size=0.5) +
        geom_vline(xintercept = -1, colour=myColors5["NPearly"], linetype="dashed", size=0.5) +
        geom_vline(xintercept = 1, colour=myColors5["GP"], linetype="dashed", size=0.5) +
        theme(plot.title = element_text(size=22, face = "italic")) +
        theme(strip.placement = "outside", strip.text.x = element_text(face="italic", size = 22), #italic
          strip.text.y = element_text(face="plain", size = 11))
p

# brewer.pal(9,"Set3")[6:9]
options(repr.plot.width=4.2, repr.plot.height=6.5)
p <- plot_cell_trajectory2(cds.sub, color_by = "lineage2", cell_size = 2, use.cells = cells) + # CellType, lineage
        facet_wrap(~stage, nrow = 2) +
        theme(legend.position = "none", legend.title=element_blank()) + #"none"
        # scale_color_discrete(name = "", labels = c("Control", "Kif7 cKO")) +
        guides(colour = guide_legend(override.aes = list(size=3))) +
        scale_color_manual(values=myColors5) +
        # add lines
        geom_vline(xintercept = -3.3, colour=myColors5["late NP"], linetype="dashed", size=0.5) +
        geom_vline(xintercept = -1, colour=myColors5["NP"], linetype="dashed", size=0.5) +
        geom_vline(xintercept = 1, colour=myColors5["GP"], linetype="dashed", size=0.5) +
        # change lengend
        theme(plot.title = element_text(size=22, face = "plain")) +
        # change strip, 17 and bold
        theme(strip.placement = "outside", strip.text.x = element_text(face="bold", size = 17), #italic
          strip.text.y = element_text(face="plain", size = 11)) +
        # change xy axis and title text size
        # axis test: 10
        # axis title: 14, with vjust or hjust
        theme(axis.text  = element_text(face="plain", angle=0, size = 10, color = "black"),
            axis.title.y =element_text(size = 14), axis.title.x =element_text(size = 14, vjust=-1)) +
        # add pct text
        geom_text(data=label.data.n1,aes(x=x, label=label), y=0, color="#ff309a", fontface="bold", size=6) +
        geom_text(data=label.data.n2,aes(x=x, label=label), y=0, color="#FCCDE5", fontface="bold", size=6) +
        geom_text(data=label.data.n3,aes(x=x, label=label), y=0, color="#B3DE69", fontface="bold", size=6) +
        geom_text(data=label.data.n4,aes(x=x, label=label), y=0, color="#FDB462", fontface="bold", size=6) +
        # add subtype text
        geom_text(data=label.data.s1,aes(x=x, label=label), y=1.1, color="#ff309a", fontface="bold", size=5) +
        geom_text(data=label.data.s2,aes(x=x, label=label), y=1.1, color="#FCCDE5", fontface="bold", size=5) +
        geom_text(data=label.data.s3,aes(x=x, label=label), y=1.1, color="#B3DE69", fontface="bold", size=5) +
        geom_text(data=label.data.s4,aes(x=x, label=label), y=1.1, color="#FDB462", fontface="bold", size=5) +
        # add x,y limit
        scale_x_continuous(limits = c(-4.3, 2.5)) +
        scale_y_continuous(limits = c(-0.8, 1.1))
p

# violin plot
options(repr.plot.width=6.5, repr.plot.height=5/6*8)
plot <- ggplot(exprM_melt, aes(x=variable, y=log2(1+value))) + 
    geom_hline(yintercept = c(2,4,6), linetype="dashed", color = "grey80", size=0.3) +
    facet_wrap( ~ group.sub, ncol=1, labeller = label_context, scales = "free_y",strip.position = "left") +
    labs(x = "", y = "log2(UMIs+1)\n") +
    theme_bw() +
    # add background color to show the significance
    geom_rect(data=subset(exprM_melt, group.sub %in% c("Control BP","Kif7 cKO BP")), xmin=0.4, xmax=1.5+1, ymin=-Inf, ymax=Inf, fill="grey90", alpha=1)+
    geom_rect(data=subset(exprM_melt, group.sub %in% c("Control BP","Kif7 cKO BP")), xmin=0.5+7, xmax=1.5+7, ymin=-Inf, ymax=Inf, fill="grey90", alpha=1)+
    geom_rect(data=subset(exprM_melt, group.sub %in% c("Control NPearly","Kif7 cKO NPearly")), xmin=0.4, xmax=1.5+3, ymin=-Inf, ymax=Inf, fill="grey90", alpha=1)+
    geom_rect(data=subset(exprM_melt, group.sub %in% c("Control NPlate","Kif7 cKO NPlate")), xmin=0.4, xmax=1.5+6, ymin=-Inf, ymax=Inf, fill="grey90", alpha=1)+
    geom_rect(data=subset(exprM_melt, group.sub %in% c("Control GP","Kif7 cKO GP")), xmin=0.5+1, xmax=1.5+1, ymin=-Inf, ymax=Inf, fill="grey90", alpha=1)+
    geom_rect(data=subset(exprM_melt, group.sub %in% c("Control GP","Kif7 cKO GP")), xmin=0.5+7, xmax=2.5+7, ymin=-Inf, ymax=Inf, fill="grey90", alpha=1)+
    #
    geom_violin(trim = F, na.rm = T, aes(fill = group.sub),colour = "black", scale="width", size=0.3, bw=0.3) + 
    scale_y_continuous(position="right", limits=c(0, 8), breaks = seq(0, 8, length.out = 3))+ #
    theme(strip.background = element_blank(),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.spacing=unit(.4, "lines"),panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
    theme(axis.text.x  = element_text(face="italic", angle=90, size = 18, color = "black", vjust=0.5),
        axis.text.y  = element_text(face="plain", size = 8, color = "black"),
        axis.title =element_text(size = 18)) +
    theme(strip.background = element_rect(fill = NA, color = NA), legend.position = "none") +
    theme(strip.placement = "outside") +
    theme(strip.text.y = element_text(angle=90, margin = margin(1,1,1,1, "mm"), size = 5, face="bold"), 
          axis.ticks.x=element_blank()) 
plot

# boxplot
options(repr.plot.width=7, repr.plot.height=20)
p <- ggplot(merge.box.data, aes(x=variable, y=pseudotime, fill=log2FC)) +
  # add background color
  geom_rect(xmin=Inf, xmax=-Inf, ymin=-0.45, ymax=0.3, fill=alpha("#B3DE69", 0.01), alpha=0.01)+
  geom_rect(xmin=Inf, xmax=-Inf, ymin=-Inf, ymax=-0.45, fill=alpha("#FDB462", 0.01), alpha=0.01)+
  geom_rect(xmin=Inf, xmax=-Inf, ymin=0.3, ymax=0.85, fill=alpha("#FCCDE5", 0.01), alpha=0.01)+
  geom_rect(xmin=Inf, xmax=-Inf, ymin=0.85, ymax=Inf, fill=alpha("#ff309a", 0.01), alpha=0.01)+
  geom_hline(yintercept = 0, color = "black", size=0.5) +
  geom_boxplot(outlier.shape = NA,  size=0.2) + # 
  labs(x = "", y = "\nGlial lineage <- BP -> Neuronal lineage\nAdjusted pseudotime", title = "") +
  scale_y_continuous(limits = c(-1, 1)) +
  # this theme will set many things to default
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.background = element_rect(fill = NA), panel.grid.major = element_line(colour = "grey90"),
        panel.ontop = F) +
  # theme(panel.background = element_blank(), panel.grid.major = element_blank(),panel.ontop = F) +
  theme(legend.title=element_text(size=13), legend.text=element_text(size=10)) +
  coord_flip() +
  scale_fill_gradient2("Expression", midpoint=0, mid="white", low="blue", high="red") +
  # theme(axis.text.x  = element_text(face="italic", angle=0, size = 10, color = "black", vjust=0.6)) +
  theme(axis.text.x  = element_text(face="plain", angle=0, size = 10, color = "black"),
        axis.text.y  = element_text(face="italic", angle=0, size = 11, colour = c.color),
        axis.title =element_text(size = 15, face="bold"))
p

# barplot
# show color in pseudotime
plot.GO.barplot.time <- function (barplot_df, time) 
{
    barplot_df$mean.time <- barplot_df[,time]
    library(Hmisc)
    library(stringr)
    library(RColorBrewer)
    for (i in 1:dim(barplot_df)[1]) {
        barplot_df[i, ]$Description <- capitalize(as.character(barplot_df[i, 
            ]$Description))
    }
    barplot_df <- barplot_df[order(barplot_df$mean.time, decreasing = F), ]
    barplot_df$Description <- factor(barplot_df$Description, levels = rev(barplot_df$Description))
    # maxpvalue <- max(-log10(barplot_df$pvalue)) + 0.5
    maxCount <- max(barplot_df$Count)
    
    g <- ggplot(data = barplot_df, aes(x = Description, y = Count)) + 
            geom_bar(stat = "identity", aes(fill = mean.time), alpha = 0.8) + 
            # add annotation upon bars, like pvalue
            # paste("P=",round(p.adjust/20, digits = 3), sep="")
            geom_text(aes(label = ifelse(p.adjust/10<0.05, "**", "*")), color = "black", vjust = 0.4, 
                      hjust = -0.5, size = 5, fontface = "plain") + 
            # ylim(0, maxpvalue * 1.1) + 
            ylim(0, maxCount * 1.2) + 
            coord_flip() + 
            labs(x = "", y = "Gene count", title = "") + 
            theme_bw() + 
            # theme(legend.position = "none") + 
            # change legend position and size
            theme(legend.title=element_text(size=13), legend.text=element_text(size=10)) +
            theme(axis.text.y = element_text(size = 18, color = "black", face = "plain"), 
                  axis.text.x = element_text(size = 10, color = "black", face = "plain"), 
                  axis.title = element_text(size = 14)) + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                plot.margin = unit(c(0, 0, 0, 0), "cm"), panel.border = element_blank()) + 
            theme(axis.line = element_line(color = "black")) + 
            scale_x_discrete(labels = function(x) str_wrap(x, width = 25)) +
            scale_fill_gradient("Activation time", low="red", high="blue", trans="reverse")
    g
}
# p-values (* , **, *** indicate statistical significance at p < 0.10, p<0.05 and p<0.01)
options(repr.plot.width=8, repr.plot.height=7)
p <- plot.GO.barplot.time(tmp.GO.core, time="mytime")
p

# Ezh2
options(repr.plot.width=5, repr.plot.height=5.2)
ggplot(data=tmp.melt, aes(x=variable, y=value, group=group, color=group, fill=group)) +
  # geom_line()+
  # geom_point() +
  geom_bar(stat="identity", width=0.5, position=position_dodge()) +
  theme_bw() + 
  # geom_hline(yintercept = 0, linetype="dashed", color = "grey20", size=0.3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "\nEmbryonic day",y = "Log2FC", title = "Ezh2 expression") + 
  theme(axis.text.x  = element_text(face="plain", angle=0, size = 10, color = "black", vjust=-0.5),
        axis.text.y  = element_text(face="plain", size = 10, color = "black", hjust=0),
        axis.title =element_text(size = 14),
        plot.title = element_text(size=22)) +
  # scale_x_continuous(breaks = c(9.5, 10.5, 11.5, 12.5, 13.5, 16.5, 19.0, 21.0), name = "\nEmbryonic day") +
  theme(legend.position = c(0.27, 0.9), legend.title=element_blank(), legend.text=element_text(size=15)) +
  scale_fill_manual(values=c("blue", "red")) +
  scale_color_manual(values=c("blue", "red"))

  



