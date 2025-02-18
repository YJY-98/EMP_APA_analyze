library(stringr)
library(ggsignif)
library(ggplot2)

# Event num diff
{
    load(file="../apaRaw.Rdata")
    APA_num_ggplot = data.frame(apply(apaRaw[,5:ncol(apaRaw)],2,function(x) sum(!is.na(x))))
    colnames(APA_num_ggplot) = c("count","group")
    ggplot(APA_num_ggplot, aes(x=group, y=count, fill=group)) + 
        geom_violin(trim=FALSE,show.legend = F,)+
        geom_boxplot(width=0.2, fill="white",show.legend = F)+
        geom_signif(comparisons=list(c("Ca","Cj")),
                    test = 'wilcox.test',
                    step_increase = 0.1,
                    map_signif_level = TRUE,
                    vjust = 0.5)+
        scale_fill_manual(values=c("#FF0000", "#0000FF"))+
}

# PDUI value diff
{
    load(file="../apaRaw.Rdata")
    PDUI_ave_ggplot = data.frame(apply(apaRaw[,5:ncol(apaRaw)],2,function(x) sum(mean(x,na.rm=T))))
    PDUI_ave_ggplot$group = str_split_fixed(rownames(PDUI_ave_ggplot),"[-]",3)[,3]
    ggplot(PDUI_ave_ggplot, aes(x=group, y=count, fill=group)) + 
        geom_violin(trim=FALSE,show.legend = F)+
        geom_boxplot(width=0.2, fill="white",show.legend = F)+
        labs(title="",x="", y = "PDUI") +
        geom_signif(comparisons=list(c("Ca","Cj")),
                    test = 'wilcox.test',
                    step_increase = 0.1,
                    map_signif_level = TRUE,
                    vjust = 0.5)+
        scale_fill_manual(values=c("#FF0000", "#0000FF"))
}

# Diff APA event
{
    library(ggpubr)
    library(ggthemes)
    library(ggrepel)
    load(file="../apaEvent.Rdata")

    shorten.genes = apaEvent[order(apaEvent$delta_PDUI)[1:5],"external_gene_name"]
    lengthen.genes = apaEvent[order(apaEvent$delta_PDUI)[(nrow(apaEvent)-4):nrow(apaEvent)],"external_gene_name"]
    apaEvent$Label = ""
    apaEvent[order(apaEvent$delta_PDUI)[1:5],"Label"]=apaEvent[order(apaEvent$delta_PDUI)[1:5],"external_gene_name"]
    apaEvent[order(apaEvent$delta_PDUI)[(nrow(apaEvent)-4):nrow(apaEvent)],"Label"]=apaEvent[order(apaEvent$delta_PDUI)[(nrow(apaEvent)-4):nrow(apaEvent)],"external_gene_name"]

    df_abline = data.frame(intercept=c(0.15,-0.15),slope=c(1,1))
    ggscatter(apaMatrix,x="cjMean",y="caMean",
            color="Group",
            palette=c("#FF0000","#BBBBBB","#0000FF"),size=2,
            alpha=1)+
    labs(x="PDUI in Normal",
            y="PDUI in Tumor")+
    geom_abline(data = df_abline, aes(intercept = intercept,slope = slope),linetype = 2)+
    xlim(c(0,1))+ylim(c(0,1))
}

# Event enrichment
{
    library(readr)
    library(dplyr)
    load(file="../enrichedPathway.Rdata")

    pathway = arrange(pathway,pathway$Group,-pathway$count)
    pathway = pathway[rev(1:nrow(pathway)),]
    rownames(pathway) = 1:nrow(pathway)
    GO_term_order = factor(as.integer(rownames(pathway)),labels = pathway$Description)

    ggplot(data=pathway,aes(x=GO_term_order,y=count,fill=Group))+
    geom_bar(stat = "identity",width = 0.8)+
    scale_fill_manual(values =c("#F1756D","#2D3091"))+
    coord_flip()+
    xlab("")+
    ylab("Gene counts")+
    labs()
}

# Heatmap plot
{
    library(RColorBrewer)
    library(pheatmap)
    load(file="../apaHeatmap.Rdata")

    pheatmap(apa_heatmap,
            color=colorRampPalette(c("#428CAE","white","#E15453"))(50),
            clustering_method="ward.D2",
            cluster_rows = F,cluster_cols = F,
            fontsize=10, fontsize_row=5,
            annotation_col = anno_name,annotation_colors = anno_col,
            show_rownames = T,fontsize_col = 10,show_colnames = F)
}
