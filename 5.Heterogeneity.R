library(stringr)
library(dplyr)
library(tidyr)
library(readr)
library(ggsignif)
library(ggplot2)
library(ggthemes)

load(file = "../apaMatrix.Rdata")
load(file = "../apaRaw.Rdata")

apaMatrix$caVar = apaMatrix$caSD ** 2
apaMatrix$cjVar = apaMatrix$cjSD ** 2
apaMatrix$delta_Var = apaMatrix$caVar - apaMatrix$cjVar

# Example
{
    event = "ENST00000227507.2|CCND1|chr11|+"
    trans = strsplit(event,"[|]")[[1]][1]
    gene = strsplit(event,"[|]")[[1]][2]
    col = c(grep("Ca",colnames(apaRaw)),grep("Cj",colnames(apaRaw)))
    density_table = apaRaw[match(event,apaRaw$Gene),col]

    density_table = density_table %>% pivot_longer(cols = contains(c("Ca", "Cj")),
                                                    names_to = "id",
                                                    values_to = "value")
    density_table$type = str_split_fixed(density_table$id,"[-]",3)[,3]

    ggplot(density_table, aes(x = value))+
        geom_density(aes(fill = type), alpha = 0.7)+
        labs(x = "PDUI", y = "Frequency density",title = gene)+
        scale_x_continuous(limits = c(0,1),breaks = c(0,0.2,0.4,0.6,0.8,1))+
        scale_fill_manual(values = c("#FF0000", "#0000FF"),
                        name = "",
                        breaks = c("Ca", "Cj"),
                        labels = c("Tumor", "Normal"))
}

# Heterogeneity
{
    library(ggpubr)
    library(ggthemes)
    library(ggrepel)

    apaMatrix$Var_group="not-significant"
    apaMatrix$Var_group[which(apaMatrix$delta_Var>0.025)]="High_variance_in_tumor"
    apaMatrix$Var_group[which(apaMatrix$delta_Var< -0.025)]="High_variance_in_normal"
    apaMatrix$number = rank(apaMatrix$delta_Var, ties.method = "random")

    apaMatrix = apaMatrix %>% mutate(
        label = case_when(
        (external_gene_name %in% intersect(apa_heatmap_genes,ICB_related)) & (Var_group=="High_variance_in_tumor") ~ external_gene_name
        )
    )

    df_abline <- data.frame(intercept=c(0.025,-0.025),slope=c(0,0))
    ggscatter(apaMatrix,x="number",y="delta_Var",
                color="Var_group",
                palette=c("#0000FF", "#FF0000", "#BBBBBB"),size=1,
                alpha=1)+
        labs(title = table(apaMatrix$Var_group) %>% str_c(names(.),sep = ":") %>% str_c(.,collapse = " "))+
        geom_text_repel(inherit.aes = F, data = apaMatrix, aes(x=number, y=delta_Var, label=label), size=3,
                        segment.size=0.5, nudge_x=100, direction = "y", hjust = 2)+
        geom_abline(data = df_abline, aes(intercept = intercept,slope = slope),linetype = 2)
}

