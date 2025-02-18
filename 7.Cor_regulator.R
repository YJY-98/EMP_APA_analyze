library(ggplot2)
library(ggpubr)
library(ggthemes)
library(ggrepel)
library(stringr)
library(dplyr)

# DEGs
{
    load(file = '../deseq.Rdata')
    final_res$Label = ""
    final_res = final_res[order(final_res$padj),]
    up.genes = head(final_res$external_gene_name[which(final_res$Group=="up-regulated")],5)
    down.genes = head(final_res$external_gene_name[which(final_res$Group=="down-regulated")],5)
    deg.top10.genes = c(as.character(up.genes),as.character(down.genes))
    final_res$Label[match(deg.top10.genes,final_res$external_gene_name)] = deg.top10.genes
    final_res$logP = -log10(final_res$padj)

    # Plot
    ggscatter(final_res,x="log2FoldChange",y="logP",
                color="Group",
                palette=c("#0000FF","#BBBBBB","#FF0000"),size=1,
                alpha=1)+
        labs(x="log2(Tumor/Normal)",
            y="-log10(Adjust P-value)")+
        geom_hline(yintercept = 1.3,linetype="dashed")+
        geom_vline(xintercept = c(-1,1),linetype="dashed")+
        xlim(-10,10)+
        coord_fixed(20/40)
}

# Correlation analysis
{
    library("ggpubr")
    load(file = "../apaRaw.Rdata")
    load(file = "../tpm.Rdata")

    log_tpm = log2(tpm+1)

    Apa_expr_cor_data=data.frame(matrix(ncol=3,nrow=dim(log_tpm)[2]))
    colnames(Apa_expr_cor_data)=c("PDUI","TPM")
    event="ENST00000227507.2|CCND1|chr11|+"

    # Example
    RBP="PPP1CB"
    Apa_expr_cor_data$PDUI=as.numeric(apaRaw[event,grep("Ca",colnames(apaRaw))])
    Apa_expr_cor_data$TPM=as.numeric(log_tpm[RBP,colnames(apaRaw)[grep("Ca",colnames(apaRaw))]])

    # Plot
    ggplot(Apa_expr_cor_data, aes(PDUI, TPM)) +
    geom_point() +
    geom_smooth(method = lm) +
    stat_cor(method = "spearman", label.x = 0)+
    labs(title="",x="PDUI of CCND1", y = "PPP1CB expression")
}
