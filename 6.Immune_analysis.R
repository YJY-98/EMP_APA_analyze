library(stringr)
library(readr)
library(tidyr)
library(ggplot2)
library(readxl)
library(ggrepel)

# Data preparation
{
    load(file = "../tpm.Rdata")
    load(file = "../deseq.Rdata")
    final_res = final_res[final_res$gene_biotype=="protein_coding",]
    tpm = tpm[intersect(rownames(final_res),rownames(tpm)),]
    rownames(tpm) = final_res[rownames(tpm),"external_gene_name"]
    log_tpm = log2(tpm+1)

    log_tpm_Ca = log_tpm[,grep("Ca",colnames(log_tpm))]
    log_tpm_Cj = log_tpm[,grep("Cj",colnames(log_tpm))]

    normalized_res = log_tpm_Ca - apply(log_tpm_Cj,1,mean)
    normalized_res = normalized_res %>% data.frame(check.names = F) %>% rownames_to_column()
}

# Immune analysis
{
    load(file = "../apaRaw.Rdata")
    load(file = "../apaEvent.Rdata")
    apaRaw_T = apaRaw[,grep("Ca",colnames(apaRaw))] 
    apaRaw_T = apaRaw_T[rownames(apaEvent),]
    
    apaRaw_T = apaRaw_T %>% group_by(Gene)
    apaRaw_T = apaRaw_T %>% summarise_all(mean,na.rm=T) %>% data.frame(check.names = F)
    rownames(apaRaw_T) = apaRaw_T$Gene
    apaRaw_T = apaRaw_T[,-1]
    
    cor_res = data.frame(matrix(nrow = nrow(apaRaw_T) , ncol = 3))
    colnames(cor_res) = c("event", "cor", "pvalue")

    # Example
    event = "ENST00000227507.2|CCND1|chr11|+"
    cor_event = cor.test(Tide_res$TIDE %>% as.numeric(), apaRaw_T[event,] %>% as.numeric())
    cor_res[i,1] = rownames(apaRaw_T)[i]
    cor_res[i,2] = cor_event$estimate
    cor_res[i,3] = cor_event$p.value

    cor_res$asb_cor = abs(cor_res$cor) 
    cor_res$log_pvalue = -log2(cor_res$pvalue)
    cor_res$group = ifelse(
        cor_res$pvalue < 0.05, "imp", "normal"
    )
    
    # Plot
    ggplot(cor_res, aes(asb_cor, log_pvalue, color = group)) +
        geom_point(size = 2.0, shape = 16, show.legend = F,alpha=1) +
        scale_color_manual(values=c("#FF0000","#000000"))+
        geom_text_repel(inherit.aes = F, data = cor_res, aes(x=asb_cor, y=log_pvalue, label=label), size=3,
                        segment.size=0.5, nudge_x=50, direction = "y", hjust = 0)+
        geom_hline(yintercept = 4.321928,linetype="dashed")+
        geom_vline(xintercept = 0.3,linetype="dashed")
}