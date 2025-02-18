library(dplyr)
library(readr)
library(readxl)

# Risk score analysis
{
    load(file = "../risk_genes.Rdata")
    load(file = "../apaRaw.Rdata")

    psd=function(x)
    {
        x=x[!is.na(x)]
        sqrt(sum((x - mean(x))^2)/length(x))
    }

    apa_score = apaRaw[risk_genes,]
    apa_score = apa_score[,c(1,grep("Ca",colnames(apa_score)))]
    apa_score = apa_score %>% group_by(Gene)
    apa_score = apa_score %>% summarise_all(mean,na.rm=T)
    rownames(apa_score) = apa_score$Gene
    apa_score = apa_score[,-1]

    for(i in 1:nrow(apa_score))
    {
        apa_score[i,]= (apa_score[i,]-mean(as.numeric(apa_score[i,]),na.rm = T))/psd(apa_score[i,])
    }

    apa_score_sample = apply(apa_score,2,function(x) mean(x,na.rm=T))
    apa_score_sample = apa_score_sample*-1
}

# Survival analysis
{
    library(survival)
    library(survminer)
    library(readxl)

    survival = read_xlsx("../data_clinical.xlsx")

    survival_table = survival
    survival_table = survival_table %>% select(
    c("ID","Months","Status")
    )
    survival_table = survival_table[!is.na(survival_table$ID),]
    survival_table = survival_table[!is.na(survival_table$Months),]

    #plot
    survival_table$apa_score = score_table[match(survival_table$ID,score_table$MM_num),"score"]
    sur_res = survival_table %>% mutate(
    group = ifelse(apa_score > quantile(apa_score,na.rm = T,0.5), "high", "low")
    ) %>% 
    filter(!is.na(apa_score))

    survobj = with(sur_res,Surv(Months,Status))
    sur_plot = survfit(survobj~group,data = sur_res)
    sur_plot = ggsurvplot(sur_plot, risk.table = TRUE,
                        pval = TRUE,
                        pval.method = TRUE,
                        linetype = "strata", # Change line type by groups
                        palette = c("#FF0000", "#000000")
    )
}
