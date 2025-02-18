library(stringr)
library(readr)
library(tidyr)
library(glmnet)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(ggpubr)
library(readxl)

# Data preparation
{
    load(file="../apaRaw.Rdata")
    load(file="../apaEvent.Rdata")
    apaRaw_rename = apaRaw[,grep("Ca",colnames(apaRaw))]
    colnames(apaRaw_rename) = colnames(apaRaw_rename) %>% substr(.,1,nchar(.)-3)
    apaEvent_PDUI = apaRaw_rename[apaEvent$Gene,]
    apaEvent_PDUI = data.matrix(apaEvent_PDUI)
    apaEvent_PDUI = t(apaEvent_PDUI)

    load(file="../mutationMatrix.Rdata")
}

# Clean
{
    common_id = intersect(rownames(apaEvent_PDUI),rownames(mutation_matrix))
    # Common samples
    mutation_matrix = mutation_matrix[common_id,]
    apaEvent_PDUI = apaEvent_PDUI[common_id,]

    # Mutation table
    mutation_matrix = data.matrix(mutation_matrix)
    mutation_matrix[which(is.na(mutation_matrix))] = 0
    mutation_matrix = mutation_matrix[,names(which(apply(mutation_matrix,2,function(x) sum(as.numeric(x)))!=0))]
    # Filter mutation
    rm_col = apply(mutation_matrix,2,sum)
    rm_col = rm_col[-which(rm_col<3)]
    mutation_matrix = mutation_matrix[,names(rm_col)]
}

# LASSO
{
    library(ipflasso)
    lasso_result = data.frame(matrix(nrow=ncol(mutation_matrix),ncol=ncol(apaEvent_PDUI)))
    rownames(lasso_result) = colnames(mutation_matrix)
    colnames(lasso_result) = colnames(apaEvent_PDUI)
    lasso_result[is.na(lasso_result)]=0

    for(i in 1:ncol(lasso_result))
    {
        print(i)
        lambdas = 10^seq(2, -3, by = -.1)
        # Setting alpha = 1 implements lasso regression
        lasso_matrix = cbind(apaEvent_PDUI[,i],mutation_matrix)
        lasso_matrix = na.omit(lasso_matrix)

        lasso_reg = cvr.glmnet(X=lasso_matrix[,2:ncol(lasso_matrix)],Y=lasso_matrix[,1],alpha = 1,
                                standardize = F, nfolds = 10, ncv=10, type.measure="mse", family = "gaussian")

        lambda_best = lasso_reg$lambda[order(lasso_reg$cvm)[1]]

        lasso_model = glmnet(x=lasso_matrix[,2:ncol(lasso_matrix)], y=lasso_matrix[,1], alpha = 1, lambda = lambda_best, standardize = F,type.measure="mse")
        lasso_coef = coef(lasso_model)
        lasso_coef = rbind(lasso_coef@Dimnames[[1]][lasso_coef@i+1][-1],lasso_coef@x[-1])
        lasso_coef = rbind(lasso_coef,abs(as.numeric(lasso_coef[2,])))
        lasso_snv = lasso_coef[1,which(lasso_coef[3,]>0)]

        for(j in lasso_snv)
        {
            lasso_result[j,i] = 1
        }
    }
}

# Mutation annotation
{
    load(file="../maf.Rdata")
    library(randomcoloR)
    SNV_frequence = data.frame(apply(lasso_result,1,sum))
    colnames(SNV_frequence) = "counts"
    SNV_frequence$refGene = maf_file[match(rownames(SNV_frequence),maf_file$SNV),"Hugo_Symbol"]
    SNV_frequence$refRegion = maf_file[match(rownames(SNV_frequence),maf_file$SNV),"Variant_Classification"]
    region_anno = data.frame(table(SNV_frequence$refRegion))
    region_anno$Per = region_anno$Freq/sum(region_anno$Freq)
    region_anno$Var1 = factor(region_anno$Var1, levels = region_anno$Var1[order(region_anno$Freq,decreasing = T)])
    ggplot(region_anno, aes(x = Var1, y = Per)) +
        geom_bar(stat = "identity", fill = brewer.pal(name="Paired", n=10)) +
        my_theme(title = 10,axis.text.x = element_text(angle = 45, hjust = 1))+
        labs(x="",y="Percentage")
}
