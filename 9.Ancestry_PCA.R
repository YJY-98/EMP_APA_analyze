library(stringr)
library(dplyr)
library(readr)
library(ggsignif)
library(ggplot2)
library(ggthemes)

# Data
{
    load(file = "../apaDiff.Rdata")
    load(file = "../group.Rdata")

    cmp_col = grep("Ca",colnames(combine_apaRaw))
    tmp_col = grep("TCGA",colnames(combine_apaRaw))

    combine_apaRaw$cmp_sd = apply(combine_apaRaw[,cmp_col],1,sd,na.rm=T)
    combine_apaRaw$cmp_mr = apply(combine_apaRaw[,cmp_col],1,function(x) sum(is.na(x))/length(x))

    combine_apaRaw$tmp_sd = apply(combine_apaRaw[,tmp_col],1,sd,na.rm=T)
    combine_apaRaw$tmp_mr = apply(combine_apaRaw[,tmp_col],1,function(x) sum(is.na(x))/length(x))

    cmp_chosen = intersect(which(combine_apaRaw$cmp_mr<0.2),which(combine_apaRaw$cmp_sd>0.05))
    tmp_chosen = intersect(which(combine_apaRaw$tmp_mr<0.2),which(combine_apaRaw$tmp_sd>0.05))
    combine_apaRaw = combine_apaRaw[intersect(cmp_chosen,tmp_chosen),]
}

# Impute
{
    library(missForest)
    library(doParallel)
    df = combine_apaRaw
    rownames(df) = df$Gene
    df = df[,-c(1,2)]
    df = t(df)
    #impute
    set.seed(111)
    registerDoParallel(cores=20)
    df_impute = missForest(df,ntree = 100, parallelize = "forests")
}

# PCA
{
    library(ggbiplot)
    library(ggthemes)
    library(patchwork)
    df_impute = df_impute$ximp
    df_impute = df_impute[c(grep("MM",rownames(df_impute)),
                            grep("TCGA",rownames(df_impute))),]

    dat.pca = prcomp(df_impute, scale. = TRUE)

    TB = data.frame(dat.pca$x)
    TB$class = group_list

    p0 = ggbiplot(dat.pca, obs.scale = 1, var.scale = 1,
                    groups = group_list, ellipse = T, circle = F ,var.axes = F) +
        #scale_color_discrete(name = '') +
        scale_color_manual(values = c("#67C2A5","#F58B63")) +
        my_theme()+theme(legend.position = 'none',
                        axis.title = element_blank())
}