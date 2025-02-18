library(stringr)
library(readr)

# Data loading
{
    apaRaw = read.table("data/homo_mm_chr.txt",header = TRUE)
    # Rename columns
    for (i in 5:ncol(apaRaw)){
    apaRaw[,i] = as.character(apaRaw[,i])
    apaRaw[,i] = as.numeric(apaRaw[,i])
    colnames(apaRaw)[i] = colnames(apaRaw)[i] %>% str_split_fixed(.,"[.]",7) %>% .[,5]
    }
    rownames(apaRaw) = apaRaw$Gene
}

# APA detection
{
    apaMatrix = apaRaw[,-c(2:4)]
    cacols = grep("CA",colnames(apaMatrix))
    cjcols = grep("CJ",colnames(apaMatrix))

    apaMatrix$caMean = apply(apaMatrix[,cacols], 1, mean,na.rm=T)
    apaMatrix$caSD = apply(apaMatrix[,cacols],1,sd,na.rm=T)
    apaMatrix$caMR = apply(apaMatrix[,cacols],1,function(x) sum(is.na(x))/length(x))
    apaMatrix$cjMean = apply(apaMatrix[,cjcols], 1, mean,na.rm=T)
    apaMatrix$cjSD = apply(apaMatrix[,cjcols],1,sd,na.rm=T)
    apaMatrix$cjMR = apply(apaMatrix[,cjcols],1,function(x) sum(is.na(x))/length(x))
    apaMatrix$delta_PDUI = apaMatrix$caMean-apaMatrix$cjMean

    caChosen = intersect(which(apaMatrix$caMR<0.2),which(apaMatrix$caSD>0.05))
    cjChosen = intersect(which(apaMatrix$cjMR<0.2),which(apaMatrix$cjSD>0.05))
    apaMatrix = apaMatrix[intersect(caChosen,cjChosen),]

    # Permutation test
    library(coin)
    permutation_test=function(x)
    {
        group = factor(c(rep("CA",length(cacols)),rep("CJ",length(cjcols))))
        x = as.numeric(x[c(cacols,cjcols)])
        data = data.frame(group,x)
        tryCatch ({p<-oneway_test(x ~ group, data = data, distribution = approximate(10000),alternative = "two.side")
        return(pvalue(p))}, 
        warning = function (w){return(NA)},
        error = function(e) {return(NA)})
    }

    p.value = apply(apaMatrix,1,permutation_test)
    apaMatrix$pvalue = p.value
    apaMatrix$fdr = p.adjust(apaMatrix$pvalue,"fdr")

    # Remove genes with no significant difference
    apaMatrix = apaMatrix[,-grep("C[AJ]",colnames(apaMatrix))]

    apaMatrix$Gene = as.character(apaMatrix$Gene)
    apaMatrix$external_gene_name = str_split_fixed(apaMatrix$Gene,'[|]',4)[,2]
    apaMatrix$transcript_ID = str_split_fixed(apaMatrix$Gene,'[|]',4)[,1]
    apaMatrix$transcript_ID = apaMatrix$transcript_ID %>% str_split_fixed(.,"[.]",2) %>% .[,1]

    apaMatrix$Group = "not-significant"
    apaMatrix$Group[intersect(which(apaMatrix$delta_PDUI>0.1),
                            which(apaMatrix$fdr<0.05))] = "Lengthen"
    apaMatrix$Group[intersect(which(apaMatrix$delta_PDUI< -0.1),
                            which(apaMatrix$fdr<0.05))] = "Shorten"
}
