library(DESeq2)
library(stringr)
library(biomaRt)

# Data loading
{
    dir = "../counts/"
    files = list.files(path = dir)
    expr = lapply(files,
                function(x){
                    expr <- read.table(file = file.path(dir,x), header = T,comment.char = "#")[,c(1:7)]
                    return(expr)
                })

    df = do.call(cbind, expr)
    df = df[,c(1:6,seq(7,ncol(df),by=7))] 

    # Rename columns
    for (i in 7:ncol(df)){
    colnames(df)[i] = str_split_fixed(colnames(df)[i],"[.]",7) %>% .[,6]
    }

    # Remove ID with _PAR_Y suffix 
    df = df[-which(substr(df$Geneid,nchar(df$Geneid)-5,nchar(df$Geneid))=="_PAR_Y"),]
    rownames(df) = str_split_fixed(df[,1],'[.]',2)[,1]
}

# DESeq2 pipeline
{
    readMatrix = df[,-c(1:6)]
    condition = factor(c(rep("MM_Ca",grep("CA",colnames(readMatrix)) %>% length()),
                       rep("MM_Cj",grep("CJ",colnames(readMatrix)) %>% length())))
    colData = data.frame(row.names=colnames(readMatrix), condition)
    dds = DESeqDataSetFromMatrix(readMatrix, colData, design= ~ condition)
    dds = DESeq(dds)

    res = results(dds,contrast=c("condition","MM_Ca","MM_Cj"))

    res = res[order(res$pvalue),]
    resNA = subset(res,pvalue >= 0)
    resNA = subset(resNA,padj >= 0)
    deMatrix = data.frame(resNA)
    deMatrix$ensembl_gene_id = rownames(deMatrix)

    # Get gene information
    my_mart = useMart("ensembl")
    my_dataset = useDataset("hsapiens_gene_ensembl",mart = my_mart)
    my_newid = getBM(attributes = c("ensembl_gene_id","external_gene_name","description","chromosome_name","gene_biotype"),filters = "ensembl_gene_id",values = deMatrix$ensembl_gene_id,mart = my_dataset)

    final_res = merge(my_newid, deMatrix, by = intersect(names(my_newid), names(deMatrix)))
    final_res = final_res[-which(final_res$external_gene_name==""),]
    
    final_res$Group="not-significant"
    final_res$Group[which((final_res$log2FoldChange < -0.5)&(final_res$padj < 0.05))]="down-regulated"
    final_res$Group[which((final_res$log2FoldChange > 0.5)&(final_res$padj < 0.05))]="up-regulated"
}