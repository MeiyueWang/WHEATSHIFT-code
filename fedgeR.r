library(edgeR)

fedgeR <- function(countTable){
  group <- c(1:2)
  dge <- DGEList(counts=countTable,group = group)
  keep <- rowSums(cpm(dge)>1) >= 1
  dge <- dge[keep,,keep.lib.sizes=F]
  dge <- calcNormFactors(dge)
  dge_bcv <- dge
  bcv <- 0.2
  et <- exactTest(dge_bcv,dispersion = bcv^2)
  results <- et$table
  results$change <- ifelse(results$PValue<0.05,ifelse(results$logFC>1,"up",ifelse(results$logFC< -1,"down","nochange")),"nochange")
  return(results)
}

args <- commandArgs(T)
CKfile <- args[1]
treatfile <- args[2]
outfile <- args[3]

CK <- read.table(CKfile,sep = "\t",stringsAsFactors = F,header = T,row.names = 1)
treat <- read.table(treatfile,sep = "\t",stringsAsFactors = F,header = T,row.names = 1)
exprSet <- data.frame(CK=CK[[6]],treat=treat[[6]])
rownames(exprSet) <- rownames(CK)
DEG <- fedgeR(exprSet[,c("CK","treat")])
DEG$id <- rownames(DEG)
DEG <- DEG[,c("id","logCPM","logFC","PValue","change")]

write.table(DEG,outfile,sep = "\t",quote = F,row.names = F,col.names = T)
