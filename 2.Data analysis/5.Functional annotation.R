## 5.Functional annotation of aging-regulated lncRNAs
# Pearson correlation>0.8 , FDR < 0.05
for (i in 1:length(organ_list)){
  organ <- organ_list[i]
  ARlncRNA <- MultiTissue[[organ]][["ARlncRNA"]]
  ARlncExpr <- MultiTissue[[organ]][["log_norm_mean"]][ARlncRNA,]
  Expr <- MultiTissue[[organ]][["log_norm_mean"]]
  mExpr <- Expr[intersect(mRNAname,rownames(Expr)),]
  corMatrix <- cor(t(ARlncExpr),t(mExpr),use = "p",method = "pearson")
  Pvalue <- corPvalueStudent(corMatrix,nrow(corMatrix)) 
  FDR <- p.adjust(Pvalue,method = 'fdr')
  dim(FDR) <- dim(Pvalue)
  rownames(FDR) <- rownames(Pvalue)
  colnames(FDR) <- colnames(Pvalue)
  MultiTissue[[organ]][["Functional Annotation"]][["total"]][["corMatrix"]] = corMatrix
  MultiTissue[[organ]][["Functional Annotation"]][["total"]][["Pvalue"]] = Pvalue
  MultiTissue[[organ]][["Functional Annotation"]][["total"]][["FDR"]] = FDR
}
library(clusterProfiler)
for (i in 1:length(organ_list)){
  organ <- organ_list[i]
  organ_corMatrix <- MultiTissue[[organ]][["Functional Annotation"]][["total"]][["corMatrix"]]
  organ_FDR <- MultiTissue[[organ]][["Functional Annotation"]][["total"]][["FDR"]] 
  organ_GO <- matrix(NA,ncol=6,nrow=0)
  colnames(organ_GO) <- c("AR-lncRNA","ID","Term","Genes","adj_pval","counts")
  for (j in 1:nrow(organ_corMatrix)){
    set1 <- colnames(organ_corMatrix[,which(abs(organ_corMatrix[j,])>0.8)])
    set2 <- colnames(organ_FDR[,which(organ_FDR[j,]<0.05)])
    set <- intersect(set1,set2)
    trans <- bitr(set, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
    list <- trans$ENTREZID 
    go <- enrichGO(list, OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'fdr',pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.1,keyType = 'ENTREZID')
    go1 <- data.frame(go)
    go1$geneID <- gsub("/",",",go1$geneID)
    go2 <- cbind(rownames(organ_corMatrix)[j],go1[,c(1,2,8,6,9)])
    organ_GO <- rbind(organ_GO,go2)
  }
  MultiTissue[[organ]][["Functional Annotation"]][["tissue"]] <- organ_GO
}
save(MultiTissue,file = "working_data/5.Functional annotation/01.Functional Annotation.R")