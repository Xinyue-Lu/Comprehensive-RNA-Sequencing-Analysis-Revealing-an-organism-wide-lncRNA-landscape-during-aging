## 6.Tissue-specificity of lncRNAs
time_list <- c("008w","026w","060w","078w","104w")
MultiAge <- vector(mode="list")
# Compare at different age point
for (i in 1:length(time_list)){
  time <- time_list[i]
  exp <- raw[,substring(colnames(raw),4,7)==time]
  exp <- exp[which(rowSums(exp)>0),]
  meanExpr <- matrix(NA,nrow = nrow(exp),ncol=11)
  rownames(meanExpr) <- rownames(exp)
  colnames(meanExpr) <- organ_list
  meanExpr <- as.data.frame(meanExpr)
  for (j in 1:length(organ_list)){
    organ <- organ_list[j]
    meanExpr[,j] <- apply(exp[,substring(colnames(exp),1,2)==organ],1,mean)
  }
  MultiAge[[time]][["MeanExpr"]] <- meanExpr
  TF <- matrix(NA,nrow = nrow(meanExpr),ncol = ncol(meanExpr))
  rownames(TF) <- rownames(meanExpr)
  colnames(TF) <- colnames(meanExpr)
  TF <- as.data.frame(TF)
  for (m in 1:nrow(TF)){
    sum <- rowSums(meanExpr[m,])
    for (n in 1:ncol(TF)){
      tf <- meanExpr[m,n]/sum
      TF[m,n] <- tf
    }
  }
  MultiAge[[time]][["Tissue specific score analysis"]][["Tissue Fraction"]] <- TF
}
save(time_list,MultiAge,file = "working_data/6.Tissue specificity analysis/01.Tissue Fraction.R")
for (i in 1:length(time_list)){
  time <- time_list[i]
  TF <- MultiAge[[time]][["Tissue specific score analysis"]][["Tissue Fraction"]]
  TS <- as.data.frame(apply(TF,1,max))
  colnames(TS) <- "Tissue specific score"
  TS_lncRNA <- as.data.frame(TS[rownames(gtf_lncRNA),])
  rownames(TS_lncRNA) <- rownames(gtf_lncRNA)
  colnames(TS_lncRNA) <- "Tissue specific score"
  TS_lncRNA <- na.omit(TS_lncRNA)
  MultiAge[[time]][["Tissue specific score analysis"]][["Tissue specific score"]] <- TS_lncRNA
}
save(MultiAge,file = "working_data/6.Tissue specificity analysis/02.Tissue specific score.R")