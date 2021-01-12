## 7.AR-lncRNA ~ AR-mRNA networks
# Combine the samples into 4 stages(008w&026w,026w&060w,060w&078w,078w&104w)
for (i in 1:length(organ_list)){
  organ <- organ_list[i]
  organ_stage <- c("Stage1","Stage2","Stage3","Stage4")
  # AR-lncRNA
  aged <- c(paste(organ,"_078w",sep = ""),paste(organ,"_104w",sep = ""))
  ARlnc <- MultiTissue[[organ]][["ARlncRNA"]]
  ARlncExpr <- MultiTissue[[organ]][["log_norm"]][ARlnc,]
  ARlncStage <- vector(mode="list")
  for (j in 1:4){
    stage <- organ_stage[j]
    rawdata <- ARlncExpr[,(5*(j-1)+1):(5*(j+1))]
    ARlncStage[[stage]][["raw"]] = rawdata
    for (k in nrow(rawdata)){
      ARlncRNA <- unlist(rawdata[k,])
      z.ARlncRNA <- (ARlncRNA-mean(ARlncRNA))/sqrt(var(ARlncRNA))
      rawdata[k,] <- z.ARlncRNA
    }
    ARlncStage[[stage]][["Zscore-norm"]] = rawdata
  }
  MultiTissue[[organ]][["Dynamic Network"]][["ARlncStage"]] <- ARlncStage
  # AR-mRNA
  ARm <- MultiTissue[[organ]][["ARmRNA"]]
  ARmExpr <- MultiTissue[[organ]][["log_norm"]][ARm,]
  ARmStage <- vector(mode="list")
  for (j in 1:4){
    stage <- organ_stage[j]
    rawdata <- ARmExpr[,(5*(j-1)+1):(5*(j+1))]
    ARmStage[[stage]][["raw"]] = rawdata
    for (k in nrow(rawdata)){
      ARmRNA <- unlist(rawdata[k,])
      z.ARmRNA <- (ARmRNA-mean(ARmRNA))/sqrt(var(ARmRNA))
      rawdata[k,] <- z.ARmRNA
    }
    ARmStage[[stage]][["Zscore-norm"]] = rawdata
  }
  MultiTissue[[organ]][["Dynamic Network"]][["ARmStage"]] <- ARmStage
}
for (i in 1:length(organ_list)){
  organ <- organ_list[i]
  organ_stage <- c("Stage1","Stage2","Stage3","Stage4")
  lnclist <- MultiTissue[[organ]][["Dynamic Network"]][["ARlncStage"]]
  mlist <- MultiTissue[[organ]][["Dynamic Network"]][["ARmStage"]]
  corInfo <- vector(mode = "list")
  for (j in 1:4){
    stage <- organ_stage[j]
    m <- t(lnclist[[stage]]$`Zscore-norm`)
    lnc <- t(mlist[[stage]]$`Zscore-norm`)
    mlncCor <- cor(m,lnc,use="p",method="pearson")
    n <- nrow(m)+nrow(lnc)
    mlncp <- corPvalueStudent(mlncCor,n)
    padj <- p.adjust(mlncp,method = "fdr")
    dim(padj) <- dim(mlncp)
    corInfo[[stage]][["cor"]] <- mlncCor
    corInfo[[stage]][["p"]] <- mlncp
    corInfo[[stage]][["fdr"]] <- padj
  }
  # keep edges >0.9 in stage1 and stage4
  stage1 <- abs(corInfo$Stage1$cor)
  stage1p <- corInfo$Stage1$fdr
  stage4 <- abs(corInfo$Stage4$cor)
  stage4p <- corInfo$Stage4$fdr
  # edges keep in Stage1
  stage1_info <- matrix(NA,ncol = ncol(stage1),nrow = nrow(stage1))
  rownames(stage1_info) <- rownames(stage1)
  colnames(stage1_info) <- colnames(stage1)
  #stage1_info <- as.matrix(stage1_info)
  for (i in 1:nrow(stage1)){
    for (j in 1:ncol(stage1)){
      if ((stage1[i,j]>0.9)&(stage1p[i,j]<0.05)){
        stage1_info[i,j] <- 1
      }
      else{
        stage1_info[i,j] <- 0
      }
    }
  }
  # edges keep in Stage4
  stage4_info <- matrix(NA,ncol = ncol(stage1),nrow = nrow(stage1))
  rownames(stage4_info) <- rownames(stage1)
  colnames(stage4_info) <- colnames(stage1)
  #stage1_info <- as.matrix(stage4_info)
  for (i in 1:nrow(stage4)){
    for (j in 1:ncol(stage4)){
      if ((stage4[i,j]>0.9)&(stage4p[i,j]<0.05)){
        stage4_info[i,j] <- 1
      }
      else{
        stage4_info[i,j] <- 0
      }
    }
  }
  stage_info <- stage1_info + stage4_info
  keep <- stage_info[which(rowSums(stage_info)>0),which(colSums(stage_info)>0)]
  graphMat <- matrix(NA,nrow = nrow(keep)*ncol(keep),ncol = 3)
  colnames(graphMat) <- c("lnc","m","weight")
  graphMat <- as.data.frame(graphMat)
  klnc <- colnames(keep)
  km <- rownames(keep)
  graphMat$lnc <- rep(klnc,each=length(km))
  graphMat$m <- rep(km,length(klnc))
  write.csv(graphMat,"stage1.csv")
  f_stage4 <- abs(corInfo$Stage1$cor[rownames(keep),colnames(keep)])
  for (k in 1:ncol(f_stage4)){
    graphMat[(nrow(f_stage4)*(k-1)+1):(nrow(f_stage4)*k),3] <- f_stage4[,k]
  }
  # Please execute the following steps on your HPC
  g4 <- graph_from_data_frame(graphMat,directed = F)
  wc <- cluster_walktrap(g4,steps = 10)
  modularity(wc)
  memb <- data.frame(EnsembleID=names(membership(wc)),Modules=as.numeric(membership(wc)),stringsAsFactors = F)
}

