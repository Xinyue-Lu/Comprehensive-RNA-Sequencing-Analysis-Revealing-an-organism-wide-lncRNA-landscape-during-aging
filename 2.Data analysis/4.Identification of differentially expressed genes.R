## 4.Identification of differentially expressed genes
# Pairwise DEG analysis
for (i in 1:length(organ_list)){
  organ <- organ_list[i]
  exp <- MultiTissue[[organ]]$raw
  # comparing the aged groups(26w,60w,78w,104w) to the younger group(8w)
  young <- paste(organ,"_008w",sep = "")
  aged <- c(paste(organ,"_026w",sep = ""),paste(organ,"_060w",sep = ""),paste(organ,"_078w",sep = ""),paste(organ,"_104w",sep = ""))
  for (j in 1:length(aged)){
    week <- aged[j]
    selected <- exp[,(substring(colnames(exp),1,7)==young)|(substring(colnames(exp),1,7)==week)]
    # ColData
    ColData <- matrix(NA,ncol = 1,nrow = ncol(selected))
    colnames(ColData) <- "Age"
    rownames(ColData) <- colnames(selected)
    ColData[,"Age"] <- rep(c("young","aged"),each = 5)
    # DEG analysis
    dds <- DESeqDataSetFromMatrix(selected, ColData, design = ~ Age)
    DEG <- DESeq(dds)
    DEG_res <- results(DEG)
    table <- cbind(DEG_res$log2FoldChange,DEG_res$pvalue,DEG_res$padj)
    rownames(table) <- rownames(selected)
    colnames(table) <- c("log2FC","p","padj")
    mRNA_total <- table[mRNAname,]
    lncRNA_total <- table[lncRNAname,]
    MultiTissue[[organ]][["Differential Expression"]][[week]][["mRNA"]][["total"]] <- mRNA_total
    MultiTissue[[organ]][["Differential Expression"]][[week]][["lncRNA"]][["total"]] <- lncRNA_total
    # tissue-specific
    fexp <- MultiTissue[[organ]]$raw_filter
    tissue_mRNA <- intersect(rownames(gtf_mRNA),rownames(fexp))
    tissue_lncRNA <- intersect(rownames(gtf_lncRNA),rownames(fexp))
    mRNA_tissue <- as.data.frame(mRNA_total[tissue_mRNA,])
    lncRNA_tissue <- as.data.frame(lncRNA_total[tissue_lncRNA,])
    mRNA_tissue <- na.omit(mRNA_tissue)
    lncRNA_tissue <- na.omit(lncRNA_tissue)
    MultiTissue[[organ]][["Differential Expression"]][[week]][["mRNA"]][["tissue"]] <- mRNA_tissue
    MultiTissue[[organ]][["Differential Expression"]][[week]][["lncRNA"]][["tissue"]] <- lncRNA_tissue
    # log2FC > 0.8 & adjust p value < 0.05
    mRNA_selected = mRNA_tissue[(abs(mRNA_tissue$log2FC)>1)&(mRNA_tissue$padj<0.05),]
    lncRNA_selected = lncRNA_tissue[(abs(lncRNA_tissue$log2FC)>1)&(lncRNA_tissue$padj<0.05),]
    MultiTissue[[organ]][["Differential Expression"]][[week]][["mRNA"]][["selected"]] = mRNA_selected
    MultiTissue[[organ]][["Differential Expression"]][[week]][["lncRNA"]][["selected"]] = lncRNA_selected
  }
}
save(MultiTissue,file = "working_data/4.DEG analysis/01.Differential Expression.R")
# load
load("working_data/4.DEG analysis/01.Differential Expression.R")

# Identify AR-mRNAs/AR-lncRNAs 
for (i in 1:length(organ_list)){
  organ <- organ_list[i]
  aged <- c(paste(organ,"_078w",sep = ""),paste(organ,"_104w",sep = ""))
  ARmRNA <- union(rownames(MultiTissue[[organ]][["Differential Expression"]][[aged[1]]][["mRNA"]][["selected"]]),
                  rownames(MultiTissue[[organ]][["Differential Expression"]][[aged[2]]][["mRNA"]][["selected"]]))
  ARlncRNA <- union(rownames(MultiTissue[[organ]][["Differential Expression"]][[aged[1]]][["lncRNA"]][["selected"]]),
                    rownames(MultiTissue[[organ]][["Differential Expression"]][[aged[2]]][["lncRNA"]][["selected"]]))
  MultiTissue[[organ]][["ARmRNA"]] <- ARmRNA
  MultiTissue[[organ]][["ARlncRNA"]] <- ARlncRNA
}
save(MultiTissue,file = "working_data/4.DEG analysis/02.AR-genes.R")
load("working_data/4.DEG analysis/02.AR-genes.R")