## 3.Normalization
library(DESeq2)
# Raw counts and processed data for different tissues
raw <- as.matrix(data[,-1])
organ_list <- c("Lv","Ki","He","Ba","iW","eW","Gm","Br","Hp","Lu","Bm")
MultiTissue = vector(mode="list")
for(i in 1:length(organ_list)){
  organ <- organ_list[i]
  # Input raw data 
  organ_raw <- raw[,substring(colnames(raw),1,2)==organ]
  MultiTissue[[organ]]$raw <- organ_raw
  ColData <- matrix(NA,ncol = 1,nrow = ncol(organ_raw))
  colnames(ColData) <- "Age"
  rownames(ColData) <- colnames(organ_raw)
  ColData[,"Age"] <- rep(1:5,each=5)
  # Create a DESeqDataSet object
  dds_norm <- DESeqDataSetFromMatrix(organ_raw, ColData, design = ~ Age)
  # Normalized and filtered expression data
  norm <- estimateSizeFactors(dds_norm)
  norm <- counts(norm,normalized=TRUE)
  log_norm <- log1p(norm)
  filt <- log_norm[rowSums(log_norm>0)>=5,]
  raw_filt <- organ_raw[rownames(filt),]
  organ_mean <- matrix(NA, nrow = nrow(filt), 5)
  rownames(organ_mean) <- rownames(filt)
  colnames(organ_mean) <- unique(substring(colnames(filt), 1, 7))
  organ_mean <- as.data.frame(organ_mean)
  for(j in 1:5){
    organ_mean[,j] = apply(filt[,(5*(j-1)+1):(5*j)],1,mean)
  }
  # Store Processed data
  MultiTissue[[organ]]$log_norm <- log_norm
  MultiTissue[[organ]]$raw_filter <- raw_filt
  MultiTissue[[organ]]$log_norm_filter <- filt
  MultiTissue[[organ]]$log_norm_mean <- organ_mean
}
save(MultiTissue,organ_list,file = "working_data/3.Normalization/01.Data for different tissues.R")
# load
load("working_data/3.Normalization/01.Data for different tissues.R")