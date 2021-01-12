## 2.Reannotate Probes
# Load gtf data
gtf <- rtracklayer::import('gencode.vM17.annotation.gtf')
gtf <- as.data.frame(gtf)
save(gtf,file = "working_data/2.Reannotate Probes/01.Load gtf data.R")
# load
load("working_data/2.Reannotate Probes/01.Load gtf data.R")
# Calculate transcript numbers
gtf$gene_id <- substring(gtf$gene_id,1,18)
gtf$transcript_id <- substring(gtf$transcript_id,1,18)
gene_id <- unique(gtf$gene_id)
transcript_number <- matrix(NA,nrow = length(gene_id),ncol = 1)
for(i in 1:length(gene_id)){
  gene <- gene_id[i]
  transcript_id <- na.omit(gtf[gtf$gene_id==gene,]$transcript_id)
  transcript_num <- length(unique(transcript_id))
  transcript_number[i] <- transcript_num
}
save(transcript_number,file = "working_data/2.Reannotate Probes/02.Calculate transcript numbers.R")
# load
load("working_data/2.Reannotate Probes/02.Calculate transcript numbers.R")
# Integrate gene information
gtf_info <- gtf[gtf$type=="gene",c("gene_id","gene_type","gene_name")]
gtf_info <- cbind(gtf_info,transcript_number=transcript_number)
rownames(gtf_info) <- gtf_info$gene_id
gtf_data <- gtf_info[rownames(data),]
gtf_data <- cbind(gtf_data,gene_length = data$Length)
gtf_mRNA <- gtf_data[gtf_data$gene_type=="protein_coding",]
gtf_lncRNA <- gtf_data[!gtf_data$gene_type=="protein_coding",]
mRNAname <- rownames(gtf_mRNA)
lncRNAname <- rownames(gtf_lncRNA)
save(gtf_data,gtf_mRNA,gtf_lncRNA,file = "working_data/2.Reannotate Probes/03.Integrate gene information.R")
# load
load("working_data/2.Reannotate Probes/03.Integrate gene information.R")
