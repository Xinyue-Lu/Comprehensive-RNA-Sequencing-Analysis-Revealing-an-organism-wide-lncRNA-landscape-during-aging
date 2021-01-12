## 1.Load Raw data
data <- read.table("featureCounts.txt",header = T,row.names = 1)
rownames(data) <- substring(rownames(data),1,18)
save(data,file = "working_data/1.Raw data/01.Raw data.R")
# load
load("working_data/1.Raw data/01.Raw data.R")