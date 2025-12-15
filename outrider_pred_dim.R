# outrider_pred.R
# input:
# 1. ctsFile
# 2. q0

# output:
# 1. norm_cts
# 2. res

args = commandArgs(trailingOnly=TRUE)
print(args)
ctsFile <- args[1]
q0 <- as.numeric(args[2])
norm_cts <- args[3]
res <- args[4]
# 打印每个参数的赋值结果
cat("ctsFile:", ctsFile, "\n")
cat("q0:", q0, "\n")
cat("norm_cts:", norm_cts, "\n")
cat("res:", res, "\n")

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(dplyr)
})

ctsTable <- read.table(ctsFile, check.names = FALSE, header = TRUE,row.names = 1)
ctsTable[1:5,1:5]
ods <- OutriderDataSet(countData=ctsTable)
ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)
ods <- estimateSizeFactors(ods)
ods <- OUTRIDER(ods, q = q0, BPPARAM = SerialParam())
# ods <- OUTRIDER(ods, q = q0, BPPARAM = SerialParam(), iterations=4)

write.table(counts(ods, normalized=TRUE) %>% round(2), file = norm_cts, sep='\t',quote = FALSE)
write.table(results(ods,all = TRUE), file = res, sep='\t',quote = FALSE)