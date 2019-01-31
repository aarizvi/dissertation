# add GSTMs
# we are going to use the same number of typed SNPs from all GSTM
# which is 25 typed SNPs
library(data.table)
typed.snps <- fread("typed.snps.txt")

setkey(typed.snps, gene)

gstm <- typed.snps[gene=="GSTM"]
gstm$gene <- "GSTM1"
gstm1 <- gstm
gstm$gene <- "GSTM2"
gstm2 <- gstm
gstm$gene <- "GSTM3"
gstm3 <- gstm
gstm$gene <- "GSTM4"
gstm4 <- gstm
gstm$gene <- "GSMT5"
gstm5 <- gstm

# remove GSTM from data.table
typed.snps <- typed.snps[!"GSTM"]

#typed.snps <- 

typed.snps <- data.table(rbind(typed.snps, gstm1, gstm2, gstm3, gstm4, gstm5))

# sort by gene name
setkey(typed.snps, gene)

# see if number is right 
# 174 * 48 = 8352
nrow(typed.snps)

write.table(typed.snps, file="typed.snps.txt", row.names=F, col.names=T, quote=F, sep="\t")


library(data.table)


v2.res <- fread("cg_v2table.txt")

typed.snps <- fread("typed.snps.txt")

setkey(typed.snps, gene, genome, disease, outcome)
setkey(v2.res, gene, genome, disease, outcome)

v2.typed.snp <- v2.res[typed.snps]
# adjust for number of typed snps
v2.typed.snp <- v2.typed.snp[,"adj.topSNP.pval":=topSNP.pval*no.typed.snps]
# adjust for number of genes
v2.typed.snps <- v2.typed.snp[,"adj.geneBased.pval":=geneBased.pval*174]
setcolorder(v2.typed.snp, c("chr", "gene", "nSNPs", "nSims","start", "stop", "Test", "geneBased.pval", "adj.geneBased.pval", "topSNP", "topSNP.pval", "adj.topSNP.pval", "genome", "cohort", "outcome", "disease", "no.typed.snps"))

write.table(v2.typed.snps, file="cg_v2table_adj.txt", col.names=T, row.names=F, quote=F, sep="\t")