library(data.table)


v2.res <- fread("cg_v2table.txt")

typed.snps <- fread("typed.snps.txt")

setkey(typed.snps, gene, genome, disease, outcome)
setkey(v2.res, gene, genome, disease, outcome)

v2.typed.snp <- v2.res[typed.snps]

v2.typed.snp <- v2.typed.snp[,"adj.topSNP.pval":=topSNP.pval/no.typed.snps]

v2.typed.snps <- v2.typed.snp[,"adj.geneBased.pval":=geneBased.pval/166]
setcolorder(v2.typed.snp, c("chr", "gene", "nSNPs", "nSims","start", "stop", "Test", "geneBased.pval", "adj.geneBased.pval", "topSNP", "topSNP.pval", "adj.topSNP.pval", "genome", "cohort", "outcome", "disease", "no.typed.snps"))

write.table(v2.typed.snps, file="cg_v2table.txt", col.names=T, row.names=F, quote=F, sep="\t")