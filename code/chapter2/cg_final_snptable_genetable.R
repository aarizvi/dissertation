library(data.table)

# going to merge 


# subsetted impute results w/ individual cohort and meta beta estimates/hr/CIs merged
path.cg <- "/projects/rpci/lsuchest/abbasriz/candidate_gene/result_files_cg/"
files.cg <- list.files(path=path.cg, pattern=".txt")[-1]
impute.res <- lapply(paste0(path.cg,files.cg), fread)

# add the columns that we want, a column for genome, outcome, and cohort
# grab that information from the file names
genome <- gsub("_.*$", "", files.cg)
disease <- sub("(.*\\_)([^.]+)(\\.[[:alnum:]]+$)", "\\2", files.cg)
cohort <- sapply(strsplit(files.cg, '_', fixed=T), "[", 2)
outcome <- sapply(strsplit(files.cg, '_', fixed=T), "[", 3)

# creating columns with that info
mapply(function(x, cat) x <- x[,'genome' := cat], impute.res, genome)
mapply(function(x, cat) x <- x[,'cohort' := cat], impute.res, cohort)
mapply(function(x, cat) x <- x[,'outcome' := cat], impute.res, outcome)
mapply(function(x, cat) x <- x[,'disease' := cat], impute.res, disease)


impute.res <- lapply(impute.res, setkey, gene)



impute.res <- do.call(rbind, impute.res)

impute.res$gene[which(impute.res$gene == "GSTM1")] <- "GSTM"

impute.res <- na.omit(impute.res)

#write.table(impute.res, file="/projects/rpci/lsuchest/abbasriz/candidate_gene/result_files_cg/cg_snptable.txt", col.names=T, row.names=F, quote=F, sep="\t")




# candidate gene list
gene.list <- unique(impute.res$gene)

gene.list <- scan("gene_list.txt", what=character())[-1]

gene.list <- c(gene.list, "GSTM1", "GSTM2", "GSTM3", "GSTM4", "GSTM5")
gene.list[gene.list=="GSTM"] <- NA

gene.list <- gene.list[complete.cases(gene.list)]

# vegas2 results
path.v2 <- "/projects/rpci/lsuchest/abbasriz/candidate_gene/result_files_cg/vegas2/results/"
files.v2 <- list.files(path=path.v2, pattern=".V2out")
v2.res <- lapply(paste0(path.v2, files.v2), fread)

v2.res <- lapply(v2.res, setkey, Gene)
v2.res <- lapply(v2.res, function(x) x[gene.list])

#v2.res <- lapply(v2.res, function(x) x[,c("Chr", "nSNPs", "nSims", "Start", "Stop", "Test"):=NULL])
v2.res <- lapply(v2.res, setNames, c("chr", "gene", "nSNPs", "nSims", "start", "stop", "Test", "geneBased.pval", "topSNP", "topSNP.pval"))

v2.res <- lapply(v2.res, setkey, gene)


genome <- sapply(strsplit(files.v2, '_', fixed=T), "[", 2)
outcome <- sapply(strsplit(files.v2, '_', fixed=T), "[", 4)
cohort <- sapply(strsplit(files.v2, '_', fixed=T), "[", 3)
disease <- sapply(strsplit(files.v2, '_', fixed=T), "[",5)
disease <- sapply(strsplit(disease, '.', fixed=T), "[",1)

# creating columns with that info
#lapply(v2.res, function(x, cat) x <- x[,'genome':=cat], genome)

#mapply(function(x, cat) x <- x[,'genome' := cat], v2.res, genome)
#mapply(function(x, cat) x <- x[,'cohort' := cat], v2.res, cohort)
#mapply(function(x, cat) x <- x[,'outcome' := cat], v2.res, outcome)
#mapply(function(x, cat) x <- x[,'disease' := cat], v2.res, disease)


for(i in 1:length(v2.res)){
	v2.res[[i]][,'genome':=genome[i]]
	v2.res[[i]][,'cohort':=cohort[i]]
	v2.res[[i]][,'outcome':=outcome[i]]
	v2.res[[i]][,'disease':=disease[i]]
}


v2.res <- do.call(rbind, v2.res)

write.table(v2.res, file="/projects/rpci/lsuchest/abbasriz/candidate_gene/result_files_cg/cg_v2_snptable.txt", col.names=T, row.names=F, quote=F, sep="\t")

## normalize top snp pvalue by number of typed snps in cohort 1
p <- c()
x <- c()
rr <- unique(impute.res[,c("gene", "genome", "disease", "outcome"),with=F])

# my code
for(i in seq(nrow(rr))){
                x <- impute.res[gene==rr$gene[i] & disease==rr$disease[i] & genome==rr$genome[i] & outcome==rr$outcome[i] & cohort=="c1"]                
                setkey(x, "impute")
                p[i] <- nrow(x[impute=="typed"])
}

rr <- rr[,no.typed.snps:=p]

write.table(file="typed_snps_table.txt", rr, quote=F, row.names=F, col.names=T, sep="\t")



