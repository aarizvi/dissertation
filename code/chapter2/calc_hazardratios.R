library(data.table)

files <- list.files(pattern=".txt")

res <- lapply(files, fread)

# 95% CI of hazard ratios
hr.ci <- function(coef.est, se){
	lb <- round(exp(coef.est-1.96*se), 5)
	ub <- round(exp(coef.est+1.96*se), 5)
	paste0("[", lb,", " , ub ,"]")
}

# make new column of HR in lists
res <- lapply(res, function(x) x <- x[,'95%-CI':= hr.ci(x$coef, x$se.coef)])


# change the imputation notation to 'imputed' or 'typed'
impute <- function(x){
        x[impute != "---", impute := "typed"]
        x[impute == "---", impute := "imputed"]
}
res <- lapply(res, impute)

res <- lapply(res, function(x) setkey(x, rsID))

# lets make alleles all upper case too
res <- lapply(res, function(x) x <- x[,'allele1':=toupper(allele1),])
res <- lapply(res, function(x) x <- x[,'allele2':=toupper(allele2),])

names(res) <- files

res.files <- paste0("/projects/rpci/lsuchest/abbasriz/candidate_gene/result_files_cg/", files)
# sort by gene
res <- lapply(res, function(x) setkey(x, gene))
# write to file
for(i in 1:length(res.files)){
	write.table(res[[i]], file=res.files[[i]], sep="\t", row.names=F, col.names=T, quote=F)
}