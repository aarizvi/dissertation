## INDEPENDENT RESULTS
library(data.table)

# read in candidate gene
info <- fread("/projects/rpci/lsuchest/abbasriz/candidate_gene/info.cg.txt", header='auto', sep="\t")
#colnames(info) <- c("chr", "position", "gene", "snp_id", "rs_id", "exp_freq_a1", "info", "certainty", "type", "info_type0", "concord_type0", "r2_type0")
setkey(info, "gene","snp_id", "rs_id")
info <- info[,c("gene", "snp_id", "rs_id")]
setkey(info, "snp_id","rs_id")


ind <- scan("/projects/rpci/lsuchest/abbasriz/candidate_gene/result_files_cg/res.directories/ind.directories.txt", what=character())

for(i in 1:length(ind)){
	ind.res <- fread(ind[i], header="auto", sep="\t", verbose=T)
	colnames(ind.res)[1:2] <- c("snp_id", "rs_id")
	setkey(ind.res, "snp_id", "rs_id")
	ind.res <- ind.res[info]
	ind.res <- na.omit(ind.res)
	ind.res[,c("z", "loglik0", "loglik"):=NULL]
	ind.res[,`95%-CI`:=NA]
	setcolorder(ind.res, c("gene", "rs_id", "CHR", "BP", "ALLELE1", "ALLELE2", "n", "coef", "se(coef)", "exp(coef)", '95%-CI', "snp_id", "Pr(>|z|)"))
	colnames(ind.res) <- c("gene", "rsID", "chr", "BP", "allele1", "allele2", "N", "coef", "se.coef",  "exp.coef", "95%-CI", "impute", "Pvalue")
	file.name <- strsplit(ind[i], split="/")[[1]][14]
	file.name <- gsub("[.]", "_", file.name)
	file.name <- strsplit(file.name, "_")[[1]]
	if(length(file.name)<6){
		file.name[6] <- "mixed"
	}	else{
			file.name <- file.name
	}
	write.table(ind.res, file=paste0(paste(file.name[1],file.name[2],file.name[4],file.name[6],sep="_"),".txt"), quote=F, sep="\t", col.names=T, row.names=F)
	rm(ind.res)
}

