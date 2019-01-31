library(data.table)

donor.files <- list.files(pattern="^D_")

for(i in 1:length(donor.files)){
	donor <- fread(donor.files[i], header="auto", sep="\t")
	setkey(donor, chr, BP, rsID)
	donor <- donor[!which(duplicated(donor$rsID))-1]
	write.table(donor, file=paste0("/projects/rpci/lsuchest/abbasriz/candidate_gene/result_files_cg/", donor.files[i]), quote=F, sep="\t", col.names=T, row.names=F)
	rm(donor)
}
