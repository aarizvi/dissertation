## CCR5 haplotypes

# chr3 46411633 46417697 6065 + CCR5

#/projects/rpci/lsuchest/lsuchest/Rserve/ImputeData/var/db/gwas/imputed_data/BMT093013_forImpute/BMT093013_forImpute.chr3-45000000-50000000.impute2

#qctool -g /projects/rpci/lsuchest/lsuchest/Rserve/ImputeData/var/db/gwas/imputed_data/BMT093013_forImpute/BMT093013_forImpute.chr3-45000000-50000000.impute2 -s /projects/rpci/lsuchest/lsuchest/Rserve/ImputeData/var/db/gwas/imputed_data/BMT093013_forImpute/BMT093013_forImpute.chr16-50000000-55000000.impute2_samples -og /projects/rpci/lsuchest/abbasriz/candidate_gene/cg_haplotypes/ccr5_rep/test/ccr5_rep_dosages_threshold.vcf -incl-rsids /projects/rpci/lsuchest/abbasriz/candidate_gene/cg_haplotypes/ccr5_rep/ccr5_snps.txt 



dosages <- read.table("ccr5_dosages.impute2", sep="\t")


library(VariantAnnotation)	
library(data.table)
library(survival)
library(dplyr)
## read in vcf file
vcf <- readVcf("/projects/rpci/lsuchest/abbasriz/candidate_gene/cg_haplotypes/ccr5_rep/ccr5_rep_dosages.vcf", genome="hg19")
vcf

## threshold genotypes
# grab genotypes of 3 SNPs for all patients
gt <- geno(vcf)$GT
gt[1:nrow(gt),1:5]

## no rs333 in our data
## proxy which rs113341849

# transpose the data.frame so we can patients as rows
gt <- data.frame(t(gt))
gt[1:5, 1:ncol(gt)]



## CCR5 rs113341849
table(gt$rs113341849)


h1h1 <- gt

# recode to homozygous major allele
h1h1 <- apply(h1h1, 2, FUN=function(x) gsub("[.]", NA, x))
h1h1 <- apply(h1h1, 2, FUN=function(x) gsub("0/0", 1, x))
h1h1 <- apply(h1h1, 2, FUN=function(x) gsub("0/1", 0, x))
h1h1 <- apply(h1h1, 2, FUN=function(x) gsub("1/1", 0, x))

library(dplyr)
library(data.table)
dosages <- fread("ccr5_dosages.impute2")
dosages <- data.table(t(dosages))
colnames(dosages) <- as.character(dosages[3,])
dosages <- dosages[-c(1:6),]
dosages <- data.frame(apply(dosages, 2, as.numeric))
dosages$ccr5 <- rowSums(dosages)
dosages <- dosages %>%
	mutate(h1h1=ifelse(ccr5<0.5,1,0)) %>%
	data.table(keep.rownames=T)


#h1h1 <- data.frame(apply(h1h1, 2, function(x) as.numeric(as.character(x))), check.names=F, row.names=row.names(h1h1))

#h1h1 <- data.table(h1h1, keep.rownames=T)

#head(h1h1)


## h1/h1 

###############################################
###### LOAD UP PATIENT ID AND COVARIATE FILES
###############################################


## grab unique identifiers of patients in different cohorts
library(data.table)
patients <- c("/projects/rpci/lsuchest/lsuchest/Rserve/BMT/genetic_data/R_c1_EA_FID_IID.txt",
              "/projects/rpci/lsuchest/lsuchest/Rserve/BMT/genetic_data/R_c2_EA_FID_IID.txt",
              "/projects/rpci/lsuchest/lsuchest/Rserve/BMT/genetic_data/D_c1_EA_FID_IID.txt",
              "/projects/rpci/lsuchest/lsuchest/Rserve/BMT/genetic_data/D_c2_EA_FID_IID.txt")

cohorts <- lapply(patients, fread)
files <- c()
for(i in 1:length(patients)) files[i] <- strsplit(patients, "/")[[i]][9]      
names(cohorts) <- files
head(cohorts)

rec.ids <- list(cohorts[[1]], cohorts[[2]])
names(rec.ids) <- files[1:2]
rec.ids[[1]]$cohort <- "c1"
rec.ids[[2]]$cohort <- "c2"


donor.ids <- list(cohorts[[3]], cohorts[[4]])
names(donor.ids) <- files[3:4]
donor.ids[[1]]$cohort <- "c1"
donor.ids[[2]]$cohort <- "c2"

## read in phenotype files files
r.pheno <- fread("/projects/rpci/lsuchest/lsuchest/Rserve/TheData/Plink_Recipient.pheno")
d.pheno <- fread("/projects/rpci/lsuchest/lsuchest/Rserve/TheData/Plink_Donor.pheno")
r.cov <- fread("/projects/rpci/lsuchest/lsuchest/Rserve/TheData/Plink_Recipient.cov")


# parse covariant/pheno files to just EA

## SUBSET EUROPEAN AMERICAN PATIENTS BASED OFF OF IID
id.subset <- function(x, pheno){pheno[pheno$IID %in% x$IID]}
# recipients
rec.pheno.ea <- lapply(rec.ids, id.subset, r.pheno)
# donors
don.pheno.ea <- lapply(donor.ids, id.subset, d.pheno)


# bring down to just FID, IID, pair_id
don.pheno.ea <- lapply(don.pheno.ea, function(x) x[,c("FID", "IID", "pair_id", "population"),with=F])
don.pheno.ea <- lapply(don.pheno.ea, setnames, c("D_FID", "D_IID", "pair_id", "population"))

# merge to recipient file

merge.pair_id <- function(recipient, donor){
        setkey(recipient, pair_id)
        setkey(donor, pair_id)
        recipient[donor]
}

merged.ea <- mapply(merge.pair_id, rec.pheno.ea, don.pheno.ea, SIMPLIFY=FALSE)

merged.ea <- data.table(do.call(rbind, merged.ea))

setnames(merged.ea, c("R_FID", "R_IID", colnames(merged.ea)[3:ncol(merged.ea)]))
# I think there are duplicated NAs because of the almost 2806 match and now 2815
merged.ea <- merged.ea[!duplicated(merged.ea$R_IID)]

ids <- merged.ea[,c("pair_id", "R_IID", "D_IID"), with=F]

# impute sample names

info.samples <- fread("/projects/rpci/lsuchest/lsuchest/Rserve/ImputeData/var/db/gwas/imputed_data/BMT093013_forImpute/BMT093013_forImpute.chr16-50000000-55000000.impute2_samples")

# top row is useless
info.samples <- info.samples[-1]

nrow(info.samples)
# 6805 samples
nrow(gt)
# nice! same number of samples...6805

# relabel the IDs so they just correspond to the sample order
# because when we converted the impute2 file to vcf, the order remains the same as seen info the impute2_samples file
# so we can just relabel the IDs as 1:6805 samples and go from there
info.samples$ID_1 <- seq(1, nrow(info.samples))


# okay now we will grab the indices from the impute_sample file and append them as columns to the ids data.table
ids$r_sample_index <- info.samples$ID_1[match(ids$R_IID, info.samples$ID_2)]
ids$d_sample_index <- info.samples$ID_1[match(ids$D_IID, info.samples$ID_2)]

# okay so now 22 donors didnt map back to that file
missing_samples <- ids[is.na(ids$d_sample_index)]
missing_samples
## weird .... but we will just omit them and move on because this is taking wayyyy too long and we'll ask lara about those samples later

ids <- na.omit(ids)
head(ids)
# now we have 2783 samples that we can go and grab all the genotype info on



##################################################
#### NOW PARSE GENOTYPE FILE BY SAMPLE ID INDICES
###################################################
## parse the unique IDs (genome/cohort) from the main genotype df into 4 dfs specific to genome/cohort
recs <- data.table(dosages[ids$r_sample_index,], keep.rownames=T)
donors <- data.table(dosages[ids$d_sample_index,], keep.rownames=T)

# remove "sample" from sample columns so we just have indices that we can map back to phenotype file
rm.sample <- function(x){ x$rn <- gsub("sample_", "", x$rn)
	class(x$rn) <- "numeric"
	x}

recs <- rm.sample(recs)
donors <- rm.sample(donors)

# match column names to ids column names
colnames(recs)[1] <- "r_sample_index"
colnames(donors)[1] <- "d_sample_index"

# grab pair ids so we can compare donor and recipient pairs
recs$pair_id <- ids$pair_id
donors$pair_id <- ids$pair_id


recs <- recs %>% na.omit()
donors <- donors %>% na.omit()

#recs$Rh1h1 <- recs %>%
#	select(grep("rs", colnames(recs))) %>%
#	mutate(Rh1h1=ifelse(rowSums(.)==4, 1, 0)) %>%
#	select(Rh1h1) %>%
#	na.omit()  

#donors$Dh1h1 <- donors %>%
#	select(grep("rs", colnames(donors))) %>%
#	mutate(Dh1h1=ifelse(rowSums(.)==4, 1, 0)) %>%
#	select(Dh1h1) %>%
#	na.omit() 

recs <- recs %>% rename(Rh1h1=h1h1)
donors <- donors %>% rename(Dh1h1=h1h1)

donors.h1h1 <- donors %>%
	select(pair_id, Dh1h1)

recs.h1h1 <- recs %>%
	select(pair_id, Rh1h1)

dr.h1h1 <- donors.h1h1 %>%
	right_join(recs.h1h1, "pair_id") %>%
	data.table()


dr.h1h1$grp2vsgrp1 <- with(dr.h1h1,
	ifelse(Dh1h1==1 & Rh1h1==0, 1, ifelse(Dh1h1==0 & Rh1h1==0, 0, NA)))
dr.h1h1$grp3vsgrp1 <- with(dr.h1h1,
	ifelse(Dh1h1==0 & Rh1h1==1, 1, ifelse(Dh1h1==0 & Rh1h1==0, 0, NA)))
dr.h1h1$grp3vsgrp2 <- with(dr.h1h1,
	ifelse(Dh1h1==0 & Rh1h1==1, 1, ifelse(Dh1h1==1 & Rh1h1==0, 0, NA)))



r.cov <- r.cov %>% right_join(dr.h1h1, "pair_id") %>% data.table()


## split into cohorts

r.cov.c1 <- r.cov[cohort1==1]

r.cov.c2 <- r.cov[cohort1==0]


#compsnps <- c("R_h1h1", "grp2vsgrp1", "grp3vsgrp1", "grp3vsgrp2")
merged.ea <- r.cov %>%
	dplyr::select(pair_id, age, distatD, PBlood, bmiOBS, bmiOVWT, MDSdummy, AMLdummy, ALLdummy, Dh1h1, Rh1h1, grp2vsgrp1, grp3vsgrp1, grp3vsgrp2) %>%
	inner_join(merged.ea, by="pair_id") %>%
	data.table()

merged.ea$population <- gsub("EA.", "", merged.ea$population)


merged.ea <- merged.ea %>% rename(cohort=population)



OScov1Y <- c("intxsurv_1Y","dead_1Y","age","distatD","PBlood")
PFScov1Y<- c("intxrel_1Y","lfs_1Y","age","distatD")
OScov.full <- c("intxsurv_1Y","dead_1Y","age","distatD","PBlood")
PFScov.full<- c("intxrel_1Y","lfs_1Y","age","distatD")




###############################################
######  SURVIVAL ANALYSIS
###############################################

########################
######## EVENTS ########
########################
## DD - disease_death_1Y
## OS  - dead_1Y
## PFS - lfs_1Y
## TRM - TRM_1Y
########################

####################################################
################## COVARIATES ######################
####################################################
## DD covariates are: age, distatD
## OS covariates are: age, distatD, Pblood
## PFS covariates: age, distatD
## TRM covariates: age, bmiOBS, bmiOVWT, PBlood
###################################################

DDcov <- c("intxsurv_1Y","disease_death_1Y","age","distatD")
OScov <- c("intxsurv_1Y","dead_1Y","age","distatD","PBlood")
PFScov<- c("intxrel_1Y","lfs_1Y","age","distatD")
TRMcov <- c("intxsurv_1Y","TRM_1Y","age","bmiOBS","bmiOVWT","PBlood")




# adjusts for 2 covariates + genotype of interest
surv_fit_cov2 <- function(geno, cov, covFile){
	outcome <- Surv(time=covFile[,cov[1]], event=covFile[,cov[2]])
	res <- coxph(outcome~covFile[,geno]+covFile[,cov[3]]+covFile[,cov[4]])
	summary(res)$coefficients[1,c(1,3,2,5)]
}

# adjusts for 3 covariates + genotype of interest
surv_fit_cov3 <- function(geno, cov, covFile){
	outcome <- Surv(time=covFile[,cov[1]], event=covFile[,cov[2]])
	res <- coxph(outcome~covFile[,geno]+covFile[,cov[3]]+covFile[,cov[4]]+covFile[,cov[5]])
	summary(res)$coefficients[1,c(1,3,2,5)]
}

# adjusts for 4 covariates + genotype of interest
surv_fit_cov4 <- function(geno, cov, covFile){
	outcome <- Surv(time=covFile[,cov[1]], event=covFile[,cov[2]])
	res <- coxph(outcome~covFile[,geno]+covFile[,cov[3]]+covFile[,cov[4]]+covFile[,cov[5]]+covFile[,cov[6]])
	summary(res)$coefficients[1,c(1,3,2,5)]
}

## WRITE TO TABLE
## NEED TO SPLIT INTO TWO TABLES FOR META-ANALYSIS
# will calculate confidence interval of hazard ratios AFTER
# set up columns so they are the same as our candidate gene table
cols <- c("gene",
          "rsID",
          "chr",
          "BP",
          "N",
          "allele1",
          "allele2",
          "coef",
          "se.coef",
          "exp.coef",
          "95%-CI",
          "Pvalue",
          "impute",
          "genome",
          "cohort",
          "outcome",
          "disease")

make.tbl <- function(outcome.vector, gene, rsID, chr, BP, N, allele1, allele2, imputed, genome, cohort, outcome, disease){
        coef <- outcome.vector[1]
        se.coef <- outcome.vector[2]
        hr <- outcome.vector[3]
        pval <- outcome.vector[4]
        res <- c(gene, rsID, chr, BP, N, allele1, allele2, coef, se.coef, hr, NA, pval, imputed, genome, cohort, outcome, disease)
        res
}





###############################
#### MASTER FILE CREATION #####
###############################


compsnps <- c("R_h1h1", "grp2vsgrp1", "grp3vsgrp1", "grp3vsgrp2")
cov.list <- c("DDcov", "PFScov", "OScov", "TRMcov")
outcomes <- c("DD", "PFS", "OS", "TRM")
cohorts <- c("c1", "c2")
diseases <- c("mixed", "AMLonly", "ALLonly", "noMDS", "noALL")
survival.functions <- c("surv_fit_cov2", "surv_fit_cov2", "surv_fit_cov3", "surv_fit_cov4")



create.master <- function(genotype, genome, outcomes, cohorts, diseases, survival.functions){
	master <- data.frame(matrix(nrow=3, ncol=2))
	#master <- list()
	for(i in 1:length(diseases)){
		master[i,] <- c(genotype, diseases[i])

	}
	colnames(master) <- c("genotype", "disease")
	master.list <- list()
	for(i in 1:length(outcomes)){
		master$outcome <- outcomes[i]
		master$covList <- cov.list[i]
		master$survivalFunc <- survival.functions[i]
		master$genome <- genome
		master.list[[i]] <- master
	}

	master.merge <- data.table(do.call(rbind, master.list))

	master.list <- list()
	for(i in 1:length(cohorts)){
		master.merge$cohort <- cohorts[i]
		master.list[[i]] <- master.merge
	}
	do.call(rbind, master.list)


}

# donors
master.dh1h1 <- create.master("Dh1h1", "D", outcomes, cohorts, diseases, survival.functions)
# recipients h1h1
master.rh1h1 <- create.master("Rh1h1", "R", outcomes, cohorts, diseases, survival.functions)

# shared
master.grp2vsgrp1 <- create.master("grp2vsgrp1", "S", outcomes, cohorts, diseases, survival.functions)
master.grp3vsgrp1 <- create.master("grp3vsgrp1", "S", outcomes, cohorts, diseases, survival.functions)
master.grp3vsgrp2 <- create.master("grp3vsgrp2", "S", outcomes, cohorts, diseases, survival.functions)



master.list <- data.table(do.call(rbind, list(master.dh1h1, master.rh1h1, master.grp2vsgrp1, master.grp3vsgrp1, master.grp3vsgrp2)))



outcome.order <- c("DD", "PFS", "OS", "TRM")
# order by outcome so we can have DD and PFS as top two outcomes b/c they both have 2 covariates
master.list.2cov <- master.list[order(match(master.list$outcome, outcome.order))][1:100]
# OS has 3 covariates so we will grab those
master.list.3cov <- master.list[order(match(master.list$outcome, outcome.order))][101:150]

# trm has 4 covs
master.list.4cov <- master.list[order(match(master.list$outcome, outcome.order))][151:200]


## subset by cohort and disease
cohort.disease <- function(pheno, master.list){
	if(master.list$cohort=="c1"){
		pheno <- pheno[cohort1==1]
	} else {
		pheno <- pheno[cohort1==0]
	}

	if(master.list$disease=="mixed"){
		data.frame(pheno)
	} else if (master.list$disease=="AMLonly"){
		data.frame(pheno[AMLdummy==1])
	} else if(master.list$disease=="ALLonly"){
		data.frame(pheno[ALLdummy==1])
	} else if (master.list$disease == "noALL"){
		data.frame(pheno[ALLdummy==0])
	} else if (master.list$disease == "noMDS"){
		data.frame(pheno[MDSdummy==0])
	}

}




# PFS 
surv.res.ddpfs <- data.frame(matrix(nrow=100, ncol=17))
colnames(surv.res.ddpfs) <- cols
for(i in 1:nrow(surv.res.ddpfs)){
	surv.res.ddpfs[i,] <- make.tbl(
		surv_fit_cov2(master.list.2cov$genotype[i],
			eval(as.name(paste(master.list.2cov$covList[i]))),
			cohort.disease(merged.ea, master.list.2cov[i,])),
		"CCR5",
		paste(master.list.2cov$genotype[i], master.list.2cov$outcome[i], master.list.2cov$disease[i], master.list.2cov$genome[i], sep="_"),
		"chr3",
		"*",
		cohort.disease(merged.ea, master.list.2cov[i,]) %>% select(cohort, eval(as.name(paste(master.list.2cov$genotype[i])))) %>% filter(cohort==master.list.2cov$cohort[i]) %>% na.omit %>% nrow,
		"*",
		"*",
		"imputed",
		master.list.2cov$genome[i],
		master.list.2cov$cohort[i],
		master.list.2cov$outcome[i],
		master.list.2cov$disease[i])
	surv.res.ddpfs
}



# OVERALL SURVIVAL 
surv.res.os <- data.frame(matrix(nrow=50, ncol=17))
colnames(surv.res.os) <- cols
for(i in 1:nrow(surv.res.os)){
	surv.res.os[i,] <- make.tbl(
		surv_fit_cov2(master.list.3cov$genotype[i],
			eval(as.name(paste(master.list.3cov$covList[i]))),
			cohort.disease(merged.ea, master.list.3cov[i,])),
		"CCR5",
		paste(master.list.3cov$genotype[i], master.list.3cov$outcome[i], master.list.3cov$disease[i], master.list.3cov$genome[i], sep="_"),
		"chr3",
		"*",
		cohort.disease(merged.ea, master.list.3cov[i,]) %>% select(cohort, eval(as.name(paste(master.list.3cov$genotype[i])))) %>% filter(cohort==master.list.3cov$cohort[i]) %>% na.omit %>% nrow,
		"*",
		"*",
		"imputed",
		master.list.3cov$genome[i],
		master.list.3cov$cohort[i],
		master.list.3cov$outcome[i],
		master.list.3cov$disease[i])
	surv.res.os
}



# TRM 
surv.res.trm <- data.frame(matrix(nrow=50, ncol=17))
colnames(surv.res.trm) <- cols
for(i in 1:nrow(surv.res.trm)){
	surv.res.trm[i,] <- make.tbl(
		surv_fit_cov2(master.list.4cov$genotype[i],
			eval(as.name(paste(master.list.4cov$covList[i]))),
			cohort.disease(merged.ea, master.list.4cov[i,])),
		"CCR5",
		paste(master.list.4cov$genotype[i], master.list.4cov$outcome[i], master.list.4cov$disease[i], master.list.4cov$genome[i], sep="_"),
		"chr3",
		"*",
		cohort.disease(merged.ea, master.list.4cov[i,]) %>% select(cohort, eval(as.name(paste(master.list.4cov$genotype[i])))) %>% filter(cohort==master.list.4cov$cohort[i]) %>% na.omit %>% nrow,
		"*",
		"*",
		"imputed",
		master.list.4cov$genome[i],
		master.list.4cov$cohort[i],
		master.list.4cov$outcome[i],
		master.list.4cov$disease[i])
	surv.res.trm
}


## save to file
surv.res <- data.table(do.call(rbind, list(surv.res.ddpfs, surv.res.os, surv.res.trm)))
# split into cohorts
surv.res.c1 <- surv.res[cohort=="c1"]
surv.res.c2 <- surv.res[cohort=="c2"]

#write.table(surv.res.c1, file="h1h1_nometa_c1.txt", sep="\t", quote=F, row.names=F, col.names=T)
#write.table(surv.res.c2, file="h1h1_nometa_c2.txt", sep="\t", quote=F, row.names=F, col.names=T)


## RAN META ANALYSIS ###

## LOAD META RESULTS

surv.res.meta <- fread("metal_ccr5_full1.tbl")
surv.res.meta <- surv.res.meta[,-7]

meta.res <- function(ccr5.meta, col.order){
	# remove direction column
	# build up columns to match how our final table is set up
	ccr5.meta$gene <- "CCR5"
	colnames(ccr5.meta)[1] <- "rsID"
	ccr5.meta$chr <- "chr3"
	ccr5.meta$exp.coef <- exp(ccr5.meta$Effect)
	ccr5.meta$genome <- sapply(strsplit(surv.res.meta$MarkerName, '_', fixed=T), "[", 4)
	ccr5.meta$cohort <- "M"
	ccr5.meta$disease <- sapply(strsplit(surv.res.meta$MarkerName, '_', fixed=T), "[", 3)
	ccr5.meta$N <- NA
	ccr5.meta$BP <- "*"
	ccr5.meta$impute <- "imputed"
	ccr5.meta$`95%-CI` <- NA
	ccr5.meta$outcome <- sapply(strsplit(surv.res.meta$MarkerName, '_', fixed=T), "[", 2)
	colnames(ccr5.meta) <- c("rsID", "allele1", "allele2", "coef", "se.coef", "Pvalue", "gene", "chr", "exp.coef",
	 "genome", "cohort", "disease", "N", "BP", "impute", "95%-CI", "outcome")
	setcolorder(ccr5.meta, cols)
	data.table(ccr5.meta)
}

surv.res.meta <- meta.res(surv.res.meta,cols)
# need to make sure in right order for N


setkey(surv.res.c1, rsID)
setkey(surv.res.c2, rsID)
setkey(surv.res.meta, rsID)
surv.res.meta$N <- as.numeric(surv.res.c1$N) + as.numeric(surv.res.c2$N)

surv.res.full <- do.call(rbind, list(surv.res.c1, surv.res.c2, surv.res.meta))

surv.res.full$rsID <- sapply(strsplit(surv.res.full$rsID, '_', fixed=T), "[", 1)


ccr5.final <- surv.res.full


# calculate confidence interval
hr.ci <- function(coef.est, se){
        lb <- round(exp(coef.est-1.96*se), 5)
        ub <- round(exp(coef.est+1.96*se), 5)
        paste0("[", lb,", " , ub ,"]")
}

for(i in 1:nrow(ccr5.final)){
        ccr5.final[i,"95%-CI"] <- hr.ci(as.numeric(ccr5.final$coef)[i], as.numeric(ccr5.final$se.coef)[i])
}

write.table(ccr5.final, file="CCR5_H1H1_FINAL_wMETA.txt", sep="\t", quote=F, row.names=F, col.names=T)