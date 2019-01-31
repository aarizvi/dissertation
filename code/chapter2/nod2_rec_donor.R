## NOD2 
library(dplyr)
library(data.table)
library(dtplyr)
library(tidyr)
library(survival)

snp13 <- fread("/projects/rpci/lsuchest/abbasriz/candidate_gene/cg_haplotypes/nod2_rep/snp13_exome/exm-rs2066847.ped")
colnames(snp13) <- c("SampleIndex", "IID", "PID", "MID", "Sex", "Affection", "allele1", "allele2")
snp13 <- within(snp13, geno <- paste(allele1, allele2, sep="/"))
table(snp13$geno)

## Grab clinical data

clin <- fread("/projects/rpci/lsuchest/lsuchest/Rserve/TheData/BMT_PHENOTYPE_FILE_QC_allclinical_5_15_15.txt")

# grab just exome ids with pair ids so we can map back to other patients
exome_ids <- na.omit(clin[,c("pair_id", "Exome.R_IID", "Exome.D_IID")])
setkey(exome_ids, pair_id)

exome.geno <- exome_ids
# ensure they are in the same order
exome.geno$R_snp13 <- snp13$geno[match(exome_ids$Exome.R_IID, snp13$IID)]
exome.geno$D_snp13 <- snp13$geno[match(exome_ids$Exome.D_IID, snp13$IID)]

table(exome.geno$R_snp13)
# D/D  I/D  I/I 
# 2236  109    3 
table(exome.geno$D_snp13)
# D/D  I/D  
# 2222  126

exome.geno <- data.frame(exome.geno)

# now we will recode these and drop the extra columns so we just have pair ids and genotypes
# if I/I = 1, I/D = 1, D/D = 0
exome.geno <- apply(exome.geno, 2, FUN=function(x) gsub("I/I", 1, x))
exome.geno <- apply(exome.geno, 2, FUN=function(x) gsub("I/D", 1, x))
exome.geno <- apply(exome.geno, 2, FUN=function(x) gsub("D/D", 0, x))

exome.geno <- data.table(exome.geno)
exome.geno <- exome.geno[,c("pair_id", "R_snp13", "D_snp13")]

recs.snp13 <- exome.geno[,c("pair_id", "R_snp13")]
donors.snp13 <- exome.geno[,c("pair_id", "D_snp13")]

recs.snp13$pair_id <- as.numeric(recs.snp13$pair_id)
donors.snp13$pair_id <- as.numeric(donors.snp13$pair_id)


#####################################
## LOAD UP THE SNP 8 and SNP 12 DATA
#####################################

library(VariantAnnotation)
## read in vcf file
vcf <- readVcf("/projects/rpci/lsuchest/abbasriz/candidate_gene/cg_haplotypes/nod2_rep/nod2_rep_dosages.vcf", genome="hg19")
vcf

##########################################
##### NOW DEAL WITH GENOTYPE FILE ########
##########################################


## threshold genotypes
# grab genotypes of 3 SNPs for all patients
gt <- geno(vcf)$GT
gt[1:nrow(gt),1:5]

# only we have both typed and imputed of rs2066844
# so we remove the first rs2066844 because it is typed and we are keeping just the imputed ones 
# we also removed the LD SNP13 because we now have the SNP13 from the EXOME CHIP
gt <- gt[-c(1:2),]
gt[1:nrow(gt),1:5]


# transpose the data.frame so we can patients as rows
gt <- data.frame(t(gt))
gt[1:5, 1:ncol(gt)]


## coding dominant model over df columns
gt <- apply(gt, 2, FUN=function(x) gsub("[.]", NA, x))
gt <- apply(gt, 2, FUN=function(x) gsub("0/0", 0, x))
gt <- apply(gt, 2, FUN=function(x) gsub("0/1", 1, x))
gt <- apply(gt, 2, FUN=function(x) gsub("1/1", 1, x))

## create a data.frame that changes the factors into numeric and keeps the row.names ...
gt <- data.frame(apply(gt, 2, function(x) as.numeric(as.character(x))), check.names=F, row.names=row.names(gt))
head(gt)


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
recs <- data.table(gt[ids$r_sample_index,], keep.rownames=T)
donors <- data.table(gt[ids$d_sample_index,], keep.rownames=T)

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


recs <- recs[,c("pair_id", "rs2066844", "rs2066845")] %>% na.omit()


donors <- donors[,c("pair_id", "rs2066844", "rs2066845")] %>% na.omit()





###########################################################
############# NOW MERGE SNP13 INTO THESE FILES AND RECODE 
###########################################################

# function for the composite snps
recode <- function(x){
#        x$snp8snp12 <- ifelse((x[,2] == 0 & x[,3] == 0), x$snp8snp12 <- 0, x$snp8snp12 <- 1)
        x$compSNPs <- ifelse((x[,2] == 0 & x[,3] == 0 & x[,4] == 0), x$compSNPs <- 0, x$compSNPs <- 1)
        data.table(x)
}

recs.nod2 <- recs %>% dtplyr:::inner_join.data.table(recs.snp13, by="pair_id") %>% na.omit() %>% recode() %>% dplyr::rename(R_compSNPs=compSNPs)
# 2048 patients
donors.nod2 <- donors %>% dtplyr:::inner_join.data.table(donors.snp13, by="pair_id") %>% na.omit() %>% recode() %>% dplyr::rename(D_compSNPs=compSNPs)
# 2056 donors
# now merge these back to merged.ea

merged.ea <- recs.nod2 %>% dplyr::select(pair_id, R_compSNPs) %>% dtplyr:::inner_join.data.table(merged.ea, by="pair_id")
merged.ea <- donors.nod2 %>% dplyr::select(pair_id, D_compSNPs) %>% dtplyr:::inner_join.data.table(merged.ea, by="pair_id")
merged.ea <- merged.ea %>% dplyr::rename(cohort=population) %>% mutate(S_compSNPs =  ifelse(R_compSNPs== 0 & D_compSNPs == 0, 0, 1))
merged.ea$cohort <- gsub("EA.", "", merged.ea$cohort)


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
## OS covariates are: age, distatD, PBlood
## PFS covariates: age, distatD
## TRM covariates: age, bmiOBS, bmiOVWT, PBlood
###################################################

DDcov <- c("intxsurv_1Y","disease_death_1Y","age","distatD")
OScov <- c("intxsurv_1Y","dead_1Y","age","distatD","PBlood")
PFScov<- c("intxrel_1Y","lfs_1Y","age","distatD")
TRMcov <- c("intxsurv_1Y","TRM_1Y","age","bmiOBS","bmiOVWT","PBlood")


# now merge covariates to merged.ea

merged.ea <- r.cov %>% dplyr::select(pair_id, age, distatD, PBlood, bmiOBS, bmiOVWT, MDSdummy, AMLdummy, ALLdummy) %>% right_join(merged.ea, by="pair_id")



# ## MIXED PATIENTS
# # 1964 patients TOTAL IN MIXED
# mixed.c1 <- data.frame(merged.ea[cohort=="c1"]) # has 1548 patients
# mixed.c2 <- data.frame(merged.ea[cohort=="c2"]) # has 416 patients

# ## ALLonly
# allonly.c1 <- data.frame(merged.ea[cohort == "c1" & ALLdummy==1]) # 357 in ALLonly Cohort 1
# allonly.c2 <- data.frame(merged.ea[cohort == "c2" & ALLdummy==1]) # 34 in ALLonly Cohort 2

# ## AMLonly
# amlonly.c1 <- data.frame(merged.ea[cohort=="c1" & AMLdummy==1]) #938 in AML only cohort 1
# amlonly.c2 <- data.frame(merged.ea[cohort=="c2" & AMLdummy==1]) #268 AML Cohort 2


# ## noMDS (AML + ALL)
# nomds.c1 <- data.frame(merged.ea[cohort=="c1" & MDSdummy==0]) #1295 in cohort 1 no MDS 
# nomds.c2 <- data.frame(merged.ea[cohort=="c2" & MDSdummy==0]) #302 in cohort 2 no MDS



# ## noALL (AML + MDS)
# noall.c1 <- data.frame(merged.ea[cohort=="c1" & ALLdummy==0]) # 1191 in cohort 1
# noall.c2 <- data.frame(merged.ea[cohort=="c2" & ALLdummy==0]) # 382 in cohort 2


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


#### MASTER FILE CREATION #####
compsnps <- c("D_compSNPs", "R_compSNPs", "DR_compSNPs")
cov.list <- c("DDcov", "PFScov", "OScov", "TRMcov")
genomes <- c("D", "R", "S")
outcomes <- c("DD", "PFS", "OS", "TRM")
cohorts <- c("c1", "c2")
diseases <- c("mixed", "AMLonly", "ALLonly", "noMDS", "noALL")
survival.functions <- c("surv_fit_cov2", "surv_fit_cov2", "surv_fit_cov3", "surv_fit_cov4")


master <- data.frame(matrix(nrow=3, ncol=2))
#master <- list()
for(i in 1:length(diseases)){
	master[i,] <- c("compSNPs", diseases[i])

}
colnames(master) <- c("genotype", "disease")
master.list <- list()
for(i in 1:length(outcomes)){
	master$outcome <- outcomes[i]
	master$covList <- cov.list[i]
	master$survivalFunc <- survival.functions[i]
	master.list[[i]] <- master
}

master.merge <- data.table(do.call(rbind, master.list))

master.list <- list()
for(i in 1:length(cohorts)){
	master.merge$cohorts <- cohorts[i]
	master.list[[i]] <- master.merge
}

## double these because needed for c1 and c2
master.list <- list(master.list[[1]],master.list[[1]],master.list[[1]],master.list[[2]],master.list[[2]],master.list[[2]])

# create columns with genomes
add.genomes <- function(x, genomes){
	x <- x[,genome:=genomes]
	x
}

master.list <- mapply(add.genomes, master.list, genomes, SIMPLIFY=F)
# make composite SNP columns proper names
master.list <- lapply(master.list, function(x){ x$genotype <- paste0(x$genome, "_", x$genotype)
x})
master.list <- data.table(do.call(rbind, master.list))


## subset by cohort and disease
cohort.disease <- function(pheno, master.list){
	if(master.list$cohorts=="c1"){
		pheno <- pheno[cohort=="c1"]
	} else {
		pheno <- pheno[cohort=="c2"]
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



#### SURVIVAL ANALYSIS #####


outcome.order <- c("DD", "PFS", "OS", "TRM")
# order by outcome so we can have DD and PFS as top two outcomes b/c they both have 2 covariates
master.list.2cov <- master.list[order(match(master.list$outcome, outcome.order))][1:60]
# OS has 3 covariates so we will grab those
master.list.3cov <- master.list[order(match(master.list$outcome, outcome.order))][61:90]

# trm has 4 covs
master.list.4cov <- master.list[order(match(master.list$outcome, outcome.order))][91:120]


surv.res.ddpfs <- data.frame(matrix(nrow=60, ncol=17))
colnames(surv.res.ddpfs) <- cols
for(i in 1:nrow(surv.res.ddpfs)){
	surv.res.ddpfs[i,] <- make.tbl(surv_fit_cov2(master.list.2cov$genotype[i], eval(as.name(paste(master.list.2cov$covList[i]))), cohort.disease(merged.ea, master.list.2cov[i,])),  "NOD2", paste("3 SNPs*", master.list.2cov$outcome[i],  master.list.2cov$disease[i], master.list.2cov$genome[i], sep="_") , "chr16", "*", nrow(cohort.disease(merged.ea, master.list.2cov[i,])), "*", "*", "imputed", master.list.2cov$genome[i], master.list.2cov$cohort[i], master.list.2cov$outcome[i], master.list.2cov$disease[i])
	surv.res.ddpfs
}

surv.res.os <- data.frame(matrix(nrow=30, ncol=17))
colnames(surv.res.os) <- cols

for(i in 1:nrow(surv.res.os)){
	surv.res.os[i,] <- make.tbl(surv_fit_cov3(master.list.3cov$genotype[i], eval(as.name(paste(master.list.3cov$covList[i]))), cohort.disease(merged.ea, master.list.3cov[i,])),  "NOD2", paste("3 SNPs*", master.list.3cov$outcome[i], master.list.3cov$disease[i], master.list.3cov$genome[i], sep="_") , "chr16", "*", nrow(cohort.disease(merged.ea, master.list.3cov[i,])), "*", "*", "imputed", master.list.3cov$genome[i], master.list.3cov$cohort[i], master.list.3cov$outcome[i], master.list.3cov$disease[i])
	surv.res.os
}


surv.res.trm <- data.frame(matrix(nrow=30, ncol=17))
colnames(surv.res.trm) <- cols

for(i in 1:nrow(surv.res.trm)){
	surv.res.trm[i,] <- make.tbl(surv_fit_cov4(master.list.4cov$genotype[i], eval(as.name(paste(master.list.4cov$covList[i]))), cohort.disease(merged.ea, master.list.4cov[i,])),  "NOD2", paste("3 SNPs*", master.list.4cov$outcome[i], master.list.4cov$disease[i], master.list.4cov$genome[i],sep="_") , "chr16", "*", nrow(cohort.disease(merged.ea, master.list.4cov[i,])), "*", "*", "imputed", master.list.4cov$genome[i], master.list.4cov$cohort[i], master.list.4cov$outcome[i], master.list.4cov$disease[i])
	surv.res.trm
}


surv.res <- data.table(do.call(rbind, list(surv.res.ddpfs, surv.res.os, surv.res.trm)))
surv.res$N <- as.numeric(surv.res$N)
surv.res$coef <- as.numeric(surv.res$coef)
surv.res$se.coef <- as.numeric(surv.res$se.coef)
surv.res$exp.coef <- as.numeric(surv.res$exp.coef)
surv.res$Pvalue <- as.numeric(surv.res$Pvalue)

## separate into cohort 1 and cohort 2 so we can do meta analysis

surv.res.c1 <- surv.res[cohort=="c1"]
surv.res.c2 <- surv.res[cohort=="c2"]


#write.table(surv.res.c1, file="NOD2_FULL_c1.txt", sep="\t", col.names=T, row.names=F, quote=F)
#write.table(surv.res.c2, file="NOD2_FULL_c2.txt", sep="\t", col.names=T, row.names=F, quote=F)

## open metal results

surv.res.meta <- fread("metal_NOD2_FULL1.tbl")

surv.res.meta <- surv.res.meta[,-7]


meta.res <- function(nod2.meta, col.order){
	# remove direction column
	# build up columns to match how our final table is set up
	nod2.meta$gene <- "NOD2"
	colnames(nod2.meta)[1] <- "rsID"
	nod2.meta$chr <- "chr16"
	nod2.meta$exp.coef <- exp(nod2.meta$Effect)
	nod2.meta$genome <- sapply(strsplit(surv.res.meta$MarkerName, '_', fixed=T), "[", 4)
	nod2.meta$cohort <- "M"
	nod2.meta$disease <- sapply(strsplit(surv.res.meta$MarkerName, '_', fixed=T), "[", 3)
	nod2.meta$N <- NA
	nod2.meta$BP <- "*"
	nod2.meta$impute <- "imputed"
	nod2.meta$`95%-CI` <- NA
	nod2.meta$outcome <- sapply(strsplit(surv.res.meta$MarkerName, '_', fixed=T), "[", 2)
	colnames(nod2.meta) <- c("rsID", "allele1", "allele2", "coef", "se.coef", "Pvalue", "gene", "chr", "exp.coef",
	 "genome", "cohort", "disease", "N", "BP", "impute", "95%-CI", "outcome")
	setcolorder(nod2.meta, cols)
	nod2.meta
}

surv.res.meta <- meta.res(surv.res.meta,cols)
# need to make sure in right order for N
setkey(surv.res.c1, rsID)
setkey(surv.res.c2, rsID)
setkey(surv.res.meta, rsID)
surv.res.meta$N <- as.numeric(surv.res.c1$N) + as.numeric(surv.res.c2$N)

surv.res.full <- do.call(rbind, list(surv.res.c1, surv.res.c2, surv.res.meta))

surv.res.full$rsID <- substr(surv.res.full$rsID, 1, 7)

nod2.final <- surv.res.full


# calculate confidence interval
hr.ci <- function(coef.est, se){
        lb <- round(exp(coef.est-1.96*se), 5)
        ub <- round(exp(coef.est+1.96*se), 5)
        paste0("[", lb,", " , ub ,"]")
}

for(i in 1:nrow(nod2.final)){
        nod2.final[i,"95%-CI"] <- hr.ci(as.numeric(nod2.final$coef)[i], as.numeric(nod2.final$se.coef)[i])
}

write.table(nod2.final, file="NOD2_DRS_ALLOUTCOMES_FINAL.txt", sep="\t", quote=F, col.names=T, row.names=F)




