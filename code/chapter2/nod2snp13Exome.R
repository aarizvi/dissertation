## Directory that Original Exome Data is in:

#/projects/rpci/lsuchest/alyssacl/Rare_Variant/Original.Exome.Files

# Use plink to parse out just that snp from the map and ped files...

# exm-rs2066847

#module load plink
#plink --file /projects/rpci/lsuchest/alyssacl/Rare_Variant/Original.Exome.Files/Copy_BMT-Exome-Final01-20130911_BC-PC_Cluster --snp exm-rs2066847 --recode --out exm-rs2066847


### FIRST WE LOAD UP EXOME rs2066847 (SNP13)
## now we load into R

library(data.table)
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

exome.geno <- data.frame(exome.snp13)

# now we will recode these and drop the extra columns so we just have pair ids and genotypes
# if I/I = 1, I/D = 1, D/D = 0
exome.geno <- apply(exome.geno, 2, FUN=function(x) gsub("I/I", 1, x))
exome.geno <- apply(exome.geno, 2, FUN=function(x) gsub("I/D", 1, x))
exome.geno <- apply(exome.geno, 2, FUN=function(x) gsub("D/D", 0, x))

exome.geno <- data.table(exome.geno)
exome.geno <- exome.geno[,c("pair_id", "R_snp13", "D_snp13")]

recs.snp13 <- exome.geno[,c("pair_id", "R_snp13")]
donors.snp13 <- exome.geno[,c("pair_id", "D_snp13")]

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


recs <- recs[,c("pair_id", "rs2066844", "rs2066845")]

donors <- donors[,c("pair_id", "rs2066844", "rs2066845")]

########################################################
############# NOW MERGE SNP13 INTO THESE FILES
######################################################


recs <- data.table(apply(recs, 2, FUN=function(x) as.numeric(x)))
donors <- data.table(apply(donors, 2, FUN=function(x) as.numeric(x)))

recs.snp13 <- data.table(apply(recs.snp13, 2, FUN=function(x) as.numeric(x)))
donors.snp13 <- data.table(apply(donors.snp13, 2, FUN=function(x) as.numeric(x)))

setkey(recs, pair_id)
setkey(donors, pair_id)
setkey(recs.snp13, pair_id)
setkey(donors.snp13, pair_id)

recs <- recs[recs.snp13]
donors <- donors[donors.snp13]
###********######
####### NOW SNP13 IS MERGED BACK TO PATIENTS
## WE HAVE 2348 DONOR AND RECIPIENT PAIRS

# ok lets merge just the snp13 stuff because we need to test just that snp in donor and recipients and then do the meta analysis, etc.
setkey(merged.ea, pair_id)
setkey(donors.snp13, pair_id)
setkey(recs.snp13, pair_id)
merged.ea <- merged.ea[recs.snp13]
setkey(merged.ea, pair_id)
merged.ea <- merged.ea[donors.snp13]

##############################################
####### COMPOSITE SNP (SNP8/SNP12/SNP13) #####
##############################################

# snp 13 = rs146528649 x[,1]
recode <- function(x){
#        x$snp8snp12 <- ifelse((x[,2] == 0 & x[,3] == 0), x$snp8snp12 <- 0, x$snp8snp12 <- 1)
        x$compSNPs <- ifelse((x[,2] == 0 & x[,2] == 0 & x[,3] == 0), x$compSNPs <- 0, x$compSNPs <- 1)
        data.table(x, keep.rownames = T)
}
recs <- recode(recs)
colnames(recs)[ncol(recs)] <- paste("R", colnames(recs)[ncol(recs)], sep="_")
donors <- recode(donors)
colnames(donors)[ncol(donors)] <- paste("D", colnames(donors)[ncol(donors)], sep="_")

dr.comp <- data.table(cbind(recs$pair_id, recs$R_compSNPs, donors$D_compSNPs))
colnames(dr.comp) <- c("pair_id", "recipient", "donor")


# function to recode 'score', again, if there is 0 in both D-R pairs the score 0, otherwise its 1
comp_score <- function(x){
        x$NOD2_score <- ifelse((x[,"recipient"] == 0 & x[,"donor"] == 0), x$NOD2_score <- 0, x$NOD2_score <- 1)
        x
}

dr.comp <- comp_score(dr.comp)


#############################################
#### MERGE NOD2 COMP SCORE BACK TO PHENO FILE
#############################################

merged.ea$NOD2_score <- dr.comp$NOD2_score

## forgot to split up by cohorts...
colnames(merged.ea)[which(colnames(merged.ea) == "population")] <- "cohort"
merged.ea$cohort <- gsub("EA.", "", merged.ea$cohort)



merged.ea.c1 <- merged.ea[cohort=="c1"] # has 1697 patients
merged.ea.c2 <- merged.ea[cohort=="c2"] # has 467 patients


# load the covariate file....
r.cov <- fread("/projects/rpci/lsuchest/lsuchest/Rserve/TheData/Plink_Recipient.cov")
r.cov <- r.cov[r.cov$pair_id %in% merged.ea$pair_id]
setkey(r.cov, pair_id)
setkey(merged.ea, pair_id)

r.cov$cohort <- merged.ea$cohort[match(r.cov$pair_id, merged.ea$pair_id)]

r.cov.c1 <- r.cov[cohort=="c1"]
r.cov.c2 <- r.cov[cohort=="c2"]

# merge the genotypes into the covariate file so we just have one covariate file...
merged.ea.c1 <- merged.ea.c1[,c("pair_id", "R_snp13", "D_snp13", "NOD2_score")]
merged.ea.c2 <- merged.ea.c2[,c("pair_id", "R_snp13", "D_snp13", "NOD2_score")]

setkey(merged.ea.c1, pair_id)
setkey(merged.ea.c2, pair_id)
setkey(r.cov.c1, pair_id)
setkey(r.cov.c2, pair_id)
r.cov.geno.c1 <- data.frame(r.cov.c1[merged.ea.c1])
r.cov.geno.c2 <- data.frame(r.cov.c2[merged.ea.c2])


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


# DD HAS 2 COVARIATES
surv_fit_cov2("NOD2_score", DDcov, r.cov.geno.c1)
surv_fit_cov2("NOD2_score", DDcov, r.cov.geno.c2)


surv_fit_cov2("R_snp13", DDcov, r.cov.geno.c1)
surv_fit_cov2("R_snp13", DDcov, r.cov.geno.c2)



surv_fit_cov2("D_snp13", DDcov, r.cov.geno.c1)
surv_fit_cov2("D_snp13", DDcov, r.cov.geno.c2)


# PFS HAS 2 COVARIATES
surv_fit_cov2("NOD2_score", PFScov, r.cov.geno.c1)
surv_fit_cov2("NOD2_score", PFScov, r.cov.geno.c2)

surv_fit_cov2("R_snp13", PFScov, r.cov.geno.c1)
surv_fit_cov2("R_snp13", PFScov, r.cov.geno.c2)

surv_fit_cov2("D_snp13", PFScov, r.cov.geno.c1)
surv_fit_cov2("D_snp13", PFScov, r.cov.geno.c2)



# OS has 3 COVARIATES
surv_fit_cov3("NOD2_score", OScov, r.cov.geno.c1)
surv_fit_cov3("NOD2_score", OScov, r.cov.geno.c2)

surv_fit_cov3("R_snp13", OScov, r.cov.geno.c1)
surv_fit_cov3("R_snp13", OScov, r.cov.geno.c2)

surv_fit_cov3("D_snp13", OScov, r.cov.geno.c1)
surv_fit_cov3("D_snp13", OScov, r.cov.geno.c2)




# TRM HAS 4 COVARIATES
surv_fit_cov4("NOD2_score", TRMcov, r.cov.geno.c1)
surv_fit_cov4("NOD2_score", TRMcov, r.cov.geno.c2)

surv_fit_cov4("R_snp13", TRMcov, r.cov.geno.c1)
surv_fit_cov4("R_snp13", TRMcov, r.cov.geno.c2)

surv_fit_cov4("D_snp13", TRMcov, r.cov.geno.c1)
surv_fit_cov4("D_snp13", TRMcov, r.cov.geno.c2)


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
## NOD2 COMPOSITE SNP RESULTS
## Cohort 1
###############################
nod2.c1 <- data.table(rbind(make.tbl(surv_fit_cov2("NOD2_score", DDcov, r.cov.geno.c1), "NOD2", "3 SNPs* DD", "chr16", "*", 1697, "*", "*", "imputed", "S", "c1", "DD", "mixed"),
                            make.tbl(surv_fit_cov2("NOD2_score", PFScov, r.cov.geno.c1), "NOD2", "3 SNPs* PFS", "chr16", "*", 1697, "*", "*","imputed", "S", "c1", "PFS", "mixed"),
                            make.tbl(surv_fit_cov3("NOD2_score", OScov, r.cov.geno.c1), "NOD2", "3 SNPs* OS", "chr16", "*", 1697, "*", "*", "imputed", "S", "c1", "OS", "mixed"),
                            make.tbl(surv_fit_cov4("NOD2_score", TRMcov, r.cov.geno.c1), "NOD2", "3 SNPs* TRM", "chr16", "*", 1697, "*", "*", "imputed", "S", "c1", "TRM", "mixed")))
colnames(nod2.c1) <- cols

write.table(nod2.c1, file="nod2.comp.c1.txt", col.names=T, row.names=F, quote=F, sep='\t')

###############################
## NOD2 COMPOSITE SNP RESULTS
## Cohort 2
###############################
nod2.c2 <- data.table(rbind(make.tbl(surv_fit_cov2("NOD2_score", DDcov, r.cov.geno.c2), "NOD2", "3 SNPs* DD", "chr16", "*", 467, "*", "*", "imputed" ,"S", "c2", "DD", "mixed"),
                            make.tbl(surv_fit_cov2("NOD2_score", PFScov, r.cov.geno.c2), "NOD2", "3 SNPs* PFS", "chr16", "*", 467, "*", "*", "imputed", "S", "c2", "PFS", "mixed"),
                            make.tbl(surv_fit_cov3("NOD2_score", OScov, r.cov.geno.c2), "NOD2", "3 SNPs* OS", "chr16", "*", 467, "*", "*","imputed", "S", "c2", "OS", "mixed"),
                            make.tbl(surv_fit_cov4("NOD2_score", TRMcov, r.cov.geno.c2), "NOD2", "3 SNPs* TRM", "chr16", "*", 467, "*", "*", "imputed", "S", "c2", "TRM", "mixed")))
colnames(nod2.c2) <- cols

write.table(nod2.c2, file="nod2.comp.c2.txt", col.names=T, row.names=F, quote=F, sep='\t')

###################################
## NOD2 SNP13 RECIPIENT SNP RESULTS
## Cohort 1
###################################
nod2.snp13.r.c1 <- data.table(rbind(make.tbl(surv_fit_cov2("R_snp13", DDcov, r.cov.geno.c1), "NOD2", "rs2066847 DD", "chr16", "*", 1697, "*", "*", "typed", "S", "c1", "DD", "mixed"),
                            make.tbl(surv_fit_cov2("R_snp13", PFScov, r.cov.geno.c1), "NOD2", "rs2066847 PFS", "chr16", "*", 1697, "*", "*", "typed", "S", "c1", "PFS", "mixed"),
                            make.tbl(surv_fit_cov3("R_snp13", OScov, r.cov.geno.c1), "NOD2", "rs2066847 OS", "chr16", "*", 1697, "*", "*", "typed", "S", "c1", "OS", "mixed"),
                            make.tbl(surv_fit_cov4("R_snp13", TRMcov, r.cov.geno.c1), "NOD2", "rs2066847 TRM", "chr16", "*", 1697, "*", "*", "typed", "S", "c1", "TRM", "mixed")))
colnames(nod2.snp13.r.c1) <- cols

write.table(nod2.snp13.r.c1, file="R_nod2.snp13.c1.txt", col.names=T, row.names=F, quote=F, sep='\t')


###################################
## NOD2 SNP13 RECIPIENT SNP RESULTS
## Cohort 2
###################################
nod2.snp13.r.c2 <- data.table(rbind(make.tbl(surv_fit_cov2("R_snp13", DDcov, r.cov.geno.c2), "NOD2", "rs2066847 DD", "chr16", "*", 467, "*", "*","typed", "R", "c2", "DD", "mixed"),
                            make.tbl(surv_fit_cov2("R_snp13", PFScov, r.cov.geno.c2), "NOD2", "rs2066847 PFS", "chr16", "*", 467, "*", "*", "typed", "R", "c2", "PFS", "mixed"),
                            make.tbl(surv_fit_cov3("R_snp13", OScov, r.cov.geno.c2), "NOD2", "rs2066847 OS", "chr16", "*", 467, "*", "*","typed", "R", "c2", "OS", "mixed"),
                            make.tbl(surv_fit_cov4("R_snp13", TRMcov, r.cov.geno.c2), "NOD2", "rs2066847 TRM", "chr16", "*", 467, "*", "*", "typed", "R", "c2", "TRM", "mixed")))
colnames(nod2.snp13.r.c2) <- cols

write.table(nod2.snp13.r.c2, file="R_nod2.snp13.c2.txt", col.names=T, row.names=F, quote=F, sep='\t')

###################################
## NOD2 SNP13 DONOR SNP RESULTS
## Cohort 1
###################################
nod2.snp13.d.c1 <- data.table(rbind(make.tbl(surv_fit_cov2("D_snp13", DDcov, r.cov.geno.c1), "NOD2", "rs2066847 DD", "chr16", "*", 1697, "*", "*", "typed", "D", "c1", "DD", "mixed"),
                            make.tbl(surv_fit_cov2("D_snp13", PFScov, r.cov.geno.c1), "NOD2", "rs2066847 PFS", "chr16", "*", 1697, "*", "*", "typed", "D", "c1", "PFS", "mixed"),
                            make.tbl(surv_fit_cov3("D_snp13", OScov, r.cov.geno.c1), "NOD2", "rs2066847 OS", "chr16", "*", 1697, "*", "*", "typed","D", "c1", "OS", "mixed"),
                            make.tbl(surv_fit_cov4("D_snp13", TRMcov, r.cov.geno.c1), "NOD2", "rs2066847 TRM", "chr16", "*", 1697, "*", "*", "typed","D", "c1", "TRM", "mixed")))
colnames(nod2.snp13.d.c1) <- cols

write.table(nod2.snp13.d.c1, file="D_nod2.snp13.c1.txt", col.names=T, row.names=F, quote=F, sep='\t')

###################################
## NOD2 SNP13 DONOR SNP RESULTS
## Cohort 2
###################################
nod2.snp13.d.c2 <- data.table(rbind(make.tbl(surv_fit_cov2("D_snp13", DDcov, r.cov.geno.c2), "NOD2", "rs2066847 DD", "chr16", "*", 467, "*", "*", "typed", "D", "c2", "DD", "mixed"),
                            make.tbl(surv_fit_cov2("D_snp13", PFScov, r.cov.geno.c2), "NOD2", "rs2066847 PFS", "chr16", "*", 467, "*", "*", "typed","D", "c2", "PFS", "mixed"),
                            make.tbl(surv_fit_cov3("D_snp13", OScov, r.cov.geno.c2), "NOD2", "rs2066847 OS", "chr16", "*", 467, "*", "*", "typed","D", "c2", "OS", "mixed"),
                            make.tbl(surv_fit_cov4("D_snp13", TRMcov, r.cov.geno.c2), "NOD2", "rs2066847 TRM", "chr16", "*", 467, "*", "*", "typed", "D", "c2", "TRM", "mixed")))
colnames(nod2.snp13.d.c2) <- cols

write.table(nod2.snp13.d.c2, file="D_nod2.snp13.c2.txt", col.names=T, row.names=F, quote=F, sep='\t')

####################
### META ANALYSIS
####################

# Function to alter meta-results

meta.res <- function(nod2.meta, genome, col.order, impute){
	# remove direction column
	nod2.meta <- nod2.meta[,-7]
	# build up columns to match how our final table is set up
	nod2.meta$gene <- "NOD2"
	colnames(nod2.meta)[1] <- "rsID"
	nod2.meta$chr <- "chr16"
	nod2.meta$exp.coef <- exp(nod2.meta$Effect)
	nod2.meta$genome <- genome
	nod2.meta$cohort <- "M"
	nod2.meta$disease <- "mixed"
	nod2.meta$N <- 2164
	nod2.meta$BP <- "*"
	nod2.meta$impute <- impute
	nod2.meta$`95%-CI` <- NA
	nod2.meta$outcome <- NA
	colnames(nod2.meta) <- c("rsID", "allele1", "allele2", "coef", "se.coef", "Pvalue", "gene", "chr", "exp.coef",
	 "genome", "cohort", "disease", "N", "BP", "impute", "95%-CI", "outcome")
	setcolorder(nod2.meta, cols)
	nod2.meta
}



###############################
## NOD2 COMPOSITE SNP RESULTS
## META-ANALYSIS
###############################
nod2.comp.meta <- fread("metal_S_M_mixed_nod21.tbl")
nod2.comp.meta <- meta.res(nod2.comp.meta, "S", cols, "imputed")

###################################
## NOD2 SNP13 DONOR SNP RESULTS
## META-ANALYSIS
###################################
donor.snp13 <- fread("metal_D_M_mixed_snp131.tbl")
donor.snp13 <- meta.res(donor.snp13, "D", cols, "typed")

###################################
## NOD2 SNP13 RECIPIENT SNP RESULTS
## META-ANALYSIS
###################################
rec.snp13 <- fread("metal_R_M_mixed_snp131.tbl")
rec.snp13 <- meta.res(rec.snp13, "R", cols, "typed")




##################################
### MERGE ALL THE TABLES TOGETHER
##################################
nod2.final <- do.call(rbind, list(nod2.c1,
				  nod2.c2,
				  nod2.comp.meta,
				  nod2.snp13.d.c1,
				  nod2.snp13.d.c2,
				  donor.snp13,
				  nod2.snp13.r.c1,
				  nod2.snp13.r.c2,
				  rec.snp13))
nod2.final$outcome[9:12] <- c("OS", "DD", "PFS", "TRM")

nod2.final$outcome[21:24] <- c("DD", "OS", "TRM", "PFS")

nod2.final$outcome[33:36] <- c("DD", "OS", "TRM", "PFS")

nod2.final$rsID <- gsub("OS", "", nod2.final$rsID)
nod2.final$rsID <- gsub("DD", "", nod2.final$rsID)
nod2.final$rsID <- gsub("PFS", "", nod2.final$rsID)
nod2.final$rsID <- gsub("TRM", "", nod2.final$rsID)

## now calculate 95% confidence interval for hazard ratio
hr.ci <- function(coef.est, se){
        lb <- round(exp(coef.est-1.96*se), 5)
        ub <- round(exp(coef.est+1.96*se), 5)
        paste0("[", lb,", " , ub ,"]")
}

for(i in 1:nrow(nod2.final)){
        nod2.final[i,"95%-CI"] <- hr.ci(as.numeric(nod2.final$coef)[i], as.numeric(nod2.final$se.coef)[i])
}



write.table(nod2.final, file="nod2_survival_pvals.txt", col.names=T, row.names=F, sep="\t", quote=F)





















