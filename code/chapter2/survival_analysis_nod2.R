library(VariantAnnotation)
## read in vcf file
vcf <- readVcf("/projects/rpci/lsuchest/abbasriz/candidate_gene/cg_haplotypes/nod2_rep/nod2_rep_dosages.vcf", genome="hg19")
vcf


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
# donor
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











##### NOW DEAL WITH GENOTYPE FILE ######
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



#### NOW PARSE GENOTYPE FILE BY SAMPLE ID INDICES
## parse the unique IDs (genome/cohort) from the main genotype df into 4 dfs specific to genome/cohort
recs <- data.table(gt[ids$r_sample_index,], keep.rownames=T)
donors <- data.table(gt[ids$d_sample_index,], keep.rownames=T)


## recode with composite coding -- 0 if snps are all 0, 1 if any of the snps are 1
# snp 8 =  rs2066844 x[,2]
# snp 12 = rs2066845 x[,3]
# snp 13 = rs146528649 x[,1]
recode <- function(x){
#        x$snp8snp12 <- ifelse((x[,2] == 0 & x[,3] == 0), x$snp8snp12 <- 0, x$snp8snp12 <- 1)
        x$compSNPs <- ifelse((x[,2] == 0 & x[,2] == 0 & x[,3] == 0), x$compSNPs <- 0, x$compSNPs <- 1)
        data.table(x, keep.rownames = T)
}
recs <- recode(recs)
donors <- recode(donors)

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

# now we have to parse down to exome chip sample size...














# function to recode 'score', again, if there is 0 in both D-R pairs the score 0, otherwise its 1
comp_score <- function(x){
        x$NOD2_score <- ifelse((x[,"recipient"] == 0 & x[,"donor"] == 0), x$NOD2_score <- 0, x$NOD2_score <- 1)
        x
}

# merge snp8/snp12 donor-recipient pair data.frame
snp8snp12.r <- recs[,c("snp8snp12", "pair_id")]
setkey(snp8snp12.r, pair_id)
snp8snp12.d <- donors[,c("snp8snp12", "pair_id")]
setkey(snp8snp12.d, pair_id)
snp8snp12 <- snp8snp12.d[snp8snp12.r]
colnames(snp8snp12) <- c("donor", "pair_id", "recipient")
snp8snp12 <- snp8snp12[,c("pair_id", "donor", "recipient")]
snp8snp12 <- comp_score(snp8snp12)

# merge snp8/snp12/snp13 donor-recipient pair data.frame
all3.r <- recs[,c("compSNPs", "pair_id")]
setkey(all3.r, pair_id)
all3.d <- donors[,c("compSNPs", "pair_id")]
setkey(all3.d, pair_id)
all3 <- all3.d[all3.r]
colnames(all3) <- c("donor", "pair_id", "recipient")
all3 <- all3[,c("pair_id", "donor", "recipient")]
all3 <- comp_score(all3)


colnames(snp8snp12)[ncol(snp8snp12)] <- "snp8snp12_score"
colnames(all3)[ncol(all3)] <- "all3_score"

snp8snp12 <- snp8snp12[,c("pair_id","snp8snp12_score")]
all3 <- all3[,c("pair_id", "all3_score")]


# merge back to phenotype file .. which is called merged.ea
setkey(merged.ea, pair_id)
setkey(snp8snp12, pair_id)
setkey(all3, pair_id)

merged.ea <- merged.ea[snp8snp12]
merged.ea <- merged.ea[all3]

# lets quickly see how many people had NA genotype scores
table(is.na(merged.ea$all3_score))
# 486 deemed NA for all 3
table(is.na(merged.ea$snp8snp12_score))
# 208 deemed NA for snp8snp12


## SO WE HAVE 2783 MATCHED DONOR-RECIPIENT PAIRS TOTAL (INCLUDING BOTH COHORT 1 AND COHORT 2)



## forgot to split up by cohorts...
colnames(merged.ea)[which(colnames(merged.ea) == "population")] <- "cohort"
merged.ea$cohort <- gsub("EA.", "", merged.ea$cohort)



merged.ea.c1 <- merged.ea[cohort=="c1"]
merged.ea.c2 <- merged.ea[cohort=="c2"]


#######################
##### SURVIVAL ANALYSIS
#######################


## do survival
library(survival)

# grab covariate files
r.cov <- fread("/projects/rpci/lsuchest/lsuchest/Rserve/TheData/Plink_Recipient.cov")

r.cov <- r.cov[r.cov$pair_id %in% ids$pair_id]
setkey(r.cov, pair_id)
r.cov$cohort <- merged.ea$cohort

r.cov.c1 <- r.cov[cohort=="c1"]
r.cov.c2 <- r.cov[cohort=="c2"]


##################################
## DD covariates are: age, distatD
##################################


## cohort 1

# snp 8/snp12
DD.c1 <- Surv(time=merged.ea.c1$intxsurv_1Y, event=merged.ea.c1$disease_death_1Y==1)
dd.snp8snp12.c1 <- coxph(DD.c1~merged.ea.c1$snp8snp12_score+r.cov.c1$distatD)
dd.snp8snp12.c1.coef <- summary(dd.snp8snp12.c1)$coef[1,1]
dd.snp8snp12.c1.se.coef <- summary(dd.snp8snp12.c1)$coef[1,3]
dd.snp8snp12.c1.hr <- summary(dd.snp8snp12.c1)$coef[1,2]
dd.snp8snp12.c1.pval <- summary(dd.snp8snp12.c1)$coef[1,5] 

# snp8/snp12/snp13
dd.all3.c1 <- coxph(DD.c1~merged.ea.c1$all3_score+r.cov.c1$distatD)


## cohort 2
# snp 8/snp12
DD.c2 <- Surv(time=merged.ea.c2$intxsurv_1Y, event=merged.ea.c2$disease_death_1Y==1)
dd.snp8snp12.c2 <- coxph(DD.c2~merged.ea.c2$snp8snp12_score+r.cov.c2$distatD)


# snp8/snp12/snp13
dd.all3.c2 <- coxph(DD.c2~merged.ea.c2$all3_score+r.cov.c2$distatD)
dd.all3.c2 <- coxph(DD.c2~merged.ea.c2$all3_score+r.cov.c2$distatD)


###########################################
## OS covariates are: age, distatD, Pblood
###########################################
        
#### OS ####
# cohort 1
## OS covariates are: age, distatD, Pblood
OS.c1 <- Surv(time=merged.ea.c1$intxsurv_1Y, event=merged.ea.c1$dead_1Y==1)
# OS SNP8/SNP12
os.snp8snp12.c1 <- coxph(OS.c1~merged.ea.c1$snp8snp12_score+r.cov.c1$age+r.cov.c1$distatD+r.cov.c1$PBlood)
summary(os.snp8snp12.c1)$coef[1,5] 

# OS SNP8/SNP12/SNP13
os.all3.c1 <- coxph(OS.c1~merged.ea.c1$all3_score+r.cov.c1$age+r.cov.c1$distatD+r.cov.c1$PBlood)
summary(os.all3.c1)$coef[1,5] 


## cohort 2
OS.c2 <- Surv(time=merged.ea.c2$intxsurv_1Y, event=merged.ea.c2$dead_1Y==1)
# OS SNP8/SNP12
os.snp8snp12.c2 <- coxph(OS.c2~merged.ea.c2$snp8snp12_score+r.cov.c2$age+r.cov.c2$distatD+r.cov.c2$PBlood)
summary(os.snp8snp12.c2)$coef[1,5] 

# OS SNP8/SNP12/SNP13
os.all3.c2 <- coxph(OS.c2~merged.ea.c2$all3_score+r.cov.c2$age+r.cov.c2$distatD+r.cov.c2$PBlood)
summary(os.all3.c2)$coef[1,5] 

####################################
## PFS covariates: age, distatD
####################################
# cohort 1
PFS.c1 <- Surv(time=merged.ea.c1$intxsurv_1Y, event=merged.ea.c1$lfs_1Y==1)
# PFS SNP8 / SNP12
pfs.snp8snp12.c1 <- coxph(PFS.c1~merged.ea.c1$snp8snp12_score+r.cov.c1$distatD+r.cov.c1$age)
summary(pfs.snp8snp12.c1)$coefficients[1,5]
# PFS SNP 8 / SNP 12 / SNP 13
pfs.all3.c1 <- coxph(PFS.c1~as.numeric(merged.ea.c1$all3_score)+r.cov.c1$distatD+as.numeric(r.cov.c1$age))
summary(pfs.all3.c1)$coefficients[1,5]


# cohort 2
PFS.c2 <- Surv(time=merged.ea.c2$intxsurv_1Y, event=merged.ea.c2$lfs_1Y==1)
# PFS SNP8 / SNP12
pfs.snp8snp12.c2 <- coxph(PFS.c2~merged.ea.c2$snp8snp12_score+r.cov.c2$distatD+r.cov.c2$age)
summary(pfs.snp8snp12.c2)$coefficients[1,5]
# PFS SNP 8 / SNP 12 / SNP 13
pfs.all3.c2 <- coxph(PFS.c2~as.numeric(merged.ea.c2$all3_score)+r.cov.c2$distatD+as.numeric(r.cov.c2$age))
summary(pfs.all3.c2)$coefficients[1,5]





###############################################
# TRM covariates: age, bmiOBS, bmiOVWT, PBlood
##############################################

# cohort 1
TRM.c1 <- Surv(time=merged.ea.c1$intxsurv_1Y, event=merged.ea.c1$TRM_1Y==1)
# TRM SNP 8 / SNP 12
trm.snp8snp12.c1 <- coxph(TRM.c1~as.numeric(merged.ea.c1$snp8snp12_score)+as.numeric(r.cov.c1$age)+r.cov.c1$PBlood+r.cov.c1$bmiOBS+r.cov.c1$bmiOVWT)
summary(trm.snp8snp12.c1)$coefficients[1,5]
# TRM SNP 8 / SNP 12 / SNP 13
trm.all3.c1 <- coxph(TRM.c1~as.numeric(merged.ea.c1$all3_score)+as.numeric(r.cov.c1$age)+r.cov.c1$PBlood+as.numeric(r.cov.c1$bmiOBS)+as.numeric(r.cov.c1$bmiOVWT))
summary(trm.all3.c1)$coefficients[1,5]



# cohort 2
TRM.c2 <- Surv(time=merged.ea.c2$intxsurv_1Y, event=merged.ea.c2$TRM_1Y==1)
# TRM SNP 8 / SNP 12
trm.snp8snp12.c2 <- coxph(TRM.c2~as.numeric(merged.ea.c2$snp8snp12_score)+as.numeric(r.cov.c2$age)+r.cov.c2$PBlood+r.cov.c2$bmiOBS+r.cov.c2$bmiOVWT)
summary(trm.snp8snp12.c2)$coefficients[1,5]
# TRM SNP 8 / SNP 12 / SNP 13
trm.all3.c2 <- coxph(TRM.c2~as.numeric(merged.ea.c2$all3_score)+as.numeric(r.cov.c2$age)+r.cov.c2$PBlood+as.numeric(r.cov.c2$bmiOBS)+as.numeric(r.cov.c2$bmiOVWT))
summary(trm.all3.c2)$coefficients[1,5]


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


## Cohort 1

make.tbl <- function(outcome.vector, gene, rsID, chr, BP, N, allele1, allele2, genome, cohort, outcome, disease){
        coef <- outcome.vector[1,1]
        se.coef <- outcome.vector[1,3]
        hr <- outcome.vector[1,2]
        pval <- outcome.vector[1,5]
        res <- c(gene, rsID, chr, BP, N, allele1, allele2, coef, se.coef, hr, NA, pval, "imputed", genome, cohort, outcome, disease)
        res
}

# cohort 1
nod2.c1 <- data.table(rbind(make.tbl(summary(dd.snp8snp12.c1)$coef, "NOD2", "2 SNPs* DD", "chr16", "*", 2033, "*", "*", "S", "c1", "DD", "mixed"),
                            make.tbl(summary(os.snp8snp12.c1)$coef, "NOD2", "2 SNPs* OS", "chr16", "*", 2033, "*", "*", "S", "c1", "OS", "mixed"),
                            make.tbl(summary(pfs.snp8snp12.c1)$coef, "NOD2", "2 SNPs* PFS", "chr16", "*", 2033, "*", "*", "S", "c1", "PFS", "mixed"),
                            make.tbl(summary(trm.snp8snp12.c1)$coef, "NOD2", "2 SNPs* TRM", "chr16", "*", 2033, "*", "*", "S", "c1", "TRM", "mixed"),
                            make.tbl(summary(dd.all3.c1)$coef, "NOD2", "3 SNPs* DD", "chr16", "*", 2033, "*", "*", "S", "c1", "DD", "mixed"),
                            make.tbl(summary(os.all3.c1)$coef, "NOD2", "3 SNPs* OS", "chr16", "*", 2033, "*", "*", "S", "c1", "OS", "mixed"),
                            make.tbl(summary(pfs.all3.c1)$coef, "NOD2", "3 SNPs* PFS", "chr16", "*", 2033, "*", "*", "S", "c1", "PFS", "mixed"),
                            make.tbl(summary(trm.all3.c1)$coef, "NOD2", "3 SNPs* TRM", "chr16", "*", 2033, "*", "*", "S", "c1", "TRM", "mixed")))
colnames(nod2.c1) <- cols

write.table(nod2.c1, file="nod2.c1.txt", col.names=T, row.names=F, quote=F, sep='\t')
## Cohort 2
nod2.c2 <- data.table(rbind(make.tbl(summary(dd.snp8snp12.c2)$coef, "NOD2", "2 SNPs* DD", "chr16", "*", 757, "*", "*", "S", "c2", "DD", "mixed"),
                            make.tbl(summary(os.snp8snp12.c2)$coef, "NOD2", "2 SNPs* OS", "chr16", "*", 757, "*", "*", "S", "c2", "OS", "mixed"),
                            make.tbl(summary(pfs.snp8snp12.c2)$coef, "NOD2", "2 SNPs* PFS", "chr16", "*", 757, "*", "*", "S", "c2", "PFS", "mixed"),
                            make.tbl(summary(trm.snp8snp12.c2)$coef, "NOD2", "2 SNPs* TRM", "chr16", "*", 757, "*", "*", "S", "c2", "TRM", "mixed"),
                            make.tbl(summary(dd.all3.c2)$coef, "NOD2", "3 SNPs* DD", "chr16", "*", 757, "*", "*", "S", "c2", "DD", "mixed"),
                            make.tbl(summary(os.all3.c2)$coef, "NOD2", "3 SNPs* OS", "chr16", "*", 757, "*", "*", "S", "c2", "OS", "mixed"),
                            make.tbl(summary(pfs.all3.c2)$coef, "NOD2", "3 SNPs* PFS", "chr16", "*", 757, "*", "*", "S", "c2", "PFS", "mixed"),
                            make.tbl(summary(trm.all3.c2)$coef, "NOD2", "3 SNPs* TRM", "chr16", "*", 757, "*", "*", "S", "c2", "TRM", "mixed")))
colnames(nod2.c2) <- cols
write.table(nod2.c2, file="nod2.c2.txt", col.names=T, row.names=F, quote=F, sep='\t')


# RUN METAL 

## append meta results
nod2.meta <- fread("metal_S_M_mixed_nod21.tbl")
# remove direction column
nod2.meta <- nod2.meta[,-7]

# build up columns to match how our final table is set up
nod2.meta$gene <- "NOD2"
colnames(nod2.meta)[1] <- "rsID"
nod2.meta$chr <- "chr16"
nod2.meta$exp.coef <- exp(nod2.meta$Effect)
nod2.meta$genome <- "S"
nod2.meta$cohort <- "M"
nod2.meta$disease <- "mixed"
nod2.meta$N <- 2790
nod2.meta$BP <- "*"
nod2.meta$impute <- "imputed"
nod2.meta$`95%-CI` <- NA
nod2.meta$outcome <- NA
colnames(nod2.meta) <- c("rsID", "allele1", "allele2", "coef", "se.coef", "Pvalue", "gene", "chr", "exp.coef", "genome", "cohort", "disease", "N", "BP", "impute", "95%-CI", "outcome")

nod2.meta$outcome[c(1,8)] <- "OS"
nod2.meta$outcome[2:3] <- "DD"
nod2.meta$outcome[c(4,6)] <- "PFS"
nod2.meta$outcome[c(5,7)] <- "TRM"

setcolorder(nod2.meta, cols)

# join all of the nod2 results (c1, c2, meta analyiss)
nod2.final <- do.call(rbind, list(nod2.c1, nod2.c2, nod2.meta))

# remove extra names in rsID column

nod2.final$rsID <- substr(nod2.final$rsID, 1, 7)

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



