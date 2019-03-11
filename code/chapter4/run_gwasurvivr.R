# run_gwasurvivr.R
library(gwasurvivr)
library(tidyverse)
library(batch)
library(magrittr)

parseCommandArgs(evaluate=TRUE)

runSurv <- function(vcf.file, ptSubset, genome, outcome, cohort, out.file){
    cohort <- as.numeric(cohort)
    switch(outcome,
           OS={
               time.to.event <- "intxsurv_1Y"
               event <- "dead_1Y"
               covariates <- c("age", "DiseaseStatusAdvanced", "PBlood")
               },
           TRM={
               time.to.event <- "intxsurv_1Y"
               event <- "TRM_1Y"
               covariates <- c("age", "bmiOBS", "bmiOVWT", "PBlood")
           },
           PFS={
               time.to.event <- "intxrel_1Y"
               event <- "lfs_1Y"
               covariates <- c("age", "DiseaseStatusAdvanced")
           },
           DRM={
               time.to.event <- "intxsurv_1Y"
               event <- "disease_death_1Y"
               covariates <- c("age", "DiseaseStatusAdvanced")
           },
           REL={
               time.to.event <- "intxrel_1Y"
               event <- "rel_1Y"
               covariates <- c("CondIntDummy", "TBIDummy")
           },
           OS_3y={
               time.to.event <- "intxsurv_3y"
               event <- "dead_3y"
               covariates <- c("age", "DiseaseStatusAdvanced", "PBlood")
           },
           REL_3y={
               time.to.event <- "intxrel_3y"
               event <- "rel_3y"
               covariates <- c("CondIntDummy", "TBIDummy")
           },
           INFX = {
               time.to.event <- "intxsurv_1Y"
               event <- "infection_1Y"
               covariates <- scan(text="age, bmiOBS, bmiOVWT, CMVpn, CMVp, CMVn",
                                  what = "character", sep=",", strip.white=TRUE,
                                  quiet = TRUE)
           },
           INFXCMB = {
               time.to.event <- "intxsurv_1Y"
               event <- "INF1yrCMB"
               covariates <- scan(text="age, bmiOBS, bmiOVWT, CMVpn, CMVp, CMVn",
                                  what = "character", sep=",", strip.white=TRUE,
                                  quiet = TRUE)
               
           },
           GVHD = {
               time.to.event <- "intxsurv_1Y"
               event <- "GVHD_death_1Y"
               covariates <- scan(text="age, dnrage, bmiOBS, bmiOVWT",
                                  what = "character", sep=",", strip.white=TRUE,
                                  quiet = TRUE)
           },
           GVHDCMB = {
               time.to.event <- "intxsurv_1Y"
               event <- "GVHD1yrCMB"
               covariates <- scan(text="age, dnrage, bmiOBS, bmiOVWT",
                                  what = "character", sep=",", strip.white=TRUE,
                                  quiet = TRUE)
           },
           OF = {
               time.to.event <- "intxsurv_1Y"
               event <- "OF_1Y"
               covariates <- scan(text="DiseaseStatusAdvanced, PBlood",
                                  what = "character", sep=",", strip.white=TRUE,
                                  quiet = TRUE)
           },
           OFCMB = {
               time.to.event <- "intxsurv_1Y"
               event <- "OF1yrCMB"
               covariates <- scan(text="DiseaseStatusAdvanced, PBlood",
                                  what = "character", sep=",", strip.white=TRUE,
                                  quiet = TRUE)
           }
           )
    sample.fun <- function(cohort, genome, ptSubset){
        if(genome=="mismatch"){
            covariate.file <- "/projects/rpci/lsuchest/lsuchest/DBMT_PhenoData/mismatch_phenotype_wide_100d_20190223.txt"
            covariate.file <- read_tsv(covariate.file)
            id.column <- "D_R"
            switch(ptSubset,
                   mixed={
                       covariates <- c(covariates, "ALLdummy", "AMLdummy")
                       sample.ids <- covariate.file %>%
                           filter(cohort==!!cohort) %>%
                           pull(D_R)
                       },
                   AMLonly={
                       covariates <- c(covariates, "AMLdummy")
                       sample.ids <- covariate.file %>%
                           filter(cohort==!!cohort,
                                  AMLdummy==1) %>%
                           pull(D_R)
                       },
                   ALLonly={
                       sample.ids <- covariate.file %>%
                           filter(cohort==!!cohort,
                                  ALLdummy==1) %>%
                           pull(D_R)
                       },
                   AMLMDS={
                       sample.ids <- covariate.file %>%
                           filter(cohort==!!cohort,
                                  ALLdummy==0) %>%
                           pull(D_R)
                       }
                   )
        } else if(genome %in% c("donor", "recipient")) {
                covariate.file <- "/projects/rpci/lsuchest/lsuchest/DBMT_PhenoData/DBMT_PhenoData_EA_long_allVar_20190223.txt"
                covariate.file <- read_tsv(covariate.file)
            id.column <- "sangerIDs"
            switch(ptSubset,
                   mixed={
                       sample.ids <- covariate.file %>%
                           filter(cohort==!!cohort,
                                  sample_type==!!genome) %>%
                           pull(sangerIDs)
                       },
                   AMLonly={
                       sample.ids <- covariate.file %>%
                           filter(cohort==!!cohort,
                                  sample_type==!!genome,
                                  AMLdummy==1) %>%
                           pull(sangerIDs)
                       },
                   ALLonly={
                       sample.ids <- covariate.file %>%
                           filter(cohort==!!cohort,
                                  sample_type==!!genome,
                                  ALLdummy==1) %>%
                           pull(sangerIDs)
                       },
                   AMLMDS={
                       sample.ids <- covariate.file %>%
                           filter(cohort==!!cohort,
                                  sample_type==!!genome,
                                  ALLdummy==0) %>%
                           pull(sangerIDs)
                       },
                   ALLsubtype_Bcell={
                       sample.ids <- covariate.file %>%
                           filter(cohort==!!cohort,
                                  sample_type==!!genome,
                                  ALLsubtype_Tcell==1) %>%
                           pull(sangerIDs)
                       },
                   ALLsubtype_Tcell={
                       sample.ids <- covariate.file %>%
                           filter(cohort==!!cohort,
                                  sample_type==!!genome,
                                  ALLsubtype_Bcell==1) %>%
                           pull(sangerIDs)
                       }
                   )
        }
        return(list(sample.ids=sample.ids, covariate.file=covariate.file, id.column=id.column))
    }
    if(cohort==1){
        sampleFun <- sample.fun(cohort=cohort, genome=genome, ptSubset=ptSubset)
        sangerCoxSurv(vcf.file=vcf.file,
                      covariate.file=sampleFun$covariate.file,
                      id.column=sampleFun$id.column,
                      sample.ids = sampleFun$sample.ids,
                      time.to.event=time.to.event,
                      event=event,
                      covariates=covariates,
                      inter.term = NULL,
                      print.covs = "only",
                      out.file=paste0(out.file, "_c1"),
                      maf.filter = 0.005,
                      info.filter = 0.8,
                      chunk.size = 10000,
                      verbose = TRUE,
                      clusterObj = NULL)
    }else if(cohort==2){
        sampleFun <- sample.fun(cohort=cohort, genome=genome, ptSubset=ptSubset)
        sangerCoxSurv(vcf.file=vcf.file,
                      covariate.file=sampleFun$covariate.file,
                      id.column=sampleFun$id.column,
                      sample.ids = sampleFun$sample.ids,
                      time.to.event=time.to.event,
                      event=event,
                      covariates=covariates,
                      inter.term = NULL,
                      print.covs = "only",
                      out.file=paste0(out.file, "_c2"),
                      maf.filter = 0.005,
                      info.filter = 0.8,
                      chunk.size = 10000,
                      verbose = TRUE,
                      clusterObj = NULL)
    }
}

options("gwasurvivr.cores"=as.numeric(ncores))

runSurv(vcf.file,
        ptSubset,
        genome,
        outcome,
        cohort,
        out.file)