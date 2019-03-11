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
               time.to.event <- "intxsurv_100d"
               event <- "dead_100d"
               covariates <- c("age", "DiseaseStatusAdvanced", "PBlood")
           },
           TRM={
               time.to.event <- "intxsurv_100d"
               event <- "trm_100d"
               covariates <- c("age", "bmiOBS", "bmiOVWT", "PBlood")
           },
           PFS={
               time.to.event <- "intxrel_100d"
               event <- "lfs_100d"
               covariates <- c("age", "DiseaseStatusAdvanced")
           },
           DRM={
               time.to.event <- "intxsurv_100d"
               event <- "drm_100d"
               covariates <- c("age", "DiseaseStatusAdvanced")
           },
           REL={
               time.to.event <- "intxrel_100d"
               event <- "rel_100d"
               covariates <- c("CondIntDummy", "TBIDummy")
           },
           INFX = {
               time.to.event <- "intxsurv_100d"
               event <- "inf_100d"
               covariates <- scan(text="age, bmiOBS, bmiOVWT, CMVpn, CMVp, CMVn",
                                  what = "character", sep=",", strip.white=TRUE,
                                  quiet = TRUE)
           },
           INFCMB = {
               time.to.event <- "intxsurv_100d"
               event <- "infcmb_100d"
               covariates <- scan(text="age, bmiOBS, bmiOVWT, CMVpn, CMVp, CMVn",
                                  what = "character", sep=",", strip.white=TRUE,
                                  quiet = TRUE)

           },
           GVHD = {
               time.to.event <- "intxsurv_100d"
               event <- "gvhd_100d"
               covariates <- scan(text="age, dnrage, bmiOBS, bmiOVWT",
                                  what = "character", sep=",", strip.white=TRUE,
                                  quiet = TRUE)
           },
           GVHDCMB = {
               time.to.event <- "intxsurv_100d"
               event <- "gvhdcmb_100d"
               covariates <- scan(text="age, dnrage, bmiOBS, bmiOVWT",
                                  what = "character", sep=",", strip.white=TRUE,
                                  quiet = TRUE)
           },
           OF = {
               time.to.event <- "intxsurv_100d"
               event <- "of_100d"
               covariates <- scan(text="DiseaseStatusAdvanced, PBlood",
                                  what = "character", sep=",", strip.white=TRUE,
                                  quiet = TRUE)
           },
           OFCMB = {
               time.to.event <- "intxsurv_100d"
               event <- "ofcmb_100d"
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
        sample.ids <- sampleFun$sample.ids
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
        sample.ids <- sampleFun$sample.ids
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
