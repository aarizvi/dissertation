library(tidyverse)
library(glue)
options(tibble.width = Inf)

# setwd("/projects/rpci/lsuchest/abbasriz/DBMT_100d/manifest")
# man <- read.table("100d_MM_INF.manifest", header=T, stringsAsFactors=F)
# man[1,]
# table(man$genome)

path1dr <- "/projects/rpci/lsuchest/lsuchest/Rserve/BMT/genetic_data/PLINK2VCF/Sanger_HRC/"
path1mm <- "/projects/rpci/lsuchest/lsuchest/Rserve/BMT/genetic_data/PLINK2VCF/Sanger_HRC/mismatch_genome/"
path2 <- ".pbwt_reference_impute.vcf.gz"
chr <- 1:22
genome <- c("donor", "recipient", "mismatch")
mem <- 24000
#outcome <- "INFX"
#ptsub <- "ALLonly"


manifest_maker <- function(man_dir, ptsub, outcome){

  inf_man <- expand.grid(path2=path2, chr=chr,
                         ptsub=ptsub, genome=genome, mem=mem,
                         outcome=outcome)
  inf_man %>%
    mutate(out1=chr,
           out2=case_when(genome == "donor" ~ "D",
                          genome == "recipient" ~ "R",
                          genome == "mismatch" ~ "MM"
           ),
           out3=ptsub,
           out4=outcome,
           path1=case_when(genome == "donor" ~ path1dr,
                           genome == "recipient" ~ path1dr,
                           genome == "mismatch" ~ path1mm
           )
    ) %>%
    unite(path_to_vcf_file, path1, chr, path2, sep = "") %>%
    unite(output_name, out1:out4) %>%
    select(path_to_vcf_file, output_name,
           patient_subset = ptsub, genome, outcome, memory=mem) %>%
    write.table(file = glue("{man_dir}/RDMM_{ptsub}_{outcome}.manifest"),
                row.names = FALSE,
                quote = FALSE,
                sep="\t")

}

#### ALLonly ####
manifest_maker(man_dir="/projects/rpci/lsuchest/abbasriz/DBMT_100d/analyses/ALLonly",
               ptsub="ALLonly",
               outcome="INFX")

manifest_maker(man_dir="/projects/rpci/lsuchest/abbasriz/DBMT_100d/analyses/ALLonly",
               ptsub="ALLonly",
               outcome="REL")

manifest_maker(man_dir="/projects/rpci/lsuchest/abbasriz/DBMT_100d/analyses/ALLonly",
               ptsub="ALLonly",
               outcome="GVHD")
#### mixed ####
manifest_maker(man_dir="/projects/rpci/lsuchest/abbasriz/DBMT_100d/analyses/mixed",
               ptsub="mixed",
               outcome="INFX")

manifest_maker(man_dir="/projects/rpci/lsuchest/abbasriz/DBMT_100d/analyses/mixed",
               ptsub="mixed",
               outcome="REL")

manifest_maker(man_dir="/projects/rpci/lsuchest/abbasriz/DBMT_100d/analyses/mixed",
               ptsub="mixed",
               outcome="GVHD")
#### AMLonly ####

manifest_maker(man_dir="/projects/rpci/lsuchest/abbasriz/DBMT_100d/analyses/AMLonly",
               ptsub="AMLonly",
               outcome="INFX")

manifest_maker(man_dir="/projects/rpci/lsuchest/abbasriz/DBMT_100d/analyses/AMLonly",
               ptsub="AMLonly",
               outcome="REL")

manifest_maker(man_dir="/projects/rpci/lsuchest/abbasriz/DBMT_100d/analyses/AMLonly",
               ptsub="AMLonly",
               outcome="GVHD")

#### AMLMDS ####
manifest_maker(man_dir="/projects/rpci/lsuchest/abbasriz/DBMT_100d/analyses/AMLMDS",
               ptsub="AMLMDS",
               outcome="INFX")

manifest_maker(man_dir="/projects/rpci/lsuchest/abbasriz/DBMT_100d/analyses/AMLMDS",
               ptsub="AMLMDS",
               outcome="REL")

manifest_maker(man_dir="/projects/rpci/lsuchest/abbasriz/DBMT_100d/analyses/AMLMDS",
               ptsub="AMLMDS",
               outcome="GVHD")

