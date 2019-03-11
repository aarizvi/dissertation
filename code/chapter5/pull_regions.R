library(tidyverse)
library(glue)
chr = 4
rsid = "rs10516806"
outcome = "R_ALLonly_DRM"
input_path = "/projects/rpci/lsuchest/abbasriz/DBMT_ALLonly/analyses/ALL_EA_results"
snipa_dir = "/projects/rpci/lsuchest/abbasriz/DBMT_ALLonly/snipa"
# input_path = "/projects/rpci/lsuchest/abbasriz/DBMT_100d/analyses/ALLonly"
# snipa_dir = "/projects/rpci/lsuchest/abbasriz/DBMT_100d/snipa"

res <- read_tsv(glue("{input_path}/out/{chr}_{outcome}.res"))

res %>%
    filter(RSID == rsid) %>%
    pull(POS) -> pos

res %>%
    filter(
        between(POS, pos - 150e3, pos + 150e3),
        HetPVal > 0.05,
        HetISq < 50,
        SAMP_MAF_c1 > 0.01,
        SAMP_MAF_c2 > 0.01,
        INFO > 0.8
   ) -> res

write_tsv(res,
          path = glue("{snipa_dir}/out/{rsid}_{chr}_{outcome}.regionres"))

write_tsv(
    select(res, RSID, PVALUE_M),
    path = glue("{snipa_dir}/out/{rsid}_{chr}_{outcome}.snipa"),
    col_names = FALSE

