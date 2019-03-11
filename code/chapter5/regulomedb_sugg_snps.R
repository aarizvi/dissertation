# read in all the data

library(tidyverse)
regdb <- data.table::fread("/projects/rpci/lsuchest/ezgi/regulomeDb/RegulomeDB.dbSNP141.lessCol.txt") 

regdb <- regdb %>%
    mutate(V1=as.double(str_replace_all(V1, "chr", "")))
colnames(regdb) <- c("CHR", "POS", "RSID", "regdb_score") 


library(tidyverse)

path <- "/projects/rpci/lsuchest/abbasriz/ALL_sugg_snps"
files <- dir(path=path, pattern=".rds", full.names=TRUE)

read_and_convert <- function(path){
    df <- readRDS(path)
    df <- df %>%
        mutate_at(vars(RefPanelAF,
                       SAMP_MAF_c1, 
                       SAMP_MAF_c2, 
                       PVALUE_M, 
                       PVALUE_c1,
                       PVALUE_c2), 
                  as.double) 
    
    pt_subset2 <- df %>% 
        pull(pt_subset) 
    
    pt_subset2 <-  strsplit(pt_subset2,"/")
    pt_subset2 <- unlist(lapply(pt_subset2, tail , 1 ))
    pt_subset2 <- str_replace_all(pt_subset2, ".res", "")
    df <- df %>%
        mutate(pt_subset=pt_subset2) 
    return(df)
}

res <- files %>%
    map_df(read_and_convert) 

res <- res %>%
    left_join({regdb %>% select(-CHR, -POS)})

rm(regdb)

saveRDS(res, file="20190306_suggestive_sig_snps_regdb.rds")