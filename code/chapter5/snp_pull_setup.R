library(tidyverse)
library(glue)
options(tibble.width=Inf)
snp_pull_setup <- function(){
    p1yr0 <<- "/projects/rpci/lsuchest/abbasriz/DBMT_ALLonly/analyses"
    p1yr1 <<- c("ALL_EA_results", "ALLsubtype_Bcell_EA_results")
    p100d <<- "/projects/rpci/lsuchest/abbasriz/DBMT_100d/analyses/ALLonly/out"
    
    ALLonly_dirs <<- dir(c(glue("{p1yr0}/{p1yr1}/out"), p100d), 
                         pattern = "*.res", full.names = TRUE)
    ALLonly_dirs <<- ALLonly_dirs[!str_detect(ALLonly_dirs, "_3y")]
    
    resCols <<- cols_only(
        RSID = col_character(),
        TYPED = col_character(),
        CHR = col_double(),
        POS = col_double(),
        REF = col_character(),
        ALT = col_character(),
        REF.O = col_character(),
        ALT.O = col_character(),
        RefPanelAF = col_double(),
        SAMP_MAF_c1 = col_double(),
        SAMP_MAF_c2 = col_double(),
        PVALUE_M = col_double(),
        HR_M = col_double(),
        HR_lower_M = col_double(),
        HR_upper_M = col_double(),
        PVALUE_c1 = col_double(),
        PVALUE_c2 = col_double(),
        INFO = col_double(),
        # HR_c1 = col_double(),
        # HR_c2 = col_double(),
        # HR_lowerCI_c1 = col_double(),
        # HR_upperCI_c1 = col_double(),
        # HR_lowerCI_c2 = col_double(),
        # HR_upperCI_c2 = col_double(),
        N_c1 = col_double(),
        Nevent_c1 = col_double(),
        N_c2 = col_double(),
        Nevent_c2 = col_double(),
        Direction = col_character(),
        HetISq = col_double(),
        HetChiSq = col_double(),
        HetDf = col_double(),
        HetPVal = col_double()
    )
}
snp_pull_setup()

snp_pull <- function(rsid, chr, gen){
    fpaths <- str_subset(ALLonly_dirs, glue("/{chr}_{gen}_"))
    res_read <- function(fpath, rsid){
        
        if(str_detect(fpath, "DBMT_100d")){
            pt_subset <- str_replace(str_split(fpath, "/",  simplify = TRUE)[1,10],
                                    ".res", "") %>%
                paste0(., "_100d")
        }else{
            pt_subset <- str_replace(str_split(fpath, "/",  simplify = TRUE)[1,10],
                                     ".res", "") %>%
                paste0(., "_1y")
        }
        
        df <- read_tsv(fpath, col_types = resCols, progress = FALSE)
        df %>%
            filter(RSID == rsid) %>%
            mutate(HR_ADJ = ifelse(REF == REF.O & ALT == ALT.O,
                                   HR_M,
                                   (1/HR_M)),
                   HR_lower_ADJ = ifelse(REF == REF.O & ALT == ALT.O,
                                         HR_lower_M,
                                         (1/HR_upper_M)),
                   HR_upper_ADJ = ifelse(REF == REF.O & ALT == ALT.O,
                                         HR_upper_M,
                                         (1/HR_lower_M)),
                   pt_subset = pt_subset
            ) %>%
            unite(c1_N_Nevent, c("N_c1", "Nevent_c1"), sep="/") %>%
            unite(c2_N_Nevent, c("N_c2", "Nevent_c2"), sep="/") %>%
            unite(REF_ALT, c("REF.O", "ALT.O"), sep="/") %>%
            select(RSID,
                   TYPED,
                   CHR,
                   POS,
                   REF_ALT,
                   RefPanelAF, 
                   SAMP_MAF_c1,
                   SAMP_MAF_c2,
                   PVALUE_M,
                   PVALUE_c1,
                   PVALUE_c2,
                   HR_ADJ,
                   HR_lower_ADJ,
                   HR_upper_ADJ,
                   c1_N_Nevent,
                   c2_N_Nevent,
                   INFO,
                   Direction,
                   pt_subset) }
    map_dfr(fpaths, res_read, rsid=rsid) %>%
        arrange(PVALUE_M)
}

system.time(rs17014016 <- snp_pull(rsid = "rs17014016", chr=4, gen="R"))
system.time(rs10022462 <- snp_pull(rsid = "rs10022462", chr=4, gen="R"))
system.time(rs2869930 <- snp_pull(rsid = "rs2869930", chr=4, gen="R"))




system.time(snp_pull(rsid = "rs138933505", chr=6, gen="R"))


system.time(rof1y <- snp_pull(rsid = "rs13362808", chr=6, gen="R"))

system.time(gata3_r <- snp_pull(rsid = "rs3824662", chr=10, gen="R"))
system.time(gata3_d <- snp_pull(rsid = "rs3824662", chr=10, gen="D"))



