library(tidyverse)
library(batch)

parseCommandArgs(evaluate=TRUE)

c1.df <- read_tsv(for_metal_cohort1)
c2.df <- read_tsv(for_metal_cohort2)

c1.df <- c1.df %>%
    mutate(COHORT=1)
c2.df <- c2.df %>%
    mutate(COHORT=2)



full.hrc <- c1.df %>%
    bind_rows(c2.df)


full.hrc <- full.hrc %>%
    unite(cols_joined,
          c(SAMP_FREQ_ALT, SAMP_MAF, PVALUE, HR, HR_lowerCI, HR_upperCI, COEF, SE.COEF, Z, N, N.EVENT),
          sep=";")
full.hrc <- full.hrc %>%
    group_by(COHORT) %>%
    spread(COHORT, cols_joined) %>%
    separate(`1`, c("SAMP_FREQ_ALT_c1",
                    "SAMP_MAF_c1",
                    "PVALUE_c1",
                    "HR_c1",
                    "HR_lowerCI_c1",
                    "HR_upperCI_c1",
                    "COEF_c1",
                    "SE.COEF_c1",
                    "Z_c1",
                    "N_c1",
                    "Nevent_c1"),
             sep=";") %>%
    separate(`2`, c("SAMP_FREQ_ALT_c2",
                    "SAMP_MAF_c2",
                    "PVALUE_c2",
                    "HR_c2",
                    "HR_lowerCI_c2",
                    "HR_upperCI_c2",
                    "COEF_c2",
                    "SE.COEF_c2",
                    "Z_c2",
                    "N_c2",
                    "Nevent_c2"),
             sep=";") %>%
    rename(MarkerName=`RSID;TYPED;CHR;POS;REF;ALT`) %>%
    select(-REF, -ALT)


#meta 1
metal.res <- read_tsv(metal_result)

# merge columns back in
metal.res <- metal.res %>%
    left_join(full.hrc)

metal.res %>% head


metal.res <- metal.res %>%
    separate(MarkerName, c("RSID", "TYPED", "CHR", "POS", "REF", "ALT"), sep=";") %>%
    rename(REF.O=REF,
           ALT.O=ALT,
           REF=Allele1,
           ALT=Allele2,
           COEF_M=Effect,
           SE.COEF_M=StdErr,
           PVALUE_M=`P-value`) %>%
    mutate(HR_M=exp(COEF_M),
           HR_lower_M=exp(COEF_M-1.96*SE.COEF_M),
          HR_upper_M=exp(COEF_M+1.96*SE.COEF_M),
           REF=toupper(REF),
           ALT=toupper(ALT))


metal.res <- metal.res[,c("RSID","TYPED","CHR","POS","REF","ALT","REF.O","ALT.O",
                          "RefPanelAF",
                          "SAMP_MAF_c1","SAMP_MAF_c2",
                          "SAMP_FREQ_ALT_c1", "SAMP_FREQ_ALT_c2",
                          "PVALUE_M","HR_M","HR_lower_M","HR_upper_M",
                          "PVALUE_c1", "PVALUE_c2",
                          "INFO",
                          "HR_c1", "HR_c2", "HR_lowerCI_c1","HR_upperCI_c1", "HR_lowerCI_c2","HR_upperCI_c2",
                          "COEF_M","SE.COEF_M","COEF_c1","SE.COEF_c1","COEF_c2","SE.COEF_c2",
                          "Z_c1","Z_c2", "N_c1", "Nevent_c1", "N_c2","Nevent_c2",
                          "Direction","HetISq","HetChiSq", "HetDf","HetPVal")]

write.table(metal.res, file=full_output, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

