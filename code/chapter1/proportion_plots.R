library(tidyverse)
library(kableExtra)
covariate.file <- read_tsv("~/DBMT_PhenoData_EA_long_allVar_20181206.txt")

dz.dist <- covariate.file %>%
    group_by(sample_type, cohort, disease) %>%
    tally() %>%
    ungroup() %>%
    mutate(cohort=ifelse(cohort==1, "cohort 1", "cohort 2")) %>% 
    arrange(cohort) %>% 
    gather(key, value, n) %>%
    unite(sample_cohort, c("sample_type", "cohort")) %>%
    spread(sample_cohort, value) %>%
    select(-key) %>%
    mutate_at(.vars=vars(contains("cohort")), .funs=funs(`per`=./sum(.))) %>%
    mutate(`Cohort 1`=paste0(`donor_cohort 1`, " (", round(`donor_cohort 1_per`, 2)*100, "%)"),
           `Cohort 2`=paste0(`donor_cohort 2`, " (", round(`donor_cohort 2_per`, 2)*100, "%)"),
           `cohort 1`=paste0(`recipient_cohort 1`, " (", round(`recipient_cohort 1_per`, 2)*100, "%)"),
           `cohort 2`=paste0(`recipient_cohort 2`, " (", round(`recipient_cohort 2_per`, 2)*100, "%)")) %>%
    select(-`donor_cohort 1`:-`recipient_cohort 2`,
           -`donor_cohort 1_per`:-`recipient_cohort 2_per`) %>%
    as.data.frame()

dz.dist <- rbind(dz.dist, c("Total",  2052, 763, 2110, 777))

rownames(dz.dist) <- dz.dist$disease
dz.dist$disease <- NULL

dz.dist %>%
    kable(format = "latex",
          caption="\\label{tab:dbmt_cohorts} Donor and Recipients Disease Proportions by Cohort in DISCOVeRY-BMT.", 
          booktabs=TRUE) %>%
    kable_styling(latex_options = "striped",
                  full_width=FALSE,
                  font_size=9) %>%
    add_header_above(c(" ", "Donor (N/Percent %)"=2, "Recipient (N/Percent %)"=2)) %>%
    row_spec(3, hline_after=TRUE) %>%
    footnote(general_title="",
             general="Shown here are disease proportions of ALL, AML and MDS in DISCOVeRY-BMT cohorts. The percentage in each cell is computed column-wise as the proportion of sample size (N) within a cohort across disease group.",
             footnote_as_chunk = TRUE,
             threeparttable = TRUE)
