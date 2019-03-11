genome <- c("Donor", "Recipient", "Mismatch")
outcomes <- c("DRM", "TRM", "OS", "PFS")
patient_subset <- c("Mixed", "AML + MDS", "AML", "ALL", "B-ALL", "T-ALL")
censor <- c("1 year", "3 years", "100 days")
full_df <- setNames(expand.grid(genome, outcomes, patient_subset, censor, stringsAsFactors = FALSE),
                    c("genome", "outcomes", "patient_subset", "cens_time"))
clean_df <- full_df %>% 
    filter(!(genome=="Mismatch" &
           patient_subset %in% c("B-ALL", "T-ALL") & 
           cens_time %in% c("1 year","3 years", "100 days")),
           !(cens_time == "3 years" &
                 patient_subset %in% c("Mixed", "AML", "AML + MDS", "B-ALL", "T-ALL") |
             cens_time == "3 years" & genome == "Mismatch" |
             cens_time == "3 years" & outcomes %in% c("DRM", "TRM", "PFS")),
           !(cens_time == "100 days" & patient_subset %in% c("B-ALL", "T-ALL")),
           !(cens_time == "1 year" & genome == "Mismatch" & patient_subset=="ALL")) %>% 
    group_by(patient_subset, cens_time) %>%
    summarize(genomes=toString(unique(genome)), 
              outcomes=toString(unique(outcomes))) %>%
    ungroup() %>%
    select(genomes, patient_subset, outcomes, cens_time) %>%
    mutate(outcomes=ifelse(cens_time=="100 days", 
                           paste0(outcomes, ", REL, INF, OF, GVHD"),
                           outcomes),
           outcomes=ifelse(patient_subset == "ALL" & cens_time %in% c("1 year"), 
                           paste0(outcomes, ", REL, OF, GVHD, INF"),
                           outcomes),
           outcomes=ifelse(patient_subset == "ALL" & cens_time %in% c("3 years"), 
                           paste0(outcomes, ", REL"),
                           outcomes),
           analyses=str_count(genomes, "\\w+") * str_count(outcomes, "\\w+")) %>%
    arrange(cens_time, desc(genomes)) 

colnames(clean_df) <- c("Genomes", "Patient Subset", "Survival Outcomes", "Censoring Time", "Analyses Run")

rbind(clean_df, c("Total", "", "", "", sum(clean_df$`Analyses Run`))) %>% 
    knitr::kable(caption="\\label{tab:jobs_run} DISCOVeRY-BMT Analyses Run Using DBMT metaPipeline.", booktabs=TRUE) %>%
    kableExtra::kable_styling(latex_options="striped", font_size=9) %>%
    kableExtra::column_spec(1, width="6em") %>%
    kableExtra::column_spec(2, width="4em") %>%
    kableExtra::column_spec(3, width="8em") %>%
    kableExtra::column_spec(4, width="6em") %>%
    kableExtra::column_spec(5, width="4em") %>%
    kableExtra::row_spec(11, hline_after=TRUE)
