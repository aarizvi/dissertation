outcome_distn <- function(disease_group){
    total <- covariate.file %>%
        group_by(sample_type, cohort, disease) %>%
        tally() %>%
        ungroup() %>%
        mutate(cohort=ifelse(cohort==1, "cohort 1", "cohort 2")) %>% 
        arrange(cohort) %>% 
        gather(key, value, n) %>%
        filter(disease %in% !!disease_group,
               sample_type=="recipient") %>% 
        group_by(cohort) %>%
        summarize(n=sum(value)) %>%
        pull(n)
    data <- covariate.file %>%
        filter(disease %in% !!disease_group,
               sample_type == "recipient") %>% 
        group_by(cohort) %>%
        summarize(OS=sum(dead_1Y),
                  DRM=sum(disease_death_1Y),
                  TRM=sum(TRM_1Y),
                  PFS=sum(lfs_1Y),
                  REL=sum(rel_1Y, na.rm=TRUE),
                  GVHD=sum(GVHD_death_1Y),
                  OF=sum(OF_1Y),
                  INF=sum(infection_1Y)) %>% 
        arrange(cohort) %>% 
        mutate(total=total) %>%
        mutate_at(.vars=vars(OS, DRM, PFS, REL, TRM, GVHD, OF, INF), .funs=funs(paste0(., " (", round((./total)*100, 2), "%)"))) %>% 
        select(-total) %>% 
        gather(key, value, -cohort) %>% 
        spread(cohort, value) %>% 
        rename(`Cohort 1`=`1`,
               `Cohort 2`=`2`) %>% 
        mutate(key=factor(key, levels=c("DRM", "TRM", "OS", "REL", "PFS", "GVHD", "INF", "OF"))) %>%
        arrange(key) %>%
        as.data.frame()
        levels(data$key) <- c("Disease related mortality (DRM)",
                              "Transplant related mortality (TRM)",
                              "Overall survival (OS)",
                              "Relapse (REL)",
                              "Progression free survival (PFS)",
                              "Graft-versus-host disease (GVHD)",
                              "Infection (INF)",
                              "Organ failure (OF)")
        data %>%
            rename(Outcome=key)
}

outcome_distributions <- rbind(
    outcome_distn(c("ALL", "AML", "MDS")),
    outcome_distn(c("AML", "MDS")),
    outcome_distn("AML"),
    outcome_distn("ALL"))

outcome_distributions %>% 
    kable(format = "latex",
          caption="\\label{tab:dbmt_props} Proportion of Events by Survival Outcome", 
          booktabs=TRUE) %>%
    kable_styling(latex_options = "striped",
                  full_width=FALSE, font_size=9) %>%
    add_header_above(c(" ", "Recipient (N/Percent %)"=2)) %>%
    group_rows("Mixed Disease", 1, 8, latex_gap_space=".5em") %>%
    group_rows("AML + MDS ", 9, 16, latex_gap_space=".5em") %>%
    group_rows("AML Only", 17, 24, latex_gap_space=".5em") %>%
    group_rows("ALL Only", 25, 32, latex_gap_space=".5em") 
