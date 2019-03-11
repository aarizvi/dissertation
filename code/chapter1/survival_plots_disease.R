library(survminer)
library(survival)
library(grid)
library(gridExtra)


covs_surv <- covariate.file %>%
    select(sample_type,
           cohort,
           disease,
           intxsurv_1Y,
           intxrel_1Y,
           dead_1Y,
           lfs_1Y,
           rel_1Y,
           disease_death_1Y,
           TRM_1Y,
           GVHD_death_1Y,
           infection_1Y,
           OF_1Y)

covs_surv_long <- covs_surv %>%
    rename(surv_time=intxsurv_1Y,
           rel_time=intxrel_1Y,
           lfs_1y=lfs_1Y,
           rel_1y=rel_1Y) %>%
    mutate_at(.vars=vars(contains("_1Y", ignore.case = FALSE)), .funs=funs(paste0(., " ", surv_time))) %>%
    mutate_at(.vars=vars(contains("_1y", ignore.case = FALSE)), .funs=funs(paste0(., " ", rel_time))) %>%
    select(-surv_time,
           -rel_time) %>%
    gather(key, value, -sample_type, -cohort, -disease) %>% 
    separate(value, c("event", "time"), sep=" ") %>%
    mutate_at(.vars=vars(event, time), .funs=as.double) %>%
    filter(sample_type=="recipient") %>%
    mutate(key=str_replace(key, "dead_1Y", "OS"),
           key=str_replace(key, "disease_death_1Y", "DRM"),
           key=str_replace(key, "TRM_1Y", "TRM"),
           key=str_replace(key, "lfs_1y", "PFS"),
           key=str_replace(key, "rel_1y", "REL"),
           key=str_replace(key, "GVHD_death_1Y", "GVHD"),
           key=str_replace(key, "infection_1Y", "INF"),
           key=str_replace(key, "OF_1Y", "OF")) %>%
    rename(outcome=key)


survivalPlot <- function(outcome_group, title){
    mixed.main <- covs_surv_long %>%
    filter(outcome %in% !!outcome_group)

    m <- survfit(Surv(time=time, event=event) ~ outcome + cohort, 
            data=mixed.main) %>%
        ggsurvplot(fit=.,
                   data=mixed.main,
                   risk.table=FALSE,
                   conf.int=FALSE,
                   title="Mixed Disease Cohorts 1 and 2 Survival Curves") 
    
    m.df <- m$data.survplot %>%
        mutate(disease="Mixed")
    
    amlmds.main <- covs_surv_long %>%
        filter(disease %in%  c("AML", "MDS"),
               outcome %in% !!outcome_group)
    amlonly <- covs_surv_long %>%
        filter(disease %in%  c("AML"),
               outcome %in% !!outcome_group)
    allonly <- covs_surv_long %>%
        filter(disease %in%  c("ALL"),
               outcome %in% !!outcome_group)
    mdsonly <- covs_surv_long %>%
        filter(disease %in%  c("MDS"),
               outcome %in% !!outcome_group)
    
    
    m2 <- survfit(Surv(time=time, event=event) ~ outcome + cohort, 
            data=amlmds.main) %>%
        ggsurvplot(fit=.,
                   data=amlmds.main,
                   risk.table=FALSE,
                   conf.int=FALSE,
                   title="AML + MDS Cohorts 1 and 2 Survival Curves") 
    
    m3 <- survfit(Surv(time=time, event=event) ~ outcome + cohort, 
            data=amlonly) %>%
        ggsurvplot(fit=.,
                   data=amlonly,
                   risk.table=FALSE,
                   conf.int=FALSE,
                   title="AML only Cohorts 1 and 2 Survival Curves")
    
    m4 <- survfit(Surv(time=time, event=event) ~ outcome + cohort, 
            data=allonly) %>%
        ggsurvplot(fit=.,
                   data=allonly,
                   risk.table=FALSE,
                   conf.int=FALSE,
                   title="ALL only Cohorts 1 and 2 Survival Curves")
    
    m5 <- survfit(Surv(time=time, event=event) ~ outcome + cohort, 
            data=mdsonly) %>%
        ggsurvplot(fit=.,
                   data=mdsonly,
                   risk.table=FALSE,
                   conf.int=FALSE,
                   title="AML + MDS Cohorts 1 and 2 Survival Curves")
    m2.df <- m2$data.survplot %>%
        mutate(disease="AML + MDS")
    m3.df <- m3$data.survplot %>%
        mutate(disease="AML")
    m4.df <- m4$data.survplot %>%
        mutate(disease="ALL")
    m5.df <- m5$data.survplot %>%
        mutate(disease="MDS")
    
    
    m.new <- rbind(m.df, m2.df, m3.df, m4.df, m5.df)
    
    m.new <- m.new %>%
        mutate(disease=factor(disease, levels=c("Mixed", "AML + MDS", "AML", "ALL", "MDS")))
    
    p <- m.new %>%
        mutate(cohort=ifelse(cohort==1, "Cohort 1", "Cohort 2")) %>%
        ggplot(aes(time, surv, color=outcome)) +
        geom_line() + 
        facet_grid(cohort~disease) +
        ylim(c(0, 1)) +
        scale_x_continuous(name="Time (months)", breaks=c(0,3, 6, 9, 12)) +
        labs(title=title,
             y="Survival Probability") +
        geom_ribbon(aes(ymin=m.new$surv-(1.96*m.new$std.err),
                        ymax=m.new$surv+(1.96*m.new$std.err)),
                    alpha=0.2,
                    linetype=0) +
        theme_bw() 
    p +
        theme(legend.position = c(0.047, 0.17),
              legend.title=element_text(size=8),
              legend.text=element_text(size=6))
}

s1 <- survivalPlot(c("TRM", "DRM", "OS"), "A. Disease, Transplant and Overall Death by Disease")
s2 <- survivalPlot(c("PFS", "REL"), "B. Progression Free Survival and Relapse by Disease")

grid.arrange(s1, s2, ncol=1)
