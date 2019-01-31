library(tidyverse)

covariate.file <- read_tsv("~/Google Drive/OSU_PHD/DBMT_100d/DBMT_PhenoData_EA_long_allVar_20181216.txt") %>%
    filter(sample_type == "recipient")


covs_surv <- covariate.file %>%
    select(sangerIDs,
           sample_type,
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
    gather(key, value, -sangerIDs,-sample_type, -cohort, -disease) %>% 
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

mixed <- covs_surv_long %>% 
    mutate(Disease="Mixed")

amlmds <- covs_surv_long %>% 
    filter(disease %in% c("AML", "MDS")) %>%
    mutate(Disease="AML + MDS")
 
amlonly <- covs_surv_long %>% 
    filter(disease %in% c("AML")) %>%
    mutate(Disease="AML only")   

allonly <- covs_surv_long %>% 
    filter(disease %in% c("ALL")) %>%
    mutate(Disease="ALL only")   

covs_surv_long <- rbind(mixed, amlmds, amlonly, allonly)
    
covs <- covariate.file %>%
    select(sangerIDs, distatD, PBlood, CondInt, bmi_cat)

covs <- covs %>% 
    right_join(covs_surv_long) %>%
    mutate(death_relapse=ifelse(outcome %in% c("REL", "PFS"), "relapse", "death"))

covs <- covs %>%
    mutate(Disease=factor(Disease, levels=c("Mixed", "AML + MDS", "AML only", "ALL only")))


distPlot <- function(dat, outcome, variable, death_relapse,  legend_title){
    if(death_relapse=="death"){
        dat <- dat %>%
            mutate(cohort=ifelse(cohort==1, "Cohort 1", "Cohort 2"),
                   PBlood=ifelse(PBlood==1, "Peripheral Blood", "Bone Marrow"))
        switch (outcome,
                OS = {dat <- dat %>% filter(outcome == "OS", event==1)},
                TRM = {dat <- dat %>% filter(outcome == "TRM", event==1)},
                DRM = {dat <- dat %>% filter(outcome == "DRM", event==1)})
        variable <- enquo(variable)
        ggplot(dat, aes(x=time, fill=as.factor(!!variable))) +   
            geom_histogram(alpha = 0.6,  binwidth = 0.25) +
            facet_grid(cohort~Disease) +
            guides(fill=guide_legend(title=legend_title)) +
            ggtitle(outcome) +
            theme_bw() +
            ylab("Number of Events") +
            xlab("Time to death after transplant (Months)") +
            #theme(legend.position = c(0.911, 0.87))
            theme(legend.position="none")
    } else if(death_relapse=="relapse"){
        dat <- dat %>%
            mutate(cohort=ifelse(cohort==1, "Cohort 1", "Cohort 2"),
                   PBlood=ifelse(PBlood==1, "Peripheral Blood", "Bone Marrow"))
        switch (outcome,
                PFS = {dat <- dat %>% filter(outcome == "PFS", event==1)},
                REL = {dat <- dat %>% filter(outcome == "REL", event==1)})
        variable <- enquo(variable)
        ggplot(dat, aes(x=time, fill=as.factor(!!variable))) +   
            geom_histogram(alpha = 0.6,  binwidth = 0.25) +
            facet_grid(cohort~Disease) +
            guides(fill=guide_legend(title=legend_title)) +
            ggtitle(outcome) +
            theme_bw() +
            ylab("Number of Events") +
            xlab("Time to relapse after transplant (Months)") +
            #theme(legend.position = c(0.911, 0.87))
            theme(legend.position="none")
    }
}

d1 <- distPlot(covs, "DRM", distatD, "death",  "Disease Status")
d2 <- distPlot(covs, "TRM", distatD, "death", "Disease Status")
d3 <- distPlot(covs, "OS", distatD, "death", "Disease Status")
d4 <- distPlot(covs, "PFS", distatD, "relapse", "Disease Status")
d5 <- distPlot(covs, "REL", distatD, "relapse",  "Disease Status")

prow <- plot_grid( d1, d2, d3, d4, d5,
                   align = 'vh',
                   labels = c("A", "B", "C", "D", "E"),
                   hjust = -1,
                   nrow = 3,
                   ncol= 2
)

legend_b <- get_legend(d1 + theme(legend.position="bottom"))

p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(5, .1))

p

library(grid)
library(gridExtra)
library(cowplot)

grid.arrange(distPlot(covs, "DRM", distatD, "death", "A. Disease Status and 1-Year Disease Related Mortality", "Disease Status"),
             distPlot(covs, "TRM", distatD, "death", "B. Disease Status and 1-Year Transplant Related Mortality", "Disease Status"),
             distPlot(covs, "OS", distatD, "death", "C. Disease Status and 1-Year Overall Survival", "Disease Status"),
             distPlot(covs, "PFS", distatD, "relapse", "D. Disease Status and 1-Year Progression-Free Survival", "Disease Status"),
             distPlot(covs, "REL", distatD, "relapse", "E. Disease Status and 1-Year Relapse", "Disease Status"),
             ncol=2)

distPlot(covs, "DRM", PBlood, "death", "Graft Source and 1-Year Disease Related Mortality", "Graft Source")
distPlot(covs, "TRM", PBlood, "death", "Graft Source and 1-Year Transplant Related Mortality", "Graft Source")
distPlot(covs, "OS", PBlood, "death", "Graft Source and 1-Year Overall Survival", "Graft Source")
distPlot(covs, "PFS", PBlood, "relapse", "Graft Source and 1-Year Progression-Free Survival", "Graft Source")
distPlot(covs, "REL", PBlood, "relapse", "Graft Source and 1-Year Relapse", "Graft Source")


distPlot(covs, "DRM", CondInt, "death", "Conditioning Intensity and 1-Year Disease Related Mortality", "Conditioning Intensity")
distPlot(covs, "TRM", CondInt, "death", "Conditioning Intensity and 1-Year Transplant Related Mortality", "Conditioning Intensity")
distPlot(covs, "OS", CondInt, "death", "Conditioning Intensity and 1-Year Overall Survival", "Conditioning Intensity")
distPlot(covs, "PFS", CondInt, "relapse", "Conditioning Intensity and 1-Year Progression-Free Survival", "Conditioning Intensity")
distPlot(covs, "REL", CondInt, "relapse", "Conditioning Intensity and 1-Year Relapse", "Conditioning Intensity")

distPlot(covs, "DRM", bmi_cat, "death", "Body Mass Index (BMI) and 1-Year Disease Related Mortality", "Body Mass Index (BMI)")
distPlot(covs, "TRM", bmi_cat, "death", "Body Mass Index (BMI) and 1-Year Transplant Related Mortality", "Body Mass Index (BMI)")
distPlot(covs, "OS", bmi_cat, "death", "Body Mass Index (BMI) and 1-Year Overall Survival", "Body Mass Index (BMI)")
distPlot(covs, "PFS", bmi_cat, "relapse", "Body Mass Index (BMI) and 1-Year Progression-Free Survival", "Body Mass Index (BMI)")
distPlot(covs, "REL", bmi_cat, "relapse", "Body Mass Index (BMI) and 1-Year Relapse", "Body Mass Index (BMI)")

distPlot(covs, "DRM", distatD, "death", "Disease Status and 1-Year Disease Related Mortality", "Disease Status")
distPlot(covs, "TRM", distatD, "death", "Disease Status and 1-Year Transplant Related Mortality", "Disease Status")
distPlot(covs, "OS", distatD, "death", "Disease Status and 1-Year Overall Survival", "Disease Status")
distPlot(covs, "PFS", distatD, "relapse", "Disease Status and 1-Year Progression-Free Survival", "Disease Status")
distPlot(covs, "REL", distatD, "relapse", "Disease Status and 1-Year Relapse", "Disease Status")

