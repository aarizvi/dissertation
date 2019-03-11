library(gridExtra)
library(grid)
mixed.age <- covariate.file %>%
    select(disease, sample_type, cohort, age, dnrage) %>% 
    gather(key, value, -disease, -sample_type, -cohort) %>%
    mutate(key=ifelse(key=="age", "recipient age", "donor age"),
           cohort=ifelse(cohort==1, "cohort 1", "cohort 2")) %>%
    filter(disease %in% c("AML", "ALL" ,"MDS")) %>%
    group_by(sample_type, cohort) %>%
    ggplot(aes(value)) +
    geom_histogram() + 
    facet_grid(cohort~key, scales = "free_y") +
    ggtitle("A. Mixed Disease") +
    xlab("age")

amlmds.age <- covariate.file %>%
    select(disease, sample_type, cohort, age, dnrage) %>% 
    gather(key, value, -disease, -sample_type, -cohort) %>%
    mutate(key=ifelse(key=="age", "recipient age", "donor age"),
           cohort=ifelse(cohort==1, "cohort 1", "cohort 2")) %>%
    filter(disease %in% c("AML", "MDS")) %>%
    group_by(sample_type, cohort) %>%
    ggplot(aes(value)) +
    geom_histogram() + 
    facet_grid(cohort~key, scales = "free_y") +
    ggtitle("B. AML + MDS") +
    xlab("age") 

aml.age <- covariate.file %>%
    select(disease, sample_type, cohort, age, dnrage) %>% 
    gather(key, value, -disease, -sample_type, -cohort) %>%
    mutate(key=ifelse(key=="age", "recipient age", "donor age"),
           cohort=ifelse(cohort==1, "cohort 1", "cohort 2")) %>%
    filter(disease %in% c("AML")) %>%
    group_by(sample_type, cohort) %>%
    ggplot(aes(value)) +
    geom_histogram() + 
    facet_grid(cohort~key,  scales = "free_y") +
    ggtitle("C. AML Only") +
    xlab("age")

all.age <- covariate.file %>%
    select(disease, sample_type, cohort, age, dnrage) %>% 
    gather(key, value, -disease, -sample_type, -cohort) %>%
    mutate(key=ifelse(key=="age", "recipient age", "donor age"),
           cohort=ifelse(cohort==1, "cohort 1", "cohort 2")) %>%
    filter(disease %in% c("ALL")) %>%
    group_by(sample_type, cohort) %>%
    ggplot(aes(value)) +
    geom_histogram() + 
    facet_grid(cohort~key, scales = "free_y") +
    ggtitle("D. ALL Only") +
    xlab("age") 

mds.age <- covariate.file %>%
    select(disease, sample_type, cohort, age, dnrage) %>% 
    gather(key, value, -disease, -sample_type, -cohort) %>%
    mutate(key=ifelse(key=="age", "recipient age", "donor age"),
           cohort=ifelse(cohort==1, "cohort 1", "cohort 2")) %>%
    filter(disease %in% c("MDS")) %>%
    group_by(sample_type, cohort) %>%
    ggplot(aes(value)) +
    geom_histogram() + 
    facet_grid(cohort~key,  scales = "free_y") +
    ggtitle("MDS Patients: Donor and Recipient Age Distribution") +
    xlab("age") +
    theme(plot.title = element_text(hjust = 0.5))

grid.arrange(mixed.age,
             amlmds.age,
             aml.age,
             all.age,
             ncol=2)
