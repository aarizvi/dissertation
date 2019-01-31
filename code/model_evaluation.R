covariate.file <- read_tsv("~/Google Drive/OSU_PHD/DBMT_100d/DBMT_PhenoData_EA_long_allVar_20181216.txt")
# cohort 1
c1 <- covariate.file %>%
    filter(sample_type=="recipient",
           cohort==1) %>%
    mutate(age_cat = ifelse(age >= 40, 1, 0),
           intxsurv_1Y_days=intxsurv_1Y*30.4167) # convert to days
library(survival)
library(survminer)

c2 <- covariate.file %>%
    filter(sample_type=="recipient",
           cohort==2) %>%
    mutate(age_cat = ifelse(age >= 40, 1, 0),
           intxsurv_1Y_days=intxsurv_1Y*30.4167) 

time <- "intxsurv_1Y"
event <- "TRM_1Y"
covariates <- c("DiseaseStatusAdvanced", "age_cat", "PBlood")


data <- c1




mod <- coxph(Surv(intxsurv_1Y, dead_1Y) ~ DiseaseStatusAdvanced + PBlood + age, data=c1)
pred.mod <- predict(mod, c1, type="expected")

# c1 <- c1 %>%
#     mutate(predicted=pred.mod)

eg <- expand.grid(intxsurv_1Y=seq(0, 12, by=0.01),
                  age=seq(0, max(c1$age), by=1))

#time.death <- c1 %>% select(intxsurv_1Y, dead_1Y)

covariate.file <- read_tsv("~/Google Drive/OSU_PHD/DBMT_100d/DBMT_PhenoData_EA_long_allVar_20181216.txt")
# cohort 1
c1 <- covariate.file %>%
    filter(sample_type=="recipient",
           cohort==1)

eg <- rbind(eg, eg, eg, eg) %>%
    mutate(PBlood=rep(c(0,1), each=nrow(eg)*2),
           DiseaseStatusAdvanced = rep(c(0,1,0,1), each=nrow(eg)),
           dead_1Y=ifelse(intxsurv_1Y==12, 0, 1))
pred.mod <- predict(mod, eg, type="expected") 



eg.pred <- eg %>% mutate(`Survival Probability`= exp(-pred.mod))

c1 <- c1 %>%
    mutate(DiseaseStatusAdvanced=ifelse(DiseaseStatusAdvanced==0, "Advanced", "Early"),
           PBlood=ifelse(PBlood==0, "Bone Marrow", "Peripheral Blood"))


eg.pred %>%
    mutate(DiseaseStatusAdvanced=ifelse(DiseaseStatusAdvanced==0, "Advanced", "Early"),
           PBlood=ifelse(PBlood==0, "Bone Marrow", "Peripheral Blood")) %>%
    ggplot(aes(x=intxsurv_1Y, y=age, fill=`Survival Probability`)) +
    facet_grid(DiseaseStatusAdvanced ~ PBlood) +
    scale_fill_gradient2(high="green",
                         mid= "red") +
    geom_tile() +
    geom_point(data=c1, fill="red") +
    scale_x_continuous(breaks=c(0, 3, 6, 9, 12))


?sample
c1.sl <- c1 %>% 
    select(intxsurv_1Y, dead_1Y, distatD, age, PBlood) %>%
    mutate(bin_age=mround(age, 20))


uniq_cov_combo <- c1.sl %>%
    select(-intxsurv_1Y, -dead_1Y, -age) %>%
    distinct() %>%
    unite(cov_combo, c(distatD, PBlood, bin_age), sep=" + ", remove=FALSE) %>% 
    arrange(bin_age)

boot_n <- 10
boot_results <- vector('list', length(boot_n))

for(boot in seq_len(length(boot_n))){
    
    for(i in seq_len(nrow(uniq_cov_combo))){
        
        x <- uniq_cov_combo %>% 
            select(-cov_combo) %>%
            slice(i) %>%
            unlist()
        
        
        survFunTotal <- exp(-cumhaz * exp(sum(x*betas)))
        sim_res[[i]] <- data.frame(sim_surv_time=tdom[colSums(outer(survFunTotal, u, `>`))])
        
        
        idx <- sample(1:nrow(c1.sl), size=nrow(c1.sl), replace=TRUE)
        c1_samp <- c1.sl[idx,]
        
        
           
    }
}





# 2d summary of a complex multidimensional space


### bootstrap at each time bin calculate a quantile .... 
# use replicate -- generate bunch a vector ... and have an apply functino that computes quantile by row
nsamps <- 2110


n_boot <- 100
boot_res <- vector('list', n_boot)    
    
for(boot in seq_len(n_boot)){
    mod <- coxph(formula(paste0("Surv(", time, ", ", event,") ~ ", paste(covariates, collapse=" + "))), data=data.frame(data))
    tdom <- basehaz(mod)$time
    cumhaz <- basehaz(mod)$hazard
    # survival probability 
    survFunTime <- exp(-cumhaz)
    u <- runif(2110)
    failtimes <- tdom[colSums(outer(survFunTime, u, `>`))]
    failtimes <- data.frame(failtimes)
    ## baseline
    # now grab covariates
    betas <- coef(mod)
    # simulation 
    sim_res <- vector("list", nrow(y))
    # total_fit <- vector('list', nrow(y))
    for (i in seq_len(nrow(y))) {
        x <- y %>% 
            select(-covariates) %>%
            slice(i) %>%
            unlist()
        survFunTotal <- exp(-cumhaz * exp(sum(x*betas)))
        sim_res[[i]] <- data.frame(sim_surv_time=tdom[colSums(outer(survFunTotal, u, `>`))])
        #total_fit[[i]] <- ggsurvplot(surv_fit(Surv(res) ~ 1, data=sim_res[[i]]))$data.survplot
    }
    # total_fit <- mapply(function(x, cat) x %>% mutate(covariates=cat), total_fit, y$covariates, SIMPLIFY = FALSE)
    # total_fit <- do.call('rbind', total_fit)
    
    sim_res <- mapply(function(x, cat, run) x %>% mutate(covariates=cat, run=run), sim_res, y$covariates, boot,SIMPLIFY = FALSE)
    
    boot_res[[boot]] <- sim_res
    
}

test <- lapply(boot_res, function(x) do.call("rbind", x))
test2 <- do.call("rbind", test)

quantiles_surv <- test2 %>%
    mutate(sim_surv_time=sim_surv_time*30.5, bin_time=round(sim_surv_time)) %>%
    group_by(bin_time, covariates, run) %>% 
    summarize(sum=length(bin_time)) %>% 
    group_by(bin_time, covariates) %>%
    summarize(q=quantile(sum, 0.975)) %>%
    ungroup()

c1_new <- c1 %>%
    mutate(surv_time=intxsurv_1Y*30.5, bin_time=round(surv_time)) %>%
    group_by(bin_time, DiseaseStatusAdvanced, age_cat, PBlood) %>%
    summarize(sum=length(bin_time)) %>%
    left_join(quantiles_surv) %>%
    mutate(qsum=(q-sum)>0) %>%
    ungroup()

mean(c1_new$qsum, na.rm = T)


c1_new %>%
    unite(covs, c(DiseaseStatusAdvanced, age_cat, PBlood), sep = " + ") %>%
    na.omit %>%
    ggplot(aes(x=bin_time, y=q-sum)) +
    geom_col(color="red", position="dodge") +
    facet_wrap(~covs) +
    coord_cartesian(xlim=c(0,360), ylim=c(-5,1800))
    


sim_res <- do.call("rbind", sim_res)
sim_res <- sim_res %>%
    mutate(covariates=ifelse(covariates==1, "no covariates", covariates),
           covariates=str_replace_all(covariates, "\\+", "+ \n")) 
eq <- paste0(time, " | ", event, " ~ \n ", paste(covariates, collapse=" + \n"))
real_data_res <- data.frame(res=data[[time]], 
                            covariates=eq) 






simFun <- function(time, event, covariates, data, nsamps){
    y <- data %>%
        select_at(vars(covariates)) %>% 
        distinct() %>%
        mutate(covariates=apply(.==1,1,function(a) {paste0(colnames(.)[a], collapse = " + ")}),
               covariates=ifelse(covariates=="", "1", covariates)) %>% 
        arrange(covariates)
    #### get real data stats ###
    models <- vector('list', nrow(y))
    for(i in seq_len(nrow(y))){
        covs <- y$covariates[i]
        models[[i]] <- ggsurvplot(surv_fit(
            formula(paste0("Surv(", time, ", ", event,") ~ ", paste(covs, collapse=" + "))),
            data=data.frame(data)
        )
        )$data.survplot
    }
    models[[1]] <- models[[1]] %>% mutate(strata=1)
    models <- mapply(function(x, cat) x %>% select_at(vars(time:strata)) %>% mutate(covariates=cat), models, y$covariates, SIMPLIFY = FALSE)
    models <- do.call("rbind", models)
    models <- models %>%
        mutate(type="observed")
    models <- models %>%
        mutate(type="observed") %>%
        filter(!str_detect(strata, pattern=regex("[0]"))) %>%
        select(-strata)
    #### for simulated data 
    mod <- coxph(formula(paste0("Surv(", time, ", ", event,") ~ ", paste(covariates, collapse=" + "))), data=data.frame(data))
    tdom <- basehaz(mod)$time
    cumhaz <- basehaz(mod)$hazard
    # survival probability 
    survFunTime <- exp(-cumhaz)
    survdf <- data.frame(idx=1:length(tdom), 
                         tdom=tdom, 
                         cumhaz=cumhaz, 
                         survFunTime=survFunTime)
    # hazard function and cumulative hazard function in cohort 1
    hazFuncumFunPlot <- survdf %>% 
        gather(key, value, -idx, -tdom) %>%
        mutate(key=ifelse(key=="cumhaz", "Cumulative Hazard Function", "Survival Function")) %>%
        ggplot(aes(tdom, value)) +
        geom_point() +
        facet_wrap(~key) +
        scale_x_continuous(breaks=c(0, 3, 6, 9, 12)) + 
        xlab("Months") +
        ylab("Probability")
    # random samples
    u <- runif(nsamps)
    failtimes <- tdom[colSums(outer(survFunTime, u, `>`))]
    failtimes <- data.frame(failtimes)
    ## baseline
    # now grab covariates
    betas <- coef(mod)
    # simulation 
    sim_res <- vector("list", nrow(y))
    total_fit <- vector('list', nrow(y))
    for (i in seq_len(nrow(y))) {
        x <- y %>% 
            select(-covariates) %>%
            slice(i) %>%
            unlist()
        survFunTotal <- exp(-cumhaz * exp(sum(x*betas)))
        sim_res[[i]] <- data.frame(res=tdom[colSums(outer(survFunTotal, u, `>=`))])
        total_fit[[i]] <- ggsurvplot(surv_fit(Surv(res) ~ 1, data=sim_res[[i]]))$data.survplot
    }
    total_fit <- mapply(function(x, cat) x %>% mutate(covariates=cat), total_fit, y$covariates, SIMPLIFY = FALSE)
    total_fit <- do.call('rbind', total_fit)
    total_fit <- total_fit %>% 
        mutate(type="predicted")
    # baseline_surv <- total_fit %>% 
    #     filter(covariates=="no covariates")
    res <- rbind(models, total_fit)
    res <- res %>%
        filter(surv!=0) %>%
        mutate(covariates=ifelse(covariates==1, "baseline", covariates),
               covariates=str_replace_all(covariates, "\\+", "+ \n")) %>%
        ggplot(aes(time, surv, color=type)) +
        geom_line() +
        facet_wrap(~covariates) +
        ylim(c(0,1)) + 
        scale_x_continuous(breaks=c(0, 3, 6, 9, 12)) +
        xlab('Months') +
        ylab('Survival Probability') +
        theme(legend.position=c(1, 0.2),
              legend.justification="right") 
    # histograms 
    sim_res <- mapply(function(x, cat) x %>% mutate(covariates=cat), sim_res, y$covariates, SIMPLIFY = FALSE)
    sim_res <- do.call("rbind", sim_res)
    sim_res <- sim_res %>%
        mutate(covariates=ifelse(covariates==1, "no covariates", covariates),
               covariates=str_replace_all(covariates, "\\+", "+ \n")) 
    eq <- paste0(time, " | ", event, " ~ \n ", paste(covariates, collapse=" + \n"))
    real_data_res <- data.frame(res=data[[time]], 
                                covariates=eq) 
    
    full_data <- rbind(real_data_res, sim_res)
    full_data <- full_data %>%
        mutate(covariates=factor(covariates, levels=c(eq, unique(sim_res$covariates)))) %>%
        ggplot(aes(x=res)) +
        geom_histogram() + 
        facet_wrap(~covariates) +
        xlab("Months")
    final <- list(hazFuncumFunPlot, res, full_data)
}





c1.os.pred <- simFun("intxsurv_1Y", "dead_1Y", c("age_cat", "PBlood", "DiseaseStatusAdvanced"), c1, 2110)
c1.os.pred[[1]]
