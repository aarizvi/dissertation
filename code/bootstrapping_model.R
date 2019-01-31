library(boot)

library(splines)


covariate.file <- read_tsv("~/Google Drive/OSU_PHD/DBMT_100d/DBMT_PhenoData_EA_long_allVar_20181216.txt")
# cohort 1
c1 <- covariate.file %>%
    filter(sample_type=="recipient",
           cohort==1)


time <- "intxsurv_1Y"
event <- "dead_1Y"
covariates <- c("distatD", "PBlood", "age")



#time.death <- c1 %>% select(intxsurv_1Y, dead_1Y)

covariate.file <- read_tsv("~/Google Drive/OSU_PHD/DBMT_100d/DBMT_PhenoData_EA_long_allVar_20181216.txt")
# cohort 1
c1 <- covariate.file %>%
    filter(sample_type=="recipient",
           cohort==1) %>%
    mutate(PBlood=ifelse(PBlood==0, "Bone Marrow", "Peripheral Blood"))

c1 <- c1 %>% select(intxsurv_1Y, dead_1Y, age, distatD, PBlood)

eg <- expand.grid(intxsurv_1Y=seq(0, 12, by=0.01),
                  age=seq(0, max(c1$age), by=1))


eg <- rbind(eg, eg, eg, eg) %>%
    mutate(PBlood=rep(c("Bone Marrow","Peripheral Blood"), each=nrow(eg)*2),
           distatD = rep(c("early","advanced","early","advanced"), each=nrow(eg)),
           dead_1Y=ifelse(intxsurv_1Y==12, 0, 1))



pred.mod <- predict(mod, eg, type="expected") 

#library(splines)# for ns
c1.cox <- coxph(Surv(intxsurv_1Y, dead_1Y) ~ age + strata(PBlood, distatD), data = c1)

#c1.surv <- survfit(c1.cox)
#agec <- cut(c1$age, c(0, 39, 49, 59, 69, 100))
c1.cens <- survfit(Surv(intxsurv_1Y - 0.001*(dead_1Y), 1-dead_1Y) ~ 1, data = c1)




c1.fun <- function(d) { 
    cox <- coxph(Surv(d$intxsurv_1Y, d$dead_1Y) ~ d$age + strata(d$distatD, d$PBlood))
    pred.mod <- predict(cox, c1, type="expected")
    exp(-pred.mod)
}
c1.str <- cbind(
    with(c1, strata(PBlood, distatD)),
    with(c1, strata(distatD, PBlood, agec))
)

c1.mod <- censboot(c1, 
                   c1.fun,
                   R = 499,
                   F.surv = c1.surv,
                   G.surv = c1.cens, 
                   cox = c1.cox, 
                   strata = c1.str, 
                   sim = "model")




mel.env <- envelope(c1.mod)$point
th <- seq(0, 100, length.out = ncol(mel.env))
plot(th, mel.env[1, ],  ylim = c(-2, 2), xlab = "age (mm)", ylab = "linear predictor", type = "n")
lines(th, c1.mod$t0, lty = 1)
matlines(th, t(mel.env), lty = 2)


#### just age
covariate.file <- read_tsv("~/Google Drive/OSU_PHD/DBMT_100d/DBMT_PhenoData_EA_long_allVar_20181216.txt")
# cohort 1
c1 <- covariate.file %>%
    filter(sample_type=="recipient",
           cohort==1) %>%
    mutate(PBlood=ifelse(PBlood==0, "Bone Marrow", "Peripheral Blood"))

c1 <- c1 %>% select(intxsurv_1Y, dead_1Y, age, distatD, PBlood)

eg <- expand.grid(intxsurv_1Y=seq(0, 12, by=0.01),
                  age=seq(0, max(c1$age), by=1))


eg <- rbind(eg, eg) %>%
    mutate(#PBlood=rep(c("Bone Marrow","Peripheral Blood"), each=nrow(eg)*2),
           distatD = rep(c("early","advanced"), each=nrow(eg)),
           dead_1Y=ifelse(intxsurv_1Y==12, 0, 1))



c1.cox <- coxph(Surv(intxsurv_1Y, dead_1Y) ~ strata(distatD), data = c1)

c1.surv <- survfit(c1.cox)
agec <- cut(c1$age, c(0, 39, 49, 59, 69, 100))
c1.cens <- survfit(Surv(intxsurv_1Y - 0.001*(dead_1Y == 1), dead_1Y != 1) ~ strata(agec), data = c1)

c1.fun <- function(d) { 
    cox <- coxph(Surv(d$intxsurv_1Y, d$dead_1Y) ~ strata(d$distatD))
    # pred.mod <- predict(cox, eg, type="expected")
    # exp(-pred.mod)
    # 
    eta <- cox$linear.predictors
    
    
    predict(cox, c1)
    

}
c1.str <- cbind(with(c1, strata(distatD),
    with(c1, strata(agec))))

c1.mod <- censboot(c1, 
                   c1.fun,
                   R = 499,
                   F.surv = c1.surv,
                   G.surv = c1.cens, 
                   cox = c1.cox, 
                   strata = strata(c1$distatD), 
                   sim = "model")




mel.env <- envelope(c1.mod)$point
th <- seq(0, 100, length.out = ncol(mel.env))
plot(th, mel.env[1, ],  ylim = c(-2, 2), xlab = "age (mm)", ylab = "linear predictor", type = "n")
lines(th, c1.mod$t0, lty = 1)
matlines(th, t(mel.env), lty = 2)







library(splines)# for ns
mel.cox <- coxph(Surv(time, status == 1) ~ ns(thickness, df=4) + strata(ulcer),
                 data = melanoma)
mel.surv <- survfit(mel.cox)
agec <- cut(melanoma$age, c(0, 39, 49, 59, 69, 100))
mel.cens <- survfit(Surv(time - 0.001*(status == 1), status != 1) ~
                        strata(agec), data = melanoma)
d <- melanoma
mel.fun <- function(d) { 
    t1 <- ns(d$thickness, df=4)
    cox <- coxph(Surv(d$time, d$status == 1) ~ t1+strata(d$ulcer))
    ind <- !duplicated(d$thickness)
    u <- d$thickness[!ind]
    eta <- cox$linear.predictors[!ind]
    sp <- smooth.spline(u, eta, df=20)
    th <- seq(from = 0.25, to = 10, by = 0.25)
    predict(sp, th)$y
}
mel.str <- cbind(melanoma$ulcer, agec)
mel.mod <- censboot(melanoma, mel.fun, R = 499, F.surv = mel.surv,
                    G.surv = mel.cens, cox = mel.cox, strata = mel.str, sim = "model")
# To plot the original predictor and a 95% pointwise envelope for it




