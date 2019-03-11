outcomes <- c("TRM", "DRM", "OS", 
              "PFS", "REL", "GVHD", 
              "OF", "INF")
time_var <- c("time to death", "time to death", "time to death", 
              "time to relapse", "time to relapse", "time to death", 
              "time to death", "time to death")
covariates <- c("recipient age, BMI, graft source", "recipient age, disease status", "age, disease status, graft source",
                "recipient age, disease status", "conditioning regimen and intensity", "recipient age, donor age, BMI",
                "disease status, graft source", "age, BMI, CMV status")

df <- setNames(data.frame(
                 outcomes, 
                 time_var, 
                 covariates,
                 stringsAsFactors = FALSE),
               c("Survival Outcomes", "Time Interval", "Covariates"))

knitr::kable(df, caption="\\label{tab:models_run} Survival Models Analyzed", booktabs=TRUE) %>%
    kableExtra::kable_styling(latex_options="striped", font_size=9) #%>%
