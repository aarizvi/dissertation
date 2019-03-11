dons <- c("Donor","age, sex, race group, blood type (ABO)")
recs <- c("Recipient",
          "age, sex, race group, blood type (ABO), graft source, causes of death (levels 1-7), conditioning regimen and intensity, TBI fractionation, GVHD prophylaxis, current disease type, current disease progression, prior disease type, cytogenetics status, time to death, time to relapse, AML type, transplant year, HLA match (8/8 or 10/10), cytomegalyvirus, Lansky Score, Karnovsky Score, Principal Components from EIGENSTRAT")

knitr::kable(setNames(data.frame(rbind(dons, recs), row.names = NULL), c("Genome", "Description")), 
             caption="\\label{tab:clin_char} Broad Overview of DISCOVeRY-BMT Clinical Characteristics", booktabs=TRUE, align='c') %>%
    kableExtra::kable_styling(latex_options = "striped", font_size=9, position='c') %>%
    kableExtra::column_spec(2, width="24em")
