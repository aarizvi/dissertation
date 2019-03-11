text_tbl <- data.frame(
    Outcomes = c("Overall Survival (OS)",
                 "Transplant Related Mortality (TRM)",
                 "Disease Related Mortality (DRM)", 
                 "Progression Free Survival (PFS)",
                 "Relapse (REL)"),
    Definitions=c(
        "defined as any patient (recipient) that died at any point within the first 12 month window post-HSCT of this observational study.",
        "defined as any cause of death except the underlying disease, pre-existing disease, accidental death or suicide, or death unrelated to the transplant.",
        "broadly defined as deaths relating to leukemia/MDS relapse/progression, including death attributed to toxicity or infection from anti-leukemic treatments post-HSCT.",
        "is defined as the time to relapse. All patients were analyzed as time to progression of disease",
        "patients who were not in CR pre-HSCT and the disease returns (relapse) after HSCT.")
)

text_tbl %>%
    kable(format = "latex",
      caption="\\label{tab:dbmt_outcomes} Definitions of Survival Outcomes", 
      booktabs=TRUE) %>%
    kable_styling(latex_options = "striped",
                  full_width=FALSE,
                  font_size=9) %>%
#    group_rows("Primary Outcomes", 1, 2, latex_gap_space=".5em") %>%
#    group_rows("Secondary Outcomes", 3, 5, latex_gap_space=".5em") %>%
    column_spec(2, width = "16em")
