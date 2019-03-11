cols <- c("path_to_vcf_file", "output_name", "patient_subset", "genome", "outcome", "memory")
desc <- c("Chromosomes 1-22 VCF files (include full directory)",
          "Output file name",
          "Patient subset - i.e. ALL only, AML only, Mixed",
          "Donor, recipient, or mismatch",
          "Outcome - i.e. DRM, TRM, OS, PFS, INF, OF, REL",
          "Random access memory (RAM) in megabytes (MB) allocated to the job")

x <- data.frame(cbind(cols, desc), row.names = NULL)

colnames(x) <- c("Column name", "Description")


knitr::kable(x, caption="Manifest file column descriptions", booktabs=TRUE, align='c') %>%
    kableExtra::kable_styling(latex_options = "striped", font_size=9, position='center') %>%
    kableExtra::column_spec(2, width="40em")
