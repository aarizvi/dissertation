
library(stringr)
library(ggplot2)
library(dplyr)

rep <- read.table("~/Google Drive/Writing paper/Tables/SignificantReplication.txt", header=T, sep="\t")
val <- read.table("~/Google Drive/Writing paper/Tables/SignificantValidation.txt", header=T, sep="\t")

cg.main <- read.table("~/Google Drive/CandidateGeneReplication/AppendingDBMT/SNP_CGwDBMT_20170413.txt", header=T, sep="\t", stringsAsFactors = F)

rep.val <- cg.main

rep.val <- data.frame(rbind(rep, val))

rep.val$Genome <- str_replace(rep.val$Genome, "D", "Donor")
rep.val$Genome <- str_replace(rep.val$Genome, "S", "Mismatch")
rep.val$Genome <- str_replace(rep.val$Genome, "R", "Recipient")

rep.val <- rep.val %>% filter(Outcome!="DD")
rep.val$Genome <- factor(rep.val$Genome, levels=c("Donor", "Recipient", "Mismatch"))

qqPlot <- function(data, ci){
        N  <- length(data$Pvalue_D.BMT)
        df <- data %>%
                select(Pvalue_D.BMT, Genome, Outcome)
        df <- df %>%
                arrange(Pvalue_D.BMT) %>%
                rename(observed=Pvalue_D.BMT, genome=Genome, outcome=Outcome)
        df <- df %>% mutate(observed = -log10(observed),
                                expected = -log10(1:N/N),
                                clower   = -log10(qbeta(ci,     1:N, N - 1:N + 1)),
                                cupper   = -log10(qbeta(1 - ci, 1:N, N - 1:N + 1)))
        log10Pe <- expression(paste("Expected -log"[10], plain(P)))
        log10Po <- expression(paste("Observed -log"[10], plain(P)))
        ggplot(df) +
                geom_point(aes(expected, observed), size = 2) +
                geom_abline(intercept = 0, slope = 1, alpha = 0.5, color="blue") +
                geom_line(aes(expected, cupper), linetype = 2, color= 'red') +
                geom_line(aes(expected, clower), linetype = 2, color = 'red') +
                xlab(log10Pe) +
                ylab(log10Po) +
                facet_grid(outcome~genome) +
                ggtitle("DISCOVeRY-BMT Replication/Validation P-values QQ Plot") +
                theme(plot.title = element_text(hjust = 0.5))
}

pdf("~/Google Drive/CandidateGeneReplication/figures/SNPbased_qqplot.pdf", width=8, height=8)
qqPlot(rep.val, 0.95)
dev.off()

