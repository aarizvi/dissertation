
##### Figure Table #####

setwd("~/Google Drive/Sucheston-Campbell Lab/Candidate Gene/Writing paper/Tables/")

# snptable <- fread("~/Google Drive/Sucheston-Campbell Lab/Candidate Gene/Tables_Figs/cg_snptable_20170315.txt")
# cg.main <- fread("~/Google Drive/Sucheston-Campbell Lab/Candidate Gene/Tables_Figs/MainCGTable/candidate_gene_main_table_20170316.txt")

library(tidyr)
library(dplyr)

CG_DBMT <- read.table("~/Google Drive/Sucheston-Campbell Lab/Candidate Gene/CandidateGeneReplication/AppendingDBMT/SNP_CGwDBMT_20170413.txt",
                      sep="\t", header = T, stringsAsFactors = F)
# CG_DBMT$Group <- ifelse((CG_DBMT$Graft == "URD" & CG_DBMT$Population == "European"), "Replication", "Validation")
 
CG_DBMT.sig <- CG_DBMT %>% filter(Significance == "Significant")

DBMT <- CG_DBMT.sig %>% select(Gene:Graft, Outcome, Genome,SNP, Pvalue_D.BMT, N_D.BMT) %>% 
  mutate(Group = "DISCOVeRY-BMT") %>%
  rename(numPvalue=Pvalue_D.BMT, N=N_D.BMT)

fig.tbl <- CG_DBMT.sig %>% select(Gene:Graft, Outcome, Genome,SNP, numPvalue, N, Group) %>% 
  bind_rows(DBMT) %>% 
  mutate(log_pval = -log(numPvalue, base = 10) )

#### plot the table #####

library(ggplot2)

#### Replication ------------------------------------ 

# now in the melted version Genes will include validation genes for D-BMT
# therefore only select replication genes
# also rs2066847 exist in validation but not 
sub.genes <- fig.tbl %>% filter( Group == "Replication") %>% select(Gene)
rep.fig.tbl <- fig.tbl %>% 
  filter(Gene %in% sub.genes$Gene & !Group == "Validation" & !Gene == "CCL2" & !SNP =="rs2066842")

# ## Fix CCR5 haplotypes before plotting -- if not fixed
# rep.fig.tbl[(rep.fig.tbl$SNP == "grp3vsgrp1" & rep.fig.tbl$Group == "Replication"),]$N <- (140+842)
# rep.fig.tbl[rep.fig.tbl$SNP == "grp3vsgrp1",]$SNP <- "D-R group 3 vs D-R group 1"
# 
# rep.fig.tbl[(rep.fig.tbl$SNP == "grp3vsgrp2" & rep.fig.tbl$Group == "Replication"),]$N <- (140+139)
# rep.fig.tbl[rep.fig.tbl$SNP == "grp3vsgrp2",]$SNP <- "D-R group 3 vs D-R group 2"


# ## Fix CCR5 haplotypes before plotting AGAIN -- if not fixed
rep.fig.tbl[rep.fig.tbl$SNP == "D-R group 3 vs D-R group 1",]$SNP <- "D-R group 3 vs\n D-R group 1"
rep.fig.tbl[rep.fig.tbl$SNP == "D-R group 3 vs D-R group 2",]$SNP <- "D-R group 3 vs\n D-R group 2"



pR <- ggplot(data=rep.fig.tbl, aes(x=SNP, y=log_pval))
pR <- pR + geom_hline(yintercept = -log(0.05, 10), colour = "red")
pR <- pR + geom_point(aes(colour=Genome, shape=Outcome, size=N), stroke = 1.5) + scale_shape(solid = FALSE)
group.facet.labs <- c(`DISCOVeRY-BMT`= "Replication in\n DISCOVeRY-BMT", Replication = "Literature P-values")
pR <- pR + facet_grid(Group ~ Gene, scales = "free_x", space = "free",
                      labeller = labeller(Group = group.facet.labs)) + guides(size=FALSE)

pR <- pR + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10),
                 axis.text.y = element_text( hjust = 1, size=10),
                 strip.text.x = element_text(size = 10, face="italic"),
                 strip.text.y = element_text(size = 13),
                 axis.title = element_text(size= 12))

pR <- pR + ylab("-log10(Pvalue)")
pR <- pR + scale_y_continuous(limits=c(0,3.2))
pR

pdf("~/Google Drive/Sucheston-Campbell Lab/Candidate Gene/Writing paper/Tables/ReplicationPlot.pdf",
    width= 11.5, height = 6.6)
pR
dev.off()

#### Validation -------------------------------------------------

# Plot genes studied more than once for main paper
m.stu.gene <- read.table("~/Google Drive/Sucheston-Campbell Lab/Candidate Gene/Writing paper/Tables/GenesReportedSummary.txt",
                         header = T, stringsAsFactors = F, sep="\t")
sig.val.ct <- m.stu.gene %>% select(Gene, SignificantValidationReports) %>% 
  arrange(desc(SignificantValidationReports)) %>% na.omit()
sub.genes2 <- sig.val.ct %>% filter(SignificantValidationReports > 1) %>% .[[1]]
val.fig.tbl <- filter(fig.tbl, Gene %in% sub.genes2 & !Group == "Replication")

pV <- ggplot(data=val.fig.tbl, aes(x=SNP, y=log_pval))
pV <- pV + geom_hline(yintercept = -log(0.05, 10), colour = "red")
pV <- pV + geom_point(aes(colour=Genome, shape=Outcome, size=N), stroke = 1.5) + scale_shape(solid = FALSE)

group.facet.labs <- c(`DISCOVeRY-BMT`= "Validation in\n DISCOVeRY-BMT", Validation = "Literature P-values")
pV <- pV + facet_grid(Group ~ Gene, scales = "free_x", space = "free",
                      labeller = labeller(Group = group.facet.labs)) + guides(size=FALSE)

pV <- pV + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10),
                 axis.text.y = element_text( hjust = 1, size=10),
                 strip.text.x = element_text(size = 10, face="italic"),
                 strip.text.y = element_text(size = 13),
                 axis.title = element_text(size= 12))
pV <- pV + ylab("-log10(Pvalue)")
pV


pdf("~/Google Drive/Sucheston-Campbell Lab/Candidate Gene/Writing paper/Tables/ValidationPlot.pdf",
    width= 11.5, height = 6.6)
pV
dev.off()

val.fig.tbl %>% filter(Gene == "IL6")
val.fig.tbl %>% filter(Gene == "TDG")

# Plot other genes for supplemantal
m.stu.gene <- read.table("~/Google Drive/Sucheston-Campbell Lab/Candidate Gene/Writing paper/Tables/GenesReportedSummary.txt",
                         header = T, stringsAsFactors = F, sep="\t")
sig.val.ct <- m.stu.gene %>% select(Gene, SignificantValidationReports) %>% 
  arrange(desc(SignificantValidationReports)) %>% na.omit()
sub.genes2 <- sig.val.ct %>% filter(SignificantValidationReports == 1) %>% .[[1]]
val.fig.tbl <- filter(fig.tbl, Gene %in% sub.genes2 & !Group == "Replication")

pV2 <- ggplot(data=val.fig.tbl, aes(x=SNP, y=log_pval))
pV2 <- pV2 + geom_hline(yintercept = -log(0.05, 10), colour = "red")
pV2 <- pV2 + geom_point(aes(colour=Genome, shape=Outcome, size=N), stroke = 1.5) + scale_shape(solid = FALSE)

group.facet.labs <- c(`DISCOVeRY-BMT`= "Validation in\n DISCOVeRY-BMT", Validation = "Literature P-values")
pV2 <- pV2 + facet_grid(Group ~ Gene, scales = "free_x", space = "free",
                      labeller = labeller(Group = group.facet.labs)) + guides(size=FALSE)

pV2 <- pV2 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10),
                 axis.text.y = element_text( hjust = 1, size=10),
                 strip.text.x = element_text(size = 10, angle = 90, face="italic"),
                 strip.text.y = element_text(size = 13),
                 axis.title = element_text(size= 12))

pV2 <- pV2 + ylab("-log10(Pvalue)")
pV2


pdf("~/Google Drive/Sucheston-Campbell Lab/Candidate Gene/Writing paper/Tables/ValidationPlotSupp.pdf",
    width= 12, height = 6)
pV2
dev.off()
