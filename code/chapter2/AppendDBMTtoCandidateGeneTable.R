
library(data.table)
library(dplyr)
library(tidyr)
library(dtplyr)

setwd("~/Google Drive/Sucheston-Campbell Lab/Candidate Gene/Tables_Figs/CGTablewDISCOVERYBMT/")

snptable <- fread("~/Google Drive/Sucheston-Campbell Lab/Candidate Gene/Tables_Figs/cg_snptable_20170315.txt")
cg.main <- fread("~/Google Drive/Sucheston-Campbell Lab/Candidate Gene/Tables_Figs/MainCGTable/candidate_gene_main_table_20170321.txt")
v2 <- fread("~/Google Drive/Sucheston-Campbell Lab/Candidate Gene/Tables_Figs/cg_v2table_adj_20170316.txt")
nod2_3snp <- fread("~/Google Drive/Sucheston-Campbell Lab/Candidate Gene/Tables_Figs/CGTablewDISCOVERYBMT/NOD2_DRS_ALLOUTCOMES_FINAL.txt")

today <- function(x, y) paste0(x, "_", gsub("-", "",Sys.Date()), y)

#### Manipulate CG Table for appropriate merge with D-BMT ####

cg.main$Disease_Paper <- cg.main$Disease  # keep the original disease group from the paper

setkey(cg.main, Gene, SNP, Genome, Outcome, Disease, Graft)

# change the disease groups to match D-BMT
cg.main[Disease=="AML"]$Disease <- "AMLonly"
cg.main[Disease=="ALL"]$Disease <-  "ALLonly"
cg.main[Disease=="Mixed"]$Disease <-  "mixed"
cg.main[Disease=="ALL-AML"]$Disease <- "noMDS"
cg.main[Disease=="Pediatric AL"]$Disease <- "mixed"

#### Manipulate D-BMT for appropriate merge with CG Table ####

snptable2 <- snptable %>% filter(!rsID == "3 SNPs*") %>% bind_rows(nod2_3snp)

colnames(snptable2) <- c("Gene", "SNP", "chr", "BP", "N_D.BMT", "allele1", "allele2",  "coef", "se.coef", "HR_D.BMT", "95%-CI_D.BMT", 
                        "Pvalue_D.BMT",
                        "impute", "Genome", "Cohort", "Outcome", "Disease")



setkey(snptable2, Gene, SNP, Genome, Outcome, Cohort, Disease)

snp.meta <- snptable2 %>% filter(Cohort=="M") %>% select(SNP, Genome, Outcome, Disease, 
                                                          Pvalue_D.BMT, HR_D.BMT, `95%-CI_D.BMT`, N_D.BMT)

#### Merge CG main table and meta DISCOVeRY-BMT results ####
CGwDBMT <- left_join(cg.main, snp.meta)

apply(CGwDBMT, 2, function(x) sum(1*is.na(x)))

CGwDBMT %>% filter(is.na(Pvalue_D.BMT)) %>% select(SNP) %>% distinct()


CGwDBMT %>% filter(PMID == 18976442) %>% data.table()

CGwDBMT %>% filter(group == "Validation" & Significance == "Significant" & Gene == "IL6") %>% arrange(Pvalue_D.BMT) %>% 
  select(-(Population:Model)) %>% data.table()

write.table(CGwDBMT, today("SNP_CGwDBMT",".txt"), sep="\t", row.names = F, quote = F)

###############################################
#### Missing observations in DISCOVeRY-BMT ####
###############################################

missing_DBMT <- filter(CGwDBMT, (is.na(Pvalue_D.BMT))) %>% select(SNP:PMID, P.Value:Statistic_95CI, Pvalue_D.BMT:N_D.BMT, Significance)
head(missing_DBMT)

haps <- c("3 SNPs*", "CYP2B6*4/*5/*6/*7",  "IL23R SNPs*", "Tagged", "deletion", "rs12953/rs17343222/rs17354990")
mis.snp <- missing_DBMT %>% filter(! ( SNP %in% haps) ) %>% select(SNP, Gene) %>% distinct


library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(MafDb.1Kgenomes.phase3.hs37d5)

snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
rsid2alleles(mis.snp$SNP, caching=TRUE)

snpsById(snps, mis.snp$SNP)

mafdb <- MafDb.1Kgenomes.phase3.hs37d5
maf_miss_snps <- mafById(mafdb, mis.snp$SNP, pop="EUR_AF")
miss.snps <- maf_miss_snps %>% filter(!is.na(EUR_AF)) %>% select(ID) %>% rename(SNP=ID) %>% distinct() %>%
   inner_join(mis.snp) %>% filter(!duplicated(SNP)) %>% distinct()

write.table(miss.snps, "./missing_snps.txt", row.names = F, quote = F)



snptable2 %>% filter(SNP == "rs10306156")


snptable[SNP %in% c(mis.snp$SNP)]

missing_DBMT %>% filter(SNP == "rs4252303" )
snptable %>% filter(SNP == "rs4252303" & Disease == "mixed" & Outcome == "OS")

125122809 125167982
125146701 

table(missing_DBMT$Outcome)

table(missing_DBMT$SNP)

#####################
##### Gene based ####

cg.main$Group <- ifelse((cg.main$Graft == "URD" & cg.main$Population == "European"), "Replication", "Validation")
geneBcg <- cg.main %>% select(Gene, Reference, PMID, Genome, Outcome, Disease, Significance) %>% distinct()

GeneBasedTable <- v2 %>% filter(Cohort == "M") %>% 
  select(Gene, chr, Genome, Outcome, Disease, geneBased.pval:adj.topSNP.pval)
 

gene.qq <- data.frame(genome= c(rep("D", 4), rep("R", 4), rep("S", 4)), 
                      outcome = rep(c("DD" , "PFS" ,"TRM" ,"OS" ), 3), 
                      title = c(rep("Donor", 4), rep("Recipient", 4), rep("Mismatch", 4)), 
                                stringsAsFactors = F )
par(mfrow=c(3,4))

for(i in seq_len(nrow(gene.qq))){
  
  my.pvalues= GeneBasedTable %>% 
                    filter(Genome == gene.qq$genome[i] & Outcome == gene.qq$outcome[i] & Disease == "mixed") %>% 
                    select(adj.geneBased.pval)
  my.pvalues=my.pvalues$adj.geneBased.pval
  
  #Calculate expectations
  exp.pvalues<-(rank(my.pvalues, ties.method="first")+.5)/(length(my.pvalues)+1)
  
  #Make plot
  plot(-log10(exp.pvalues), -log10(my.pvalues), asp=1, main= paste(gene.qq$title[i], gene.qq$outcome[i], sep=" "))
  abline(0,1)
  
}
  
# Read data from the web
url = "https://gist.githubusercontent.com/stephenturner/806e31fce55a8b7175af/raw/1a507c4c3f9f1baaa3a69187223ff3d3050628d4/results.txt"

results = read.table(url, header=TRUE)
results = mutate(results, sig=ifelse(results$padj<0.05, "FDR<0.05", "Not Sig"))

p = ggplot(results, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col=sig)) +
  scale_color_manual(values=c("red", "black"))
p
p+geom_text(data=filter(results, padj<0.05), aes(label=Gene))  

write.table(CGwDBMT, today("Gene_CGwD-BMT",".txt"), sep="\t")
