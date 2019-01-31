##### Figure Table #####

setwd("~/Google Drive/OSU_PHD/dissertation/code/chapter2/")
# snptable <- fread("~/Google Drive/Sucheston-Campbell Lab/Candidate Gene/Tables_Figs/cg_snptable_20170315.txt")
# cg.main <- fread("~/Google Drive/Sucheston-Campbell Lab/Candidate Gene/Tables_Figs/MainCGTable/candidate_gene_main_table_20170316.txt")

library(tidyverse)
library(stringr)
library(grid)
library(gtable)
library(gridExtra)



#### plot the table #####


#### Replication ------------------------------------

# now in the melted version Genes will include validation genes for D-BMT
# therefore only select replication genes
# also rs2066847 exist in validation but not

# 
# rep.fig.tbl <- fig.tbl %>%
#     filter(Gene %in% sub.genes$Gene,
#            !Group == "Validation",
#            !Gene == "CCL2" & !SNP =="rs2066842")

# ## Fix CCR5 haplotypes before plotting -- if not fixed
# rep.fig.tbl[(rep.fig.tbl$SNP == "grp3vsgrp1" & rep.fig.tbl$Group =="Replication"),]$N <- (140+842)
# rep.fig.tbl[rep.fig.tbl$SNP == "grp3vsgrp1",]$SNP <- "D-R group 3 vs D-R group 1"
#
# rep.fig.tbl[(rep.fig.tbl$SNP == "grp3vsgrp2" & rep.fig.tbl$Group == "Replication"),]$N <- (140+139)
# rep.fig.tbl[rep.fig.tbl$SNP == "grp3vsgrp2",]$SNP <- "D-R group 3 vs D-R group 2"


# ## Fix CCR5 haplotypes before plotting AGAIN -- if not fixed
# rep.fig.tbl[rep.fig.tbl$SNP == "D-R group 3 vs D-R group 1",]$SNP <- "D-R group 3 vs\n D-R group 1"
# rep.fig.tbl[rep.fig.tbl$SNP == "D-R group 3 vs D-R group 2",]$SNP <- "D-R group 3 vs\n D-R group 2"
# rep.fig.tbl$Group <- as.factor(rep.fig.tbl$Group)
# levels(rep.fig.tbl$Group) <- c("Replication in DISCOVeRY-BMT", "Literature")
# rep.fig.tbl$Group <- as.character(rep.fig.tbl$Group)

library(stringr)
library(grid)
library(gtable)
library(gridExtra)





CG_DBMT <- read.table(
    "~/Google Drive/OSU_PHD/dissertation/code/chapter2/SNP_CGwDBMT_20170522.txt",
    sep = "\t",
    header = T,
    stringsAsFactors = F
)
# CG_DBMT$Group <- ifelse((CG_DBMT$Graft == "URD" & CG_DBMT$Population == "European"), "Replication", "Validation")

CG_DBMT.sig <- CG_DBMT %>%
    filter(Significance == "Significant")

DBMT <- CG_DBMT.sig %>%
    select(Gene:Graft, Outcome, Genome,SNP, Pvalue_D.BMT, N_D.BMT) %>%
    mutate(Group = "DISCOVeRY-BMT") %>%
    rename(numPvalue=Pvalue_D.BMT, N=N_D.BMT)

fig.tbl <- CG_DBMT.sig %>%
    select(Gene:Graft, Outcome, Genome,SNP,
           numPvalue, N, Group) %>%
    bind_rows(DBMT) %>%
    mutate(log_pval = -log(numPvalue, base = 10) )

mac <- read.table("~/Google Drive/OSU_PHD/dissertation/code/chapter2/MacMillan2003.txt", header=T, sep="\t", stringsAsFactors = F)
mac <- mac %>%
    mutate(log_pval=-log10(numPvalue),
           Group=ifelse(Group=="Replication in DISCOVeRY-BMT", "DISCOVeRY-BMT", "Replication"))
fig.tbl <- rbind(fig.tbl, mac)


plotFun <- function(repval){
    
    sub.genes <- fig.tbl %>%
        filter(Group == repval) %>%
        pull(Gene)
    
    fig.tbl <- fig.tbl %>%
        filter(Gene %in% c(sub.genes, "CCL2", "MIF"),
               Group %in% c(repval, "DISCOVeRY-BMT"),
               SNP != "rs2066842") 
    
    fig.tbl$Group <- as.factor(fig.tbl$Group)
    levels(fig.tbl$Group) <- c(paste0(repval, " in DISCOVeRY-BMT"), "Literature")
    fig.tbl$Group <- as.character(fig.tbl$Group)
    
    if(repval=="Replication"){
        group.facet.labs <- c(Literature="Literature", 
                              `Replication in DISCOVeRY-BMT`="Replication in \n DISCOVeRY-BMT")
    }else if(repval=="Validation"){
        group.facet.labs <- c(Literature="Literature", 
                              `Validation in DISCOVeRY-BMT`="Validation in \n DISCOVeRY-BMT")
    }
    
    pR <- fig.tbl %>%
        mutate(SNP=str_replace(SNP, "D-R group 3 vs D-R group 1", "D-R group 3 vs\n D-R group 1"),
               SNP=str_replace(SNP, "D-R group 3 vs D-R group 2", "D-R group 3 vs\n D-R group 2"),
               Outcome=str_replace(Outcome, "DD", "DRM"),
               Genome=str_replace(Genome, "^R$", "Recipient"),
               Genome=str_replace(Genome, "^D$", "Donor"),
               Genome=str_replace(Genome, "^S$", "Mismatch"), 
               Gene=str_replace(Gene, "NOD2", "NOD2/\nCARD15")) %>% 
        ggplot(aes(SNP, log_pval)) +
        geom_hline(yintercept = -log10(0.05), color="red") +
        geom_point(aes(color=Genome, shape=Outcome, size=N), stroke=1.5) +
        facet_grid(Group~Gene, 
                   scales="free_x", 
                   space="free", 
                   labeller = labeller(Group=group.facet.labs)) +
        guides(size=FALSE) +
        theme(rect = element_rect(fill = "white"),
              legend.key = element_rect(fill = "white"),
              panel.border = element_rect(colour = "gray85", fill=NA),
              panel.background = element_rect(fill = "white"),
              panel.grid.major = element_line(colour = "gray85"),
              plot.background = element_rect(fill = "white"),
              axis.text.x = element_text(angle = 90, hjust = 1, size=18),
              axis.text.y = element_text(hjust = 1, size=18, angle=90),
              strip.text.x = element_text(size = 18, face="italic", angle=90),
              strip.text.y = element_text(size = 18, angle=90),
              axis.title = element_text(size=18),
              axis.title.x = element_text(size=18, angle=90),
              panel.spacing.x=unit(0.1, "lines"),
              panel.spacing.y=unit(0.1, "lines"),
              legend.position=c(0.9, 0.5)) +
        ylab("-log10(Pvalue)") +
        scale_shape_discrete(solid=FALSE, 
                             guide = guide_legend(direction = "horizontal", 
                                                  title.position = "left",
                                                  title.theme=element_text(angle=90, size=18),
                                                  label.position="bottom",
                                                  label.hjust = 0.5, 
                                                  label.vjust = 0.5,
                                                  label.theme = element_text(angle = 90, size=18))) +
        scale_color_discrete(guide = guide_legend(direction = "horizontal",
                                                  title.position = "left",
                                                  title.theme=element_text(angle=90, size=18),
                                                  label.position="bottom", 
                                                  label.hjust = 0.5, 
                                                  label.vjust = 0.5,
                                                  label.theme = element_text(angle = 90, size=18))) 
    
    
    xr <- ggplotGrob(pR)
    #labelR <- "Pvalues"
    labelT <- "Genes"
    # get the positions of strips in the gtable: t = top, l=left
    #posR <- subset(xr$layout, grepl("strip-r", name), select=t:r)
    posT <- subset(xr$layout, grepl("strip-t", name), select=t:r)
    # add a new column to the right of current right strips,
    # and a new row on top of current top stips
    #width <- xr$widths[max(posR$r)]
    height <- xr$heights[min(posT$t)]
    #xr <- gtable_add_cols(xr, width, max(posR$r))
    xr <- gtable_add_rows(xr, height, min(posT$t)-1)
    # Construct the new strip grobs
    # stripR <- gTree(name = "Strip_right", children = gList(
    #     rectGrob(gp = gpar(col = NA,fill = "grey85")),
    #     textGrob(labelR, rot = -90, gp = gpar(fontsize = 13, col = "grey10"))))
    stripT <- gTree(name = "Strip_top",
                    children = gList(
                        rectGrob(gp = gpar(col = NA, fill = "grey85")),
                        textGrob(labelT, rot = 90,gp = gpar(fontsize = 18, col = "grey10"))))
    # Position the grobs in the gtable
    #xr <- gtable_add_grob(xr, stripR, t = min(posR$t)+1, l = max(posR$r) +1, b = max(posR$b)+1, name = "strip-right")
    xr <- gtable_add_grob(xr,stripT,t = min(posT$t),l = min(posT$l),r = max(posT$r),name = "strip-top")
    # Add small gaps between strips
    #xr <- gtable_add_cols(xr, unit(1/5, "line"), max(posR$r))
    xr <- gtable_add_rows(xr,
                          unit(1/5, "line"),
                          min(posT$t))
    # Draw it
    grid.newpage()
    grid.draw(xr)
} 

plotFun("Replication")

plotFun("Validation")


# fig.tbl %>%
#     filter(Gene %in% c(sub.genes, "CCL2"),
#            Group %in% c(!!repval, "DISCOVeRY-BMT"),
#            SNP != "rs2066842") %>%
#     mutate(SNP=str_replace(SNP, "D-R group 3 vs D-R group 1", "D-R group 3 vs\n D-R group 1"),
#            SNP=str_replace(SNP, "D-R group 3 vs D-R group 2", "D-R group 3 vs\n D-R group 2"),
#            Group=as.character(as.factor(Group)),
#            Outcome=str_replace(Outcome, "DD", "DRM"),
#            Genome=str_replace(Genome, "R", "Recipient"),
#            Genome=str_replace(Genome, "D", "Donor"),
#            Genome=str_replace(Genome, "S", "Mismatch"),
#            Gene=str_replace(Gene, "NOD2", "NOD2/CARD15")) %>%
#     filter(Gene != "MIF",
#            Gene %in% sub.genes$Gene,
#            Group != "Validation",
#            Gene != "CCL2",
#            SNP !="rs2066842") %>%
#     ggplot(aes(SNP, log_pval)) +
#     geom_hline(yintercept = -log10(0.05), color="red") +
#     geom_point(aes(color=Genome, shape=Outcome, size=N), stroke=1.5) +
#     facet_grid(Group ~ Gene,
#                scales="free_x",
#                space="free",
#                labeller = labeller(Group=group.facet.labs)) +
#     guides(size=FALSE) +
#     theme(rect = element_rect(fill = "white"),
#          legend.key = element_rect(fill = "white"),
#          panel.border = element_rect(colour = "gray85", fill=NA),
#          panel.background = element_rect(fill = "white"),
#          panel.grid.major = element_line(colour = "gray85"),
#          plot.background = element_rect(fill = "white"),
#          axis.text.x = element_text(angle = 90, hjust = 1, size=12),
#          axis.text.y = element_text( hjust = 1, size=12),
#          strip.text.x = element_text(size = 12, face="italic"),
#          strip.text.y = element_text(size = 12),
#          axis.title = element_text(size=12),
#          panel.spacing.x=unit(0.0001, "lines"),
#          panel.spacing.y=unit(0.0001, "lines")) +
#          #,strip.background.x =  =  element_rect(color = "black", size = 1)) +
#     ylab("-log10(Pvalue)") +
#     scale_y_continuous(limits=c(0,3.2))  +
#     scale_shape_discrete(solid=FALSE)

# pR

##### TO ADD GENE LABEL ON TOP OF FACETS #######
# grab  ggplot object info
xr <- ggplotGrob(pR)
#labelR <- "Pvalues"
labelT <- "Genes"
# get the positions of strips in the gtable: t = top, l=left
#posR <- subset(xr$layout, grepl("strip-r", name), select=t:r)
posT <- subset(xr$layout, grepl("strip-t", name), select=t:r)
# add a new column to the right of current right strips,
# and a new row on top of current top stips
#width <- xr$widths[max(posR$r)]
height <- xr$heights[min(posT$t)]
#xr <- gtable_add_cols(xr, width, max(posR$r))
xr <- gtable_add_rows(xr, height, min(posT$t)-1)
# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
rectGrob(gp = gpar(col = NA,fill = "grey85")),
textGrob(labelR, rot = -90, gp = gpar(fontsize = 13, col = "grey10"))))
stripT <- gTree(name = "Strip_top",
children = gList(
rectGrob(gp = gpar(col = NA, fill = "grey85")),
textGrob(labelT, gp = gpar(fontsize = 13, col = "grey10"))))
# Position the grobs in the gtable
#xr <- gtable_add_grob(xr, stripR, t = min(posR$t)+1, l = max(posR$r) +1, b = max(posR$b)+1, name = "strip-right")
xr <- gtable_add_grob(xr,stripT,t = min(posT$t),l = min(posT$l),r = max(posT$r),name = "strip-top")
# Add small gaps between strips
#xr <- gtable_add_cols(xr, unit(1/5, "line"), max(posR$r))
xr <- gtable_add_rows(xr,
unit(1/5, "line"),
min(posT$t))
# Draw it
grid.newpage()
#pdf("~/Google Drive/CandidateGeneReplication/figures/ReplicateNewColwLabels.pdf",
#width= 11.5,
#height = 6.6)
grid.draw(xr)
dev.off()















#### Validation -------------------------------------------------

# Plot genes studied more than once for main paper
m.stu.gene <- read.table("~/Google Drive/Writing paper/Tables/GenesReportedSummary.txt",
header = T, stringsAsFactors = F, sep="\t")
sig.val.ct <- m.stu.gene %>% select(Gene, SignificantValidationReports) %>%
        arrange(desc(SignificantValidationReports)) %>%
        na.omit()
sub.genes2 <- sig.val.ct %>% filter(SignificantValidationReports > 1) %>% .[[1]]
val.fig.tbl <- filter(fig.tbl, Gene %in% sub.genes2 & !Group == "Replication")


val.fig.tbl$Group <- as.factor(val.fig.tbl$Group)

val.fig.tbl$SNP == "rs"



levels(val.fig.tbl$Group) <- c("Validation in DISCOVeRY-BMT", "Literature")
val.fig.tbl$Group <- as.character(val.fig.tbl$Group)
# change labels for genomes
library(stringr)
val.fig.tbl$Outcome <- str_replace(val.fig.tbl$Outcome, "DD", "DRM")
val.fig.tbl$Outcome <- as.factor(val.fig.tbl$Outcome)
val.fig.tbl$Genome <- str_replace(val.fig.tbl$Genome, "R", "Recipient")
val.fig.tbl$Genome <- str_replace(val.fig.tbl$Genome, "D", "Donor")
val.fig.tbl$Genome <- str_replace(val.fig.tbl$Genome, "S", "Mismatch")

val.fig.tbl$Gene <- str_replace(val.fig.tbl$Gene, "NOD2", "NOD2/CARD15")


levels(val.fig.tbl$Outcome) <- c("OS", "PFS", "TRM", "DRM")

pV <- ggplot(data=val.fig.tbl, aes(x=SNP, y=log_pval))
pV <- pV + geom_hline(yintercept = -log(0.05, 10), colour = "red")
pV <- pV + geom_point(aes(colour=Genome, shape=Outcome, size=N), stroke
= 1.5) + scale_shape(solid = FALSE)

group.facet.labs <- c(`Validation in DISCOVeRY-BMT`= "Validation in\n
DISCOVeRY-BMT", Literature = "Literature")
pV <- pV + facet_grid(Group ~ Gene, scales = "free_x", space = "free",
labeller = labeller(Group = group.facet.labs)) +
guides(size=FALSE)

pV <- pV +  theme(rect = element_rect(fill = "white"),
legend.key = element_rect(fill = "white"),
panel.border = element_rect(colour = "gray85", fill=NA),
panel.background = element_rect(fill = "white"),
panel.grid.major = element_line(colour = "gray85"),
plot.background = element_rect(fill = "white"),
axis.text.x = element_text(angle = 90, hjust = 1, size=10),
axis.text.y = element_text( hjust = 1, size=10),
strip.text.x = element_text(size = 10, face="italic"),
strip.text.y = element_text(size = 10),
axis.title = element_text(size= 12))
pV <- pV + ylab("-log10(Pvalue)")
pV


### ADD GENE LABEL ABOVE FACETS ###
# grab  ggplot object info
xv <- ggplotGrob(pV)

#labelR <- "Pvalues"
labelT <- "Genes"

# get the positions of strips in the gtable: t = top, l=left
#posR <- subset(xr$layout, grepl("strip-r", name), select=t:r)
posT <- subset(xv$layout, grepl("strip-t", name), select=t:r)


# add a new column to the right of current right strips,
# and a new row on top of current top stips

#width <- xr$widths[max(posR$r)]
height <- xv$heights[min(posT$t)]


#xr <- gtable_add_cols(xr, width, max(posR$r))
xv <- gtable_add_rows(xv, height, min(posT$t)-1)

# Construct the new strip grobs
#stripR <- gTree(name = "Strip_right", children = gList(
#rectGrob(gp = gpar(col = NA, fill = "grey85")),
#textGrob(labelR, rot = -90, gp = gpar(fontsize = 13, col = "grey10"))))

stripT <- gTree(name = "Strip_top", children = gList(
rectGrob(gp = gpar(col = NA, fill = "grey85")),
textGrob(labelT, gp = gpar(fontsize = 13, col = "grey10"))))

# Position the grobs in the gtable
#xr <- gtable_add_grob(xr, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
xv <- gtable_add_grob(xv, stripT, t = min(posT$t), l = min(posT$l), r =
max(posT$r), name = "strip-top")

# Add small gaps between strips
#xr <- gtable_add_cols(xr, unit(1/5, "line"), max(posR$r))
xv <- gtable_add_rows(xv, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()

pdf("~/Google Drive/CandidateGeneReplication/figures/ValidationNewColwLabels.pdf",
width= 11.5, height = 6.6)
grid.draw(xv)
dev.off()
























val.fig.tbl %>% filter(Gene == "IL6")
val.fig.tbl %>% filter(Gene == "TDG")

# Plot other genes for supplemantal
m.stu.gene <- read.table("~/Desktop/OSU_PHD/candidate_gene/AppendingDBMT/GenesReportedSummary.txt",
                         header = T,
                         stringsAsFactors = F,
                         sep="\t")
sig.val.ct <- m.stu.gene %>%
        select(Gene, SignificantValidationReports) %>%
        arrange(desc(SignificantValidationReports)) %>%
        na.omit()
sub.genes2 <- sig.val.ct %>%
        filter(SignificantValidationReports == 1) %>%
        .[[1]]
val.fig.tbl <- filter(fig.tbl, Gene %in% sub.genes2 & !Group =="Replication")



val.fig.tbl$Group <- as.factor(val.fig.tbl$Group)


levels(val.fig.tbl$Group) <- c("Validation in DISCOVeRY-BMT", "Literature")
val.fig.tbl$Group <- as.character(val.fig.tbl$Group)
# change labels for genomes
library(stringr)

val.fig.tbl$Genome <- str_replace(val.fig.tbl$Genome, "R", "Recipient")
val.fig.tbl$Genome <- str_replace(val.fig.tbl$Genome, "D", "Donor")
val.fig.tbl$Genome <- str_replace(val.fig.tbl$Genome, "S", "Mismatch")

val.fig.tbl$Gene <- str_replace(val.fig.tbl$Gene, "NOD2", "NOD2/CARD15")

val.fig.tbl <- val.fig.tbl %>% filter(!Gene %in% c("IL1A", "IL1B", "PTPN22"))



pV2 <- ggplot(data=val.fig.tbl, aes(x=SNP, y=log_pval))
pV2 <- pV2 + geom_hline(yintercept = -log(0.05, 10), colour = "red")
pV2 <- pV2 + geom_point(aes(colour=Genome, shape=Outcome, size=N),
stroke = 1.5) + scale_shape(solid = FALSE)

group.facet.labs <- c(`Validation in DISCOVeRY-BMT`= "Validation in\n
DISCOVeRY-BMT", Literature = "Literature")
pV2 <- pV2 + facet_grid(Group ~ Gene, scales = "free_x", space = "free",
labeller = labeller(Group = group.facet.labs)) +
guides(size=FALSE)

pV2 <- pV2 + theme(rect = element_rect(fill = "white"),
legend.key = element_rect(fill = "white"),
panel.border = element_rect(colour = "gray85", fill=NA),
panel.background = element_rect(fill = "white"),
panel.grid.major = element_line(colour = "gray85"),
plot.background = element_rect(fill = "white"),
axis.text.x = element_text(angle = 90, hjust = 1,
size=10),
axis.text.y = element_text( hjust = 1, size=10),
strip.text.x = element_text(size = 10, face="italic"),
strip.text.y = element_text(size = 10),
axis.title = element_text(size= 12))

pV2 <- pV2 + ylab("-log10(Pvalue)")
pV2


### ADD GENE LABEL ABOVE FACETS ###
# grab  ggplot object info
xv2 <- ggplotGrob(pV2)

#labelR <- "Pvalues"
labelT <- "Genes"

# get the positions of strips in the gtable: t = top, l=left
#posR <- subset(xr$layout, grepl("strip-r", name), select=t:r)
posT <- subset(xv2$layout, grepl("strip-t", name), select=t:r)


# add a new column to the right of current right strips,
# and a new row on top of current top stips

#width <- xr$widths[max(posR$r)]
height <- xv$heights[min(posT$t)]


#xr <- gtable_add_cols(xr, width, max(posR$r))
xv2 <- gtable_add_rows(xv2, height, min(posT$t)-1)

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(
rectGrob(gp = gpar(col = NA, fill = "grey85")),
textGrob(labelR, rot = -90, gp = gpar(fontsize = 13, col =
"grey10"))))

stripT <- gTree(name = "Strip_top", children = gList(
rectGrob(gp = gpar(col = NA, fill = "grey85")),
textGrob(labelT, gp = gpar(fontsize = 13, col = "grey10"))))

# Position the grobs in the gtable
#xr <- gtable_add_grob(xr, stripR, t = min(posR$t)+1, l = max(posR$r) +
1, b = max(posR$b)+1, name = "strip-right")
xv2 <- gtable_add_grob(xv2, stripT, t = min(posT$t), l = min(posT$l), r
= max(posT$r), name = "strip-top")

# Add small gaps between strips
#xr <- gtable_add_cols(xr, unit(1/5, "line"), max(posR$r))
xv2 <- gtable_add_rows(xv2, unit(1/5, "line"), min(posT$t))

# Draw it
grid.newpage()

pdf("~/Google Drive/CandidateGeneReplication/figures/SuppValidationNewColwLabels.pdf",
width=15.9, height = 6.6)
grid.draw(xv2)
dev.off()
