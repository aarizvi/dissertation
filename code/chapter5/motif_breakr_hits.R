library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP142.GRCh37) # dbSNP142 in hg19
library(BSgenome.Hsapiens.UCSC.hg19)  

snps.mb <- snps.from.rsid(rsid = "rs180962157",
                          dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37,
                          search.genome = BSgenome.Hsapiens.UCSC.hg19)



motifbreakr.results <- motifbreakR(snpList = snps.mb,
                                   pwmList = MotifDb, 
                                   threshold = 0.9)


motifbreakr.results


t <- motifbreakr.results[motifbreakr.results$geneSymbol %in% c("BRCA1", "CEBPB"),]

plotMB(results = t, rsid = "rs150712913", effect = "strong")





# library(MotIV)
# path <- system.file(package="MotIV")
# jaspar <- readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
# jaspar.scores <- readDBScores(paste(path,"/extdata/jaspar2010_PCC_SWU.scores",sep=""))
# example.motifs <- readPWMfile(paste(path,"/extdata/example_motifs.txt",sep=""))
# example.jaspar <- motifMatch(inputPWM=example.motifs,align="SWU",cc="PCC",
#                              database=jaspar,
#                              DBscores=jaspar.scores
#                              ,top=5)
# example.motifs
# viewAlignments(example.jaspar)[[1]]
# plot(example.jaspar[1:4],ncol=2,top=5, cex=0.8)
