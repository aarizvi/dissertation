library(tidyverse)
library(ggrepel)
library(glue)
library(batch)

batch::parseCommandArgs(evaluate=TRUE)

files <- dir(path=input_path,
             pattern=glue::glue("*_{patt}.res"),
             full.names=TRUE)

data <- map_dfr(files, read_tsv)

title_name <- files %>%
    str_split("/") %>%
    map_chr(c(10,1)) %>%
    str_replace(".res", "") %>%
    str_replace("[0-9]+_", "") %>%
    unique() %>%
    str_replace_all("_", " ")

if(str_detect(title_name, "^D")){
    title_name <- str_replace(title_name, "^D", "Donor")
}else if(str_detect(title_name, "^R")){
    title_name <- str_replace(title_name, "^R", "Recipient")
}else{
    title_name <- str_replace(title_name, "^MM", "Mismatch")
}

manPlot <- function(data, plot.name, title_name, output_path){


    non_hits <- data %>%
	filter(PVALUE_M < 5e-08,
	      (PVALUE_c1 > 5e-05 | PVALUE_c2 > 1e-01 | SAMP_MAF_c1 < 0.01 | SAMP_MAF_c2 < 0.008)) %>%
	pull(RSID)

    data <- data %>%
        filter(Direction=="--" | Direction == "++",
               HR_M < 15 & HR_M > (1/15),
               SAMP_MAF_c1 > 0.005,
               SAMP_MAF_c2 > 0.005,
               !RSID %in% non_hits,
               str_detect(RSID, "rs"),
               HetPVal > 0.05,
               HetISq < 50,
	       INFO > 0.8) %>%
        select(RSID, CHR, POS, PVALUE_M)

    filtered.res <- data %>%
        filter((is.character(RSID) & is.numeric(CHR) & is.numeric(POS) & is.numeric(PVALUE_M))) %>%
        arrange(CHR, POS) %>%
        mutate(logp=-log10(PVALUE_M),
               pos=NA,
               index=NA)

    ind <- 0
    for (i in unique(filtered.res$CHR)){
        ind <- ind + 1
        filtered.res[filtered.res$CHR==i,]$index <- ind
    }

    nchr <- length(unique(filtered.res$CHR))
    lastbase <- 0
    ticks <- NULL

    for (i in unique(filtered.res$index)) {
        if (i==1) {
            filtered.res[filtered.res$index==i, ]$pos <- filtered.res[filtered.res$index==i, ]$POS
        } else {
            lastbase <- lastbase + tail(subset(filtered.res,index==i-1)$POS, 1)
            filtered.res[filtered.res$index==i, ]$pos <- filtered.res[filtered.res$index==i, ]$POS+lastbase
        }
        ticks <- c(ticks, (min(filtered.res[filtered.res$CHR == i,]$pos) + max(filtered.res[filtered.res$CHR == i,]$pos))/2 + 1)
    }

    xlabel <- 'Chromosome'
    labs <- unique(filtered.res$CHR)
    filtered.res <- filtered.res %>%
        na.omit
    xmax <- ceiling(max(filtered.res$POS) * 1.03)
    xmin <- floor(max(filtered.res$POS) * -0.03)
    ymin <- floor(min(filtered.res$logp))
    ymax <- ceiling(max(filtered.res$logp))
    # Set ymax to even number.
    if ((ymax %% 2) != 0 ) {
        ymax <- ymax + 1
    }
    cols <- c("gray10", "gray60")
    mycols <- rep(cols, nchr/2+1)
    suggestiveline <- -log10(5e-5)
    genomewideline <- -log10(5e-8)

    filtered.res <- filtered.res %>%
        mutate(col_pos = pos, col_logp = logp)

    g <- filtered.res %>%
        ggplot(aes(x=pos, y=logp, label=as.character(RSID))) +
        geom_point(aes(color=factor(CHR)), size=4) +
        geom_text_repel(data = filter(filtered.res, logp > -log10(5e-08)),
                        nudge_y = 36 - {filter(filtered.res, logp > -log10(5e-08)) %>% pull(logp)},
                        segment.size = 0.2,
                        segment.color = "grey50",
                        direction = "x",
                        size=5,
			angle=90,
			force=10) +
        ylab(expression(-log[10](italic(p)))) +
        scale_x_continuous(name=xlabel, breaks=ticks, labels=labs) +
        scale_y_continuous(breaks=seq(2, ymax, 2), labels=seq(2, ymax, 2),
                           limits = c(ymin-0.5, ymax), expand=c(0, 0)) +
        scale_colour_manual(values=mycols) +
        # geom_point(aes(x=col_pos, y=col_logp), color="red", size=5) +
        theme(panel.background=element_blank(),
              panel.grid.minor=element_blank(),
              axis.line=element_line(),
              text=element_text(size=22),
              plot.title=element_text(hjust = 0.5),
              panel.grid = element_blank(),
              panel.border = element_blank()) +
        geom_hline(yintercept=suggestiveline, color="blue", size=2.5) +
        geom_hline(yintercept=genomewideline, color="red", size=2.5, alpha=.5) +
        ggtitle(glue::glue("100 Days {title_name} Manhattan Plot")) +
        guides(shape = guide_legend(override.aes = list(color = "black"))) +
        guides(color=FALSE)
    ggsave(filename=glue::glue("{output_path}/{plot.name}_100d.jpeg"),
           device = "png",
           plot=g,
           width = 20,
           height = 12)
}

manPlot(data, patt, title_name, output_path)
