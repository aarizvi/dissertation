# running statistics

library(tidyverse)
x <- read_tsv("~/Desktop/run_stats.txt")


plot_maker <- function(data, mapping=aes(), plot.title, ...){
    
    snps_rm <- data %>%
        group_by(chr) %>%
        summarize(avg_snps_removed=mean(snps_removed, na.rm=TRUE)) %>%
        pull(avg_snps_removed) %>%
        sum() / 1e06 
    
    snps_rm <- round(snps_rm, 3)
    
    data_size <- sum(data$size)/1000/2
    
    
    snps_run <- data %>%
        pull(snps_analyzed) %>%
        unique() %>%
        sum() /1e06 
    snps_run <- round(snps_run, 3)
    
    total_snps = snps_rm + snps_run
    
    ggplot(data, mapping) + 
        geom_col(...) +
        geom_text(position=position_dodge(width=0.9),
                  hjust=1.2, angle=0, color="white", size=5.5, fontface="bold")  +
        theme_bw() +
        theme(plot.title = element_text(size=24),
              axis.text.x = element_text(size=20),
              axis.text.y = element_text(size=20),
              axis.title.x = element_text(size=22),
              axis.title.y = element_text(size=22),
              legend.text = element_text(size=18),
              legend.key.size = unit(5, "mm"),
              strip.text.x = element_text(size=20),
              legend.position="bottom",
              legend.title = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        ggtitle(plot.title) +
        ylab("Runtime (hours)") +
        xlab("Chromosome") +
        coord_flip(ylim = c(.1,4.5)) +
        annotate("text",
                 y=3.4,
                 x=4, 
                 label=glue::glue("Average SNPs Total: {total_snps} million", 
                              "\n",
                              "Average SNPs removed by QC: {snps_rm} million",
                              "\n",
                              "SNPs after meta-analysis: {snps_run} million",
                              "\n",
                              "Filtered genome size: {data_size} GB"),
                 size=8) +
        scale_fill_manual(values=c("grey50", "black"))
} 

x %>%
    unite(cohort_1, c(cohort_1, snps_removed_c1), sep="_") %>%
    unite(cohort_2, c(cohort_2, snps_removed_c2), sep="_") %>%
    gather(key, value, cohort_1:cohort_2) %>% 
    separate(value, c("time", "snps_removed"), sep="_") %>% 
    mutate_at(vars(time, snps_removed), as.double) %>% 
    filter(genome=="D", outcome=="TRM") %>%
    mutate(chr=as.factor(chr),
           key=ifelse(key=="cohort_1", "Cohort 1", "Cohort 2")) %>% 
    plot_maker(aes(fct_rev(chr), time/60, label=paste0(size, " MB"), fill=key),
               plot.title="100 Days: Mixed Donor TRM All Chromosomes Runtime (Filtered MAF > 0.005, INFO > 0.8)",
               position="dodge")

