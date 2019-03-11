# grabbing file size from CCR
library(gdata)
library(tidyverse)
trm_mixed_100d <- data.frame(chromosome=dir(), 
                             size=humanReadable(file.size(dir())), 
                             stringsAsFactors = FALSE)

trm_mixed_100d <- trm_mixed_100d %>%
    mutate(chromosome=str_replace(chromosome, ".res", "")) %>%
    separate(chromosome, c("chr", 
                           "genome", 
                           "disease", 
                           "outcome"), sep="_") %>%
    arrange(genome, chr) %>% 
    mutate(chr=as.double(chr),
           size=trimws(size)) %>%
    separate(size, c("size", "size_unit"), sep=" ") %>%
    mutate(size=as.double(size)) %>%
    as_tibble()


log_files <- dir(path="/projects/rpci/lsuchest/abbasriz/DBMT_100d/analyses/mixed/log",
                 pattern=".err", full.names=TRUE)

pull_times <- function(log_files){
    pulled_times <- system(sprintf("awk '$1 == \"Analysis\" {print $0}' %s", log_files), intern=TRUE)
    dates <- pulled_times %>%
        str_extract('[0-9]+-[0-9]+-[0-9]+') 
    times <- pulled_times %>%
        str_extract('[0-9]+:[0-9]+:[0-9]+') 
    start_finish <- as.POSIXct(paste(dates, times),
                               format="%Y-%m-%d %H:%M:%S")
    as.double(difftime(start_finish[2], start_finish[1],
                       units="min"))
}



run_times <- sapply(log_files, pull_times)

res <- data.frame(files=names(run_times),
                  mins=as.double(run_times),
                  stringsAsFactors = FALSE) %>%
    mutate(time_units="min")

res$files <- res$files %>%
    str_split(pattern="/") %>%
    map_chr(c(10, 1)) %>%
    str_replace(".err", "") 


res <- res %>%
    na.omit() %>%
    separate(files,
             c("jobid", "analysis"),
             sep="_",
             extra = "merge") %>%
    separate(analysis, c("chr",
                         "genome",
                         "disease",
                         "outcome"), 
             sep="_") %>%
    mutate(chr=as.double(chr),
           cohort=rep(1:2, nrow(.)/2)) %>%
    select(-jobid) %>%
    gather(key, value, cohort, -mins) %>% 
    unite(value, c(key, value), sep="_") %>%
    spread(value, mins) %>%
    as_tibble()


times_collected <- res %>%
    inner_join(trm_mixed_100d) %>%
    arrange(genome, chr)



# snps removed
library(tidyverse)
path <- "/projects/rpci/lsuchest/abbasriz/DBMT_100d/analyses/mixed/temp"
rm_files <- dir(path=path,
                pattern="snps_removed")

pull_snps_removed <- function(rm_files, path){
    as.double(str_extract(system(sprintf("wc -l %s/%s", path, rm_files), intern=TRUE), "\\S+"))
}


snps_removed <- sapply(rm_files, pull_snps_removed, path)
snps_removed_df <- data.frame(analysis=names(snps_removed), snps_removed=snps_removed, row.names=NULL, stringsAsFactors = FALSE)
snps_removed_df <- snps_removed_df %>%
    mutate(analysis=str_replace(analysis, ".snps_removed", "")) %>%
    separate(analysis, c("chr", "genome", "disease", "outcome", "cohort"), sep="_") %>%
    gather(key, value, cohort, -chr:-outcome) %>%
    spread(value, snps_removed) %>%
    select(-key) %>%
    rename(snps_removed_c1=c1,
           snps_removed_c2=c2)


path_analyzed <- "/projects/rpci/lsuchest/abbasriz/DBMT_100d/analyses/mixed/out"
res_files <- dir(path=path_analyzed,
                 pattern=".res")
grab_snps_analyzed <- function(res_files, path_analyzed){
    as.double(str_extract(system(sprintf("wc -l %s/%s", path_analyzed, res_files), intern=TRUE), "\\S+"))
}

snps_analyzed <- sapply(res_files, grab_snps_analyzed, path_analyzed)
snps_analyzed_df <- data.frame(analysis=names(snps_analyzed), snps_analyzed=snps_analyzed, row.names=NULL, stringsAsFactors = FALSE)
snps_analyzed_df <- snps_analyzed_df %>%
    mutate(analysis=str_replace(analysis, ".res", "")) %>%
    separate(analysis, c("chr", "genome", "disease", "outcome"))

snps_final <- snps_analyzed_df %>%
    left_join(snps_removed_df) %>%
    mutate(chr=as.double(chr)) %>%
    as_tibble()

final <- times_collected %>% 
    inner_join(snps_final)



