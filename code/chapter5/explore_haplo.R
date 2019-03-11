motifs <- haplo_sig %>% 
    select(query_snp_rsid, rsid, Motifs)

library(tidytext)
# motifs
# haplo_sig %>%
#     select(query_snp_rsid, rsid, Proteins) %>%  
#     unnest_tokens(proteins, Proteins, token="regex", pattern = ";") %>% 
#     na.omit() %>%
#     unnest_tokens(Proteins, proteins) %>%
#     count(Proteins) %>%
#     arrange(desc(n)) %>%
#     print(n=50)



library(glue)
library(tidyverse)
library(tidytext)

joint_res <- joint_res %>% full_join(regdb %>% rename(RSID=rsid))


joint_res_split <- joint_res %>% 
    filter(PVALUE_M < 5e-06, PVALUE_c1 < 0.05 , PVALUE_c2 < 0.05) %>%
    mutate(pt_subset=str_replace(pt_subset, "ALLsubtype_Bcell", "BALL")) %>%  
    separate(pt_subset, c("donor", "disease", "outcome", "interval"), sep="_")

e_codes <- read_tsv("~/Google Drive/OSU_PHD/DBMT_ALLonly/hits/ecode_epigenome.txt")
e_codes <- e_codes %>%
    mutate(`Standardised epigenome name`=str_replace_all(`Standardised epigenome name`, "\xca", " ")) %>%
    mutate_all(tolower) %>%
    mutate_all(trimws)
e_codes_vec <- e_codes$`Standardised epigenome name`
names(e_codes_vec) <- e_codes$EID

explore_haplo <- function(donor, disease, outcome, interval, haplo_col){
    df <- joint_res_split %>%
        filter(donor==!!donor,
               disease==!!disease,
               outcome==!!outcome,
               interval==!!interval)
    haplo_df <- haplo_sig %>% 
        filter(query_snp_rsid %in% df$RSID)
    
    switch(haplo_col,
       motifs={
           haplo_df %>%
               filter(Motifs != ".") %>%
               unnest_tokens(motifs, Motifs, token="regex", pattern = ";") %>% 
               na.omit() %>% 
               mutate_all(trimws) %>%
               count(motifs) %>%
               arrange(desc(n))
        },
        proteins={
            haplo_df %>%
                filter(Proteins != ".") %>%
                select(query_snp_rsid, rsid, Proteins) %>%
                unnest_tokens(proteins, Proteins, token="regex", pattern = ";") %>%
                na.omit() %>% 
                separate(proteins, c('cell lines', 'proteins', 'group', 'treatment'), sep=",") %>%
                unite(cell_line_protein, c('cell lines', "proteins"), sep=";") %>%
                mutate_all(trimws) %>%
                count(cell_line_protein) %>%
                arrange(desc(n)) %>%
                separate(cell_line_protein, c("cell_line", "protein"), sep=";")
        },
        chromatin_marks={
            haplo_df %>%
                select(query_snp_rsid, rsid, Chromatin_Marks) %>%
                unnest_tokens(chromatin_marks, Chromatin_Marks, token="regex", pattern = ";") %>%
                na.omit() %>%
                mutate_all(trimws) %>%
                count(chromatin_marks) %>%
                filter(!(chromatin_marks == "none"),
                       (chromatin_marks!= ".")) %>%
                arrange(desc(n)) %>%
                separate(chromatin_marks, c("eid", "chromatin_mark"), sep=",") %>%
                mutate(eid=e_codes_vec[eid])
        },
        promoter={
            haplo_df %>%
                select(query_snp_rsid, rsid, Promoter_histone_marks) %>%
                unnest_tokens(promoter_histone_marks, Promoter_histone_marks, token="regex", pattern = ";") %>%
                na.omit() %>%
                unnest_tokens(Promoter_histone_marks, promoter_histone_marks, token="regex", pattern = ",|\\n") %>%
                mutate_all(trimws) %>%
                count(Promoter_histone_marks) %>%
                filter(!(Promoter_histone_marks == "none"),
                       (Promoter_histone_marks != ".")) %>%
                arrange(desc(n))
        },
        enhancer={
            haplo_df %>%
                select(query_snp_rsid, rsid, Enhancer_histone_marks) %>%
                unnest_tokens(enhancer_histone_marks, Enhancer_histone_marks, token="regex", pattern = ";") %>%
                na.omit() %>%
                unnest_tokens(Enhancer_histone_marks, enhancer_histone_marks, token="regex", pattern = ",|\\n") %>%
                mutate_all(trimws) %>%
                count(Enhancer_histone_marks) %>%
                filter(!(Enhancer_histone_marks == "none"),
                       (Enhancer_histone_marks != ".")) %>%
                arrange(desc(n)) 
        },
       eqtl={
           haplo_df %>%
               select(query_snp_rsid, rsid, eQTL) %>%
               unnest_tokens(eqtl, eQTL, token="regex", pattern = ";") %>%
               filter(eqtl != ".") %>%
               na.omit() %>%
               separate(eqtl, c("gtex", "cell_lines", "gene", "pvalue"), sep=",") %>%
               select(-rsid, -pvalue, -gtex) %>%
               distinct() %>%
               group_by(cell_lines, gene) %>%
               summarize(n=n()) %>%
               ungroup() %>%
               arrange(desc(n))
       },
       gwas={
               haplo_df %>%
                   filter(gwas != ".") %>%
                   select(query_snp_rsid, r2, rsid, gwas)     
           },
       grasp={
           haplo_df %>%
               filter(grasp != ".") %>%
               select(query_snp_rsid, r2, rsid, grasp) %>% 
               mutate(grasp=ifelse(str_detect(grasp, ";"), grasp, paste0(grasp, ";"))) %>%
               unnest_tokens(Grasp, grasp, token="regex", pattern=";") %>%
               separate(Grasp, c("pmid", "grasp_trait", "grasp_pvalue"), sep=",") %>%
               mutate_at(vars(grasp_pvalue), as.double) %>%
               group_by(rsid, pmid, grasp_trait, grasp_pvalue) %>%
               summarize(query_snp_rsids=toString(query_snp_rsid)) %>%
               ungroup() %>%
               arrange(grasp_pvalue)
       }
    )
}
explore_haplo("R", "ALLonly", "DRM", "1y", "motifs") %>% 
    print(n=35)
explore_haplo("R", "ALLonly", "DRM", "1y", "proteins") %>% 
    print(n=20)
explore_haplo("R", "ALLonly", "DRM", "1y", "chromatin_marks")
explore_haplo("R", "ALLonly", "DRM", "1y", "promoter")
explore_haplo("R", "ALLonly", "DRM", "1y", "enhancer")

explore_haplo("R", "ALLonly", "DRM", "1y", "grasp") %>% 
    arrange(grasp_pvalue) %>%
    print(n=122)

explore_haplo("R", "ALLonly", "DRM", "1y", "eqtl") %>% 
    print(n=122)



all_opts <- expand.grid(genome=c("R","MM", "D"), 
            outcome=c("TRM", "DRM", "OS", "PFS", "GVHD", "OF", "INFX", "REL"),
            interval=c("100d", "1y", "3y"),
            stringsAsFactors = FALSE) %>% 
    mutate(disease="ALLonly") %>%
    filter(!(interval=="3y" & outcome %in% c("DRM", "GVHD", "INFX", "PFS", "TRM")),
           !(interval=="3y" & genome == "MM"),
           genome != "MM") %>%
    as_tibble()
    
    

    
explore_haplo("D", "ALLonly", "GVHD", "1y", "eqtl") %>% print(n=92)
# explore_haplo("R", "ALLonly", "TRM", "1y", "gwas") %>% print(n=30)

res <- vector('list', nrow(all_opts))
for(i in seq_len(nrow(all_opts))){

    
    res[[i]] <- explore_haplo(all_opts[i,]$genome, 
                              all_opts[i,]$disease,
                              all_opts[i,]$outcome,
                              all_opts[i,]$interval,
                              "gwas") %>%
        mutate(analysis=paste(c(all_opts[i,]$genome, 
                                all_opts[i,]$disease,
                                all_opts[i,]$outcome,
                                all_opts[i,]$interval),
                              collapse = "_"))

}

res <- do.call("rbind", res)

gwas_hits <- res %>% 
    na.omit() %>%
    group_by(rsid, gwas, analysis) %>%
    summarize(query_snp_rsid=toString(query_snp_rsid)) %>%
    ungroup() %>% 
    unnest_tokens(gwas, gwas, token="regex", pattern=";") %>%
    separate(gwas, c("pmid", "gwas_trait", "gwas_pvalue"), sep=",") %>%
    mutate_at(vars(gwas_pvalue), as.double) %>%
    select(rsid, pmid, gwas_trait, gwas_pvalue, analysis, query_snp_rsid) %>%
    filter(!str_detect(analysis, "3y"))

#saveRDS(gwas_hits, file="~/Google Drive/OSU_PHD/DBMT_ALLonly/hits/20190228_gwas_hits.rds")    
 
