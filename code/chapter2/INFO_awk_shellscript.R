#awk '{if ( $1  ~ /chr7-/ && $4 >87123179 && $4< 87352639) print "chr7" "\t" $4 "\t" "ABCB1" "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11}' Impute2.INFO  >> /projects/rpci/lsuchest/ezgi/candidate_gene/CG_test/info.test.txt



genes <- read.table("gene_locations.txt", header=T, stringsAsFactors = F)

#head(genes)

#awk <- c()
for(i in 1:nrow(genes)){
        cat("awk '{if ($1 ~ /",
            genes$seqnames[i],
            "-/ && $4 > ",
            genes$start_position[i],
            " && $4 < ",
            genes$end_position[i],
            ") print ",
            "\"",
            genes$seqnames[i],
            "\" \"\\t\" $4 \"\\t\" \"",
            genes$gene_symbol[i],
            "\"",
            " \"\\t\" $2",
            " \"\\t\" $3",
            " \"\\t\" $5",
            " \"\\t\" $6",
            " \"\\t\" $7",
            " \"\\t\" $8",
            " \"\\t\" $9",
            " \"\\t\" $10",
            " \"\\t\" $11}' ",
            "/projects/rpci/lsuchest/lsuchest/Rserve/ImputeData/var/db/gwas/imputed_data/BMT093013_forImpute/Impute2_summary/Impute2.INFO >> /projects/rpci/lsuchest/abbasriz/candidate_gene/info.test.txt",
            "\n",
            sep="")        
}


