---
title: "appendix"
author: "Abbas Rizvi"
date: "9/14/2018"
output: pdf_document
---





# Chapter 2

## Candidate Gene Analyses

We conducted an exhaustive literature search looking at papers from 71 candidate gene studies that were found via Pubmed search. The candidate gene studies focused on any blood disorders as long as one of the disorders was either ALL or AML. The candidate gene studies typically looked at several candidate SNPs in candidate genes and determined if there was an association with the SNP and survival outcomes (DRM, OS, PFS, or TRM). We gathered all of the SNPs that were studied in these candidate gene studies (even if they were not significant in the original study). 

 The final table can be found on UB CCR at:    
 `/projects/rpci/lsuchest/abbasriz/candidate_gene/ \`  
`result_files_cg/final_table/cg_snptable.txt`.  

 The ultimate goal of this study was to test all of the candidate SNPs in DISCOVeRY-BMT to replicate or validate the original studies' findings. Replication means that the phenotype and population are the same. Validation means that the phenotype is the same but the study population is not necessarily the same (i.e. different ethnic group). And also, because these studies essentially had gene based hypotheses, we sought to test DISCOVeRY-BMT survival data using gene-based statistical approaches (VEGAS2 software) to test whether the aggregation of the entire gene locus was significant or not. 

 Here we present our data analysis.


\scriptsize


```bash
# grab gene list in unix shell from final table and unique it
[abbasriz@rush:]$ </projects/rpci/lsuchest/abbasriz/candidate_gene/ \
    result_files_cg/final_table/cg_snptable.txt \
    tr -s ' ' ' ' | \
    awk '{print $1}' | \
    uniq > gene_list.txt
```

 \normalsize



\scriptsize


```r
# load gene list
# skip column header "gene" 
# this list has genes from chrX already removed
genes.w.locs <- read.table("/projects/rpci/lsuchest/lsuchest/
                            CandidateGeneReplication/
                           final_gene_list/
                           gene_locations_20170222.txt",
                           header=TRUE,
                           stringsAsFactors = FALSE)
# order by chromosome and sort by numeric part of the 
# character vector so they are ordered correctly
vals <- as.numeric(gsub("chr", "", genes.w.locs$seqnames))
genes.w.locs <- genes.w.locs[order(vals),]
head(genes.w.locs)
```

 \normalsize


## Final Candidate Gene List

\scriptsize


```r
gene.list <- as.character(sort(genes.w.locs$gene_symbol))
write.table(gene.list, 
            file="final_cg_genelist.txt", 
            quote=FALSE, sep="\t",
            row.names=FALSE,
            col.names=FALSE)
```

 \normalsize



 We need to grab the typed and imputed SNPs from DISCOVeRY-BMT. These are located in:    
 `/projects/rpci/lsuchest/lsuchest/Rserve/ImputeData/var/db/gwas/ \`  
`imputed_data/BMT093013_forImpute/Impute2_summary/Impute2.INFO`   

 file. 


\scriptsize


```bash
[abbasriz@rush:]$ head /projects/rpci/lsuchest/lsuchest/Rserve/ImputeData/var/db/gwas/ \
imputed_data/BMT093013_forImpute/Impute2_summary/Impute2.INFO
region snp_id rs_id position exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0
chr10-0-5000000 --- rs148087467 60523 0.002 0.686 0.998 0 -1 -1 -1
chr10-0-5000000 --- rs187110906 60969 0.092 0.499 0.891 0 -1 -1 -1
chr10-0-5000000 --- rs192025213 61005 0.003 0.371 0.996 0 -1 -1 -1
chr10-0-5000000 --- rs115033199 61020 0.000 0.341 0.999 0 -1 -1 -1
chr10-0-5000000 --- rs183305313 61334 0.005 0.302 0.990 0 -1 -1 -1
chr10-0-5000000 --- rs186558141 65978 0.000 0.114 1.000 0 -1 -1 -1
chr10-0-5000000 --- rs190079063 66269 0.000 0.137 1.000 0 -1 -1 -1
chr10-0-5000000 --- rs12260013 66326 0.030 0.649 0.972 0 -1 -1 -1
chr10-0-5000000 --- chr10:66627:D 66627 0.548 0.590 0.734 0 -1 -1 -1
```

 \normalsize


 Now that we have the gene locations, we are going to create a shell script that will be able to grab the `Impute2.INFO` file regions and let us collect all of the SNPs that we have from DISCOVeRY-BMT that are typed and imputed. The shell script can be found at:  

 `/projects/rpci/lsuchest/lsuchest/CandidateGeneReplication/ \`  
`parse_Impute2INFO/awk_Impute2.INFO_by_geneLoc_cg.sh`

This shell script rearranges the columns to the following: `chr`, `position`, `gene`, `snp_id`, `rs_id`, `exp_freq_q1`, `info`, `certainty`, `type`, `info_type0`, `concord_type0`, `r2_type0`. 


\scriptsize


```r
capture.output(file="/projects/rpci/lsuchest/lsuchest/CandidateGeneReplication/
               parse_Impute2INFO/awk_Impute2.INFO_by_geneLoc_cg.sh",{
for(i in 1:nrow(genes.w.locs)){
        cat("awk '{if ($1 ~ /",
            genes.w.locs$seqnames[i],
            "-/ && $4 > ",
            genes.w.locs$start_position[i],
            " && $4 < ",
            genes.w.locs$end_position[i],
            ") print ",
            "\"",
            genes.w.locs$seqnames[i],
            "\" \"\\t\" $4 \"\\t\" \"",
            genes.w.locs$gene_symbol[i],
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
            "/projects/rpci/lsuchest/lsuchest/Rserve/ImputeData/var/db/
            gwas/imputed_data/BMT093013_forImpute/Impute2_summary/Impute2.INFO >>
            /projects/rpci/lsuchest/abbasriz/candidate_gene/info.cg.txt",
            "\n",
            sep="")        
}
})
```

 \normalsize



 We edited this to make this into a SLURM script would and after running we have a file `info.cg.txt`. 


\scriptsize


```bash
[abbasriz@rush:/projects/rpci/lsuchest/abbasriz/candidate_gene]$ head -3 info.cg.txt
chr1	70866967	CTH	---	rs77482262	0.069	0.984	0.997	0	-1	-1	-1
chr1	70866988	CTH	---	rs187366946	0.002	0.548	0.997	0	-1	-1	-1
chr1	70867130	CTH	---	rs190937437	0.001	0.208	0.999	0	-1	-1	-1
```

 \normalsize


 `info.cg.txt` does not have a header so are going to add the header in R and re-write the file (alternatively, we could have just assigned this in each of the survival parsing R scripts that we are about to make.)  


\scriptsize


```r
library(data.table)
info <- fread("/projects/rpci/lsuchest/abbasriz/candidate_gene/info.cg.txt")
colnames(info) <- c("chr", "position", "gene", "snp_id", "rs_id", 
                    "exp_freq_a1", "info", "certainty", "type", 
                    "info_type0", "concord_type0", "r2_type0")
write.table(info, file="/projects/rpci/lsuchest/abbasriz/candidate_gene/info.cg.txt", 
            quote=FALSE, 
            col.names=TRUE,
            row.names=FALSE, 
            sep="\t")
```

 \normalsize

 

## Survival Results Directories
Individual (donor/recipient genotyping cohorts 1 and 2) and shared (mismatch between donor-recipient pairs from cohorts 1 and 2), done for 4 different outcomes (DRM, PFS, OS, TRM) and 4 different disease groups (AMLonly, ALLonly, mixed, noALL), and their corresponding meta analyses (combining cohort 1 and 2 using a fixed effect model from METAL software) are located in different directories on UB CCR. 3 genomes (donor, recipient, shared) x 3 cohorts (c1, c2, meta) x 4 outcomes (DRM, PFS, OS, TRM) x 4 diseases (AML, ALL, mixed, noALL) = 144 analyses.

 Individual directory:     
 `/projects/rpci/lsuchest/lsuchest/Rserve/ImputeData/var/db/gwas/ \`  
`imputed_data/BMT093013_forImpute/analyses/`  
 Meta-individual directory:    
 `/projects/rpci/lsuchest/lsuchest/Rserve/ImputeData/var/db/gwas/ \`  
`imputed_data/BMT093013_forImpute/analyses/METAL.results/`  
 Shared directory:    
 `/projects/rpci/lsuchest/lsuchest/Rserve/ImputeData/var/db/gwas/ \`  
`imputed_data/SHARED/analyses/METAL.RESULTS.SHARED/`  
 Meta-shared directory:    
 `/projects/rpci/lsuchest/lsuchest/Rserve/ImputeData/var/db/gwas/ \`    
`imputed_data/BMT093013_forImpute/analyses/METAL.results/`  

 As a collective, these directories contain files of 144 analyses. However, the result files are not so well organized. Using some unix commands to capture a clean amount of directories and files, as well as *ad hoc* manual curation, final directory lists (`/projects/rpci/lsuchest/lsuchest/CandidateGeneReplication/ \`  
 `survival_results_directories/`) can be found in the following files: `ind.directories.txt`, `shared.directories.txt`, `meta.ind.directories.txt`, and `meta.shared.directories.txt`.

### Parsing Result Files for High Quality Candidate Gene SNPs
We have survival results located in these directories, now we want to subset each of these survival results for just all of the high quality SNPs in our candidate genes. Here I will show only 1 of 4 (`independent_results.R`) parsing results. The others are `meta_independent.results.R`, `shared_results.R`, and `meta_results.R`. I split these up to decrease computational time and to have some safety checks at a smaller scale. 

 The result files go into the directory: `/projects/rpci/lsuchest/abbasriz/\`
`candidate_gene/result_files_cg/impute.results.w.typedsnps`  

 The result files will be in the following format:  
`genome_cohort_outcome_disease.txt`, e.g. (`D_c1_DRM_ALLonly.txt`, for donor, cohort 1, death to due to disease, ALL only subset)



\scriptsize


```r
## INDEPENDENT RESULTS
library(data.table)
# read in candidate gene
info <- fread("/projects/rpci/lsuchest/abbasriz/candidate_gene/info.cg.txt",
              header='auto', sep="\t")
#colnames(info) <- c("chr", "position", "gene", "snp_id", "rs_id", 
# "exp_freq_a1", "info", "certainty", "type", "info_type0", 
# "concord_type0", "r2_type0")
setkey(info, "gene","snp_id", "rs_id")
info <- info[,c("gene", "snp_id", "rs_id")]
setkey(info, "snp_id","rs_id")
ind <- scan("/projects/rpci/lsuchest/abbasriz/candidate_gene/
            result_files_cg/res.directories/ind.directories.txt",
            what=character())
for(i in 1:length(ind)){
	ind.res <- fread(ind[i], header="auto", sep="\t", verbose=TRUE)
	colnames(ind.res)[1:2] <- c("snp_id", "rs_id")
	setkey(ind.res, "snp_id", "rs_id")
	ind.res <- ind.res[info]
	ind.res <- na.omit(ind.res)
	ind.res[,c("z", "loglik0", "loglik"):=NULL]
	ind.res[,`95%-CI`:=NA]
	setcolorder(ind.res, c("gene", "rs_id", "CHR", "BP", "ALLELE1",
	                       "ALLELE2", "n", "coef", "se(coef)", "exp(coef)",
	                       '95%-CI', "snp_id", "Pr(>|z|)"))
	colnames(ind.res) <- c("gene", "rsID", "chr", "BP", "allele1", "allele2", "N",
	                       "coef", "se.coef",  "exp.coef", "95%-CI",
	                       "impute", "Pvalue")
	file.name <- strsplit(ind[i], split="/")[[1]][14]
	file.name <- gsub("[.]", "_", file.name)
	file.name <- strsplit(file.name, "_")[[1]]
	if(length(file.name)<6){
		file.name[6] <- "mixed"
	}	else{
			file.name <- file.name
	}
	write.table(ind.res, file=paste0("/projects/rpci/lsuchest/abbasriz/candidate_gene/
	                                 results_files_cg/impute.results.w.typedsnps/",
	                                 paste(file.name[1],
	                                       file.name[2],
	                                       file.name[4],
	                                       file.name[6],
	                                       sep="_"),
	                                 ".txt"),
	            quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
	rm(ind.res)
}
```

 \normalsize

 
 \textbf{Note}: Don't run the following shell script without editing the R script above (and similar ones) on UB CCR, as the directory may not be correct. The following shell script (`/projects/rpci/lsuchest/abbasriz/candidate_gene/ \`  
 `result_files_cg/res.directories/independent_results.sh`) was written to run this command using SLURM:  



\scriptsize


```bash
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=5000
#SBATCH --job-name=myjobname
#SBATCH --output=myjob.out
#SBATCH --error=myjoberr.err
#SBATCH --partition=general-compute
#SBATCH --mail-user=abbasriz@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END

#Get date and time
tstart=`date`
echo "###### start time:"$tstart
cd /projects/rpci/lsuchest/abbasriz/candidate_gene/\
    result_files_cg/res.directories/
echo "Run program"
module load R 
R CMD BATCH independent_results.R
echo "program finished"
echo "All Done!"
tend=`date`
echo "###### end time: $tend" 
```

 \normalsize



 Again, the 144 result files can be found in: `/projects/rpci/lsuchest/abbasriz/\`  
`candidate_gene/result_files_cg/impute.results.w.typedsnps/`

## Remove Duplicates 
METAL and VEGAS2 ignore duplicates rsids. In DISCOVeRY-BMT often has duplicates rsids due to some SNPs being typed AND imputed. We filtered these files to removed the typed SNPs, as the imputed SNP may be more reliable as it is calculated using the reference genome. Again, we separated this into 3 jobs, divvying up the results by donor, recipient, and shared. Instead of using two columns to remove duplicates, we removed duplicates by searching for duplicates and removing the one with lower base pair (i.e. position 111110 would be typed and position 111111 would be imputed, so we would remove position 111110 in this case). Here we demonstrate this using only the donor files (donor_remove_dups.R and donor_remove_dups.sh)



\scriptsize


```r
# donor_remove_dups.R
library(data.table)
donor.files <- list.files(pattern="^D_")
for(i in 1:length(donor.files)){
	donor <- fread(donor.files[i], header="auto", sep="\t")
	setkey(donor, chr, BP, rsID)
	donor <- donor[!which(duplicated(donor$rsID))-1]
	write.table(donor, file=paste0("/projects/rpci/lsuchest/abbasriz/
	                               candidate_gene/result_files_cg/", 
	                               donor.files[i]),
	            quote=F, sep="\t", col.names=T, row.names=F)
	rm(donor)
}
```

 \normalsize


\scriptsize


```bash
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem=5000
#SBATCH --job-name=myjobname
#SBATCH --output=myjob.out
#SBATCH --error=myjoberr.err
#SBATCH --partition=general-compute
#SBATCH --mail-user=abbasriz@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END

#Get date and time
tstart=`date`
echo "###### start time:"$tstart
cd /projects/rpci/lsuchest/abbasriz/candidate_gene/\
    result_files_cg/impute.results.w.typedsnps/
echo "Run program"
module load R 
R CMD BATCH donor_remove_dups.R
echo "program finished"
echo "All Done!"
tend=`date`
echo "###### end time: $tend"  
```

 \normalsize


We did this for recipient (`recipient_remove_dups.R` and `recipient_remove_`  
 `dups.sh`) and shared (`shared_remove_dups.R` and `shared_remove_dups.sh`) genomes as well. Again, the donor, recipient, and shared result files and their corresponding shell scripts are in `/projects/rpci/lsuchest/ \`
 `abbasriz/candidate_gene/result_files_cg/impute.results.w.typedsnps`. The output for this analysis can be found in: `/projects/rpci/lsuchest/ \`
 `abbasriz/candidate_gene/result_files_cg/`, however, the files in these directories are the most current files, and they have been updated from this point in the analysis and include meta results, hazard ratios, and 95% confidence intervals for the HRs. 

## Meta-Analysis
So since the meta analysis files didn't have hazard ratio coefficients and confidence intervals, we had to re-run the metal analysis using `STDERR` option to get these results. The files from the meta-analysis can be found in `/projects/rpci/lsuchest/abbasriz/candidate_gene/\`  
 `result_files_cg/metal_results/` These analyses are run pairwise (cohort 1 and cohort 2), e.g. `/projects/rpci/lsuchest/ \`   
 `abbasriz/candidate_gene/result_files_cg/D_c1_DRM_ALLonly.txt` and `/projects/rpci/ \`  
 `lsuchest/abbasriz/candidate_gene/result_files_cg/D_c2_DRM_ALLonly.txt`.

A metal script can be run using METAL software:  

## Example of METAL script (e.g. metal.txt)

\scriptsize


```bash
module load metal 
cat metal.txt
# META ANALYSIS FOR COHORT 1 and COHORT 2 FROM DISCOVERY-BMT 
# CANDIDATE GENE REPLICATION/VALIDATION SUBSET
# THE RESULTS ARE STORED IN FILES metal_D_M_DRM_ALLonly.tbl 
# and metal_D_M_DRM_ALLonly.tbl.info
SCHEME STDERR
# LOAD COHORT 1 and COHORT 2 FILES
# === DESCRIBE AND PROCESS THE FIRST INPUT FILE ===
MARKER rsID
ALLELE allele1 allele2
EFFECT coef
STDERR se.coef
PVALUE Pvalue
WEIGHT N
PROCESS /projects/rpci/lsuchest/abbasriz/candidate_gene/\
    result_files_cg/D_c1_DRM_ALLonly.txt
# === THE SECOND INPUT FILE HAS THE SAME FORMAT AND CAN BE PROCESSED IMMEDIATELY ===
PROCESS /projects/rpci/lsuchest/abbasriz/candidate_gene/\
    result_files_cg/D_c2_DRM_ALLonly.txt

OUTFILE /projects/rpci/lsuchest/abbasriz/candidate_gene/\
    result_files_cg/metal_results/metal_D_M_DRM_ALLonly .tbl
MINWEIGHT 10
```

 \normalsize



## Python script to generate metal files

To automate this for our cohorts we wrote a Python script to generate all of these different pairwise metal scripts. 


\scriptsize


```python
#!/usr/bin/python
import glob
import io
# grab file names
files_c1=sorted(glob.glob("/projects/rpci/lsuchest/abbasriz/candidate_gene/result_files_cg/*c1*.txt"))
files_c2=sorted(glob.glob("/projects/rpci/lsuchest/abbasriz/candidate_gene/result_files_cg/*c2*.txt"))
# set up file names
meta_file_names = []
for i in range(len(files_c1)):
	last_item = files_c1[i].split("/")[-1]
	last_item = last_item.replace("c1", "M")
	last_item = last_item.replace(".txt", "")
	meta_file_names.append(last_item)
# generate metal scripts
#for i in meta_file_names:
for i in range(len(meta_file_names)):
    with open("run_metal_{}.txt".format(meta_file_names[i]), "w") as f:
    	f.write("\nSCHEME STDERR")
    	f.write("\n# LOAD COHORT 1 and COHORT 2 FILES")
    	f.write("\n# === DESCRIBE AND PROCESS THE FIRST INPUT FILE ===")
	f.write("\nMARKER rsID")
	f.write("\nALLELE allele1 allele2")
	f.write("\nEFFECT coef")
	f.write("\nSTDERR se.coef")
	f.write("\nPVALUE Pvalue")
	f.write("\nWEIGHT N")
	f.write("\nPROCESS "+ files_c1[i])
	f.write("\n# === THE SECOND INPUT FILE HAS THE SAME FORMAT AND CAN BE PROCESSED IMMEDIATELY ===")
	f.write("\nPROCESS " + files_c2[i])
	f.write("\n")
	f.write("\nOUTFILE /projects/rpci/lsuchest/abbasriz/",
	"candidate_gene/result_files_cg/metal_results/metal_"+meta_file_names[i] + " .tbl")
	f.write("\nMINWEIGHT 10")
	f.write("\nANALYZE")
```

 \normalsize


 Another python script was written to create the shell script to run these METAL scripts


\scriptsize


```python
#!/usr/bin/python
import glob
import io
# grab file names
files = sorted(glob.glob("/projects/rpci/lsuchest/abbasriz/",
            "candidate_gene/result_files_cg/metal_results/run_*.txt"))
f = open('run_metal.sh', 'w')
f.write("#!/bin/bash")
f.write("\n#SBATCH --time=24:00:00")
f.write("\n#SBATCH --nodes=14")
f.write("\n#SBATCH --ntasks-per-node=4")
f.write("\n#SBATCH --mem-per-cpu=5000")
f.write("\n#SBATCH --job-name=myjobname")
f.write("\n#SBATCH --output=myjob.out")
f.write("\n#SBATCH --error=myjoberr.err")
f.write("\n#SBATCH --partition=general-compute")
f.write("\n#SBATCH --mail-user=abbariz@buffalo.edu")
f.write("\n#SBATCH --mail-type=ALL")
f.write("\n#SBATCH --mail-type=END")
f.write("\n")
f.write("\n#Get date and time")
f.write("\ntstart=`date`")
f.write("\necho \"###### start time:\"$tstart")
f.write("\nmodule load python")
f.write('\nmodule load metal')
f.write("\n")
for i in range(len(files)):
	f.write("metal " + files[i] + "\n")
f.write("\n")
f.write("program finished")
f.write("All done!")
f.write("tend=`date`")
f.write("echo \"###### end time:\"$tend")
f.write("exit")
f.close()
```

 \normalsize


 The results are stored in `.tbl` files in `/projects/rpci/lsuchest/abbasriz/ \`   
`candidate_gene/metal_results/`.

## Merge meta-analys results back into original files

Now that we had the meta-analysis results, we merged the beta coefficients and standard errors, back into the `._M_$'` files so that we could have the same format. The script that does this can be found at:  
`/projects/rpci/lsuchest/abbasriz/candidate_gene/ \`  
 `metal_results/merge_metal_res_w_impute_res.R`


\scriptsize


```r
# merge_metal_res_w_impute_res.R
library(data.table)
path <- "/projects/rpci/lsuchest/abbasriz/candidate_gene/metal_results/"
metal.res.files <- list.files(path=path, pattern=".tbl$")
metal.files <- lapply(metal.res.files, fread)
cnames <- c("rsID", "Allele1", "Allele2", 
            "coef", "se.coef", "Pvalue", "Direction")
metal.files <- lapply(metal.files, setNames, cnames)
metal.res.files <- list.files(pattern=".tbl$")
res.files <- list.files(pattern="._M_", 
                        path="/projects/rpci/lsuchest/abbasriz/candidate_gene/
                        result_files_cg/")
res.files <- paste0("/projects/rpci/lsuchest/abbasriz/
                    candidate_gene/result_files_cg/", res.files)
impute.res <- lapply(res.files, fread)
metal.files <- lapply(metal.files, 
                      function(x) setkey(x, rsID))
metal.files <- lapply(metal.files,
                      function(x) x[,c("Allele1","Allele2","Direction"):=NULL,])
metal.files <- lapply(metal.files, 
                      function(x) x[,exp.coef:=exp(coef),])
impute.res <- lapply(impute.res, function(x) setkey(x, rsID))
impute.res <- lapply(impute.res,  function(x) x[,c("coef", 
                                                   "se.coef",
                                                   "exp.coef", 
                                                   "Pvalue"):=NULL,])

for(i in 1:length(impute.res)){
	impute.res[[i]] <- impute.res[[i]][metal.files[[i]]]
}
impute.res <- lapply(impute.res, setcolorder, c("gene", "rsID", "chr", "BP",
                                                "allele1", "allele2", "N", 
                                                "coef", "se.coef", "exp.coef",
                                                "95%-CI", "impute", "Pvalue"))
res.files <- list.files(pattern="._M_",
                        path="/projects/rpci/lsuchest/
                        abbasriz/candidate_gene/result_files_cg/")
names(impute.res) <- res.files
res.files <- paste0("/projects/rpci/lsuchest/abbasriz/
                    candidate_gene/result_files_cg/", res.files)
for(i in 1:length(res.files)){
	write.table(impute.res[[i]], 
	            file=res.files[[i]], sep="\t",
	            row.names=F, col.names=T, quote=F)
}
```

 \normalsize

## Calculate Hazard Ratios and 95% CI


Now that all the files are in the same format, we can automate the way that we calculate the hazard ratios and HR 95 CI, as all files should be the same for all analyses. Here we used the R script: `/projects/rpci/lsuchest/abbasriz/ \`   
 `candidate_gene/result_files_cg/calc_hazardratios.R`. This R script also ensures both allele1 and allele2 are uppercase and reassigns the rows in the `impute` column to `typed` (repeated rsid) and `imputed` (---).


\scriptsize


```r
library(data.table)
path <- "/projects/rpci/lsuchest/abbasriz/candidate_gene/result_files_cg/"
files <- list.files(path=path, pattern=".txt")
res <- lapply(files, fread)
# 95% CI of hazard ratios
hr.ci <- function(coef.est, se){
	lb <- round(exp(coef.est-1.96*se), 5)
	ub <- round(exp(coef.est+1.96*se), 5)
	paste0("[", lb,", " , ub ,"]")
}
# make new column of HR in lists
res <- lapply(res, function(x) x <- x[,'95%-CI':= hr.ci(x$coef, x$se.coef)])
# change the imputation notation to 'imputed' or 'typed'
impute <- function(x){
        setkey(x, impute)
        x[impute != "---", impute := "typed"]
        x[impute == "---", impute := "imputed"]
}
res <- lapply(res, impute)
# lets make alleles all upper case too
res <- lapply(res, function(x) x <- x[,'allele1':=toupper(allele1),])
res <- lapply(res, function(x) x <- x[,'allele2':=toupper(allele2),])
res <- lapply(res, setkey, gene)
names(res) <- files
res.files <- paste0("/projects/rpci/lsuchest/abbasriz/candidate_gene/result_files_cg/", files)
# sort by gene
res <- lapply(res, setkey, gene)
# write to file
for(i in 1:length(res.files)){
	write.table(res[[i]], file=res.files[[i]],
	            sep="\t", row.names=F, col.names=T, quote=F)
}
```

 \normalsize

## Create genome, cohort, outcome, disease columns

The script can be found at `/projects/rpci/lsuchest/abbasriz/ \`   
 `candidate_gene/result_files_cg/descript_columns.R`



\scriptsize


```r
library(data.table)
# subsetted impute results w/ individual cohort and meta beta estimates/hr/CIs merged
path.cg <- "/projects/rpci/lsuchest/abbasriz/candidate_gene/result_files_cg/"
files.cg <- list.files(path=path.cg, pattern=".txt")[-1]
impute.res <- lapply(paste0(path.cg,files.cg), fread)

# add the columns that we want, a column for genome, outcome, and cohort
# grab that information from the file names
genome <- gsub("_.*$", "", files.cg)
disease <- sub("(.*\\_)([^.]+)(\\.[[:alnum:]]+$)", "\\2", files.cg)
cohort <- sapply(strsplit(files.cg, '_', fixed=T), "[", 2)
outcome <- sapply(strsplit(files.cg, '_', fixed=T), "[", 3)

# creating columns with that info
mapply(function(x, cat) x <- x[,'genome' := cat], impute.res, genome)
mapply(function(x, cat) x <- x[,'cohort' := cat], impute.res, cohort)
mapply(function(x, cat) x <- x[,'outcome' := cat], impute.res, outcome)
mapply(function(x, cat) x <- x[,'disease' := cat], impute.res, disease)

impute.res <- lapply(impute.res, setkey, gene)
impute.res <- do.call(rbind, impute.res)
impute.res$gene[which(impute.res$gene == "GSTM1")] <- "GSTM"
impute.res <- na.omit(impute.res)

write.table(impute.res, file="/projects/rpci/lsuchest/abbasriz/
            candidate_gene/result_files_cg/final_table/
            cg_snptable.txt", col.names=T,
            row.names=F, quote=F, sep="\t")
```

 \normalsize



## VEGAS2
To run VEGAS2 you need two unlabeled columns: 1.) RSID and 2.) P-values. We created a shell script through R to do this, which grabs the 2nd and 13th column from our phenotype subsets.


\scriptsize


```r
path <- "/projects/rpci/lsuchest/abbasriz/candidate_gene/result_files_cg/"
x <- list.files(path=path, pattern=".txt")
capture.output(file="v2_parse.sh", {
for (i in 1:length(x)){
	cat("awk '{if (NR!=1) {print $2 \"\\t\" $13}}'",
	    x[i],
	    ">",
	    paste0("/projects/rpci/lsuchest/abbasriz/candidate_gene/
	           result_files_cg/vegas2/v2_",
	           x[i]),
	           "\n",
	           sep=" ")
}
})
```

 \normalsize


 The shell script can be found at: `/projects/rpci/lsuchest/abbasriz/ \`  
 `candidate_gene/result_files_cg/vegas2/v2_parse.sh` The results look like the following:  


\scriptsize


```bash
[abbasriz@rush:/projects/rpci/lsuchest/abbasriz/\
candidate_gene/result_files_cg/vegas2]$ head v2_D_c1_DRM_ALLonly.txt 
rs6696752	0.827182998293298
rs72638700	0.874653327370856
rs6659287	0.526921155628948
rs141567582	0.874653327370856
rs6699669	0.262962565390368
rs6659541	0.627228150206333
rs12132479	0.890715829741203
rs6699881	0.829409938995696
rs41275458	0.874758581607631
rs7538516	0.756659623949985
```

 \normalsize


 The vegas2 results were split into 3 shell scripts (`run_vegas2_donor.sh`, `run_vegas2_recipient.sh`, `run_vegas2_shared.sh`). To create these files we used perl.


\scriptsize


```perl
ls *.txt | grep '_D_' | perl \
    -lane '$new=$_; $new =~ s/\.txt$//; print "vegas2 $_ \
    -pop 1000GEURO \
    -subpop EURO \
    -genesize 10kbloc \
    -top 100 \
    -sex BothMnF -max 1000000 \
    -out ./results/$new.V2out "'  > run_vegas2.donor.sh
ls *.txt | grep '_R_' | perl \
    -lane '$new=$_; $new =~ s/\.txt$//; print "vegas2 $_ \
    -pop 1000GEURO \
    -subpop EURO \
    -genesize 10kbloc \
    -top 100 \
    -sex BothMnF \
    -max 1000000 \
    -out ./results/$new.V2out "' > run_vegas2_recipient.sh
ls *.txt | grep '_S_' | perl \
    -lane '$new=$_; $new =~ s/\.txt$//; print "vegas2 $_ \
    -pop 1000GEURO \
    -subpop EURO \
    -genesize 10kbloc \
    -top 100 \
    -sex BothMnF \
    -max 1000000 \
    -out ./results/$new.V2out "' > run_vegas2_shared.sh
```

 \normalsize


 These were put into SLURM scripts and the following lines were run in the programs



\scriptsize


```bash
module load R
module load plink
module load vegas2
```

 \normalsize


 An example of the arguments that were used for vegas2 produced from the perl script can be seen here:  


\scriptsize


```bash
[abbasriz@rush:/projects/rpci/lsuchest/abbasriz/candidate_gene/vegas2]$ vegas2 v2_D_DRM_ALLonly.txt /
    -pop 1000GEURO /
    -subpop EURO / 
    -genesize 10kbloc /
    -top 100 / 
    -sex BothMnF / 
    -max 1000000 / 
    -out /projects/rpci/lsuchest/abbasriz/candidate_gene/vegas2/results/v2_D_DRM_ALLonly.V2out
```

 \normalsize


 VEGAS2 results are stored in: `/projects/rpci/lsuchest/abbasriz/ \`  
 `candidate_gene/result_files_cg/vegas2/results` with the extension `.V2out`.

## Finalizing VEGAS2 results


\scriptsize


```r
# candidate gene list
gene.list <- scan("gene_list.txt", what=character())
gene.list <- c(gene.list, "GSTM1", "GSTM2", "GSTM3", "GSTM4", "GSTM5")
gene.list[gene.list=="GSTM"] <- NA
gene.list <- gene.list[complete.cases(gene.list)]
# vegas2 results
path.v2 <- "/projects/rpci/lsuchest/abbasriz/candidate_gene/result_files_cg/vegas2/results/"
files.v2 <- list.files(path=path.v2, pattern=".V2out")
v2.res <- lapply(paste0(path.v2, files.v2), fread)
v2.res <- lapply(v2.res, setkey, Gene)
v2.res <- lapply(v2.res, function(x) x[gene.list])
v2.res <- lapply(v2.res, setNames, c("chr", "gene",
                                     "nSNPs", "nSims", "start", "stop",
                                     "Test", "geneBased.pval", "topSNP", 
                                     "topSNP.pval"))
v2.res <- lapply(v2.res, setkey, gene)
# parse the file names to get the info that we want for the columns
genome <- sapply(strsplit(files.v2, '_', fixed=T), "[", 2)
outcome <- sapply(strsplit(files.v2, '_', fixed=T), "[", 4)
cohort <- sapply(strsplit(files.v2, '_', fixed=T), "[", 3)
disease <- sapply(strsplit(files.v2, '_', fixed=T), "[",5)
disease <- sapply(strsplit(disease, '.', fixed=T), "[",1)

# create descriptive columns
for(i in 1:length(v2.res)){
	v2.res[[i]][,'genome':=genome[i]]
	v2.res[[i]][,'cohort':=cohort[i]]
	v2.res[[i]][,'outcome':=outcome[i]]
	v2.res[[i]][,'disease':=disease[i]]
}

v2.res <- do.call(rbind, v2.res)
head(v2.res)
write.table(v2.res, file="/projects/rpci/lsuchest/abbasriz/candidate_gene/
            result_files_cg/final_table/cg_v2table.txt", 
            col.names=T, row.names=F,
            quote=F, sep="\t")
```

 \normalsize


\scriptsize



 \normalsize

## Top SNP and gene Based Pvalue adjustment

Grab number of typed SNPs from just cohort 1 (because it's the largest cohort). The script to count the typed SNPs is called `typed_snps.R` and `typed_snps.sh`. 


\scriptsize


```r
#typed snps
library(data.table)
impute.res <- fread("/projects/rpci/lsuchest/abbasriz/candidate_gene/
                    result_files_cg/final_table/cg_snptable.txt")

## normalize top snp pvalue by number of typed snps in cohort 1
p <- c()
x <- c()
typed.snps <- unique(impute.res[,c("gene", "genome", 
                                   "disease", "outcome"),with=F])

# my code
for(i in seq(nrow(typed.snps))){
                x <- impute.res[gene==rr$gene[i] & disease==rr$disease[i] & 
                                    genome==rr$genome[i] &
                                    outcome==rr$outcome[i] & cohort=="c1"]                
                setkey(x, "impute")
                p[i] <- nrow(x[impute=="typed"])
}

typed.snps <- typed.snps[,no.typed.snps:=p]
head(typed.snps)
write.table(typed.snps,
            file="typed.snps.txt", 
            col.names=T, row.names=F,
            sep="\t", quote=F)
```

 \normalsize

\scriptsize



 \normalsize


 We adjusted the topSNP pvalues by the number of typed SNPs in that gene (in cohort 1). We also adjusted the gene-based pvalue by the total number of genes that we tested in the candidate gene replication study. 


\scriptsize


```r
library(data.table)
# vegas 2 results
v2.res <- fread("/projects/rpci/lsuchest/abbasriz/candidate_gene/
                result_files_cg/final_table/cg_v2table.txt")
# number of typed snps in each gene
typed.snps <- fread("/projects/rpci/lsuchest/abbasriz/candidate_gene/
                    result_files_cg/final_table/typed.snps.txt")

setkey(typed.snps, gene, genome, disease, outcome)
setkey(v2.res, gene, genome, disease, outcome)
v2.typed.snp <- v2.res[typed.snps]
# adjust for number of typed snps
v2.typed.snp <- v2.typed.snp[,"adj.topSNP.pval":=topSNP.pval*no.typed.snps]
# adjust for number of genes
v2.typed.snps <- v2.typed.snp[,"adj.geneBased.pval":=geneBased.pval*174]
setcolorder(v2.typed.snp, c("chr", "gene", "nSNPs", "nSims","start", "stop",
                            "Test", "geneBased.pval", "adj.geneBased.pval", 
                            "topSNP", "topSNP.pval", "adj.topSNP.pval", 
                            "genome", "cohort", "outcome", "disease",
                            "no.typed.snps"))

write.table(v2.typed.snps, file="cg_v2table_adj.txt",
            col.names=T, row.names=F, 
            quote=F, sep="\t")
```

 \normalsize



\scriptsize



 \normalsize



## NOD2 analysis
NOD2 is a gene that has been frequently been studied in candidate gene studies looking at genetic variants from patients who have been treated with blood and marrow transplants (BMT). Oftentimes when NOD2 is studied, three SNPs frequently appear which have been labeled as SNP8 (rs2066844), SNP12 (rs2066845), and SNP13 (rs2066847). SNP8 and SNP12 were genotyped in in the DISCOVeRY-BMT GWAS and SNP13 was not genotyped in DISCOVeRY-BMT. In order to still consider this SNP, we searched for a SNP that may be in LD with rs2066847, and the best LD pair found in DISCOVeRY-BMT was rs146528649 ($r^2$ = 0.7422). All of the groups that studied NOD2, also chose the same design in the way they analyzed this SNP, such that only the patient-donor pairs that were homozygous wildtype for all 3 SNPs were deemed “wild type”, and on the contrary, if any one variant was present in any of the 3 SNPs in either of the donor-recipient pair, they were deemed to have a variant. Survival analysis was conducted using this aforementioned classification method.

 Here we will discuss how we pre-processed the DISCOVeRY-BMT data to grab these SNPs of interest and how we implemented a similar composite scoring design and subsequent survival analysis. In order to replicate these studies we will subset rs2066844, rs2066845 and rs146528649 from the imputed data. Each of these SNPs lie on  chr16 between the ranges of 50000000-55000000. The imputed data was found on the UB supercomputer at the following location:  


\scriptsize


```r
## File location:
"/projects/rpci/lsuchest/lsuchest/Rserve/
ImputeData/var/db/gwas/imputed_data/
BMT093013_forImpute/BMT093013_forImpute.chr16-50000000-55000000.impute2"
```

 \normalsize


 The imputed data was extracted and loaded into R so that the SNPs of interest and their corresponding genotype probabilities for each sample can be pulled. The impute2 data was coverted into a VCF file. Also we need a file listing the SNPs in an unlabeled column vector (`nod2_snps.txt`). We will use QCTOOL to convert `.impute2` to `.vcf`.


\scriptsize


```bash
module load qctool
qctool \
    -g BMT093013_forImpute.chr16-50000000-55000000.impute2 \
    -og /projects/rpci/lsuchest/abbasriz/\
    candidate_gene/nod2_rep/nod2_rep_dosages.vcf \
    -incl-rsids /projects/rpci/lsuchest/abbasriz/\
    candidate_gene/nod2_rep/nod2_snps.txt
```

 \normalsize


 Now that we have the vcf file, we can use `VariantAnnotation` library to easily pull this data into R.


\scriptsize


```r
library(VariantAnnotation)
## read in vcf file
vcf <- readVcf("/projects/rpci/lsuchest/abbasriz/
               candidate_gene/cg_haplotypes/nod2_rep/
               nod2_rep_dosages.vcf",
               genome="hg19")
```

 \normalsize


 Now we will grab the donor and recipients FID and IID. These correspond to the indices for the samples that we will subset.


\scriptsize


```r
# grab unique identifiers of patients in different cohorts
library(data.table)
patients <- c("/projects/rpci/lsuchest/lsuchest/Rserve/BMT/
              genetic_data/R_c1_EA_FID_IID.txt",
              "/projects/rpci/lsuchest/lsuchest/Rserve/BMT/
              genetic_data/R_c2_EA_FID_IID.txt",
              "/projects/rpci/lsuchest/lsuchest/Rserve/BMT/
              genetic_data/D_c1_EA_FID_IID.txt",
              "/projects/rpci/lsuchest/lsuchest/Rserve/BMT/
              genetic_data/D_c2_EA_FID_IID.txt")
cohorts <- lapply(patients, fread)
files <- c()
for(i in 1:length(patients)) files[i] <- strsplit(patients, "/")[[i]][9]      
names(cohorts) <- files
head(cohorts)
rec.ids <- list(cohorts[[1]], cohorts[[2]])
names(rec.ids) <- files[1:2]
donor.ids <- list(cohorts[[3]], cohorts[[4]])
names(donor.ids) <- files[3:4]
## read in phenotype files files
r.pheno <- fread("/projects/rpci/lsuchest/lsuchest/Rserve/
                 TheData/Plink_Recipient.pheno")
d.pheno <- fread("/projects/rpci/lsuchest/lsuchest/Rserve/
                 TheData/Plink_Donor.pheno")
r.cov <- fread("/projects/rpci/lsuchest/lsuchest/Rserve/
               TheData/Plink_Recipient.cov")
## SUBSET EUROPEAN AMERICAN PATIENTS BASED OFF OF IID
id.subset <- function(x, pheno){pheno[pheno$IID %in% x$IID]}
# recipients
rec.pheno.ea <- lapply(rec.ids, id.subset, r.pheno)
# donor
don.pheno.ea <- lapply(donor.ids, id.subset, d.pheno)
# bring down to just FID, IID, pair_id
don.pheno.ea <- lapply(don.pheno.ea, function(x) x[,c("FID",
                                                      "IID", 
                                                      "pair_id"),
                                                   with=F])
don.pheno.ea <- lapply(don.pheno.ea, setnames, c("D_FID", 
                                                 "D_IID", "pair_id"))
# merge to recipient file
merge.pair_id <- function(recipient, donor){
        setkey(recipient, pair_id)
        setkey(donor, pair_id)
        recipient[donor]
}
merged.ea <- mapply(merge.pair_id, rec.pheno.ea,
                    don.pheno.ea, SIMPLIFY=FALSE)
merged.ea <- data.table(do.call(rbind, merged.ea))
setnames(merged.ea, c("R_FID", "R_IID",
                      colnames(merged.ea)[3:ncol(merged.ea)]))
# I think there are duplicated NAs because of the almost 
# 2806 match and now 2815
merged.ea <- merged.ea[!duplicated(merged.ea$R_IID)]
ids <- merged.ea[,c("pair_id", "R_IID", "D_IID"), with=F]
# impute sample names
info.samples <- fread("/projects/rpci/lsuchest/lsuchest/Rserve/ImputeData/var/db/
                      gwas/imputed_data/BMT093013_forImpute/
                      BMT093013_forImpute.chr16-50000000-55000000.impute2_samples")
info.samples <- info.samples[-1]
nrow(info.samples)
# 6805 samples
nrow(gt)
# same number of samples...6805
# relabel the IDs so they just correspond to the sample order
# because when we converted the impute2 file to vcf, 
# the order remains the same as seen info the impute2_samples file
# so we can just relabel the IDs as 1:6805 samples and go from there
info.samples$ID_1 <- seq(1, nrow(info.samples))
# okay now we will grab the indices from the impute_sample 
# file and append them as columns to the ids data.table
ids$r_sample_index <- info.samples$ID_1[match(ids$R_IID,
                                              info.samples$ID_2)]
ids$d_sample_index <- info.samples$ID_1[match(ids$D_IID,
                                              info.samples$ID_2)]
# okay so now 22 donors didnt map back to that file
missing_samples <- ids[is.na(ids$d_sample_index)]
missing_samples
ids <- na.omit(ids)
head(ids)
# now we have 2783 samples that we can go 
# and grab all the genotype info on
```

 \normalsize


 Now we we will grab the genotype threshold values


\scriptsize


```r
## threshold genotypes
# grab genotypes of 3 SNPs for all patients
gt <- geno(vcf)$GT
gt[1:nrow(gt),1:5]
```

 \normalsize


 We can see that we have both typed and imputed for rs2066844. We will keep only the imputed ones (which is the second rs2066844)


\scriptsize


```r
# only we have both typed and imputed of rs2066844
# so we remove the first rs2066844 because it 
# is typed and we are keeping just the imputed ones
gt <- gt[-2,]
gt[1:nrow(gt),1:5]
```

 \normalsize


 Now we are going to transpose this data frame so we can have the samples as the rows and the SNPs as the columns


\scriptsize


```r
# transpose the data.frame so we can patients as rows
gt <- data.frame(t(gt))
gt[1:5, 1:ncol(gt)]
```

 \normalsize


 Now we code the alleles in a dominant model. So if there is any 1 alleles, we will just relabel it as 1 and if there are homozygous 0, we will relabel as 0. We will also relabel the . as NA and then remove the NAs


\scriptsize


```r
## coding dominant model over df columns
gt <- apply(gt, 2, FUN=function(x) gsub("[.]", NA, x))
gt <- apply(gt, 2, FUN=function(x) gsub("0/0", 0, x))
gt <- apply(gt, 2, FUN=function(x) gsub("0/1", 1, x))
gt <- apply(gt, 2, FUN=function(x) gsub("1/1", 1, x))
## create a data.frame that changes the factors into numeric
# and keeps the row.names ...
gt <- data.frame(apply(gt, 2, function(x) as.numeric(as.character(x))),
                 check.names=F, row.names=row.names(gt))
head(gt)
```

 \normalsize


 Now we will subset based off of the genome and cohort we are interested and create variables to do the composite SNP testing. So if SNP 8 and SNP 12 had 0 for both SNPs, they will be recoded as 0, and if any had 1 in it, it will be recoded as a 1. This will also done for the 3 SNPs (SNP8/SNP12/SNP13).


\scriptsize


```r
## parse the unique IDs (genome/cohort) from the main 
# genotype df into 4 dfs specific to genome/cohort
recs <- gt[ids$r_sample_index,]
donors <- gt[ids$d_sample_index,]
## recode with composite coding -- 0 if snps are all 0, 1 
# if any of the snps are 1
# snp 8 =  rs2066844 x[,2]
# snp 12 = rs2066845 x[,3]
# snp 13 = rs146528649 x[,1]
recode <- function(x){
        x$snp8snp12 <- ifelse((x[,2] == 0 & x[,3] == 0), x$snp8snp12 <- 0, x$snp8snp12 <- 1)
        x$compSNPs <- ifelse((x[,1] == 0 & x[,2] == 0 & x[,3] == 0), x$compSNPs <- 0, x$compSNPs <- 1)
        data.table(x, keep.rownames = T)
}
recs <- recode(recs)
donors <- recode(donors)
# remove "sample" from sample columns so we just have indices 
# that we can map back to phenotype file
rm.sample <- function(x){ x$rn <- gsub("sample_", "", x$rn)
class(x$rn) <- "numeric"
x}
recs <- rm.sample(recs)
donors <- rm.sample(donors)
# match column names to ids column names
colnames(recs)[1] <- "r_sample_index"
colnames(donors)[1] <- "d_sample_index"
# grab pair ids so we can compare donor and recipient pairs
recs$pair_id <- ids$pair_id
donors$pair_id <- ids$pair_id
# function to recode 'score', again, if there is 
# 0 in both D-R pairs the score 0, otherwise its 1
comp_score <- function(x){
    x$NOD2_score <- ifelse((x[,"recipient"] == 0 & x[,"donor"]
                            == 0), x$NOD2_score <- 0, 
                           x$NOD2_score <- 1)
    x
}
# merge snp8/snp12 donor-recipient pair data.frame
snp8snp12.r <- recs[,c("snp8snp12", "pair_id")]
setkey(snp8snp12.r, pair_id)
snp8snp12.d <- donors[,c("snp8snp12", "pair_id")]
setkey(snp8snp12.d, pair_id)
snp8snp12 <- snp8snp12.d[snp8snp12.r]
colnames(snp8snp12) <- c("donor", "pair_id", "recipient")
snp8snp12 <- snp8snp12[,c("pair_id", "donor", "recipient")]
snp8snp12 <- comp_score(snp8snp12)
# merge snp8/snp12/snp13 donor-recipient pair data.frame
all3.r <- recs[,c("compSNPs", "pair_id")]
setkey(all3.r, pair_id)
all3.d <- donors[,c("compSNPs", "pair_id")]
setkey(all3.d, pair_id)
all3 <- all3.d[all3.r]
colnames(all3) <- c("donor", "pair_id", "recipient")
all3 <- all3[,c("pair_id", "donor", "recipient")]
all3 <- comp_score(all3)
colnames(snp8snp12)[ncol(snp8snp12)] <- "snp8snp12_score"
colnames(all3)[ncol(all3)] <- "all3_score"
snp8snp12 <- snp8snp12[,c("pair_id","snp8snp12_score")]
all3 <- all3[,c("pair_id", "all3_score")]
```

 \normalsize


 Now we have the recoded SNP8 and SNP12. We need to append these back to the phenotype files so we can run survival analyses on these. Here we grab the recipient and donor phenotype files and merge them based off of pair_id column.


\scriptsize


```r
# merge back to phenotype file .. which is called merged.ea
setkey(merged.ea, pair_id)
setkey(snp8snp12, pair_id)
setkey(all3, pair_id)
merged.ea <- merged.ea[snp8snp12]
merged.ea <- merged.ea[all3]
# lets quickly see how many people had NA genotype scores
table(is.na(merged.ea$all3_score))
# 486 deemed NA for all 3
table(is.na(merged.ea$snp8snp12_score))
# 208 deemed NA for snp8snp12
# write.table(merged.ea, file="pheno_merged_nod2haps.txt",
# sep="\t", row.names=F, col.names=T, quote=F)
```

 \normalsize


 Survival analysis


\scriptsize


```r
#######################
##### SURVIVAL ANALYSIS
######################w
## do survival
library(survival)
# grab covariate files
r.cov <- fread("/projects/rpci/lsuchest/lsuchest/Rserve/TheData/
               Plink_Recipient.cov")

r.cov <- r.cov[r.cov$pair_id %in% ids$pair_id]
setkey(r.cov, pair_id)
r.cov$cohort <- merged.ea$cohort
r.cov.c1 <- r.cov[cohort=="c1"]
r.cov.c2 <- r.cov[cohort=="c2"]
##################################
## DRM covariates are: age, distatD
#################################
## cohort 1
# snp 8/snp12
DRM.c1 <- Surv(time=merged.ea.c1$intxsurv_1Y, 
              event=merged.ea.c1$disease_death_1Y==1)
DRM.snp8snp12.c1 <- coxph(DRM.c1~merged.ea.c1$snp8snp12_score+
                             r.cov.c1$distatD)
DRM.snp8snp12.c1.coef <- summary(DRM.snp8snp12.c1)$coef[1,1]
DRM.snp8snp12.c1.se.coef <- summary(DRM.snp8snp12.c1)$coef[1,3]
DRM.snp8snp12.c1.hr <- summary(DRM.snp8snp12.c1)$coef[1,2]
DRM.snp8snp12.c1.pval <- summary(DRM.snp8snp12.c1)$coef[1,5] 
# snp8/snp12/snp13
DRM.all3.c1 <- coxph(DRM.c1~merged.ea.c1$all3_score+
                        r.cov.c1$distatD)
## cohort 2
# snp 8/snp12
DRM.c2 <- Surv(time=merged.ea.c2$intxsurv_1Y, 
              event=merged.ea.c2$disease_death_1Y==1)
DRM.snp8snp12.c2 <- coxph(DRM.c2~merged.ea.c2$snp8snp12_score+
                             r.cov.c2$distatD)
# snp8/snp12/snp13
DRM.all3.c2 <- coxph(DRM.c2~merged.ea.c2$all3_score+
                        r.cov.c2$distatD)
DRM.all3.c2 <- coxph(DRM.c2~merged.ea.c2$all3_score+
                        r.cov.c2$distatD)
###########################################
## OS covariates are: age, distatD, Pblood
###########################################
#### OS ####
# cohort 1
## OS covariates are: age, distatD, Pblood
OS.c1 <- Surv(time=merged.ea.c1$intxsurv_1Y,
              event=merged.ea.c1$dead_1Y==1)
# OS SNP8/SNP12
os.snp8snp12.c1 <- coxph(OS.c1~merged.ea.c1$snp8snp12_score+
                             r.cov.c1$age+
                             r.cov.c1$distatD+
                             r.cov.c1$PBlood)
summary(os.snp8snp12.c1)$coef[1,5] 
# OS SNP8/SNP12/SNP13
os.all3.c1 <- coxph(OS.c1~merged.ea.c1$all3_score+
                        r.cov.c1$age+
                        r.cov.c1$distatD+
                        r.cov.c1$PBlood)
summary(os.all3.c1)$coef[1,5] 
## cohort 2
OS.c2 <- Surv(time=merged.ea.c2$intxsurv_1Y, 
              event=merged.ea.c2$dead_1Y==1)
# OS SNP8/SNP12
os.snp8snp12.c2 <- coxph(OS.c2~merged.ea.c2$snp8snp12_score+
                             r.cov.c2$age+
                             r.cov.c2$distatD+
                             r.cov.c2$PBlood)
summary(os.snp8snp12.c2)$coef[1,5] 

# OS SNP8/SNP12/SNP13
os.all3.c2 <- coxph(OS.c2~merged.ea.c2$all3_score+
                        r.cov.c2$age+
                        r.cov.c2$distatD+
                        r.cov.c2$PBlood)
summary(os.all3.c2)$coef[1,5] 
####################################
## PFS covariates: age, distatD
####################################
# cohort 1
PFS.c1 <- Surv(time=merged.ea.c1$intxsurv_1Y, 
               event=merged.ea.c1$lfs_1Y==1)
# PFS SNP8 / SNP12
pfs.snp8snp12.c1 <- coxph(PFS.c1~merged.ea.c1$snp8snp12_score+
                              r.cov.c1$distatD+r.cov.c1$age)
summary(pfs.snp8snp12.c1)$coefficients[1,5]
# PFS SNP 8 / SNP 12 / SNP 13
pfs.all3.c1 <- coxph(PFS.c1~as.numeric(merged.ea.c1$all3_score)+
                         r.cov.c1$distatD+
                         as.numeric(r.cov.c1$age))
summary(pfs.all3.c1)$coefficients[1,5]
# cohort 2
PFS.c2 <- Surv(time=merged.ea.c2$intxsurv_1Y,
               event=merged.ea.c2$lfs_1Y==1)
# PFS SNP8 / SNP12
pfs.snp8snp12.c2 <- coxph(PFS.c2~merged.ea.c2$snp8snp12_score+
                              r.cov.c2$distatD+
                              r.cov.c2$age)
summary(pfs.snp8snp12.c2)$coefficients[1,5]
# PFS SNP 8 / SNP 12 / SNP 13
pfs.all3.c2 <- coxph(PFS.c2~as.numeric(merged.ea.c2$all3_score)+
                         r.cov.c2$distatD+
                         as.numeric(r.cov.c2$age))
summary(pfs.all3.c2)$coefficients[1,5]
###############################################
# TRM covariates: age, bmiOBS, bmiOVWT, PBlood
##############################################
# cohort 1
TRM.c1 <- Surv(time=merged.ea.c1$intxsurv_1Y, event=merged.ea.c1$TRM_1Y==1)
# TRM SNP 8 / SNP 12
trm.snp8snp12.c1 <- coxph(TRM.c1~as.numeric(merged.ea.c1$snp8snp12_score)+
                              as.numeric(r.cov.c1$age)+
                              r.cov.c1$PBlood+r.cov.c1$bmiOBS+
                              r.cov.c1$bmiOVWT)
summary(trm.snp8snp12.c1)$coefficients[1,5]
# TRM SNP 8 / SNP 12 / SNP 13
trm.all3.c1 <- coxph(TRM.c1~as.numeric(merged.ea.c1$all3_score)+
                         as.numeric(r.cov.c1$age)+
                         r.cov.c1$PBlood+
                         as.numeric(r.cov.c1$bmiOBS)+
                         as.numeric(r.cov.c1$bmiOVWT))
summary(trm.all3.c1)$coefficients[1,5]
# cohort 2
TRM.c2 <- Surv(time=merged.ea.c2$intxsurv_1Y, 
               event=merged.ea.c2$TRM_1Y==1)
# TRM SNP 8 / SNP 12
trm.snp8snp12.c2 <- coxph(TRM.c2~as.numeric(merged.ea.c2$snp8snp12_score)+
                              as.numeric(r.cov.c2$age)+
                              r.cov.c2$PBlood+
                              r.cov.c2$bmiOBS+
                              r.cov.c2$bmiOVWT)
summary(trm.snp8snp12.c2)$coefficients[1,5]
# TRM SNP 8 / SNP 12 / SNP 13
trm.all3.c2 <- coxph(TRM.c2~as.numeric(merged.ea.c2$all3_score)+
                         as.numeric(r.cov.c2$age)+
                         r.cov.c2$PBlood+
                         as.numeric(r.cov.c2$bmiOBS)+
                         as.numeric(r.cov.c2$bmiOVWT))
summary(trm.all3.c2)$coefficients[1,5]


## WRITE TO TABLE
## NEED TO SPLIT INTO TWO TABLES FOR META-ANALYSIS
# will calculate confidence interval of hazard ratios AFTER
# set up columns so they are the same as our candidate
# gene table
cols <- c("gene",
          "rsID",
          "chr",
          "BP",
          "N",
          "allele1",
          "allele2",
          "coef",
          "se.coef",
          "exp.coef",
          "95%-CI",
          "Pvalue",
          "impute",
          "genome",
          "cohort",
          "outcome",
          "disease")
## Cohort 1
make.tbl <- function(outcome.vector,
                     gene,
                     rsID,
                     chr,
                     BP,
                     N,
                     allele1,
                     allele2,
                     genome, 
                     cohort,
                     outcome,
                     disease){
        coef <- outcome.vector[1,1]
        se.coef <- outcome.vector[1,3]
        hr <- outcome.vector[1,2]
        pval <- outcome.vector[1,5]
        res <- c(gene, rsID, chr, BP, N, allele1, 
                 allele2, coef, se.coef, hr, NA, pval,
                 "imputed", genome, cohort, outcome, disease)
        res
}

# cohort 1
nod2.c1 <- data.table(rbind(
    make.tbl(summary(DRM.snp8snp12.c1)$coef, "NOD2", "2 SNPs* DRM", "chr16",
             "*", 2033, "*", "*", "S", "c1", "DRM", "mixed"),
    make.tbl(summary(os.snp8snp12.c1)$coef, "NOD2", "2 SNPs* OS", "chr16",
             "*", 2033, "*", "*", "S", "c1", "OS", "mixed"),
    make.tbl(summary(trm.snp8snp12.c1)$coef, "NOD2", "2 SNPs* TRM", "chr16",
             "*", 2033, "*", "*", "S", "c1", "TRM", "mixed"),
    make.tbl(summary(pfs.snp8snp12.c1)$coef, "NOD2", "2 SNPs* PFS", "chr16",
             "*", 2033, "*", "*", "S", "c1", "PFS", "mixed"),
    make.tbl(summary(DRM.all3.c1)$coef, "NOD2", "3 SNPs* DRM", "chr16", "*",
             2033, "*", "*", "S", "c1", "DRM", "mixed"),
    make.tbl(summary(os.all3.c1)$coef, "NOD2", "3 SNPs* OS", "chr16", "*", 
             2033, "*", "*", "S", "c1", "OS", "mixed"),
    make.tbl(summary(pfs.all3.c1)$coef, "NOD2", "3 SNPs* PFS", "chr16", "*",
             2033, "*", "*", "S", "c1", "PFS", "mixed"),
    make.tbl(summary(trm.all3.c1)$coef, "NOD2", "3 SNPs* TRM", "chr16", "*",
             2033, "*", "*", "S", "c1", "TRM", "mixed")))
colnames(nod2.c1) <- cols
write.table(nod2.c1, file="nod2.c1.txt", col.names=T, row.names=F, quote=F, sep='\t')
## Cohort 2
nod2.c2 <- data.table(rbind(
    make.tbl(summary(DRM.snp8snp12.c2)$coef, "NOD2", "2 SNPs* DRM", "chr16", "*",
             757, "*", "*", "S", "c2", "DRM", "mixed"),
    make.tbl(summary(os.snp8snp12.c2)$coef, "NOD2", "2 SNPs* OS", "chr16", "*",
             757, "*", "*", "S", "c2", "OS", "mixed"),
    make.tbl(summary(pfs.snp8snp12.c2)$coef, "NOD2", "2 SNPs* PFS", "chr16", "*", 
             757, "*", "*", "S", "c2", "PFS", "mixed"),
    make.tbl(summary(trm.snp8snp12.c2)$coef, "NOD2", "2 SNPs* TRM", "chr16", "*",
             757, "*", "*", "S", "c2", "TRM", "mixed"),
    make.tbl(summary(DRM.all3.c2)$coef, "NOD2", "3 SNPs* DRM", "chr16", "*",
             757, "*", "*", "S", "c2", "DRM", "mixed"),
    make.tbl(summary(os.all3.c2)$coef, "NOD2", "3 SNPs* OS", "chr16", "*", 
             757, "*", "*", "S", "c2", "OS", "mixed"),
    make.tbl(summary(pfs.all3.c2)$coef, "NOD2", "3 SNPs* PFS", "chr16", "*",
             757, "*", "*", "S", "c2", "PFS", "mixed"),
    make.tbl(summary(trm.all3.c2)$coef, "NOD2", "3 SNPs* TRM", "chr16", "*",
             757, "*", "*", "S", "c2", "TRM", "mixed")))
colnames(nod2.c2) <- cols
write.table(nod2.c2, file="nod2.c2.txt", col.names=T, 
            row.names=F, quote=F, sep='\t')
```

 \normalsize


 Meta-analysis was performed using METAL software. The hazard ratios and 95\% confidence intervals were computed using the standard errors from the METAL output. 


\scriptsize


```r
# RUN METAL 
## append meta results
nod2.meta <- fread("metal_S_M_mixed_nod21.tbl")
# remove direction column
nod2.meta <- nod2.meta[,-7]
# build up columns to match how our final table is set up
nod2.meta$gene <- "NOD2"
colnames(nod2.meta)[1] <- "rsID"
nod2.meta$chr <- "chr16"
nod2.meta$exp.coef <- exp(nod2.meta$Effect)
nod2.meta$genome <- "S"
nod2.meta$cohort <- "M"
nod2.meta$disease <- "mixed"
nod2.meta$N <- 2790
nod2.meta$BP <- "*"
nod2.meta$impute <- "imputed"
nod2.meta$`95%-CI` <- NA
nod2.meta$outcome <- NA
colnames(nod2.meta) <- c("rsID", "allele1", "allele2", "coef", "se.coef",
                         "Pvalue","gene", "chr", "exp.coef", "genome", 
                         "cohort",
                         "disease", "N", "BP", "impute", "95%-CI", "outcome")

nod2.meta$outcome[c(1,8)] <- "OS"
nod2.meta$outcome[2:3] <- "DRM"
nod2.meta$outcome[c(4,6)] <- "PFS"
nod2.meta$outcome[c(5,7)] <- "TRM"
setcolorder(nod2.meta, cols)
# join all of the nod2 results (c1, c2, meta analyiss)
nod2.final <- do.call(rbind, list(nod2.c1, nod2.c2, nod2.meta))
# remove extra names in rsID column
nod2.final$rsID <- substr(nod2.final$rsID, 1, 7)
## now calculate 95% confidence interval for hazard ratio
hr.ci <- function(coef.est, se){
        lb <- round(exp(coef.est-1.96*se), 5)
        ub <- round(exp(coef.est+1.96*se), 5)
        paste0("[", lb,", " , ub ,"]")
}
for(i in 1:nrow(nod2.final)){
        nod2.final[i,"95%-CI"] <-
            hr.ci(as.numeric(nod2.final$coef)[i], 
                  as.numeric(nod2.final$se.coef)[i])
}
write.table(nod2.final, file="nod2_survival_pvals.txt",
            col.names=T, row.names=F, sep="\t", quote=F)
```

 \normalsize


### CCR5 analysis



Creating a file with the CCR SNPs in it:  

\scriptsize


```bash
printf "rs1799987\nrs333\nrs1800023\nrs1800024\nrs113341849" > ccr5_snps_tt.txt
```

 \normalsize

 
\scriptsize


```r
library(VariantAnnotation)	
library(data.table)
library(survival)
library(dplyr)
## read in vcf file
vcf <- readVcf("/projects/rpci/lsuchest/abbasriz/
           candidate_gene/cg_haplotypes/ccr5_rep/
           ccr5_rep_dosages.vcf",
           genome="hg19")
## threshold genotypes
# grab genotypes of 3 SNPs for all patients
gt <- geno(vcf)$GT
gt[1:nrow(gt),1:5]
## no rs333 in our data
## proxy which rs113341849
# transpose the data.frame so we can patients as rows
gt <- data.frame(t(gt))
gt[1:5, 1:ncol(gt)]
## CCR5 rs113341849
table(gt$rs113341849)
h1h1 <- gt
# recode to homozygous major allele
h1h1 <- apply(h1h1, 2, FUN=function(x) gsub("[.]", NA, x))
h1h1 <- apply(h1h1, 2, FUN=function(x) gsub("0/0", 1, x))
h1h1 <- apply(h1h1, 2, FUN=function(x) gsub("0/1", 0, x))
h1h1 <- apply(h1h1, 2, FUN=function(x) gsub("1/1", 0, x))

dosages <- fread("ccr5_dosages.impute2")
dosages <- data.table(t(dosages))
colnames(dosages) <- as.character(dosages[3,])
dosages <- dosages[-c(1:6),]
dosages <- data.frame(apply(dosages, 2, as.numeric))
dosages$ccr5 <- rowSums(dosages)
dosages <- dosages %>%
	mutate(h1h1=ifelse(ccr5<0.5,1,0)) %>%
	data.table(keep.rownames=T)
###############################################
###### LOAD UP PATIENT ID AND COVARIATE FILES
###############################################
## grab unique identifiers of patients in different cohorts
library(data.table)
patients <-
c(
"/projects/rpci/lsuchest/lsuchest/
Rserve/BMT/genetic_data/R_c1_EA_FID_IID.txt",
"/projects/rpci/lsuchest/lsuchest
/Rserve/BMT/genetic_data/R_c2_EA_FID_IID.txt",
"/projects/rpci/lsuchest/lsuchest/
Rserve/BMT/genetic_data/D_c1_EA_FID_IID.txt",
"/projects/rpci/lsuchest/lsuchest/
Rserve/BMT/genetic_data/D_c2_EA_FID_IID.txt"
)
cohorts <- lapply(patients, fread)
files <- c()
for(i in 1:length(patients)) files[i] <- 
    strsplit(patients, "/")[[i]][9]      
names(cohorts) <- files
head(cohorts)
rec.ids <- list(cohorts[[1]], cohorts[[2]])
names(rec.ids) <- files[1:2]
rec.ids[[1]]$cohort <- "c1"
rec.ids[[2]]$cohort <- "c2"
donor.ids <- list(cohorts[[3]], cohorts[[4]])
names(donor.ids) <- files[3:4]
donor.ids[[1]]$cohort <- "c1"
donor.ids[[2]]$cohort <- "c2"
## read in phenotype files files
r.pheno <- fread("/projects/rpci/lsuchest/lsuchest/Rserve/
                 TheData/Plink_Recipient.pheno")
d.pheno <- fread("/projects/rpci/lsuchest/lsuchest/Rserve/
                 TheData/Plink_Donor.pheno")
r.cov <- fread("/projects/rpci/lsuchest/lsuchest/Rserve/
               TheData/Plink_Recipient.cov")
# parse covariant/pheno files to just EA
## SUBSET EUROPEAN AMERICAN PATIENTS BASED OFF OF IID
id.subset <- function(x, pheno){pheno[pheno$IID %in% x$IID]}
# recipients
rec.pheno.ea <- lapply(rec.ids, id.subset, r.pheno)
# donors
don.pheno.ea <- lapply(donor.ids, id.subset, d.pheno)
# bring down to just FID, IID, pair_id
don.pheno.ea <- 
    lapply(don.pheno.ea,
           function(x) x[,c("FID", 
                            "IID", 
                            "pair_id", 
                            "population"),with=F])
don.pheno.ea <- 
    lapply(don.pheno.ea,
           setnames, c("D_FID",
                       "D_IID", "pair_id", "population"))
# merge to recipient file
merge.pair_id <- function(recipient, donor){
        setkey(recipient, pair_id)
        setkey(donor, pair_id)
        recipient[donor]
}
merged.ea <- mapply(merge.pair_id,
                    rec.pheno.ea,
                    don.pheno.ea, 
                    SIMPLIFY=FALSE)
merged.ea <- data.table(do.call(rbind, merged.ea))
setnames(merged.ea,
         c("R_FID", "R_IID", 
           colnames(merged.ea)[3:ncol(merged.ea)]))
# I think there are duplicated NAs because
# of the almost 2806 match and now 2815
merged.ea <- merged.ea[!duplicated(merged.ea$R_IID)]
ids <- merged.ea[,c("pair_id", "R_IID", "D_IID"), with=F]
# impute sample names
info.samples <- fread("/projects/rpci/
                      lsuchest/lsuchest/Rserve/ImputeData/
                      var/db/gwas/imputed_data/
                      BMT093013_forImpute/
                      BMT093013_forImpute.chr16-
                      50000000-55000000.impute2_samples")
# top row is useless
info.samples <- info.samples[-1]
nrow(info.samples)
# 6805 samples
nrow(gt)
# nice! same number of samples...6805
# relabel the IDs so 
#they just correspond to the sample order
# because when we converted the impute2 file to vcf, 
# the order remains the same as seen info the 
# impute2_samples file
# so we can just relabel the IDs as 1:6805 samples and go 
# from there
info.samples$ID_1 <- seq(1, nrow(info.samples))
# okay now we will grab the indices from the impute_sample 
# file and append them as columns to the ids data.table
ids$r_sample_index <-
    info.samples$ID_1[match(ids$R_IID, info.samples$ID_2)]
ids$d_sample_index <- 
    info.samples$ID_1[match(ids$D_IID, info.samples$ID_2)]
# okay so now 22 donors didnt map back to that file
missing_samples <- ids[is.na(ids$d_sample_index)]
missing_samples
ids <- na.omit(ids)
# now we have 2783 samples that we can go and grab 
# all the genotype info on
##################################################
#### NOW PARSE GENOTYPE FILE BY SAMPLE ID INDICES
###################################################
## parse the unique IDs (genome/cohort) from the main 
# genotype df into 4 dfs specific to genome/cohort
recs <- data.table(dosages[ids$r_sample_index,],
                   keep.rownames=T)
donors <- data.table(dosages[ids$d_sample_index,],
                     keep.rownames=T)
# remove "sample" from sample columns so we just have 
# indices that we can map back to phenotype file
rm.sample <- function(x){ x$rn <- gsub("sample_", "", x$rn)
	class(x$rn) <- "numeric"
	x}
recs <- rm.sample(recs)
donors <- rm.sample(donors)
# match column names to ids column names
colnames(recs)[1] <- "r_sample_index"
colnames(donors)[1] <- "d_sample_index"
# grab pair ids so we can compare donor and recipient pairs
recs$pair_id <- ids$pair_id
donors$pair_id <- ids$pair_id
recs <- recs %>% na.omit()
donors <- donors %>% na.omit()
recs <- recs %>% rename(Rh1h1=h1h1)
donors <- donors %>% rename(Dh1h1=h1h1)
donors.h1h1 <- donors %>%
	select(pair_id, Dh1h1)
recs.h1h1 <- recs %>%
	select(pair_id, Rh1h1)
dr.h1h1 <- donors.h1h1 %>%
	right_join(recs.h1h1, "pair_id") %>%
	data.table()
dr.h1h1$grp2vsgrp1 <- with(dr.h1h1,
	ifelse(Dh1h1==1 & Rh1h1==0, 1,
	       ifelse(Dh1h1==0 & Rh1h1==0, 0, NA)))
dr.h1h1$grp3vsgrp1 <- with(dr.h1h1,
	ifelse(Dh1h1==0 & Rh1h1==1, 1,
	       ifelse(Dh1h1==0 & Rh1h1==0, 0, NA)))
dr.h1h1$grp3vsgrp2 <- with(dr.h1h1,
	ifelse(Dh1h1==0 & Rh1h1==1, 1, 
	       ifelse(Dh1h1==1 & Rh1h1==0, 0, NA)))
r.cov <- r.cov %>%
    right_join(dr.h1h1, "pair_id") %>%
    data.table()
## split into cohorts
r.cov.c1 <- r.cov[cohort1==1]
r.cov.c2 <- r.cov[cohort1==0]
merged.ea <- r.cov %>%
dplyr::select(
    pair_id,
    age,
    distatD,
    PBlood,
    bmiOBS,
    bmiOVWT,
    MDSdummy,
    AMLdummy,
    ALLdummy,
    Dh1h1,
    Rh1h1,
    grp2vsgrp1,
    grp3vsgrp1,
    grp3vsgrp2
) %>%
inner_join(merged.ea, by = "pair_id") %>%
data.table()
merged.ea$population <- gsub("EA.", "",
                             merged.ea$population)
OScov1Y <- c("intxsurv_1Y", 
             "dead_1Y", "age", "distatD", "PBlood")
PFScov1Y <- c("intxrel_1Y", "lfs_1Y", "age", "distatD")
OScov.full <- c("intxsurv_1Y",
                "dead_1Y", "age", "distatD", "PBlood")
PFScov.full <- c("intxrel_1Y", "lfs_1Y", "age", "distatD")
###############################################
######  SURVIVAL ANALYSIS
###############################################
########################
######## EVENTS ########
########################
## DRM - disease_death_1Y
## OS  - dead_1Y
## PFS - lfs_1Y
## TRM - TRM_1Y
########################
####################################################
################## COVARIATES ######################
####################################################
## DRM covariates are: age, distatD
## OS covariates are: age, distatD, Pblood
## PFS covariates: age, distatD
## TRM covariates: age, bmiOBS, bmiOVWT, PBlood
###################################################
DRMcov <- c("intxsurv_1Y", 
           "disease_death_1Y", "age", "distatD")
OScov <- c("intxsurv_1Y", 
           "dead_1Y", "age", "distatD", "PBlood")
PFScov <- c("intxrel_1Y",
            "lfs_1Y", "age", "distatD")
TRMcov <-c("intxsurv_1Y",
           "TRM_1Y", "age", "bmiOBS", "bmiOVWT", "PBlood")
# adjusts for 2 covariates + genotype of interest
surv_fit_cov2 <- function(geno, cov, covFile) {
outcome <- Surv(time = covFile[, cov[1]], 
                event = covFile[, cov[2]])
res <-
coxph(outcome ~ covFile[, geno] + 
          covFile[, cov[3]] +
          covFile[, cov[4]])
summary(res)$coefficients[1, c(1, 3, 2, 5)]
}
# adjusts for 3 covariates + genotype of interest
surv_fit_cov3 <- function(geno, cov, covFile) {
outcome <- Surv(time = covFile[, cov[1]],
                event = covFile[, cov[2]])
res <-
coxph(outcome ~ covFile[, geno] + 
          covFile[, cov[3]] + 
          covFile[, cov[4]] +
covFile[, cov[5]])
summary(res)$coefficients[1, c(1, 3, 2, 5)]
}
# adjusts for 4 covariates + genotype of interest
surv_fit_cov4 <- function(geno, cov, covFile) {
    outcome <- Surv(time = covFile[, cov[1]], 
                    event = covFile[, cov[2]])
    res <-
    coxph(outcome ~ covFile[, geno] + 
              covFile[, cov[3]] +
              covFile[, cov[4]] +
    covFile[, cov[5]] + 
        covFile[, cov[6]])
    summary(res)$coefficients[1, c(1, 3, 2, 5)]
}
## WRITE TO TABLE
## NEED TO SPLIT INTO TWO TABLES FOR META-ANALYSIS
# will calculate confidence interval of HRs AFTER
# set up columns so they are the same as our CG table
cols <- c("gene",
          "rsID",
          "chr",
          "BP",
          "N",
          "allele1",
          "allele2",
          "coef",
          "se.coef",
          "exp.coef",
          "95%-CI",
          "Pvalue",
          "impute",
          "genome",
          "cohort",
          "outcome",
          "disease")

make.tbl <-
    function(outcome.vector,
    gene,
    rsID,
    chr,
    BP,
    N,
    allele1,
    allele2,
    imputed,
    genome,
    cohort,
    outcome,
    disease) {
    coef <- outcome.vector[1]
    se.coef <- outcome.vector[2]
    hr <- outcome.vector[3]
    pval <- outcome.vector[4]
    res <-
    c(
    gene,
    rsID,
    chr,
    BP,
    N,
    allele1,
    allele2,
    coef,
    se.coef,
    hr,
    NA,
    pval,
    imputed,
    genome,
    cohort,
    outcome,
    disease
    )
    res
    }
###############################
#### MASTER FILE CREATION #####
###############################
compsnps <- c("R_h1h1", 
              "grp2vsgrp1",
              "grp3vsgrp1", 
              "grp3vsgrp2")
cov.list <- c("DRMcov", 
              "PFScov", 
              "OScov",
              "TRMcov")
outcomes <- c("DRM",
              "PFS",
              "OS",
              "TRM")
cohorts <- c("c1",
             "c2")
diseases <- c("mixed",
              "AMLonly",
              "ALLonly",
              "noMDS",
              "noALL")
survival.functions <-
    c("surv_fit_cov2",
    "surv_fit_cov2",
    "surv_fit_cov3",
    "surv_fit_cov4")
create.master <- function(genotype,
                          genome, 
                          outcomes, 
                          cohorts, 
                          diseases, 
                          survival.functions){
	master <- data.frame(matrix(nrow=3, ncol=2))
	#master <- list()
	for(i in 1:length(diseases)){
		master[i,] <- c(genotype, diseases[i])

	}
	colnames(master) <- c("genotype", "disease")
	master.list <- list()
	for(i in 1:length(outcomes)){
		master$outcome <- outcomes[i]
		master$covList <- cov.list[i]
		master$survivalFunc <- survival.functions[i]
		master$genome <- genome
		master.list[[i]] <- master
	}
	master.merge <- data.table(do.call(rbind, master.list))
	master.list <- list()
	for(i in 1:length(cohorts)){
		master.merge$cohort <- cohorts[i]
		master.list[[i]] <- master.merge
	}
	do.call(rbind, master.list)
}
# donors
master.dh1h1 <- create.master("Dh1h1",
                              "D",
                              outcomes,
                              cohorts, 
                              diseases, 
                              survival.functions)
# recipients h1h1
master.rh1h1 <- create.master("Rh1h1", 
                              "R", 
                              outcomes,
                              cohorts,
                              diseases,
                              survival.functions)
# shared
master.grp2vsgrp1 <-
    create.master("grp2vsgrp1",
    "S",
    outcomes,
    cohorts,
    diseases,
    survival.functions)
master.grp3vsgrp1 <-
    create.master("grp3vsgrp1",
    "S",
    outcomes,
    cohorts,
    diseases,
    survival.functions)
master.grp3vsgrp2 <-
    create.master("grp3vsgrp2",
    "S",
    outcomes,
    cohorts,
    diseases,
    survival.functions)
master.list <-
    data.table(do.call(
    rbind,
    list(
    master.dh1h1,
    master.rh1h1,
    master.grp2vsgrp1,
    master.grp3vsgrp1,
    master.grp3vsgrp2
    )
    ))
outcome.order <- c("DRM", "PFS", "OS", "TRM")
# order by outcome so we can have DRM and PFS as \
# top two outcomes b/c they both have 2 covariates
master.list.2cov <- master.list[order(
    match(master.list$outcome, outcome.order))][1:100]
# OS has 3 covariates so we will grab those
master.list.3cov <- master.list[order(
    match(master.list$outcome, 
          outcome.order))][101:150]
# trm has 4 covs
master.list.4cov <- master.list[order(
    match(master.list$outcome,
          outcome.order))][151:200]
## subset by cohort and disease
cohort.disease <- function(pheno, master.list){
	if(master.list$cohort=="c1"){
		pheno <- pheno[cohort1==1]
	} else {
		pheno <- pheno[cohort1==0]
	}
	if(master.list$disease=="mixed"){
		data.frame(pheno)
	} else if (master.list$disease=="AMLonly"){
		data.frame(pheno[AMLdummy==1])
	} else if(master.list$disease=="ALLonly"){
		data.frame(pheno[ALLdummy==1])
	} else if (master.list$disease == "noALL"){
		data.frame(pheno[ALLdummy==0])
	} else if (master.list$disease == "noMDS"){
		data.frame(pheno[MDSdummy==0])
	}

}
# PFS 
surv.res.DRMpfs <- data.frame(matrix(nrow=100, ncol=17))
colnames(surv.res.DRMpfs) <- cols
for(i in 1:nrow(surv.res.DRMpfs)){
	surv.res.DRMpfs[i,] <- make.tbl(
		surv_fit_cov2(master.list.2cov$genotype[i],
			eval(as.name(
			    paste(master.list.2cov$covList[i])
			    )),
			cohort.disease(merged.ea, master.list.2cov[i,])),
		"CCR5",
		paste(master.list.2cov$genotype[i],
		      master.list.2cov$outcome[i], 
		      master.list.2cov$disease[i],
		      master.list.2cov$genome[i], sep="_"),
		"chr3",
		"*",
		cohort.disease(merged.ea,
		               master.list.2cov[i,]) %>% 
		    select(cohort, eval(as.name(
		        paste(master.list.2cov$genotype[i])))) %>%
		    filter(cohort==master.list.2cov$cohort[i]) %>%
		    na.omit %>%
		    nrow,
		"*",
		"*",
		"imputed",
		master.list.2cov$genome[i],
		master.list.2cov$cohort[i],
		master.list.2cov$outcome[i],
		master.list.2cov$disease[i])
	surv.res.DRMpfs
}
# OVERALL SURVIVAL 
surv.res.os <- data.frame(matrix(nrow=50, ncol=17))
colnames(surv.res.os) <- cols
for(i in 1:nrow(surv.res.os)){
	surv.res.os[i,] <- make.tbl(
		surv_fit_cov2(master.list.3cov$genotype[i],
			eval(as.name(
			    paste(master.list.3cov$covList[i]))),
			cohort.disease(merged.ea,
			               master.list.3cov[i,])),
		"CCR5",
		paste(master.list.3cov$genotype[i],
		      master.list.3cov$outcome[i], 
		      master.list.3cov$disease[i], 
		      master.list.3cov$genome[i],
		      sep="_"),
		"chr3",
		"*",
		cohort.disease(merged.ea, 
		               master.list.3cov[i,]) %>%
		    select(cohort, eval(
		        as.name(paste(master.list.3cov$genotype[i]))
		        )) %>%
		    filter(cohort==master.list.3cov$cohort[i]) %>%
		    na.omit %>%
		    nrow,
		"*",
		"*",
		"imputed",
		master.list.3cov$genome[i],
		master.list.3cov$cohort[i],
		master.list.3cov$outcome[i],
		master.list.3cov$disease[i])
	surv.res.os
}
# TRM 
surv.res.trm <- data.frame(matrix(nrow=50, ncol=17))
colnames(surv.res.trm) <- cols
for(i in 1:nrow(surv.res.trm)){
	surv.res.trm[i,] <- make.tbl(
		surv_fit_cov2(master.list.4cov$genotype[i],
			eval(as.name(
			    paste(master.list.4cov$covList[i]))),
			cohort.disease(merged.ea,
			               master.list.4cov[i,])),
		"CCR5",
		paste(master.list.4cov$genotype[i],
		      master.list.4cov$outcome[i], 
		      master.list.4cov$disease[i], 
		      master.list.4cov$genome[i], sep="_"),
		"chr3",
		"*",
		cohort.disease(merged.ea, 
		               master.list.4cov[i,]) %>%
		    select(cohort, eval(as.name(paste(
		        master.list.4cov$genotype[i])))) %>%
		    filter(cohort==master.list.4cov$cohort[i]) %>%
		    na.omit %>% 
		    nrow,
		"*",
		"*",
		"imputed",
		master.list.4cov$genome[i],
		master.list.4cov$cohort[i],
		master.list.4cov$outcome[i],
		master.list.4cov$disease[i])
	surv.res.trm
}
## save to file
surv.res <- data.table(do.call(
    rbind, list(surv.res.DRMpfs,
                surv.res.os,
                surv.res.trm)
    ))
# split into cohorts
surv.res.c1 <- surv.res[cohort=="c1"]
surv.res.c2 <- surv.res[cohort=="c2"]
write.table(
    surv.res.c1,
    file = "h1h1_nometa_c1.txt",
    sep = "\t",
    quote = F,
    row.names = F,
    col.names = T
    )
write.table(
    surv.res.c2,
    file = "h1h1_nometa_c2.txt",
    sep = "\t",
    quote = F,
    row.names = F,
    col.names = T
    )
## RAN META ANALYSIS ###
```

 \normalsize


 CCR5 meta analysis was run and the results are loaded in below:   


\scriptsize


```r
## LOAD META RESULTS
surv.res.meta <- fread("metal_ccr5_full1.tbl")
surv.res.meta <- surv.res.meta[,-7]
meta.res <- function(ccr5.meta, col.order){
# remove direction column
# build up columns to match how our final table is set up
	ccr5.meta$gene <- "CCR5"
	colnames(ccr5.meta)[1] <- "rsID"
	ccr5.meta$chr <- "chr3"
	ccr5.meta$exp.coef <- exp(ccr5.meta$Effect)
	ccr5.meta$genome <-
	    sapply(strsplit(surv.res.meta$MarkerName, '_',
	                    fixed = T), "[", 4)
	ccr5.meta$cohort <- "M"
	ccr5.meta$disease <-
	    sapply(strsplit(surv.res.meta$MarkerName, '_',
	                    fixed = T), "[", 3)
	ccr5.meta$N <- NA
	ccr5.meta$BP <- "*"
	ccr5.meta$impute <- "imputed"
	ccr5.meta$`95%-CI` <- NA
	ccr5.meta$outcome <-
	    sapply(strsplit(surv.res.meta$MarkerName, '_', 
	                    fixed = T), "[", 2)
	colnames(ccr5.meta) <-
	    c(
	    "rsID",
	    "allele1",
	    "allele2",
	    "coef",
	    "se.coef",
	    "Pvalue",
	    "gene",
	    "chr",
	    "exp.coef",
	    "genome",
	    "cohort",
	    "disease",
	    "N",
	    "BP",
	    "impute",
	    "95%-CI",
	    "outcome"
	    )
	setcolorder(ccr5.meta, cols)
	data.table(ccr5.meta)
}

surv.res.meta <- meta.res(surv.res.meta,cols)
# need to make sure in right order for N
setkey(surv.res.c1, rsID)
setkey(surv.res.c2, rsID)
setkey(surv.res.meta, rsID)
surv.res.meta$N <- as.numeric(surv.res.c1$N) +
    as.numeric(surv.res.c2$N)
surv.res.full <- do.call(rbind,
                         list(surv.res.c1,
                              surv.res.c2,
                              surv.res.meta)
                         )
surv.res.full$rsID <-
    sapply(strsplit(surv.res.full$rsID,
                    '_', 
                    fixed = T), "[", 1)
ccr5.final <- surv.res.full
# calculate confidence interval
hr.ci <- function(coef.est, se){
        lb <- round(exp(coef.est-1.96*se), 5)
        ub <- round(exp(coef.est+1.96*se), 5)
        paste0("[", lb,", " , ub ,"]")
}
for(i in 1:nrow(ccr5.final)){
    ccr5.final[i, "95%-CI"] <-
        hr.ci(as.numeric(ccr5.final$coef)[i],
        as.numeric(ccr5.final$se.coef)[i])
}
write.table(
    ccr5.final,
    file = "CCR5_H1H1_FINAL_wMETA.txt",
    sep = "\t",
    quote = F,
    row.names = F,
    col.names = T
    )
```

 \normalsize

## Candidate Replication and Validation Plots

\scriptsize


```r
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
mac <- read.table("~/Google Drive/OSU_PHD/dissertation/code/chapter2/MacMillan2003.txt",
                  header=T, 
                  sep="\t", 
                  stringsAsFactors = F)
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
    if(repval=="Replication"){
        group.facet.labs <- c(`Replication in DISCOVeRY-BMT`="Replication in \n DISCOVeRY-BMT",
                              Literature="Literature")
    }else if(repval=="Validation"){
        group.facet.labs <- c(Literature="Literature",
                              `Validation in DISCOVeRY-BMT`="Validation in \n DISCOVeRY-BMT")
    }
    pR <- fig.tbl %>%
        mutate(SNP=str_replace(SNP, "D-R group 3 vs D-R group 1", 
                               "D-R group 3 vs\n D-R group 1"),
               SNP=str_replace(SNP, "D-R group 3 vs D-R group 2",
                               "D-R group 3 vs\n D-R group 2"),
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
              axis.text.x = element_text(hjust = 1, size=18, angle = -90),
              axis.text.y = element_text(hjust = 1, size=18, angle = -90),
              strip.text.x = element_text(size = 18, face="italic", angle=-90),
              strip.text.y = element_text(size = 18, angle=-90),
              axis.title.y = element_text(size=18, angle=-90),
              axis.title.x = element_text(size=18, angle=-90),
              panel..x=unit(0.1, "lines"),
              panel..y=unit(0.1, "lines"),
              legend.position=c(0.05, 0.565)) + 
        scale_y_reverse() +
        ylab("-log10(Pvalue)") +
        scale_shape_discrete(solid=FALSE, 
                             guide = guide_legend(direction = "horizontal", 
                                                  title.position = "right",
                                                  title.theme=element_text(angle = -90,
                                                                           size=18),
                                                  label.position="bottom",
                                                  label.hjust = 0.5, 
                                                  label.vjust = 0.5,
                                                  label.theme = element_text(angle = -90,
                                                                             size=18)
                                                  )
                             ) +
        scale_color_discrete(guide = guide_legend(direction = "horizontal",
                                                  title.position = "right",
                                                  title.theme=element_text(angle=-90,
                                                                           size=18),
                                                  label.position="bottom", 
                                                  label.hjust = 0.5, 
                                                  label.vjust = 0.5,
                                                  label.theme = element_text(angle = -90,
                                                                             size=18)
                                                  )
                             ) 
    
    
    xr <- ggplotGrob(pR)
    #labelR <- "Pvalues"
    labelT <- "Genes"
    posT <- subset(xr$layout,
                   grepl("strip-t", name),
                   select=t:r)
    height <- xr$heights[min(posT$t)]
    xr <- gtable_add_rows(xr, height, min(posT$t)-1)
    # Construct the new strip grobs
    stripT <- gTree(name = "Strip_top",
                    children = gList(
                        rectGrob(gp = gpar(col = NA, 
                                           fill = "grey85")),
                        textGrob(labelT, rot = -90,
                                 gp = gpar(fontsize = 18,
                                           col = "grey10"))))
    # Position the grobs in the gtable
    xr <- gtable_add_grob(xr,
                          stripT,
                          t = min(posT$t),l = min(posT$l)
                          ,r = max(posT$r),
                          name = "strip-top")
    # Add small gaps between strips
    xr <- gtable_add_rows(xr,
                          unit(1/5, "line"),
                          min(posT$t))
    # Draw it
    grid.newpage()
    grid.draw(xr)
} 

plotFun("Replication")
plotFun("Validation")
```

 \normalsize

## RegulomeDB of Candidate SNPs
\scriptsize


```r
library(tidyverse)
x <- readxl::read_xlsx("~/Google Drive/OSU_PHD/dissertation/code/chapter2/repval_supps.xlsx",
                       sheet=1,
                       na = "NS") ## this is supp table 1 from CG paper
x <- x %>%
    mutate(P.Value = as.numeric(P.Value),
           pval_cat = case_when(
               is.na(P.Value) ~ "Not Significant",
               P.Value >= 0.05 ~ "Not Significant",
               P.Value < 0.05 ~ "Significant")) %>%
    select(SNP, pval_cat) %>%
    distinct()
regdb <- read_delim("~/Downloads/RegulomeDB.dbSNP141.lessCol.txt", delim=" ", col_names = FALSE)
colnames(regdb) <- c("chr", "pos", "SNP", "RegDB_Score")
# big file takes up too much memory
rm(regdb)
cg.regDB <- regdb %>%
    right_join(x) %>%
    distinct() %>%
    group_by(RegDB_Score, pval_cat) %>%
    summarise(Count=n()) %>%
    ungroup() %>%
    rename(RegulomeDB_Score = RegDB_Score) %>%
    mutate(Freq=Count/sum(Count))
write.table(cg.regDB,
            file="~/Google Drive/OSU_PHD/dissertation/code/chapter2/cg_regulomedb.txt",
            sep="\t", 
            quote=FALSE, 
            row.names=FALSE, 
            col.names=TRUE)
cg.regDB <- read.table("~/Google Drive/OSU_PHD/dissertation/code/chapter2/cg_regulomedb.txt",
                       header = TRUE,
                       sep="\t")
p <- cg.regDB %>%
    na.omit() %>%
    ggplot(aes(x=fct_rev(RegulomeDB_Score), y=Count, fill=pval_cat) ) +
    geom_col() +
    scale_fill_manual(values=c("black", "gray70")) +
    guides(fill=guide_legend(title=NULL), cex=1) + 
    xlab("RegulomeDB Score") +
    ylab("Number of previously reported SNPs") + 
    labs(fill="Reported Association to Survival") +
    coord_flip() + 
    theme_bw() +
    theme(legend.position = c(0.7, 0.9), 
          legend.text=element_text(size=14),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16))
p
```

 \normalsize
