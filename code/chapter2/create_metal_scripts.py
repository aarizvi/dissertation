#!/usr/bin/python

import glob
import io

# grab file names
files_c1 = sorted(glob.glob("/projects/rpci/lsuchest/abbasriz/candidate_gene/result_files_cg/*c1*.txt"))

files_c2 = sorted(glob.glob("/projects/rpci/lsuchest/abbasriz/candidate_gene/result_files_cg/*c2*.txt"))


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
    	f.write("# META ANALYSIS FOR COHORT 1 and COHORT 2 FROM DISCOVERY-BMT CANDIDATE GENE REPLICATION/VALIDATION SUBSET")
    	f.write("\n")
    	f.write("\n# THE RESULTS ARE STORED IN FILES metal_"+meta_file_names[i]+".tbl and metal_"+meta_file_names[i]+".tbl.info")
    	f.write("\n")
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
	f.write("\nOUTFILE /projects/rpci/lsuchest/abbasriz/candidate_gene/result_files_cg/metal_results/metal_"+meta_file_names[i] + " .tbl")
	f.write("\nMINWEIGHT 10")
	f.write("\nANALYZE")


