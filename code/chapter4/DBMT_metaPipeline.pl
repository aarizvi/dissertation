#! /usr/bin/perl
use warnings;
use strict;
use Cwd;
use Getopt::Long;

my $manifest;
my $walltime;
my $email;
my $rscript;

GetOptions ("m=s" => \$manifest, "w=s" => \$walltime, "e=s" => \$email, "rscript=s" => \$rscript) or die$!;

my $wd = getcwd;

if (!-d "$wd/job_scripts") { mkdir "$wd/job_scripts"; }
if (!-d "$wd/log") { mkdir "$wd/log"; }
if (!-d "$wd/job_scripts/metal_scripts") { mkdir "$wd/job_scripts/metal_scripts"; }
if (!-d "$wd/temp") { mkdir "$wd/temp"; }
if (!-d "$wd/out") { mkdir "$wd/out"; }

open MANIFEST, "$manifest" or die$!;

my $line = 1;

while(<MANIFEST>) {

	if ($line == 1) {
		$line++;
		next;
	}

	chomp($_);
	my @temparray = split("\t", $_);

	my $vcf = "$temparray[0]";
	my $out = $temparray[1];
	my $ptSubset = $temparray[2];
	my $genome = $temparray[3];
	my $outcome = $temparray[4];
	my $memory = $temparray[5];

	open JOB1, ">${out}_1.sh" or die$!;

	print JOB1 "\#!/bin/bash\n\#SBATCH --time=$walltime\n\#SBATCH --nodes=1\n\#SBATCH --mem=$memory\n\#SBATCH --mail-user=$email\n\#SBATCH --ntasks-per-node=4";
	print JOB1 "\n\#SBATCH --mail-type=END\n\#SBATCH --partition=general-compute --qos general-compute\n\#SBATCH--job-name=$out\_c1\n\#SBATCH --output=$wd/log/\%j\_$out.out\n\#SBATCH --error=$wd/log/\%j\_$out.err\n\n";
	print JOB1 "\#Get date and time\ntstart=\$(date +\%s)\necho \"\#\#\#\#\#\# start time:\"`date`\n\n";
	print JOB1 "echo \"SLURM_JOB_ID\"=\$SLURM_JOB_ID\necho \"SLURM_JOB_NODELIST\"=\$SLURM_JOB_NODELIST\necho \"SLURM_NNODES\"=\$SLURM_NNODES\necho \"SLURM_NTASKS\"=\$SLURM_NTASKS\n";
	print JOB1 "echo \"SLURMTMPDIR\"=\$SLURMTMPDIR\necho \"working directory\"=\$SLURM_SUBMIT_DIR\n\necho \"************************\"\n\n";
	print JOB1 "NPROCs=`srun --nodes=\${SLURM_NNODES} bash -c 'hostname' | wc -l`\necho NPROCS=\$NPROCS\n";
	print JOB1 "\nmodule load R/3.4.1\n\nR --file=$rscript -q --args cohort 1 vcf.file $vcf out.file $wd/$out ptSubset $ptSubset genome $genome outcome $outcome ncores \$SLURM_NTASKS\n\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n";
	close JOB1;

	system "mv ${out}_1.sh $wd/job_scripts/\n";
	my $job1 =`sbatch $wd/job_scripts/${out}_1.sh`;
	chomp $job1;
	my @jobid = split(" ", $job1);
	$job1 = $jobid[3];
	print "$out Cohort 1 job submitted $job1 \n";


	open JOB2, ">${out}_2.sh" or die$!;

	print JOB2 "\#!/bin/bash\n\#SBATCH --time=$walltime\n\#SBATCH --nodes=1\n\#SBATCH --mem=$memory\n\#SBATCH --mail-user=$email\n\#SBATCH --ntasks-per-node=4";
	print JOB2 "\n\#SBATCH --mail-type=END\n\#SBATCH --partition=general-compute --qos=general-compute\n\#SBATCH --job-name=$out\_c2\n\#SBATCH --output=$wd/log/\%j\_$out.out\n\#SBATCH --error=$wd/log/\%j\_$out.err\n\n";
	print JOB2 "\#Get date and time\ntstart=\$(date +\%s)\necho \"\#\#\#\#\#\# start time:\"`date`\n\n";
	print JOB2 "echo \"SLURM_JOB_ID\"=\$SLURM_JOB_ID\necho \"SLURM_JOB_NODELIST\"=\$SLURM_JOB_NODELIST\necho \"SLURM_NNODES\"=\$SLURM_NNODES\necho \"SLURM_NTASKS\"=\$SLURM_NTASKS\n";
	print JOB2 "echo \"SLURMTMPDIR\"=\$SLURMTMPDIR\necho \"working directory\"=\$SLURM_SUBMIT_DIR\n\necho \"************************\"\n\n";
	print JOB2 "NPROCs=`srun --nodes=\${SLURM_NNODES} bash -c 'hostname' | wc -l`\necho NPROCS=\$NPROCS\n";
	print JOB2 "\nmodule load R/3.4.1\n\nR --file=$rscript -q --args cohort 2 vcf.file $vcf out.file $wd/$out ptSubset $ptSubset genome $genome outcome $outcome ncores \$SLURM_NTASKS\n\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n";
	close JOB2;

	system "mv ${out}_2.sh $wd/job_scripts/\n";
	my $job2 =`sbatch $wd/job_scripts/${out}_2.sh`;
	chomp $job2;
	@jobid = split(" ", $job2);
	$job2 = $jobid[3];
	print "$out Cohort 2 job submitted $job2 \n";

	open JOB3, ">${out}_3.sh" or die$!;

	print JOB3 "\#!/bin/bash\n\#SBATCH --time=$walltime\n\#SBATCH --nodes=1\n\#SBATCH --mem=8000\n\#SBATCH --mail-user=$email";
	print JOB3 "\n\#SBATCH --mail-type=END\n\#SBATCH --partition=general-compute --qos=general-compute\n\#SBATCH --job-name=$out\_meta\n\#SBATCH --output=$wd/log/\%j\_$out.out\n\#SBATCH --error=$wd/log/\%j\_$out.err\n\n";
	print JOB3 "\#Get date and time\ntstart=\$(date +\%s)\necho \"\#\#\#\#\#\# start time:\"`date`\n\n";
	print JOB3 "echo \"SLURM_JOB_ID\"=\$SLURM_JOB_ID\necho \"SLURM_JOB_NODELIST\"=\$SLURM_JOB_NODELIST\necho \"SLURM_NNODES\"=\$SLURM_NNODES\n";
	print JOB3 "echo \"SLURMTMPDIR\"=\$SLURMTMPDIR\necho \"working directory\"=\$SLURM_SUBMIT_DIR\n\necho \"************************\"\n\n";
	print JOB3 "NPROCs=`srun --nodes=\${SLURM_NNODES} bash -c 'hostname' | wc -l`\necho NPROCS=\$NPROCS\n";


	print JOB3 "awk -F \"\\t\" 'OFS=\";\" { print \$1,\$2,\$3,\$4,\$5,\$6}' $wd/$out\_c1.coxph > $wd/$out\_first6.txt\n";
	print JOB3 "awk -F \"\\t\" 'OFS=\"\\t\" {for (i=5; i<NF; i++) printf \$i \"\\t\"; print \$NF}' $wd/$out\_c1.coxph > $wd/$out\_last.txt\n";
	print JOB3 "paste -d \"\\t\" $wd/$out\_first6.txt $wd/$out\_last.txt > $wd/$out\_c1.formetal\n";

	print JOB3 "awk -F \"\\t\" 'OFS=\";\" { print \$1,\$2,\$3,\$4,\$5,\$6}' $wd/$out\_c2.coxph > $wd/$out\_first6.txt\n";
	print JOB3 "awk -F \"\\t\" 'OFS=\"\\t\" {for (i=5; i<NF; i++) printf \$i \"\\t\"; print \$NF}' $wd/$out\_c2.coxph > $wd/$out\_last.txt\n";
	print JOB3 "paste -d \"\\t\" $wd/$out\_first6.txt $wd/$out\_last.txt > $wd/$out\_c2.formetal\n";

	#create the metal_arguments.txt file for this job - these are distinguished by the variable $out, so that needs to be unique for each run.
	open METAL, ">$wd/job_scripts/metal_scripts/metal_arguments_$out.txt" or die$!;
	print METAL "\#Meta analysis for cohort 1 and cohort 2\nSCHEME STDERR\nSEPARATOR TAB\nMARKER RSID;TYPED;CHR;POS;REF;ALT\n";
	print METAL "ALLELE REF ALT\nEFFECT COEF\nSTDERR SE.COEF\nPVALUE PVALUE\nWEIGHT N\nPROCESS $wd/$out\_c1.formetal\n";
	print METAL "PROCESS $wd/$out\_c2.formetal\nOUTFILE $wd/$out\_metal_out. .tbl\nMINWEIGHT 10\nANALYZE HETEROGENEITY\n\nEXIT\n";
	close METAL;

	#lines to run metal for this $out
	print JOB3 "module load metal\n\n";
	print JOB3 "metal $wd/job_scripts/metal_scripts/metal_arguments_$out.txt\n";

	#clean up the name of the metal outfile
	print JOB3 "cp $wd/$out\_metal_out*tbl $wd/$out.tbl\n";

	#end up with $out.tbl, $out_c1.formetal, and $out_c2.formetal - put these together in R
	print JOB3 "module load R/3.4.1\n";
	print JOB3 "R --file=/projects/rpci/lsuchest/lsuchest/DBMT_metaPipeline/spread_metal.R -q --args metal_result $wd/$out.tbl for_metal_cohort1 $wd/$out\_c1.formetal for_metal_cohort2 $wd/$out\_c2.formetal full_output $wd/$out.res\n\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n";

	print JOB3 "rm $wd/$out\_c1_temp.formetal\n";
	print JOB3 "rm $wd/$out\_c2_temp.formetal\n";
	print JOB3 "rm $wd/$out\_first6.txt\n";
	print JOB3 "rm $wd/$out\_last.txt\n";
	print JOB3 "mv $wd/$out\_c1* $wd/temp\n";
	print JOB3 "mv $wd/$out\_c2* $wd/temp\n";
	print JOB3 "mv $wd/$out\_metal* $wd/temp\n";
	print JOB3 "mv $wd/$out.tbl $wd/temp\n";
	print JOB3 "mv $wd/$out.res $wd/out\n";
	close JOB3;

	system "mv ${out}_3.sh $wd/job_scripts/\n";
	my $job3 =`sbatch --dependency=afterok:${job1}:${job2} $wd/job_scripts/${out}_3.sh`;
	chomp $job3;
        @jobid = split(" ", $job3);
        $job3 = $jobid[3];
	print "$out METAL script for cohorts 1 and 2 submitted $job3 \n";
}
