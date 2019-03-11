declare -a arr=("D_ALLonly_DRM"
                "D_ALLonly_TRM"
                "D_ALLonly_PFS"
                "D_ALLonly_OS"
                "D_ALLonly_OS_3y"
                "D_ALLonly_REL"
                "D_ALLonly_REL_3y"
                "R_ALLonly_DRM"
                "R_ALLonly_TRM"
                "R_ALLonly_PFS"
                "R_ALLonly_OS"
                "R_ALLonly_OS_3y"
                "R_ALLonly_REL"
                "R_ALLonly_REL_3y")

for pattern in "${arr[@]}"
do

cat <<EOM > ${pattern}.sh
#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --mem=32000
#SBATCH --mail-user=rizvi.33@osu.edu
#SBATCH --mail-type=END
#SBATCH --partition=general-compute --qos=general-compute
#SBATCH --job-name=ALL_only_Plots
#SBATCH --output=/projects/rpci/lsuchest/abbasriz/DBMT_ALLonly/analyses/ALL_EA_results/plots/log/%j_ALLonly.out
#SBATCH --error=/projects/rpci/lsuchest/abbasriz/DBMT_ALLonly/analyses/ALL_EA_results/plots/log/%j_ALLonly.err

DIR=/projects/rpci/lsuchest/abbasriz/DBMT_ALLonly/analyses/ALL_EA_results

module load R

R --file=\$DIR/plots/manPlot.R -q --args \\
    input.path \$DIR/out \\
    pattern ${pattern} \\
    output.path \$DIR/plots/manPlots_Feb2019

EOM

echo -e "\tsubmitting file: ${pattern}\n";
sbatch ${pattern}.sh;

done