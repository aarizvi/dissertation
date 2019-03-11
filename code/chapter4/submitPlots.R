declare -a arr=("D_ALLonly_DRM"
                "D_ALLonly_TRM"
                "D_ALLonly_PFS"
                "D_ALLonly_OS"
                "R_ALLonly_DRM"
                "R_ALLonly_TRM"
                "R_ALLonly_PFS"
                "R_ALLonly_OS")

for patt in "${arr[@]}";
do

cat <<EOM > ${patt}.sh
#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --mem=32000
#SBATCH --mail-user=rizvi.33@osu.edu
#SBATCH --mail-type=END
#SBATCH --partition=general-compute --qos=general-compute
#SBATCH --job-name=${patt}_100d_Plots
#SBATCH --output=/projects/rpci/lsuchest/abbasriz/DBMT_100d/analyses/ALLonly/plots/log/%j_${patt}.out
#SBATCH --error=/projects/rpci/lsuchest/abbasriz/DBMT_100d/analyses/ALLonly/plots/log/%j_${patt}.err

DIR=/projects/rpci/lsuchest/abbasriz/DBMT_100d/analyses/ALLonly

module load R/3.4.1

R --file=\$DIR/plots/manPlot.R -q --args input_path \$DIR/out patt ${patt} output_path \$DIR/plots/manPlots_Feb2019

exit
EOM

echo -e "\tsubmitting file: ${patt}\n";
sbatch ${patt}.sh;

done
