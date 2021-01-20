#PBS -N AML
#PBS -l nodes=1:ppn=1,walltime=6:00:00
#PBS -M dinglei@indiana.edu


module load nlopt
module load curl
module load java/1.7.0_51
module load r/3.6.0
cd Desktop/Lei/AIMER2/ISMB/real/AML/
R CMD BATCH --no-save --no-restore analysis.R 
