#!/bin/bash
n=3
for xy in $(eval echo {0..$(($n-1))}{0..$(($n-1))}
do
cat <<EOF>xy$xy
#!/bin/bash
#SBATCH --job-name xy$xy  
#SBATCH -N1 --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-socket=1
#SBATCH --time=24:00:00
#SBATCH --account=IscrC_PNAsens
#SBATCH --partition=m100_usr_prod
###SBATCH --qos=m100_qos_dbg
##SBATCH --partition=m100_all_serial
#load needed modules
module purge 
module load profile/chem-phys
module load anaconda
conda init bash
source /m100/home/userexternal/gsacco00/.conda/envs/racfenv/bin/activate
conda run -n racfenv python3 ../mfield.py --dir . --rcut 200 -n 10000 20000 -ng $n -xy $xy 
EOF
sbatch xy$xy
mv xy$xy jobs
done
