#!/bin/bash
for i in 0 5000 10000 15000
do
declare -i i=$i
declare -i j=0
j=$(( i + 5000 ))
cat <<EOF>mf$i
#!/bin/bash
#SBATCH --job-name mf$i  
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
conda run -n racfenv python3 ../mfield.py --dir . --rcut 200 -n $i $j
EOF
sbatch mf$i
mv mf$i jobs
done
