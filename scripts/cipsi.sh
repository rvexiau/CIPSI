#!/bin/sh
#SBATCH --job-name=Cipsi
# #SBATCH --partition=300m			# submission queue
#SBATCH --time=7-20:00:00			# 1-1 means one day and one hour
#SBATCH --mail-type=END	
#SBATCH --mail-user=romain.vexiau@u-psud.fr	
#SBATCH --output=./lumat.stdout		# if --error is absent, includes also the errors
#SBATCH --mem=28G 				# T-tera, G-giga, M-mega
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1


echo "------------------------------------------------------------------------------"
echo "hostname                     =   $(hostname)"
echo "SLURM_JOB_NAME               =   $SLURM_JOB_NAME"
echo "SLURM_SUBMIT_DIR             =   $SLURM_SUBMIT_DIR"
echo "SLURM_JOBID                  =   $SLURM_JOBID"
echo "SLURM_JOB_ID                 =   $SLURM_JOB_ID"
echo "SLURM_NODELIST               =   $SLURM_NODELIST"
echo "SLURM_JOB_NODELIST           =   $SLURM_JOB_NODELIST"
echo "SLURM_TASKS_PER_NODE         =   $SLURM_TASKS_PER_NODE"
echo "SLURM_JOB_CPUS_PER_NODE      =   $SLURM_JOB_CPUS_PER_NODE"
echo "SLURM_TOPOLOGY_ADDR_PATTERN  =   $SLURM_TOPOLOGY_ADDR_PATTERN"
echo "SLURM_TOPOLOGY_ADDR          =   $SLURM_TOPOLOGY_ADDR"
echo "SLURM_CPUS_ON_NODE           =   $SLURM_CPUS_ON_NODE"
echo "SLURM_NNODES                 =   $SLURM_NNODES"
echo "SLURM_JOB_NUM_NODES          =   $SLURM_JOB_NUM_NODES"
echo "SLURMD_NODENAME              =   $SLURMD_NODENAME"
echo "SLURM_NTASKS                 =   $SLURM_NTASKS"
echo "SLURM_NPROCS                 =   $SLURM_NPROCS"
echo "SLURM_MEM_PER_NODE           =   $SLURM_MEM_PER_NODE"
echo "SLURM_PRIO_PROCESS           =   $SLURM_PRIO_PROCESS"
echo "------------------------------------------------------------------------------"

# USER Commands
# special commands for openmpi/intel
module load intel/14.0.5
export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE

cp *.dat $TMPDIR/
cp Autocip13.exe $TMPDIR/
cp auto.in $TMPDIR/
cd $TMPDIR

# export CIPSI_WORK=$TMPDIR 
./Autocip13.exe<auto.in

find . -name "*_mat*" -type f -delete
find . -name "*fort*" -type f -delete
find . -name "*BD_SCRATCH*" -type f -delete
find . -name "COMB_VCPP" -type f -delete
find . -name "*_ijkl*" -type f -delete
find . -name "*_f50*" -type f -delete
find . -name "*_scr*" -type f -delete
find . -name "*_pqrs*" -type f -delete   
find . -name "*_f25*" -type f -delete   
find . -name "*_det*" -type f -delete
find . -name "*_bda*" -type f -delete
find . -name "F10_*" -type f -delete
find . -name "*_dcla*" -type f -delete
find . -name "*_dspina*" -type f -delete      

result=`sed -n 4p auto.in`
tar -zcf lumat.tar.gz *
    
mv lumat.tar.gz /home/romain_vexiau/results/Cipsi/$result

# end of the USER commands
