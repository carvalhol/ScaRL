#!/bin/bash
#SBATCH --nodes=#NNODES                      #Numero de Nós
#SBATCH --ntasks-per-node=#NTASKSNODE        #Numero de tarefas por Nó
#SBATCH --ntasks=#NMPI                     #Numero total de tarefas MPI
#SBATCH -p cpu_dev                      #Fila (partition) a ser utilizada
#SBATCH -J ScaRL          #Nome job
#SBATCH --exclusive                    #Utilização exclusiva dos nós durante a execução do job
#SBATCH --time=40:00
#SBATCH -o hostname_%j.out      # File to which STDOUT will be written
#SBATCH -e hostname_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=lucio.a.c@gmail.com # Email to which notifications will be sent

#Exibe os nós alocados para o Job
echo $SLURM_JOB_NODELIST
nodeset -e $SLURM_JOB_NODELIST

cd $SLURM_SUBMIT_DIR
echo $SLURM_SUBMIT_DIR
#ou
#intel MPI
source /scratch/app/modulos/intel-psxe-2017.sh 
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
#module load hdf5/1.8.14_intel
#module load fftw/3.3.5_intel
export FFTW3_ROOT=/scratch/msgeo/lucio.correa/LOCAL/FFT

#export I_MPI_EXTRA_FILESYSTEM=on
#export I_MPI_EXTRA_FILESYSTEM_LIST=lustre

export LD_LIBRARY_PATH=/scratch/msgeo/lucio.correa/LOCAL/FFTW/lib:$LD_LIBRARY_PATH

#Configura o executavel
EXEC=/scratch/msgeo/lucio.correa/ScaRL/build/ScaRL.exe

#exibe informações sobre o executável
/usr/bin/ldd $EXEC

#printenv 

pwd
#ls
echo "$EXEC"
#ls "$EXEC"

#ls -lArt > list_of_files.txt
srun -n "$SLURM_NTASKS" "$EXEC"
