#!/bin/bash

echo "Bash version ${BASH_VERSION}"

model_folder="./MODEL"
this_folder="/workdir/carvalho/RF_Test/BIG_TEST/WEAK_Scaling"
#declare -a nPROCS=(9 8 24 24 24 24)
#declare -a nSELECT=(3 8 6 9 15 22)
#declare -a nPROCSTOTAL=(27 64 125 216 343 512)
declare -a nPROCS=(1 8 9 8 24 24 24 24)
declare -a nSELECT=(1 1 3 8 6 9 15 22)
declare -a nPROCSTOTAL=(1 8 27 64 125 216 343 512)
declare -a INITIAL_SIZE=('63.75' '117.25' '180.75' '239.25' '297.75' '356.25' '414.75' '473.25') #Sizes OK
size_decrease=10
mem_per_proc_mb=2500
tests_folder="TESTS"
older_prefix="WEAK"
count_max=10 #Number of runs for each case

#READING ARGUMENTS
if [ $# -gt 0  ]; then
    it_procs=$1
else
    it_procs=0
fi
printf "it_procs     =$it_procs\n"

if [ $# -gt 1  ]; then
    it_size=$2
else
    it_size=$it_procs
fi
printf "it_size     =$it_size\n"

if [ $# -gt 2  ]; then
    count=$3
else
    count=1
fi
printf "count     =$count\n"

if [ $count -gt  $count_max ]; then
    printf "count=count_max NO MORE ITERATIONS\n"
    exit 2
fi

#DEFINING PARAMETERS
np=${nPROCS[$it_procs]}
ns=${nSELECT[$it_procs]}
npt=${nPROCSTOTAL[$it_procs]}
nmem=$(expr $mem_per_proc_mb \* $np)
sz=${INITIAL_SIZE[$it_size]}
printf "nPROCS     =$np\n"
printf "nSELECT    =$ns\n"
printf "nPROCSTOTAL=$npt\n"
printf "nMEMORY    =$nmem\n"
printf "size       =$sz\n"
printf "==========================\n"

#CREATING FOLDERS
#echo "$np"
folder_name=$(printf "$folder_prefix_P%03d" $npt)
#echo "$folder_name"
cd $this_folder
mkdir -p $tests_folder/$folder_name
cp $model_folder/* $tests_folder/$folder_name/

#REPLACING VARIABLES INSIDE EACH FOLDER
echo "$tests_folder/$folder_name"
cd $tests_folder/$folder_name
    to_find=#NPROCS#
    replace_by=$np
    in_file="run.pbs"
    sed -i -e "s/$to_find/$replace_by/g" $in_file
    to_find=#NSELECT#
    replace_by=$ns
    in_file="run.pbs"
    sed -i -e "s/$to_find/$replace_by/g" $in_file
    to_find=#MEMmb#
    replace_by=$nmem
    in_file="run.pbs"
    sed -i -e "s/$to_find/$replace_by/g" $in_file
    to_find=#NPROCSTOTAL#
    replace_by=$npt
    in_file="run.pbs"
    sed -i -e "s/$to_find/$replace_by/g" $in_file
    to_find=#SIZE#
    replace_by=$sz
    in_file="input_ScaRL.txt"
    sed -i -e "s/$to_find/$replace_by/g" $in_file
#cd ~-

#RUNNING CASE
    job_id=$(qsub run.pbs)
    qsub -v \
         it_procs=$it_procs,it_size=$it_size,count=$count \
         -W depend=afterany:$job_id verify.pbs
cd ~-
