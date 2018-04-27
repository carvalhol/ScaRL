#!/bin/bash

echo "Bash version ${BASH_VERSION}"

test_type=1


if [ $test_type -eq 1 ]; then
#FINDING MAXIMUM SIZE WE CAN RUN WHEN FIXING THE MEMORY BY PROC
declare -a nPROCS=(1 2 4 8 9 8 24 24 24 24)
declare -a nSELECT=(1 1 1 1 3 8 6 9 15 22)
declare -a nPROCSTOTAL=(1     2   4   8   27  64 125 216 343 512)
declare -a INITIAL_SIZE=(90.0 130.0 150.0 190.0 310.0 380.0 540.0 550.0 600.0 700.0) #DIMENSIONED FOR ONE OK THEN  THEN SECOND FAIL
folder_prefix="MAX"
queue="haswellq"
mem_per_proc_mb=2500
walltime="01:00:00"
max_count=20 #MAXIMAL NUMBER OF ITERATION IF IT DOES NOT FAIL
size_increment=10.0


elif [ $test_type -eq 2 ]; then
#FINDING MAXIMUM SIZE WE CAN RUN MONOPROC
declare -a nPROCS=(1)
declare -a nSELECT=(1)
declare -a nPROCSTOTAL=(1)
declare -a INITIAL_SIZE=(63.75)
folder_prefix="monoMAX"
queue="uvq"
mem_per_proc_mb=250000
walltime="24:00:00"
max_count=50 #MAXIMAL NUMBER OF ITERATION IF IT DOES NOT FAIL
size_increment=53.5 #SO WE HAVE THE SAME VALUES THAT FOR WEAK SCALING


elif [ $test_type -eq 3 ]; then
#WEAK SCALING MPI
declare -a nPROCS=(1 8 9 8 24 24 24 24)
declare -a nSELECT=(1 1 3 8 6 9 15 22)
declare -a nPROCSTOTAL=(1 8 27 64 125 216 343 512)
declare -a INITIAL_SIZE=(63.75 117.25 180.75 239.25 297.75 356.25 414.75 473.25) #Sizes OK and power of two
folder_prefix="WEAK"
queue="haswellq"
mem_per_proc_mb=2500
walltime="01:00:00"
max_count=20 #MAXIMAL NUMBER OF ITERATION IF IT DOES NOT FAIL
size_increment=0.0 #WE WILL RUN 20 TIMES EACH WITH SAME SIZE


elif [ $test_type -eq 4 ]; then
#WEAK SCALING MONOPROC
declare -a nPROCS=(1 1 1 1 1 1 1 1)
declare -a nSELECT=(1 1 1 1 1 1 1 1)
declare -a nPROCSTOTAL=(1 1 1 1 1 1 1 1)
declare -a INITIAL_SIZE=(63.75 117.25 180.75 239.25 297.75 356.25 414.75 473.25) #Sizes OK and power of two
folder_prefix="monoWEAK"
queue="uvq"
mem_per_proc_mb=250000
walltime="24:00:00"
max_count=20 #MAXIMAL NUMBER OF ITERATION IF IT DOES NOT FAIL
size_increment=0.0 #WE WILL RUN 20 TIMES EACH WITH SAME SIZE

else
    echo "test_type not implemented"
    exit -1
fi

N=${#nPROCS[@]} #NUMBER OF CASES TO CALCULATE
N=1 #TEST


#######################################################
#######################################################

this_folder=$(pwd)
model_folder=$(printf "$this_folder/MODEL")
tests_folder=$(printf "$this_folder/RESULTS")
inputs_folder=$(printf "$this_folder/INPUTS")
exec_path="/home/carvalho/ScaRL/build/ScaRL.exe"
out_file="out_run.txt"
exit_when_ok=0 #1 for true /0 for false


#CREATING FOLDERS
#echo "$np"
mkdir $inputs_folder
mkdir $tests_folder

##MAKING CASES
run_f=$(printf $(pwd)"/run_cases_$folder_prefix.sh")
rm $run_f
touch $run_f
chmod u+x $run_f

for ((pp = 0; pp < N; pp++)); do
    #DEFINING PARAMETERS
    np=${nPROCS[$pp]}
    ns=${nSELECT[$pp]}
    npt=${nPROCSTOTAL[$pp]}
    init_size=${INITIAL_SIZE[$pp]}
    int_init_size=${init_size%.*} 
    nmem=$(expr $mem_per_proc_mb \* $np)
    printf "nPROCS     =$np\n"
    printf "nSELECT    =$ns\n"
    printf "nPROCSTOTAL=$npt\n"
    printf "nMEMORY    =$nmem\n"
    printf "INITIAL_SIZE    =$init_size\n"
    printf "==========================\n"
    model_proc=$(printf "$inputs_folder/$folder_prefix-P%03d-SZ%03d" $npt $int_init_size)
    mkdir -p $model_proc
    cp $model_folder/* $model_proc/

    printf "echo nPROCS     =$np\n" >> $run_f

    #REPLACING VARIABLES INSIDE EACH FOLDER
    echo "$model_proc"
    cd $model_proc 
        to_find=#NPROCS#
        replace_by=$np
        in_file="run.pbs"
        sed -i -e "s/$to_find/$replace_by/g" $in_file
        in_file="run_all.pbs"
        sed -i -e "s/$to_find/$replace_by/g" $in_file
        
        to_find=#QUEUE#
        replace_by=$queue
        in_file="run.pbs"
        sed -i -e "s/$to_find/$replace_by/g" $in_file
        in_file="run_all.pbs"
        sed -i -e "s/$to_find/$replace_by/g" $in_file
        
        to_find=#WALLTIME#
        replace_by=$walltime
        in_file="run.pbs"
        sed -i -e "s/$to_find/$replace_by/g" $in_file
        in_file="run_all.pbs"
        sed -i -e "s/$to_find/$replace_by/g" $in_file
        
        to_find=#NSELECT#
        replace_by=$ns
        in_file="run.pbs"
        sed -i -e "s/$to_find/$replace_by/g" $in_file
        in_file="run_all.pbs"
        sed -i -e "s/$to_find/$replace_by/g" $in_file
        
        to_find=#MEMmb#
        replace_by=$nmem
        in_file="run.pbs"
        sed -i -e "s/$to_find/$replace_by/g" $in_file
        in_file="run_all.pbs"
        sed -i -e "s/$to_find/$replace_by/g" $in_file
       
        to_find=#NPROCSTOTAL#
        replace_by=$npt
        in_file="run.pbs"
        sed -i -e "s/$to_find/$replace_by/g" $in_file
        in_file="run_all.pbs"
        sed -i -e "s/$to_find/$replace_by/g" $in_file
        
        to_find=#JOBNAME#
        run_proc_name=$(printf "$folder_prefix-P%03d" $npt)
        replace_by=$run_proc_name
        in_file="run_all.pbs"
        sed -i -e "s/$to_find/$replace_by/g" $in_file
        
        to_find=#FOLDERPREFIX#
        replace_by=$folder_prefix
        in_file="run_all.pbs"
        sed -i -e "s/$to_find/$replace_by/g" $in_file
        
        to_find=#INIT_SIZE#
        replace_by=$init_size
        in_file="run_all.pbs"
        sed -i -e "s/$to_find/$replace_by/g" $in_file

        to_find=#OUTPUTFILE#
        run_proc_name=$(printf "out_run$folder_prefix-P%03d-iSZ%03d.txt" $npt $int_init_size)
        replace_by=$run_proc_name
        in_file="run_all.pbs"
        sed -i -e "s/$to_find/$replace_by/g" $in_file

        to_find=#INCREMENT#
        replace_by=$size_increment
        in_file="run_all.pbs"
        sed -i -e "s/$to_find/$replace_by/g" $in_file

        to_find=#MAX_COUNT#
        replace_by=$max_count
        in_file="run_all.pbs"
        sed -i -e "s/$to_find/$replace_by/g" $in_file
        
        to_find=#EXIT_WHEN_OK#
        replace_by=$exit_when_ok
        in_file="run_all.pbs"
        sed -i -e "s/$to_find/$replace_by/g" $in_file
        
        to_find=#MODEL_FOLDER#
        replace_by=$model_proc
        in_file="run_all.pbs"
        sed -i -e "s+$to_find+$replace_by+g" $in_file
        
        to_find=#TESTS_FOLDER#
        replace_by=$tests_folder
        in_file="run_all.pbs"
        sed -i -e "s+$to_find+$replace_by+g" $in_file
        
        to_find=#EXEC_PATH#
        replace_by=$exec_path
        in_file="run_all.pbs"
        sed -i -e "s+$to_find+$replace_by+g" $in_file
        
        
        to_find=#OUTPUTFILE#
        replace_by=$out_file
        in_file="run_all.pbs"
        sed -i -e "s/$to_find/$replace_by/g" $in_file
       
        run_proc_name=$(printf "run_all_$folder_prefix-P%03d-iSZ%03d.pbs" $npt $int_init_size)
        mv run_all.pbs $inputs_folder/$run_proc_name
        printf "(cd $inputs_folder; qsub $run_proc_name)\n" >> $run_f

    cd ~-
done
