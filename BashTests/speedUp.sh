#!/bin/bash

echo "Bash version ${BASH_VERSION}"

model_folder="./MODEL"
this_folder="/workdir/carvalhol/RF_Test/BIG_TEST/AUTOTESTS"
#declare -a nPROCS=(9 8 24 24 24 24)
#declare -a nSELECT=(3 8 6 9 15 22)
#declare -a nPROCSTOTAL=(27 64 125 216 343 512)
declare -a nPROCS=(1 2 4 8 9 8 24 24 24 24)
declare -a nSELECT=(1 1 1 1 3 8 6 9 15 22)
declare -a nPROCSTOTAL=(1 2 4 8 27 64 125 216 343 512)
declare -a INITIAL_SIZE=(120 200 200 300 350)
size_decrease=10
mem_per_proc_mb=2500
tests_folder="TESTS"
folder_prefix="BIG"

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
    it_size=0
fi
printf "it_size     =$it_size\n"

#DEFINING PARAMETERS
np=${nPROCS[$it_procs]}
ns=${nSELECT[$it_procs]}
npt=${nPROCSTOTAL[$it_procs]}
nmem=$(expr $mem_per_proc_mb \* $np)
sz=$(expr ${INITIAL_SIZE[$it_procs]} - $it_size \* $size_decrease)
printf "nPROCS     =$np\n"
printf "nSELECT    =$ns\n"
printf "nPROCSTOTAL=$npt\n"
printf "nMEMORY    =$nmem\n"
printf "size       =$sz\n"
printf "==========================\n"

#CREATING FOLDERS
#echo "$np"
folder_name=$(printf "$folder_prefix-P%03d-SZ%03d" $npt $sz)
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
#cd $this_folder
#cd $tests_folder/$folder_name; 
    job_id=$(qsub run.pbs)
    qsub -v \
         it_procs=$it_procs,it_size=$it_size \
         -W depend=afterany:$job_id verify.pbs
cd ~-


##for ((sz = 500; sz > 190 ; sz=sz-10));
##do
##    sizes[$ii]=$sz
##    ii=$ii+1
##done
#ii=0;
#
##echo ${sizes[*]}
#
##START BASH FILE
#bash_run=run_all_pbs.sh
#rm $bash_run
#touch $bash_run
#chmod u+x $bash_run
#printf "#!/bin/bash\n\n" >> $bash_run
#
#bash_report=make_report.sh
#rm $bash_report
#touch $bash_report
#chmod u+x $bash_report
#printf "#!/bin/bash\n\n" >> $bash_report
#printf "instruction='tail -n2 out_run.txt'\n" >> $bash_report
#
##MAKING CASES
#N=${#nPROCS[@]}
#mkdir $tests_folder
##for np in "${nPROCS[@]}"
#for ((ii = 0; ii < N; ii++));
#do
#    printf "echo; echo ----------------------------\n" >> $bash_report
#    printf "echo; echo ----------------------------\n" >> $bash_report
#    printf "echo $ii\) nPROCS=$np, nSELECT=$ns, nPROCSTOTAL=$npt, nMEMORY=$nmem\n\n" >> $bash_run
#    np=${nPROCS[$ii]}
#    ns=${nSELECT[$ii]}
#    npt=${nPROCSTOTAL[$ii]}
#    nmem=$(expr $mem_per_proc_mb \* $np)
#    first=1
#    printf "echo nPROCS     =$np\n"  >> $bash_report
#    printf "echo nSELECT    =$ns\n"  >> $bash_report
#    printf "echo nPROCSTOTAL=$npt\n" >> $bash_report
#    printf "echo nMEMORY    =$nmem\n" >> $bash_report
#    printf "echo; echo ----------------------------\n" >> $bash_report
#
#    echo; echo "-------------------------------------"    
#    echo "$ii) nPROCS=$np, nSELECT=$ns, nPROCSTOTAL=$npt, nMEMORY=$nmem"
#    
#    for sz in "${sizes[@]}"
#    do
#
#        #CREATING FOLDERS
#        #echo "$np"
#        folder_name=$(printf "$folder_prefix-P%03d-SZ%03d" $npt $sz)
#        #echo "$folder_name"
#        #rm -r $folder_name
#        cp -r $model_folder $tests_folder/$folder_name
#        
#        #REPLACING VARIABLES INSIDE EACH FOLDER
#        echo "$tests_folder/$folder_name"
#        cd $tests_folder/$folder_name
#            to_find=#NPROCS#
#            replace_by=$np
#            in_file="run.pbs"
#            sed -i -e "s/$to_find/$replace_by/g" $in_file
#            to_find=#NSELECT#
#            replace_by=$ns
#            in_file="run.pbs"
#            sed -i -e "s/$to_find/$replace_by/g" $in_file
#            to_find=#MEMmb#
#            replace_by=$nmem
#            in_file="run.pbs"
#            sed -i -e "s/$to_find/$replace_by/g" $in_file
#            to_find=#NPROCSTOTAL#
#            replace_by=$npt
#            in_file="run.pbs"
#            sed -i -e "s/$to_find/$replace_by/g" $in_file
#            to_find=#SIZE#
#            replace_by=$sz
#            in_file="input_ScaRL.txt"
#            sed -i -e "s/$to_find/$replace_by/g" $in_file
#        cd ~-
#    
#        if [ $first -eq 1  ]; then
#            echo FIRST!!
#            first=0
#            printf "OR=\$(pwd)\n" >> $bash_run
#            printf "cd \$OR; cd $tests_folder/$folder_name;job_id=\$(qsub run.pbs)\n" >> $bash_run
#            printf "         qsub -v \\ \n" >> $bash_run
#            printf "         NEXT=$next_folder,END_FILE=\$END_FILE,COUNTER=\$COUNTER,PATH_FOLDER=\$PATH_FOLDER \\ \n" >> $bash_run
#            printf "         -W depend=afterany:$job_id verify.pbs\n" >> $bash_run
#        else
#            echo OTHER!
#        fi
#
#        #MAKE BASH TO RUN ALL PBS FILES AND CREATE REPORT
#        printf "(cd $tests_folder/$folder_name; qsub run.pbs)\n" >> $bash_run
#        printf "echo ------------ size = $sz ------------\n" >> $bash_report
#        printf "echo \(cd $tests_folder/$folder_name\; qsub run.pbs\)\n" >> $bash_report
#        printf "(cd $tests_folder/$folder_name; \$instruction )\n" >> $bash_report
#    done
#done
