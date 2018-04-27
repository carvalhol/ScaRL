#!/bin/bash

Delete_Results=0
Open_Output=1
output_Name="out_run.txt"
PBS_Name="run.pbs"
time_to_find="20"
Build_Path="/home/carvalho/ScaRL/build"
Result_Folder="SAMPLES"
arg=$1 
User="carvalhol"
echo "Argument " $arg 
 
if [[ -n "$arg"  ]]; then 
    echo "I've got an argument" 
else 
    echo "No argument, i'll use build" 
    arg="build"
    echo "OPTIONS: "
    echo "build - builds the program" 
    echo "roda - runs the program" 
    echo "all - build + roda" 
fi 
 
if [[ "$arg" == "roda" ]]; then 
     echo "And it is roda" 
     Build=0
     Run=1 
elif [[ "$arg" == "build" ]]; then 
     Build=1
     Run=0 
elif [[ "$arg" == "all" ]]; then 
     Build=1
     Run=1 
fi 



if [ "$Build" -eq "1" ]; then
    echo " "
    echo " "
    echo "Building--------"
    echo " "
    echo " "
    (cd $Build_Path; make all) 
fi

if [ "$Run" -eq "1" ]; then
    echo " "
    echo " "
    echo "Running--------" 
    echo " "
    echo " "
    rm $output_Name
    rm -r $Result_Folder
    qsub $PBS_Name
#    qsub -v NP=$NP,Run=$Run \
#    -S /bin/bash \
#    -N $list_Name \
#    -o $output_Name \
#    -j oe \
#    -l walltime=$W_TIME \
#    -l select=$N_SELECT:ncpus=$N_CPU:mpiprocs=$N_CPU:mem=$mem \
#    -q $Queue \
#    -P $project_name \
#    -M $Mail \
#    $PBS_Name
    qstat -u $User
fi

if [ "$Run" -eq "1" ]; then
if [ "$Open_Output" -eq "1" ]; then
    COUNTER=$time_to_find
    until [ $COUNTER -lt 1 ]; do
        let COUNTER-=1
        sleep 1
        if [ -f $output_Name ]; then
            let COUNTER=0
            echo ":) FILE FOUND!!"
            more $output_Name
        fi
    done
    if [ "$COUNTER" -lt 0 ]; then
        echo ":( I'VE BEEN SEARCHING FOR" $time_to_find "s AND I THE FILE ISN'T HERE YET"
    fi
fi
fi

# INSTRUCIONS TO CREATE A BASH FILE
# create a file like this "FILETEST.sh"
# change the mode of the file, so you can execute this: "chmod u+x FILETEST.sh"
# it can be useful to make things you do every day automatically by executing this file
# more about it on http://tldp.org/LDP/Bash-Beginners-Guide/html/
