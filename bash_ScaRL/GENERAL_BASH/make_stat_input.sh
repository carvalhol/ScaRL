#!/bin/bash

output_Name="input_stat.txt"
arg=$1 

echo "Argument " $arg 
 
if [[ -n "$arg"  ]]; then 
    echo "MAKING STAT INPUT FROM: $arg" 
else 
    echo "STOP - No argument, select the folder" 
    echo "SYNTAX: ./make_stat_input FOLDER"
    exit -1
fi

folder=$arg

cd $folder
cases=($(ls -d OVER*.h5))
path=$(pwd)
cd ~-

printf '%s\n' "${cases[@]}"
nCases=${#cases[@]}

rm $output_Name
touch $output_Name

printf "$nCases #Number of Samples\n" >> $output_Name
printf "1 #Calculate Correlation Legth\n" >> $output_Name
printf "0 #Delete Sample In The End\n" >> $output_Name
for ((cc = 0; cc < nCases; cc++)); do
    printf "\"$path/${cases[$cc]}\"\n" >> $output_Name  
    printf "\"samples\"\n" >> $output_Name  
done
