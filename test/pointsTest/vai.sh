#!/bin/bash

root_folder="."
src_dir="/workdir/carvalhol/SEM_Test/eurodyn2017/SRC_lcp" #Only important when copying pbs files
#src_dir="/workdir/carvalhol/SEM_Test/eurodyn2017/SRC_pipo" #Only important when copying pbs files

arg=$1
echo "Argument: " $arg

if [[ -n "$arg"  ]]; then
    echo "I've got an argument :-)"
else
    echo "No argument (it should be list/clean/copy/prep/mvmesh/gen/run/cleanRun), the syntax is: './vai.sh argument'"
fi


or_path="$(pwd)"
cd $root_folder
abs_path="$(pwd)"

#PRE-PROCESSING
if [[ "$arg" == "copy" ]]; then
    echo "Copying pbs files from: $src_dir"
elif [[ "$arg" == "prep" ]]; then
    rm $or_path/out_PBS_prep.txt.txt
elif [[ "$arg" == "gen" ]]; then
    rm $or_path/out_PBS_gen.txt
elif [[ "$arg" == "run" ]]; then
    rm $or_path/out_PBS_run.txt
elif [[ "$arg" == "cleanRun" ]]; then
    counter=0
    while IFS='' read -r line || [[ -n "$line" ]]; do
        isOdd=$(expr $counter % 2)
        if [[ $isOdd == 0 ]]; then
            echo "Folder: $line"
            echo "(cd $line; rm -r mat prot res traces out_run.txt ro.txt)"
            (cd $line; rm -r mat prot res traces fin_sem out_run.txt ro.txt)
        else
            echo "PBS id: $line"
            echo "qdel -W force $line"
            qdel -W force $line
        fi
        (( counter ++ ))
    done < $or_path/out_PBS_run.txt
    rm $or_path/out_PBS_run.txt
elif [[ "$arg" == "list" ]]; then
    echo "Listing files: "
elif [[ "$arg" == "clean" ]]; then
    echo "Cleaning folders"
elif [[ "$arg" == "mvmesh" ]]; then
    echo "Moving meshes to ./sem folder"
else
    echo "Argument: $arg undentified (it should be list/clean/copy/prep/mvmesh/gen/run/cleanRun)"
fi

#ACTIONS ON ALL FOLDERS THAT HAVE A run.pbs FILE
for f in $(find . -print | grep -i "run.pbs"); do
    #Do something, the file is accessible with $f:
    folder="$(dirname $f)"

if [[ "$arg" == "list" ]]; then
    echo "-------------"
    echo "ABS: " $abs_path/$folder
    echo "REL: " $root_folder/$folder
    (cd $folder; ls)
elif [[ "$arg" == "clean" ]]; then
    (cd $folder; rm -r SAMPLES out_*)
elif [[ "$arg" == "copy" ]]; then
    cp $src_dir/prepro.pbs $src_dir/run.pbs $src_dir/mesh.input $folder/
else
    (cd /home/carvalhol/Projects/ScaRL/build; make all)
    (cd $folder; rm -r *SAMPLES out_*)
    (cd $folder; ./vai.sh $arg)
fi
done
#done >output_file
