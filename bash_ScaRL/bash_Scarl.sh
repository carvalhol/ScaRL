#!/bin/bash

#Made by LUCIANO DE CARVALHO PALUDO (lucianopaludoecp@gmail.com)
src_folder="`(pwd)`/src_folder" #absolute path of the files to copy
dest_folder="./ScaRL_WEAK"
run_system=2 #1=PBS, 2=SLURM


init_sz=31.75
iter_sz=26.75

sizes=()
sizes+=($(echo "$init_sz+0.0*$iter_sz" | bc -l))
sizes+=($(echo "$init_sz+1.0*$iter_sz" | bc -l))
sizes+=($(echo "$init_sz+2.0*$iter_sz" | bc -l))
sizes+=($(echo "$init_sz+3.0*$iter_sz" | bc -l))
sizes+=($(echo "$init_sz+4.0*$iter_sz" | bc -l))
sizes+=($(echo "$init_sz+5.0*$iter_sz" | bc -l))
sizes+=($(echo "$init_sz+6.0*$iter_sz" | bc -l))
sizes+=($(echo "$init_sz+7.0*$iter_sz" | bc -l))
sizes+=($(echo "$init_sz+8.0*$iter_sz" | bc -l))
sizes+=($(echo "$init_sz+9.0*$iter_sz" | bc -l))
sizes+=($(echo "$init_sz+10.0*$iter_sz" | bc -l))
sizes+=($(echo "$init_sz+11.0*$iter_sz" | bc -l))
sizes+=($(echo "$init_sz+12.0*$iter_sz" | bc -l))
sizes+=($(echo "$init_sz+13.0*$iter_sz" | bc -l))

procs=()
procs+=($((1**3)))
procs+=($((2**3)))
procs+=($((3**3)))
procs+=($((4**3)))
procs+=($((5**3)))
procs+=($((6**3)))
procs+=($((7**3)))
procs+=($((8**3)))
procs+=($((9**3)))
procs+=($((10**3)))
procs+=($((11**3)))
procs+=($((12**3)))
procs+=($((13**3)))
procs+=($((14**3)))


printf "\nProcs: \n"
printf '%d\n' "${procs[@]}"
printf "\nSize : \n"
printf '%f\n' "${sizes[@]}"

sizes_sz=${#sizes[@]}
procs_sz=${#procs[@]}

if [[ "$sizes_sz" -ne "$procs_sz" ]]; then
    printf "ERROR sizes and procs don't have the same size\n"
    printf "sizes_sz=$sizes_sz; procs_sz=$procs_sz\n"
    exit 1
fi

mkdir $dest_folder
cd $dest_folder
vai=`(pwd)`"/vai.sh"
touch $vai
echo "#!/bin/bash" >> $vai

ii=0
while [ $ii -lt $sizes_sz ];
do
    printf "\n\n"
    echo "Procs: " ${procs[$ii]}
    echo "Size : " ${sizes[$ii]}
    folder=`printf "%05g_procs" ${procs[$ii]}`

    ret_folder=`(pwd)`
    mkdir $folder
    cd $folder
    echo "Folder: " `(pwd)`

    #####COPY FILES AND REPLACE VALUES
    cp $src_folder/input_ScaRL.txt ./ 
    sed -i -e "s/#MAX_X/${sizes[$ii]}/" ./input_ScaRL.txt
    sed -i -e "s/#MAX_Y/${sizes[$ii]}/" ./input_ScaRL.txt
    sed -i -e "s/#MAX_Z/${sizes[$ii]}/" ./input_ScaRL.txt

    if [ "$run_system" -eq "1" ]; then
        nselect=$((1+(${procs[$ii]}-1)/24))
         
        cp $src_folder/run.pbs ./ 
        sed -i -e "s/#NSELECT/$nselect/" ./run.pbs
        sed -i -e "s/#NCPU/24/" ./run.pbs
        sed -i -e "s/#NMPI/24/" ./run.pbs
        sed -i -e "s/#NPROC/${procs[$ii]}/" ./run.pbs
        echo "(cd `(pwd)`; qsub run.pbs)" >> $vai
    elif [ "$run_system" -eq "2" ]; then
        cp $src_folder/run.slurm ./ 
        sed -i -e "s/#NNODES/${procs[$ii]}/" ./run.slurm
        sed -i -e "s/#NTASKSNODE/1/" ./run.slurm
        sed -i -e "s/#NMPI/${procs[$ii]}/" ./run.slurm
        echo "(cd `(pwd)`; sbatch run.slurm)" >> $vai
    fi
    ###### RETURNING
    cd $ret_folder
    (( ii ++ ))
done

chmod u+x $vai
chmod a+r $vai

exit 1
