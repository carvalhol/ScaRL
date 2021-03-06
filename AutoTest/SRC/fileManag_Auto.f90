module fileManag_Auto

    use charFunctions
    use constants_Auto

    implicit none

contains

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine makeCase(nDim, Nmc, corrMod, margiFirst, corrL, fieldAvg, fieldVar, method, &
                        seedStart, independent, overlap, &
                        xMinGlob, xMaxGlob, pointsPerCorrL, &
                        nProcsTotal, nProcsPerChunk, &
                        localizationLevel, nFields, &
                        nChunks, memPerChunk, queue, wallTime, cluster, &
                        folderPath, runPath, iter)
        implicit none
        !INPUT
        integer, intent(in) :: nDim, Nmc, corrMod, margiFirst, method, seedStart, independent, iter
        double precision, dimension(:), intent(in) :: corrL, overlap
        double precision, dimension(:), intent(in) :: xMinGlob, xMaxGlob
        integer, dimension(:), intent(in) :: pointsPerCorrL
        double precision, intent(in) :: fieldAvg, fieldVar
        integer, intent(in) :: nProcsTotal, nProcsPerChunk, nChunks, memPerChunk, cluster
        character(len=*), intent(in) :: wallTime
        character(len=*), intent(in) :: queue
        character(len=tSize) :: folderPath
        integer, intent(in) :: localizationLevel
        integer, dimension(:), intent(in) :: nFields
        !OUTPUT
        character(len=*), intent(out), optional :: runPath
        !LOCAL
        character(len=tSize) :: QManagerFile_path
        character(len=tSize) :: gen_path, gen_name
        character(len=tSize) :: mesh_path, mesh_name
        character(len=tSize) :: main_input_path
        character(len=tSize) :: command_path, vai_path
        character(len=tSize) :: jobName
        integer :: i
        character :: indepChar

        write(*,*) " "
        write(*,*) "------------------------------------------------------"
        write(*,*) "Making case on: ", trim(adjustL(folderPath))
        write(*,*) "        xRange: ", xMaxGlob-xMinGlob
        write(*,*) "   nProcsTotal: ", nProcsTotal
        write(*,*) "    memPerProc: ", dble(memPerChunk)/dble(nProcsPerChunk)
        write(*,*) "       nFields: ", nFields


        gen_name = "./gen_input"
        mesh_name = "./mesh_input"
        command_path = string_join_many(folderPath,"/","run.command")
        vai_path = string_join_many(folderPath,"/","vai.sh")
        gen_path  = string_join_many(folderPath,"/",gen_name)
        mesh_path = string_join_many(folderPath,"/",mesh_name)
        main_input_path = string_join_many(folderPath,"/","input_ScaRL.txt")
        
        call write_Scarl_input_file(main_input_path, &
                                      xMinGlob,             &
                                      xMaxGlob,             &
                                      corrL,             &
                                      overlap,   &
                                      pointsPerCorrL,             &
                                      corrMod,             &
                                      margiFirst,             &
                                      fieldAvg,             &
                                      fieldVar,             &
                                      seedStart  &
                                      )

        !call write_RF_main_input_file(main_input_path, gen_name, mesh_name)

        !call write_mesh_file(nDim, xMinGlob, xMaxGlob, pointsPerCorrL, mesh_path)

        !call write_gen_file(nDim, Nmc, corrMod, margiFirst, corrL, fieldAvg, fieldVar, method, &
        !                    seedStart, independent, overlap, gen_path, &
        !                    localizationLevel, nFields)

        indepChar = "g"
        if(independent == 1) indepChar = "l"
        jobName = string_join_many("M",numb2String(method),"-",indepChar,"_",numb2String(iter,2))

        select case (cluster) 
            case(IGLOO, FUSION)
                QManagerFile_path  = string_join_many(folderPath,"/","run.pbs")
                call writePBSfile(nDim, nProcsTotal, nProcsPerChunk, nChunks, &
                              memPerChunk, wallTime, queue, QManagerFile_path, jobName, cluster, vai_path)
            case(OCCYGEN, S_DUMONT)
                QManagerFile_path  = string_join_many(folderPath,"/","run.slurm")
                call writeSlurmfile(nDim, nProcsTotal, nProcsPerChunk, nChunks, &
                              memPerChunk, wallTime, queue, QManagerFile_path, jobName, cluster)
            case(LOCAL_MAC)
                call write_command_file(nProcsTotal, command_path)
        end select

        if(present(runPath)) runPath = QManagerFile_path

    end subroutine makeCase
    
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    subroutine write_Scarl_input_file(main_input_path, &
                                      xMinGlob,        &
                                      xMaxGlob,        &
                                      corrL,           &
                                      overlap,         &
                                      pointsPerCorrL,  &
                                      corrMod,         &
                                      margiFirst,      &
                                      fieldAvg,        &
                                      fieldVar,        &
                                      seedStart        &
                                      )

        implicit none
        !INPUT
        character(len=*), intent(in) :: main_input_path
        double precision, dimension(:), intent(in) :: xMinGlob, xMaxGlob
        double precision, dimension(:), intent(in) :: corrL
        integer, dimension(:), intent(in) :: pointsPerCorrL
        integer, intent(in) :: corrMod, margiFirst
        double precision, intent(in) :: fieldAvg, fieldVar
        integer, intent(in) :: seedStart
        double precision, dimension(:), intent(in) :: overlap
        !LOCAL
        integer :: fileId
        character(len=1024) :: subfolder = "SAMPLES"

        fileID = 18

        open (unit = fileId , file = main_input_path, action = 'write')

        write(fileId,*) "#Output Folder <string, in-between quotes>"
        write(fileId,*) '"',trim(adjustL(subfolder)),'"'
        write(fileId,*) "#Number of Samples <int>"
        write(fileId,*) 1
        write(fileId,*) "#--------------------------------"
        write(fileId,*) "#Sample name <string, in-between quotes>"
        write(fileId,*) '"Sample_1"' 
        write(fileId,*) "#Coodinate Min (X,Y,Z) <double>"
        write(fileId,*) xMinGlob 
        write(fileId,*) "#Coordinate Max (X, Y, Z) <double>"
        write(fileId,*) xMaxGlob
        write(fileId,*) "#Correlation Lenghts (X,Y,Z) <double>"
        write(fileId,*) corrL
        write(fileId,*) "#Overlap (in correlation lengths)(X,Y,Z) <double>"
        write(fileId,*) overlap
        write(fileId,*) "#Points per Correlation Length (X,Y,Z) <integer>"
        write(fileId,*) pointsPerCorrL
        write(fileId,*) "#Correlation Model: 1-Gaussian <int> "
        write(fileId,*) corrMod
        write(fileId,*) "#First-order marginal density: 1-Gaussian/2-Lognormal <int>"
        write(fileId,*) margiFirst
        write(fileId,*) "#Average <double>"
        write(fileId,*) fieldAvg 
        write(fileId,*) "#Standard Deviation <double>"
        write(fileId,*) fieldVar 
        write(fileId,*) "#Random seed: <0-random (changes at every run)/>=0 deterministic"
        write(fileId,*) seedStart     
        write(fileId,*) "#--------------------------------"

        close(fileId)

        call system("chmod a+r "//trim(main_input_path))

    end subroutine write_ScaRL_input_file
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_RF_main_input_file(main_input_path, gen_path, mesh_path)

        implicit none
        !INPUT
        character(len=*), intent(in) :: main_input_path, gen_path, mesh_path
        !LOCAL
        integer :: nSamples = 1
        character(len=tsize) :: out_folder = "./results/res", out_name = "sample_1"
        integer :: fileID

        fileID = 18

        open (unit = fileId , file = main_input_path, action = 'write')

        write(fileId,"(A)") "$application 1"
        write(fileId,"(A)") "$nSamples "//numb2String(nSamples)
        write(fileId,"(A)") "$timeFolder 1"
        write(fileId,*) " "
        write(fileId,"(A)") "$calculateCorrL 0"
        write(fileId,"(A)") "$deleteSampleAfterStatistics 1"
        write(fileId,*) " "
        write(fileId,"(A)") '$mesh_input_1 '//trim(string_join_many('"',mesh_path,'"'))
        write(fileId,"(A)") '$gen_input_1  '//trim(string_join_many('"',gen_path,'"'))
        write(fileId,"(A)") '$out_folder_1 '//trim(string_join_many('"',out_folder,'"'))
        write(fileId,"(A)") '$out_name_1   '//trim(string_join_many('"',out_name,'"'))

        close(fileId)

        !call system("chmod a+r "//trim(mesh_path))

    end subroutine write_RF_main_input_file

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_mesh_file(nDim, xMinGlob, xMaxGlob, pointsPerCorrL, mesh_path)

        implicit none
        !INPUT
        integer, intent(in) :: nDim
        double precision, dimension(:), intent(in) :: xMinGlob, xMaxGlob
        integer, dimension(:), intent(in) :: pointsPerCorrL
        character(len=*), intent(in) :: mesh_path
        !LOCAL
        integer :: i
        integer :: fileId

        fileID = 18

        open (unit = fileId , file = mesh_path, action = 'write')

        write(fileId,*) "$$nDim ", nDim
        write(fileId,*) "$$meshMod 1"
        write(fileId,*) "          $Min            $Max           $pointsPerCorrL"
        do i = 1, nDim
            write(fileId, "(2(F15.5, A), (I15))") xMinGlob(i), " ", xMaxGlob(i), " ", pointsPerCorrL(i)
        end do

        close(fileId)

        call system("chmod a+r "//trim(mesh_path))

    end subroutine write_mesh_file

    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    subroutine write_gen_file(nDim, Nmc, corrMod, margiFirst, corrL, fieldAvg, fieldVar, method, &
                              seedStart, independent, overlap, gen_path,                         &
                              localizationLevel, nFields)

        implicit none
        !INPUT
        integer, intent(in) :: nDim, Nmc, corrMod, margiFirst, method, seedStart, independent
        double precision, dimension(:), intent(in) :: corrL, overlap
        double precision, intent(in) :: fieldAvg, fieldVar
        character(len=*), intent(in) :: gen_path
        integer, intent(in) :: localizationLevel
        integer, dimension(:), intent(in) :: nFields
        !LOCAL
        integer :: fileId

        fileID = 18

        open (unit = fileId , file = gen_path, action = 'write')

        write(fileId,*) "$$nDim ", nDim
        write(fileId,*) "$$Nmc ", Nmc
        write(fileId,*) "$$corrMod ", corrMod
        write(fileId,*) "$$margiFirst ", margiFirst
        write(fileId,*) "$$localizationLevel ", localizationLevel
        write(fileId,*) "$nFields "
        write(fileId,*) nFields
        write(fileId,*) "$corrL "
        write(fileId,*) corrL
        write(fileId,*) "$$fieldAvg "
        write(fileId,*) fieldAvg
        write(fileId,*) "$$fieldVar "
        write(fileId,*) fieldVar
        write(fileId,*) "$$method "
        write(fileId,*) method
        write(fileId,*) "$$seedStart"
        write(fileId,*) seedStart
        !write(fileId,*) "$$independent"
        !write(fileId,*) independent
        write(fileId,*) "$overlap"
        write(fileId,*) overlap

        close(fileId)

        call system("chmod a+r "//trim(gen_path))

    end subroutine write_gen_file

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_command_file(nProcsTotal, folderPath)

        implicit none
        !INPUT
        integer, intent(in) :: nProcsTotal
        character(len=*), intent(in) :: folderPath
        !LOCAL
        integer :: i
        integer :: fileId
        character(len=20) :: NP

        fileID = 18

        open (unit = fileId , file = folderPath, action = 'write')

        NP = string_join_many("NP=",numb2String(nProcsTotal))

        write(fileId,"(A)") "#!/bin/bash"
        write(fileId,"(A)") NP
        !write(fileId,"(A)") "(cd /Users/carvalhol/Desktop/GITs/RANDOM_FIELD/build; make all)"
        write(fileId,"(A)") 'echo ""'
        write(fileId,"(A)") 'echo "---------------------------------"'
        write(fileId,"(A)") 'echo ""'
        !write(fileId,"(A)") 'make all'
        !write(fileId,"(A)") 'cd '//folderPath
        write(fileId,"(A)") 'rm  log*'
        write(fileId,"(A)") 'rm  out_.*'
        !write(fileId,"(A)") 'rm  -r results'
        !write(fileId,"(A)") '#sleep 1'
        write(fileId,"(A)") 'mpirun --allow-run-as-root -np $NP /Users/carvalhol/Desktop/GITs/RANDOM_FIELD/build/randomField.exe'
        !write(fileId,"(A)") &
        !'mpirun --allow-run-as-root -np $NP /Users/carvalhol/Desktop/GITs/RANDOM_FIELD/build/statistics.exe<stat_input'
         write(fileId,"(A)") &
        'mpirun --allow-run-as-root -np 1 /Users/carvalhol/Desktop/GITs/RANDOM_FIELD/build/statistics.exe<stat_input'
        
        write(fileId,"(A)") ' '
        write(fileId,"(A)") 'ls'

        close(fileId)

        call system("chmod u+x "//trim(folderPath))
        call system("chmod a+r "//trim(folderPath))

    end subroutine write_command_file

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine writePBSfile(nDim, nProcsTotal, nProcsPerChunk, nChunks, memPerChunk, &
                            wallTime, queue, QManagerFile_path, jobName, cluster, vai_path)

        implicit none
        !INPUT
        integer, intent(in) :: nDim, nProcsTotal, nProcsPerChunk, nChunks, memPerChunk
        character(len=*), intent(in) :: wallTime
        character(len=*), intent(in) :: QManagerFile_path, vai_path
        character(len=*), intent(in) :: queue
        !character(len=50) :: name
        character(len=*), intent(in) :: jobName
        integer, intent(in) :: cluster

        !LOCAL
        integer :: nProcsPerChunk_chSz, nProcsTotal_chSz
        integer :: nChunks_chSz
        integer :: memPerChunk_chSz
        integer :: nDim_chSz
        integer :: fileId, fid2
        character(len=200) :: format
        character(len=50) :: outName
        integer :: i

        fileID = 28
        fid2 = 29
        outName = "output_run"
        nDim_chSz = findCharSize(nDim)
        nProcsPerChunk_chSz = findCharSize(nProcsPerChunk)
        nChunks_chSz = findCharSize(nChunks)
        memPerChunk_chSz = findCharSize(memPerChunk)

        open (unit = fileId , file = trim(adjustL(QManagerFile_path)), action = 'write')

        write(fileId,"(A)") "#!/bin/bash"
        write(fileId,"(A)") ""
        write(fileId,"(A)") "#PBS -S /bin/bash"
        write(fileId,"(A8,A50)") "#PBS -N ", jobName
        write(fileId,"(A8,A50)") "#PBS -o ", outName
        write(fileId,"(A)") "#PBS -j oe"
        write(fileId,"(A17,A8)") "#PBS -l walltime=", wallTime
        format = string_join_many("(A15,A",numb2String(nChunks_chSz),",A7,A",numb2String(nProcsPerChunk_chSz), &
                                   ", A10, A", numb2String(nProcsPerChunk_chSz),", A5, A", &
                                   numb2String(memPerChunk_chSz),", A2  )")
        write(fileId,format) "#PBS -l select=", numb2String(nChunks), ":ncpus=",numb2String(nProcsPerChunk),&
                            ":mpiprocs=",numb2String(nProcsPerChunk),":mem=", numb2String(memPerChunk), "mb"
        write(fileId,"(A)") "#PBS -q "//queue
        write(fileId,"(A)") "#PBS -M lucianopaludoecp@gmail.com"
        if(cluster == FUSION) write(fileId,"(A)") "#PBS -P omaha"
        write(fileId,"(A)") ""
        write(fileId,"(A)") 'if [ $NP ]'
        write(fileId,"(A)") 'then'
        write(fileId,"(A)") '    echo "NP = " $NP'
        write(fileId,"(A)") 'else'
        write(fileId,"(A)") '    NP='//trim(numb2String(nProcsTotal))
        write(fileId,"(A)") 'fi'
        write(fileId,"(A)") ' '
        write(fileId,"(A)") 'if [ $nRuns ]'
        write(fileId,"(A)") 'then'
        write(fileId,"(A)") '    echo "nRuns_in = " $nRuns'
        write(fileId,"(A)") 'else'
        write(fileId,"(A)") '    nRuns=1'
        write(fileId,"(A)") 'fi'
        write(fileId,"(A)") ' '
        write(fileId,"(A)") 'if [ $Run_Stat ]'
        write(fileId,"(A)") 'then'
        write(fileId,"(A)") '    echo "Run_Stat_in = " $Run_Stat'
        write(fileId,"(A)") 'else'
        write(fileId,"(A)") '    Run_Stat=0'
        write(fileId,"(A)") 'fi'
        write(fileId,"(A)") ' '
        write(fileId,"(A)") 'if [ $Run_RF ]'
        write(fileId,"(A)") 'then'
        write(fileId,"(A)") '    echo "Run_RF_in = " $Run_RF'
        write(fileId,"(A)") 'else'
        write(fileId,"(A)") '    Run_RF=1'
        write(fileId,"(A)") 'fi'
        write(fileId,"(A)") ''
        write(fileId,"(A)") 'echo "NP       = " $NP'
        write(fileId,"(A)") 'echo "Run_RF   = " $Run_RF'
        write(fileId,"(A)") 'echo "Run_Stat = " $Run_Stat'
        write(fileId,"(A)") 'echo "   nRuns = " $nRuns'

        if(cluster == IGLOO) then
            write(fileId,"(A)") "# chargement des modules"
            write(fileId,"(A)") "module purge" 
            write(fileId,"(A)") "module load intel-compiler/15.0.1"
            write(fileId,"(A)") "module load intel-mkl/11.2.1"
            write(fileId,"(A)") "module load intel-mpi/5.0.2"
            write(fileId,"(A)") "module load phdf5/1.8.15"
            write(fileId,"(A)") "module load fftw/3.3.4-intelmpi5.0.2"
        end if

        if(cluster == FUSION) then
            write(fileId,"(A)") "# chargement des modules"
            write(fileId,"(A)") "module purge"
            write(fileId,"(A)") "module load intel-compilers/16.0.3"
            write(fileId,"(A)") "module load intel-mkl/11.3.3"
            write(fileId,"(A)") "module load intel-mpi/5.1.2"
            write(fileId,"(A)") "module load phdf5/1.8.17-IntelMPI"
            write(fileId,"(A)") "module load fftw/3.3.4"
            !write(fileId,"(A)") "export LD_LIBRARY_PATH=~/LOCAL/lib:/gpfs/opt/compilers/intel/compilers_and_libraries_2016.3.210/linux/compiler/lib/intel64"
            !write(fileId,"(A)") "export LIBRARY_PATH=~/LOCAL/lib:/gpfs/opt/compilers/intel/compilers_and_libraries_2016.3.210/linux/compiler/lib/intel64"
        end if

        write(fileId,"(A)") ""
        write(fileId,"(A)") "# On se place dans le repertoire depuis lequel le job a ete soumis"
        write(fileId,"(A)") "cd $PBS_O_WORKDIR"
        write(fileId,"(A)") ""
        write(fileId,"(A)") "cat $PBS_NODEFILE | uniq > mpd.hosts"
        write(fileId,"(A)") "nb_nodes=`cat mpd.hosts|wc -l`"
        write(fileId,"(A)") ""
        write(fileId,"(A)") "startTime=`date`"
        write(fileId,"(A)") "for ((i=1;  i<=$nRuns; i++))"
        write(fileId,"(A)") "do"
        write(fileId,"(A)") 'echo "Running " $i'
        write(fileId,"(A)") 'if [ "$Run_RF" -eq "1"  ]'
        write(fileId,"(A)") 'then'
        write(fileId,"(A)") '    mpirun -np $NP '//trim(execPath)
        write(fileId,"(A)") 'fi'
        write(fileId,"(A)") ""
        write(fileId,"(A)") 'if [ "$Run_Stat" -eq "1"  ]'
        write(fileId,"(A)") 'then'
        write(fileId,"(A)") '    mpirun -np 1 '//trim(exec2Path)//"<stat_input"
        write(fileId,"(A)") 'fi'
        write(fileId,"(A)") 'done'
        write(fileId,"(A)") "endTime=`date`"
        write(fileId,"(A)") ""
        write(fileId,"(A)") 'echo "     nRuns=$nRuns"'
        write(fileId,"(A)") 'echo "Start time $startTime"'
        write(fileId,"(A)") 'echo "  End time $endTime"'

        close(fileId)

        call system("chmod a+r "//trim(QManagerFile_path))


        !WRITING VAI FILE 
        open (unit = fid2 , file = trim(adjustL(vai_path)), action = 'write')
        write(fid2,"(A)") '#!/bin/bash'
        write(fid2,"(A)") 'N_SELECT='//trim(numb2String(nChunks))
        write(fid2,"(A)") 'N_CPU='//trim(numb2String(nProcsPerChunk))
        write(fid2,"(A)") '# NP=$(($N_SELECT*$N_CPU))'
        write(fid2,"(A)") 'NP='//trim(numb2String(nProcsTotal))
        write(fid2,"(A)") 'Only_Build=0'
        write(fid2,"(A)") 'Delete_Results=0'
        write(fid2,"(A)") 'Run=0'
        write(fid2,"(A)") 'Open_Output=0'
        write(fid2,"(A)") ''
        write(fid2,"(A)") 'Build_Path="/home/carvalhol/ScaRL/build"'
        write(fid2,"(A)") 'Queue="'//trim(queue)//'"'
        write(fid2,"(A)") 'Mail="lucianopaludoecp@gmail.com"'
        write(fid2,"(A)") 'User="carvalhol"'
        write(fid2,"(A)") 'W_TIME="'//trim(wallTime)//'"'
        write(fid2,"(A)") 'PBS_Name="run.pbs"'
        write(fid2,"(A)") 'list_Name="VAI_Test"'
        write(fid2,"(A)") 'output_Name="output_run"'
        write(fid2,"(A)") 'mem="'//trim(numb2String(memPerChunk))//'mb"'
        write(fid2,"(A)") 'time_to_find="20"'
        write(fid2,"(A)") ''
        write(fid2,"(A)") 'arg=$1'
        write(fid2,"(A)") 'echo "Argument " $arg'
        write(fid2,"(A)") ''
        write(fid2,"(A)") 'if [[ -n "$arg"  ]]; then'
        write(fid2,"(A)") '   echo "Ive got an argument"'
        write(fid2,"(A)") 'else'
        write(fid2,"(A)") '   echo "No argument, ill use build"'
        write(fid2,"(A)") '   arg="build"'
        write(fid2,"(A)") 'fi'
        write(fid2,"(A)") ''
        write(fid2,"(A)") 'if [[ "$arg" == "roda" ]]; then'
        write(fid2,"(A)") '   echo "And it is roda"'
        write(fid2,"(A)") '   Run=1'
        write(fid2,"(A)") '   Only_Build=0'
        write(fid2,"(A)") '   Delete_Results=1'
        write(fid2,"(A)") '   Open_Output=1'
        write(fid2,"(A)") '   (cd $Build_Path; make all)'
        write(fid2,"(A)") 'elif [[ "$arg" == "build" ]]; then'
        write(fid2,"(A)") '   Run=1'
        write(fid2,"(A)") '   Only_Build=1'
        write(fid2,"(A)") '   Delete_Results=0'
        write(fid2,"(A)") '   Open_Output=0'
        write(fid2,"(A)") '   (cd $Build_Path; make all)'
        write(fid2,"(A)") 'elif [[ "$arg" == "all" ]]; then'
        write(fid2,"(A)") '   Run=1'
        write(fid2,"(A)") '   Only_Build=0'
        write(fid2,"(A)") '   Delete_Results=1'
        write(fid2,"(A)") '   Open_Output=1'
        write(fid2,"(A)") '   (cd $Build_Path; make all)'
        write(fid2,"(A)") 'fi'
        write(fid2,"(A)") ''
        write(fid2,"(A)") 'if [[ "$Delete_Results" -eq "1" ]]; then'
        write(fid2,"(A)") '   rm -r results logs'
        write(fid2,"(A)") '   echo "Results Deleted"'
        write(fid2,"(A)") 'fi'
        write(fid2,"(A)") ''
        write(fid2,"(A)") 'if [[ "$Only_Build" -eq "0" ]]; then'
        write(fid2,"(A)") '   rm $output_Name'
        write(fid2,"(A)") '   #qsub -v NP=$NP,Run_Stat=$Run_Stat,Run_RF=$Run_RF \'
        write(fid2,"(A)") '   #-S /bin/bash \'
        write(fid2,"(A)") '   #-N $list_Name \'
        write(fid2,"(A)") '   #-o $output_Name \'
        write(fid2,"(A)") '   #-j oe \'
        write(fid2,"(A)") '   #-l walltime=$W_TIME \'
        write(fid2,"(A)") '   #-l select=$N_SELECT:ncpus=$N_CPU:mpiprocs=$N_CPU:mem=$mem \'
        write(fid2,"(A)") '   #-q $Queue \'
        write(fid2,"(A)") '   #-M $Mail \'
        if(cluster == FUSION) write(fid2,"(A)") "#    -P omaha \"
        write(fid2,"(A)") '   #$PBS_Name'
        write(fid2,"(A)") '   qsub -v NP=$NP,Run_Stat=$Run_Stat,Run_RF=$Run_RF $PBS_Name'
        write(fid2,"(A)") '   qstat -u $User'
        write(fid2,"(A)") 'fi'
        write(fid2,"(A)") '  '
        write(fid2,"(A)") 'if [[ "$Open_Output" -eq "1" ]]; then'
        write(fid2,"(A)") '   COUNTER=$time_to_find'
        write(fid2,"(A)") '   until [ $COUNTER -lt 1 ]; do'
        write(fid2,"(A)") '       let COUNTER-=1'
        write(fid2,"(A)") '       sleep 1'
        write(fid2,"(A)") '       if [[ -f $output_Name ]]; then'
        write(fid2,"(A)") '           let COUNTER=0'
        write(fid2,"(A)") '           echo ":) FILE FOUND!!"'
        write(fid2,"(A)") '           more $output_Name'
        write(fid2,"(A)") '       fi'
        write(fid2,"(A)") '   done'
        write(fid2,"(A)") ' '
        write(fid2,"(A)") '   if [[ "$COUNTER" -lt 0 ]]; then'
        write(fid2,"(A)") '       echo ":( IVE BEEN SEARCHING FOR"$time_to_find "s AND THE FILE ISNT HERE YET"'
        write(fid2,"(A)") '   fi'
        write(fid2,"(A)") 'fi'

        close(fid2)

        call system("chmod a+r "//trim(vai_path))
        call system("chmod u+x "//trim(vai_path))

    end subroutine writePBSfile


    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine writeSlurmfile(nDim, nProcsTotal, nProcsPerChunk, nChunks, memPerChunk, wallTime, &
                              queue, QManagerFile_path, jobName, cluster)

        implicit none
        !INPUT
        integer, intent(in) :: nDim, nProcsTotal, nProcsPerChunk, nChunks, memPerChunk
        character(len=8), intent(in) :: wallTime
        character(len=200) :: QManagerFile_path
        character(len=*), intent(in) :: queue
        !character(len=50) :: name
        character(len=50), intent(in) :: jobName
        integer, intent(in) :: cluster

        !LOCAL
        integer :: nProcsPerChunk_chSz, nProcsTotal_chSz
        integer :: nChunks_chSz
        integer :: memPerChunk_chSz
        integer :: nDim_chSz
        integer :: fileId
        integer :: memTot, memTot_chSz
        character(len=200) :: format
        character(len=50) :: outName
        integer :: i

        fileID = 28
        outName = "output_RF"
        !memTot = memPerChunk*nChunks/1000
        memTot = 64
        nDim_chSz = findCharSize(nDim)
        nProcsPerChunk_chSz = findCharSize(nProcsPerChunk)
        nChunks_chSz = findCharSize(nChunks)
        memPerChunk_chSz = findCharSize(memPerChunk)
        nProcsTotal_chSz = findCharSize(nProcsTotal)
        memTot_chSz = findCharSize(memTot)

        open (unit = fileId , file = QManagerFile_path, action = 'write')

    !#!/bin/bash
    !
    !#SBATCH -J Mesh_SEM
    !#SBATCH --nodes=1
    !#SBATCH --ntasks=1
    !#SBATCH --ntasks-per-node=1
    !#SBATCH --threads-per-core=1
    !#SBATCH --time=00:01:00
    !#SBATCH --output output.txt
    !#SBATCH --mail-type=ALL
    !#SBATCH --mail-user=lucio.a.c@gmail.com
    !
    !module purge
    !module load intel/15.0.0.090
    !module load bullxmpi/1.2.8.4
    !module load hdf5/1.8.14
    !srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS /panfs/panasas/cnt0025/mss7417/abreul/SEM/build/MESH/mesher<mesh.input
            !,",A7,A",numb2String(nProcsPerChunk_chSz),", A10, A", numb2String(nProcsPerChunk_chSz),", A5, A", numb2String(memPerChunk_chSz),", A2  )"])


!#!/bin/bash
!
!#SBATCH -J M4-l_01
!#SBATCH --nodes=1
!#SBATCH --ntasks=24
!#SBATCH --ntasks-per-node=24
!#SBATCH --threads-per-core=1
!#SBATCH --time=04:00:00
!#SBATCH --mem=64GB
!#SBATCH --output out_RF
!#SBATCH --mail-type=ALL
!#SBATCH --mail-user=lucio.a.c@gmail.com
!
!module load intel-compiler/15.0.3.187
!module load intelmpi/5.0.3.048
!module load hdf5/1.8.14
!module load fftw3/3.3.4
!srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS /home/abreul/RF/Novo/build/randomField.exe
!##srun --mpi=pmi2 -K1 --resv-ports -n 1 /home/abreul/RF/Novo/build/statistics.exe<stat_input


        write(fileId,"(A)") "#!/bin/bash"
        write(fileId,"(A)") ""
        write(fileId,"(A11,A50)") "#SBATCH -J ", jobName
        format = string_join_many("(A16,A",numb2String(nChunks_chSz),")")
        write(fileId,format) "#SBATCH --nodes=", numb2String(nChunks)
        format = string_join_many("(A17,A",numb2String(nProcsTotal_chSz),")")
        write(fileId,format) "#SBATCH --ntasks=", trim(numb2String(nProcsTotal))
        format = string_join_many("(A26,A",numb2String(nProcsPerChunk_chSz),")")
        write(fileId,format) "#SBATCH --ntasks-per-node=", numb2String(nProcsPerChunk)
        format = string_join_many("(A14,A",numb2String(memTot_chSz),",A2)")
        write(fileId,format) "#SBATCH --mem=", numb2String(memTot),"GB"
        write(fileId,"(A)") "#SBATCH --threads-per-core=1"
        write(fileId,"(A15,A8)") "#SBATCH --time=", wallTime
        write(fileId,"(A17,A)") "#SBATCH --output ", trim(outName)
        !write(fileId,"(A)") "#SBATCH --mail-type=ALL"
        !write(fileId,"(A)") "#SBATCH --mail-user=lucio.a.c@gmail.com"
        write(fileId,"(A)") ""
        !format = string_join_many("(A15,A",numb2String(nChunks_chSz),",A7,A",numb2String(nProcsPerChunk_chSz),", A10, A", numb2String(nProcsPerChunk_chSz),", A5, A", numb2String(memPerChunk_chSz),", A2  )")

        write(fileId,"(A)") "module load intel-compiler/15.0.3.187"
        write(fileId,"(A)") "module load intelmpi/5.0.3.048"
        write(fileId,"(A)") "module load hdf5/1.8.14"
        write(fileId,"(A)") "module load fftw3/3.3.4"
        write(fileId,"(A)") "srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS "//trim(execPath)
        write(fileId,"(A)") "#srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS "//trim(exec2Path)//"<stat_input"
        !write(fileId,"(A)") "mpirun --rsh=ssh -n $nb_nodes -f mpd.hosts -np "//trim(numb2String(nProcsTotal))//" "//trim(execPath)

        close(fileId)

        call system("chmod a+r "//trim(QManagerFile_path))

    end subroutine writeSlurmfile

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function findCharSize(number) result(nSize)
        !INPUT
        integer, intent(in) :: number
        !LOCAL
        integer :: nSize
        integer :: comp

        comp = 9
        nSize = 1

        do while(comp < abs(number))
            comp = comp+9*(10**nSize)
            nSize = nSize + 1
        end do

        if(number<0) nSize = nSize + 1

    end function findCharSize

end module fileManag_Auto
