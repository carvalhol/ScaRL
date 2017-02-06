program main_ScaRL

    use mpi

    implicit none

    !INPUTS
    integer :: comm, rang, nb_procs
    !LOCAL VARIABLES

    !INPUT VARIABLES

    !Initializing MPI
    call init_communication(MPI_COMM_WORLD, comm, rang, nb_procs)

    print*, "I'm processor number ", rang, " of ", nb_procs

    !Finalizing MPI
    call end_communication()





        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
    contains

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine init_communication(comm_local, comm, rang, nb_procs)
            implicit none
            !INPUT
            integer, intent(in) :: comm_local
            !OUTPUT
            integer, intent(out) :: comm, rang, nb_procs
            !LOCAL
            integer :: code

            !call MPI_Init_thread(MPI_THREAD_MULTIPLE, &provided
            call MPI_INIT(code)
            call MPI_COMM_RANK(comm_local, rang, code)
            call MPI_COMM_SIZE(comm_local, nb_procs, code)

            comm = comm_local

        end subroutine init_communication

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine end_communication()
            implicit none
            !LOCAL
            integer :: code

            !call finalize_log_file()
            call MPI_FINALIZE(code)

        end subroutine end_communication

end program main_ScaRL
