module constants_ScaRL

    implicit none

    double precision, parameter :: PI = 3.1415926535898d0;
    double precision, parameter :: SQRT_2PI = 2.50662827463d0;
    double precision, parameter :: ZERO = 0d0;
    double precision, parameter :: TOLERANCE = 0.000000001
    double precision, parameter :: MIN_DOUBLE = -2.0D+307
    double precision, parameter :: MAX_DOUBLE = 2.0D+307
    !METHOD (gen_meth)
    integer, parameter :: FFT = 1, &
                          FFT_MPI = 2
    !Correlation Model
    integer, parameter :: cm_GAUSSIAN = 1, &
                          cm_VON_KARMAN = 2
    !First-order Marginal Density
    integer, parameter :: fom_GAUSSIAN = 1, &
                          fom_LOGNORMAL = 2
    !Mesh Mode
    integer, parameter :: msh_AUTO = 1, msh_UNV = 2

    integer, parameter :: SCREEN=6

    !APPLICATION
    integer, parameter :: NATIVE = 1
    integer, parameter :: SEM = 2

end module constants_ScaRL
