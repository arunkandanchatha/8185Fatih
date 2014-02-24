MODULE nrtype
    INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
    INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
    INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
    INTEGER, PARAMETER :: SP = KIND(1.0)
    INTEGER, PARAMETER :: DP = KIND(1.0D0)
    INTEGER, PARAMETER :: SPC = KIND((1.0_sp,1.0_sp))
    INTEGER, PARAMETER :: DPC = KIND((1.0_dp,1.0_dp))
    INTEGER, PARAMETER :: LGT = KIND(.true.)
    REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
    REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
    REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
    REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
    REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
    REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
    REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
    REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
    TYPE sprs2_sp
        INTEGER(I4B) :: n,len
        REAL(SP), DIMENSION(:), POINTER :: val
        INTEGER(I4B), DIMENSION(:), POINTER :: irow
        INTEGER(I4B), DIMENSION(:), POINTER :: jcol
    END TYPE sprs2_sp
    TYPE sprs2_dp
        INTEGER(I4B) :: n,len
        REAL(DP), DIMENSION(:), POINTER :: val
        INTEGER(I4B), DIMENSION(:), POINTER :: irow
        INTEGER(I4B), DIMENSION(:), POINTER :: jcol
    END TYPE sprs2_dp
    ABSTRACT INTERFACE
        FUNCTION template_function(point) RESULT(z)
            !
            ! the function we are trying to minimize
            ! INPUTS: point - the value at which we are evaluating the function
            ! OUTPUTS: z - the value of the function at that point
            REAL(KIND=8), INTENT(IN) :: point
            REAL(KIND=8) :: z
        END FUNCTION template_function
    END INTERFACE
    ABSTRACT INTERFACE
        PURE FUNCTION template_function2(point1, point2, point3) RESULT(z)
            !
            ! the function we are trying to minimize
            ! INPUTS: point1, point2 - the value at which we are evaluating the function
            ! OUTPUTS: z - the value of the function at that point
            REAL(KIND=8), INTENT(IN) :: point1, point2, point3
            REAL(KIND=8) :: z
        END FUNCTION template_function2
    END INTERFACE
    ABSTRACT INTERFACE
        FUNCTION template_function3(interest,state,employment) RESULT(z)
            REAL(KIND=8), INTENT(IN) :: interest, state, employment
            REAL(KIND=8) :: z
        END FUNCTION template_function3
    END INTERFACE
    type household
        integer :: employmentState
        real(DP) :: capital
    end type household

END MODULE nrtype
