MODULE nr
    use nrtype
    implicit none
contains
    FUNCTION linear(func,evalPoint,gridPoints) RESULT(z)
        !only makes sense in two dimensions
        PROCEDURE(template_function), POINTER, INTENT(in) :: func
        REAL(DP),DIMENSION(:),INTENT(in) :: gridPoints
        REAL(DP), INTENT(in) :: evalPoint
        REAL(DP) :: z

        !local variables
        REAL(DP), DIMENSION(1) :: baseGrid, nextGrid, fBase, fNext
        REAL(DP) :: slope
        integer, DIMENSION(1) :: indexInGridTemp
        integer :: indexInGrid

        indexInGridTemp=minloc(abs(gridPoints-evalPoint))
        indexInGrid = indexInGridTemp(1)

        baseGrid(1)=gridPoints(indexInGrid)
        nextGrid(1)=gridPoints(indexInGrid+floor(sign(1.0D0,evalPoint-baseGrid(1))))

        fBase = func(baseGrid)
        fNext = func(nextGrid)
        slope = (fNext(1)-fBase(1))/(nextGrid(1)-baseGrid(1))

        z=fBase(1)+slope*(evalPoint-baseGrid(1))
    END FUNCTION linear
end module nr

MODULE  hmwk2
    use nrtype
    use nr
    implicit none
    INTEGER(I8B), parameter :: NUM_ITER=100000000

contains

    !-------------------------------------
    subroutine sub_mystop(calling)
        !-------------------------------------

        ! a personal stop subroutine. Makes it easier to edit behaviour of stop. All
        ! functions and subroutines should call this.
        !
        ! INPUTS: calling - a string indicating where this subroutine was called from

        CHARACTER (LEN=*), intent(in) :: calling
        print *, "STOP: ", calling
        STOP 0
    end subroutine sub_mystop


    !-------------------------------------
    SUBROUTINE q2a(func1,func2,func3)
        !-------------------------------------
        PROCEDURE(template_function), POINTER, INTENT(in) :: func1
        PROCEDURE(template_function), POINTER, INTENT(in) :: func2
        PROCEDURE(template_function), POINTER, INTENT(in) :: func3
        INTEGER, PARAMETER :: numSteps=10
        INTEGER, PARAMETER :: numSteps2=100
        REAL(DP), DIMENSION(numSteps+1) :: gridPoints
        REAL(DP), DIMENSION(numSteps2+1) :: gridPoints2
        REAL(DP) :: result1, result11, result2, pointToEval
        REAL(DP), DIMENSION(1) :: pointToEvalTemp
        INTEGER :: i

        do i=1,numSteps+1
            gridPoints(i)=.05D0+(2.0D0-0.05)/numSteps*(i-1)
        end do
        do i=1,numSteps2+1
            gridPoints2(i)=.05D0+(2.0D0-0.05)/numSteps2*(i-1)
        end do

        !minimizing evaluations of function.

        !note: linear only requires two calls to the utility function, regardless
        !      of the number of grid points. So the finer, the better. (Assuming
        !      we only want to make on comparison). If we want to make many
        !      interpolations, then it is linear in the number, as we could simply
        !      evaluate on all grid points.
        print *,"Linear interpolation"
        print *,"Point                  ",numSteps,"          ",numSteps2
        do i=1,floor(2.0D0/.05D0)
            pointToEval = 0.05D0*i
            result1=linear(func1,pointToEval,gridPoints)
            result11=linear(func1,pointToEval,gridPoints2)
            pointToEvalTemp(1)=pointToEval
            result2=func1(pointToEvalTemp)
            print *, pointToEval,abs((result2-result1)/result2),abs((result2-result11)/result2)
        end do
    END SUBROUTINE q2a
END MODULE hmwk2

program q2
    use nrtype
    use hmwk2
    implicit none

    PROCEDURE(template_function), POINTER :: func1, func2, func3

    func1 => u1
    func2 => u2
    func3 => u3

    CALL q2a(func1, func2, func3)

CONTAINS

    !-----------------
    FUNCTION u1(point) RESULT(z)
        !
        ! the function we are trying to minimize
        ! INPUTS: point - the value at which we are evaluating the function
        ! OUTPUTS: z - the value of the function at that point
        REAL(DP), DIMENSION(:), INTENT(IN) :: point
        REAL(DP) :: x,z
        x=point(1)
        z = log(x)
    END FUNCTION u1

    !-----------------
    FUNCTION u2(point) RESULT(z)
        !
        ! the function we are trying to minimize
        ! INPUTS: point - the value at which we are evaluating the function
        ! OUTPUTS: z - the value of the function at that point
        REAL(DP), DIMENSION(:), INTENT(IN) :: point
        REAL(DP) :: x,z

        x=point(1)
        z = sqrt(x)

    END FUNCTION u2

    !-----------------
    FUNCTION u3(point) RESULT(z)
        !
        ! the function we are trying to minimize
        ! INPUTS: point - the value at which we are evaluating the function
        ! OUTPUTS: z - the value of the function at that point
        REAL(DP), DIMENSION(:), INTENT(IN) :: point
        REAL(DP) :: z
        REAL(DP) :: x,alpha

        x=point(1)
        alpha=point(2)
        z = x**(1-alpha)/(1-alpha)

    END FUNCTION u3

end program q2
