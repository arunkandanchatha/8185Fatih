MODULE nr
    use nrtype
    implicit none
contains
    SUBROUTINE tridag(a,b,c,r,u)
        USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
        IMPLICIT NONE
        REAL(dp), DIMENSION(:), INTENT(IN) :: a,b,c,r
        REAL(dp), DIMENSION(:), INTENT(OUT) :: u
        REAL(dp), DIMENSION(size(b)) :: gam
        INTEGER(I4B) :: n,j
        REAL(dp) :: bet
        n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_ser')
        bet=b(1)
        if (bet == 0.0) call nrerror('tridag_ser: Error at code stage 1')
        u(1)=r(1)/bet
        do j=2,n
            gam(j)=c(j-1)/bet
            bet=b(j)-a(j-1)*gam(j)
            if (bet == 0.0) &
                call nrerror('tridag_ser: Error at code stage 2')
            u(j)=(r(j)-a(j-1)*u(j-1))/bet
        end do
        do j=n-1,1,-1
            u(j)=u(j)-gam(j+1)*u(j+1)
        end do
    END SUBROUTINE tridag

    FUNCTION splint(xa,ya,y2a,x)
        ! Use: call spline once, then use splint to get interpolated values
        !   xa: tabulated function (could use values passed to spline, I think)
        !   ya: function values at each xa
        !   y2a: output from spline
        !   x: value on which we want to interpolate
        !   splint: retuns the interpolated value
        USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
        IMPLICIT NONE
        REAL(dp), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
        REAL(dp), INTENT(IN) :: x
        REAL(dp) :: splint
        INTEGER(I4B) :: khi,klo,n
        INTEGER(I4B),DIMENSION(1)::myLoc
        REAL(dp) :: a,b,h
        n=assert_eq(size(xa),size(ya),size(y2a),'splint')
        myLoc=minloc(abs(xa-x))
        if((xa(myLoc(1))-x)>0.0D0) then
            myLoc=myLoc-1
        end if
        klo=max(min(myLoc(1),n-1),1)
        khi=klo+1
        h=xa(khi)-xa(klo)
        if (h == 0.0) call nrerror('bad xa input in splint')
        a=(xa(khi)-x)/h
        b=(x-xa(klo))/h
        splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_dp
    END FUNCTION splint

    SUBROUTINE spline(x,y,yp1,ypn,y2)
        ! Use: call spline once, then use splint to get interpolated values
        ! x: grid of points at which function is evaluated
        ! y: values of function at each grid point, x
        ! yp1: First derivative of function at x(1)
        ! ypn: First derivative of function at x(n)
        ! y2: second derivative of interpolating function at
        !     each point in x
        USE nrtype; USE nrutil, ONLY : assert_eq
        IMPLICIT NONE
        REAL(dp), DIMENSION(:), INTENT(IN) :: x,y
        REAL(dp), INTENT(IN) :: yp1,ypn
        REAL(dp), DIMENSION(:), INTENT(OUT) :: y2
        INTEGER(I4B) :: n
        REAL(dp), DIMENSION(size(x)) :: a,b,c,r
        n=assert_eq(size(x),size(y),size(y2),'spline')
        c(1:n-1)=x(2:n)-x(1:n-1)
        r(1:n-1)=6.0_dp*((y(2:n)-y(1:n-1))/c(1:n-1))
        r(2:n-1)=r(2:n-1)-r(1:n-2)
        a(2:n-1)=c(1:n-2)
        b(2:n-1)=2.0_dp*(c(2:n-1)+a(2:n-1))
        b(1)=1.0
        b(n)=1.0
        if (yp1 > 0.99e30_dp) then
            r(1)=0.0
            c(1)=0.0
        else
            r(1)=(3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
            c(1)=0.5
        end if
        if (ypn > 0.99e30_dp) then
            r(n)=0.0
            a(n)=0.0
        else
            r(n)=(-3.0_dp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
            a(n)=0.5
        end if
        call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
    END SUBROUTINE spline

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
    SUBROUTINE q2a(func1,func2,func3,Dfunc1,Dfunc2,Dfunc3)
        !-------------------------------------
        PROCEDURE(template_function), POINTER, INTENT(in) :: func1
        PROCEDURE(template_function), POINTER, INTENT(in) :: func2
        PROCEDURE(template_function), POINTER, INTENT(in) :: func3
        PROCEDURE(template_derivative), POINTER, INTENT(in) :: Dfunc1
        PROCEDURE(template_derivative), POINTER, INTENT(in) :: Dfunc2
        PROCEDURE(template_derivative), POINTER, INTENT(in) :: Dfunc3

        INTEGER, PARAMETER :: numSteps=10
        INTEGER, PARAMETER :: numSteps2=20
        REAL(DP), DIMENSION(numSteps+1) :: gridPoints, funcAtGrid, splinePoints
        REAL(DP), DIMENSION(numSteps2+1) :: gridPoints2, funcAtGrid2, splinePoints2
        REAL(DP) :: result1, result11, result2, pointToEval
        REAL(DP), DIMENSION(1) :: pointToEvalTemp
        REAL(DP), DIMENSION(1) :: dFx1, dFxn, dFx12, dFxn2
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

        !note: spline requires as many calls to the utility function, as
        !      the number of grid points. This is the same as linear, if
        !      we are making a lot of interpolations. More if we are making
        !      just one
        !      It also requires two calls for the function derivative
        print *,"Cubic Spline interpolation"
        do i=1,numSteps+1
            pointToEvalTemp(1)=gridPoints(i)
            funcAtGrid(i)=func1(pointToEvalTemp)
            if(i == 1)then
                dFx1 = Dfunc1(pointToEvalTemp)
            else if (i==(numSteps+1)) then
                dFxn = Dfunc1(pointToEvalTemp)
            end if
        end do
        do i=1,numSteps2+1
            pointToEvalTemp(1)=gridPoints2(i)
            funcAtGrid2(i)=func1(pointToEvalTemp)
            if(i == 1)then
                dFx12 = Dfunc1(pointToEvalTemp)
            else if (i==(numSteps2+1)) then
                dFxn2 = Dfunc1(pointToEvalTemp)
            end if
        end do
        call spline(gridPoints,funcAtGrid,dFx1(1),dFxn(1),splinePoints)
        call spline(gridPoints2,funcAtGrid2,dFx12(1),dFxn2(1),splinePoints2)
        print *,"Point                  ",numSteps,"          ",numSteps2
        do i=1,floor(2.0D0/.05D0)
            pointToEval = 0.05D0*i
            result1=splint(gridPoints,funcAtGrid,splinePoints,pointToEval)
            result11=splint(gridPoints2,funcAtGrid2,splinePoints2,pointToEval)
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
    PROCEDURE(template_derivative), POINTER :: Dfunc1, Dfunc2, Dfunc3

    func1 => u1
    func2 => u2
    func3 => u3

    Dfunc1 => Du1
    Dfunc2 => Du2
    Dfunc3 => Du3

    CALL q2a(func1, func2, func3, Dfunc1, Dfunc2, Dfunc3)

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

    FUNCTION Du1(point) RESULT(z)
        !
        ! the function we are trying to minimize
        ! INPUTS: point - the value at which we are evaluating the function
        ! OUTPUTS: z - the value of the function at that point
        REAL(DP), DIMENSION(:), INTENT(IN) :: point
        REAL(DP), DIMENSION(size(point)) :: z
        REAL(DP) :: x
        x=point(1)
        z = 1.0/x
    END FUNCTION Du1

    !-----------------
    FUNCTION Du2(point) RESULT(z)
        !
        ! the function we are trying to minimize
        ! INPUTS: point - the value at which we are evaluating the function
        ! OUTPUTS: z - the value of the function at that point
        REAL(DP), DIMENSION(:), INTENT(IN) :: point
        REAL(DP), DIMENSION(size(point)) :: z
        REAL(DP) :: x

        x=point(1)
        z = 0.5D0*x**(-0.5D0)

    END FUNCTION Du2

    !-----------------
    FUNCTION Du3(point) RESULT(z)
        !
        ! the function we are trying to minimize
        ! INPUTS: point - the value at which we are evaluating the function
        ! OUTPUTS: z - the value of the function at that point
        REAL(DP), DIMENSION(:), INTENT(IN) :: point
        REAL(DP), DIMENSION(size(point)) :: z
        REAL(DP) :: x,alpha

        x=point(1)
        alpha=point(2)
        z = x**(-alpha)

    END FUNCTION Du3

end program q2
