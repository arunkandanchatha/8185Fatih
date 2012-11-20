MODULE nr
    use nrtype
    implicit none
contains
    FUNCTION chebft(a,b,n,func)
        ! Called prior to using chebev.
        ! INPUTS: a,b - bounds of the interval
        !         n - the maximum degree for the estimate
        USE nrtype; USE nrutil, ONLY : arth,outerprod
        IMPLICIT NONE
        REAL(dp), INTENT(IN) :: a,b
        INTEGER(I4B), INTENT(IN) :: n
        REAL(dp), DIMENSION(n) :: chebft
        PROCEDURE(template_function), POINTER, INTENT(in) :: func
        REAL(DP) :: bma,bpa
        REAL(DP), DIMENSION(n) :: theta,funcEvals
        REAL(DP),DIMENSION(1) :: pointToEval
        INTEGER(I4B) :: i

        bma=0.5_dp*(b-a)
        bpa=0.5_dp*(b+a)
        theta(:)=PI_D*arth(0.5_dp,1.0_dp,n)/n
        do i=1,n
            pointToEval(1)=real(cos(theta(i))*bma+bpa,dp)
            funcEvals(i)=func(pointToEval)
        end do
        chebft(:)=matmul(cos(outerprod(arth(0.0_dp,1.0_dp,n),theta)), &
            &funcEvals)*2.0_dp/n
    END FUNCTION chebft

    FUNCTION chebev(a,b,c,x)
        ! Called after using chebft.
        ! INPUTS: a,b - bounds of the interval (must be same as used
        !               in chebft call)
        !         c - the first m (<n) coefficients returned from
        !             chebft
        !         x - the point where we want to estimate the value
    USE nrtype; USE nrutil, ONLY : nrerror
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: a,b,x
    REAL(dp), DIMENSION(:), INTENT(IN) :: c
    REAL(dp) :: chebev
    INTEGER(I4B) :: j,m
    REAL(dp) :: d,dd,sv,y,y2
    if ((x-a)*(x-b) > 0.0) call nrerror('x not in range in chebev')
    m=size(c)
    d=0.0
    dd=0.0
    y=(2.0_dp*x-a-b)/(b-a)
    y2=2.0_dp*y
    do j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
    end do
    chebev=y*d-dd+0.5_dp*c(1)
    END FUNCTION chebev

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

        SUBROUTINE amoeba(p,y,ftol,func,iter)
        USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,iminloc,nrerror,swap
        IMPLICIT NONE
        INTEGER(I4B), INTENT(OUT) :: iter
        REAL(dp), INTENT(IN) :: ftol
        REAL(dp), DIMENSION(:), INTENT(INOUT) :: y   ! "func" evaluated at the n vertices provided in "p"
        REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: p ! vertices. If we have n vertices, then we must be
                                                         ! in n-1 dimensional space (we need one extra vertex
                                                         ! than dimensions. For each row, the n-1 vector
                                                         ! specifies the vertex
        PROCEDURE(template_function), POINTER, INTENT(in) :: func
        INTEGER(I4B), PARAMETER :: ITMAX=5000
        REAL(dp), PARAMETER :: TINY=1.0e-10
        INTEGER(I4B) :: ihi,ndim
        REAL(dp), DIMENSION(size(p,2)) :: psum
        call amoeba_private
    CONTAINS
        !BL
        SUBROUTINE amoeba_private
            IMPLICIT NONE
            INTEGER(I4B) :: i,ilo,inhi
            REAL(dp) :: rtol,ysave,ytry,ytmp

            ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,'amoeba')
            iter=0
            psum(:)=sum(p(:,:),dim=1)
            do
                ilo=iminloc(y(:))
                ihi=imaxloc(y(:))
                ytmp=y(ihi)
                y(ihi)=y(ilo)
                inhi=imaxloc(y(:))
                y(ihi)=ytmp
                rtol=2.0_dp*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
                if (rtol < ftol) then
                    call swap(y(1),y(ilo))
                    call swap(p(1,:),p(ilo,:))
                    RETURN
                end if
                if (iter >= ITMAX) call nrerror('ITMAX exceeded in amoeba')
                ytry=amotry(-1.0_dp)
                iter=iter+1
                if (ytry <= y(ilo)) then
                    ytry=amotry(2.0_dp)
                    iter=iter+1
                else if (ytry >= y(inhi)) then
                    ysave=y(ihi)
                    ytry=amotry(0.5_dp)
                    iter=iter+1
                    if (ytry >= ysave) then
                        p(:,:)=0.5_dp*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
                        do i=1,ndim+1
                            if (i /= ilo) y(i)=func(p(i,:))
                        end do
                        iter=iter+ndim
                        psum(:)=sum(p(:,:),dim=1)
                    end if
                end if
            end do
        END SUBROUTINE amoeba_private
        !BL
        FUNCTION amotry(fac)
            IMPLICIT NONE
            REAL(dp), INTENT(IN) :: fac
            REAL(dp) :: amotry
            REAL(dp) :: fac1,fac2,ytry
            REAL(dp), DIMENSION(size(p,2)) :: ptry
            fac1=(1.0_dp-fac)/ndim
            fac2=fac1-fac
            ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
            ytry=func(ptry)
            if (ytry < y(ihi)) then
                y(ihi)=ytry
                psum(:)=psum(:)-p(ihi,:)+ptry(:)
                p(ihi,:)=ptry(:)
            end if
            amotry=ytry
        END FUNCTION amotry
    END SUBROUTINE amoeba

        SUBROUTINE sobseq(x,init)
        USE nrtype; USE nrutil, ONLY : nrerror
        IMPLICIT NONE
        REAL(dp), DIMENSION(2), INTENT(OUT) :: x
        INTEGER(I4B), OPTIONAL, INTENT(IN) :: init
        INTEGER(I4B), PARAMETER :: MAXBIT=30,MAXDIM=6
        REAL(dp), SAVE :: fac
        INTEGER(I4B) :: i,im,ipp,j,k,l
        INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE:: iu
        INTEGER(I4B), SAVE :: in
        INTEGER(I4B), DIMENSION(MAXDIM), SAVE :: ip,ix,mdeg
        INTEGER(I4B), DIMENSION(MAXDIM*MAXBIT), SAVE :: iv
        DATA ip /0,1,1,2,1,4/, mdeg /1,2,3,3,4,4/, ix /6*0/
        DATA iv /6*1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9,156*0/
        if (present(init)) then
            ix=0
            in=0
            if (iv(1) /= 1) RETURN
            fac=1.0_dp/2.0_dp**MAXBIT
            allocate(iu(MAXDIM,MAXBIT))
            iu=reshape(iv,shape(iu))
            do k=1,MAXDIM
                do j=1,mdeg(k)
                    iu(k,j)=iu(k,j)*2**(MAXBIT-j)
                end do
                do j=mdeg(k)+1,MAXBIT
                    ipp=ip(k)
                    i=iu(k,j-mdeg(k))
                    i=ieor(i,i/2**mdeg(k))
                    do l=mdeg(k)-1,1,-1
                        if (btest(ipp,0)) i=ieor(i,iu(k,j-l))
                        ipp=ipp/2
                    end do
                    iu(k,j)=i
                end do
            end do
            iv=reshape(iu,shape(iv))
            deallocate(iu)
        else
            im=in
            do j=1,MAXBIT
                if (.not. btest(im,0)) exit
                im=im/2
            end do
            if (j > MAXBIT) call nrerror('MAXBIT too small in sobseq')
            im=(j-1)*MAXDIM
            j=min(size(x),MAXDIM)
            ix(1:j)=ieor(ix(1:j),iv(1+im:j+im))
            x(1:j)=ix(1:j)*fac
            in=in+1
        end if
    END SUBROUTINE sobseq

end module nr

MODULE splineParams
    use nrtype
    use nr
    implicit none

        PROCEDURE(template_function), SAVE, POINTER :: func
        PROCEDURE(template_derivative), SAVE, POINTER :: Dfunc
        INTEGER, SAVE :: steps
        REAL(DP), SAVE :: g_min, g_max
        REAL(DP), ALLOCATABLE, DIMENSION(:) :: evaluationPoints

contains
    !-------------------------------------
    SUBROUTINE splineInit(func1,Dfunc1,numSteps,pointsToEvaluate,a_min,a_max)
    !-------------------------------------
        PROCEDURE(template_function), POINTER, INTENT(in) :: func1
        PROCEDURE(template_derivative), POINTER, INTENT(in) :: Dfunc1
        INTEGER, INTENT(IN) :: numSteps
        REAL(DP), INTENT(IN) :: a_min, a_max
        REAL(DP), DIMENSION(:), INTENT(IN) :: pointsToEvaluate

        func => func1
        Dfunc => Dfunc1
        steps = numSteps
        g_min = a_min
        g_max = a_max
        allocate(evaluationPoints(size(pointsToEvaluate)))
        evaluationPoints = pointsToEvaluate
    END SUBROUTINE splineInit

    !-------------------------------------
    SUBROUTINE splineDestroy()
    !-------------------------------------
        deallocate(evaluationPoints)
    end SUBROUTINE splineDestroy


    !-------------------------------------
    FUNCTION splineErr(x) RESULT (z)
    !-------------------------------------
        REAL(DP), DIMENSION(:), INTENT(IN) :: x
        REAL(DP) :: z
        REAL(DP), DIMENSION(steps) :: grid, funcAtGrid, splinePoints
        REAL(DP), DIMENSION(size(evaluationPoints)) :: err
        REAL(DP), DIMENSION(1) :: pointToEvalTemp, dFx1, dFxn
        REAL(DP) :: result1, result2, pointToEval
        INTEGER(I4B) :: i

        call sub_grid_generation(grid, g_min, g_max, x(1),1)

        do i=1,steps
            pointToEvalTemp(1)=grid(i)
            funcAtGrid(i)=func(pointToEvalTemp)
            if(i == 1)then
                dFx1 = Dfunc(pointToEvalTemp)
            else if (i==(steps)) then
                dFxn = Dfunc(pointToEvalTemp)
            end if
        end do
        call spline(grid,funcAtGrid,dFx1(1),dFxn(1),splinePoints)
        do i=1,size(evaluationPoints)
            pointToEval = evaluationPoints(i)
            result1=splint(grid,funcAtGrid,splinePoints,pointToEval)
            pointToEvalTemp(1)=pointToEval
            result2=func(pointToEvalTemp)
            err(i) = abs((result2-result1)/result2)
            if(ISNAN(err(i))) then
                flush(6)
                print *,grid
                print *,splinePoints
                print *,funcAtGrid
                print *,pointToEval,result1,result2,dFx1,dFxn
                flush(6)
                stop 0
            end if
        end do

        z=maxval(err)
    END FUNCTION splineErr

    SUBROUTINE sub_grid_generation(x,xcentre,xbounds,s, gentype)
        ! Purpose: Generate grid x on [xcentre*(1-xbounds),xcentre*(1+xbounds)] with "n" steps using spacing parameter s set as follows:
        ! s=1       linear spacing
        ! s>1       left skewed grid spacing with power s
        ! 0<s<1     right skewed grid spacing with power s
        ! s<0       geometric spacing with distances changing by a factor -s^(1/(n-1)), (>1 grow, <1 shrink)
        ! s=-1      logarithmic spacing with distances changing by a factor (xmax-xmin+1)^(1/(n-1))
        ! s=0       logarithmic spacing with distances changing by a factor (xmax/xmin)^(1/(n-1)), only if xmax,xmin>0
        !
        ! INPUTS: Inputs depend on gentype. If gentype is missing:
        !               xcentre - The centre of the grid
        !               xbounds - how far away from the centre the bounds go.
        !               s - skewness of grid - see above
        !         If gentype is provided and equals 1
        !               xcentre - The min point of the grid
        !               xbounds - The max point of the grid
        !               s - skewness of grid - see above
        ! OUTPUT: x - the generated grid, size n

        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(OUT) :: x
        REAL(DP), INTENT(IN) :: xcentre,xbounds, s
        INTEGER, OPTIONAL, INTENT(IN) :: gentype
        REAL(DP) :: c ! growth rate of grid subintervals for logarithmic spacing
        REAL(DP) :: xmax, xmin
        INTEGER :: n,i

        if(present(gentype)) then
            if(gentype==1) then
                xmin=xcentre
                xmax=xbounds
            else
                PRINT *, 'grid_generation: unsuported gentype: ',gentype
                stop 0
            endif
        else
            xmax=xcentre*(1+xbounds);
            xmin=xcentre*(1-xbounds);
        endif
        n=size(x)

        FORALL(i=1:n) x(i)=(i-1)/real(n-1,DP)
        IF (s>0.0_DP) THEN
            x=x**s*(xmax-xmin)+xmin
        ELSE
            IF (s==-1.0_DP) THEN
                c=xmax-xmin+1
            ELSE
                c=-s
            END IF
            x=((xmax-xmin)/(c-1))*(c**x)-((xmax-c*xmin)/(c-1))
        END IF

        !don't want two points being the same
        do i=1,n-1
            if(x(n+1-i) .le. x((n+1)-(i+1))) then
                x((n+1)-(i+1))=x(n+1-i)-epsilon(1.0D0)
            endif
        end do

    END SUBROUTINE sub_grid_generation

END MODULE splineParams

MODULE hwutil
    use nrtype
    use splineParams
contains


    FUNCTION iterateMin(func,dfunc,numSteps,pointsToEvaluate,a_min,a_max,tol) RESULT(y)
        !-------------------------------------

        !use nelder-meade
        !
        ! INPUTS:
        !   1) func - the function we are trying to find the min value of
        !   2) startPoint - the co-ordinates of the starting point simplex
        !   3) startVals - the value of the function at the starting points
        PROCEDURE(template_function), POINTER, INTENT(in) :: func
        PROCEDURE(template_derivative), POINTER, INTENT(in) :: dfunc
        INTEGER(I4B), INTENT(IN) :: numSteps
        REAL(DP), INTENT(IN) :: a_min, a_max, tol
        REAL(DP), DIMENSION(:), INTENT(IN) :: pointsToEvaluate
        REAL(DP), DIMENSION(2) :: y

        REAL(DP), DIMENSION(2,1) :: startPoints
        REAL(DP), DIMENSION(2) :: startVals
        PROCEDURE(template_function), POINTER :: splinePointer

        splinePointer => splineErr
        call splineInit(func,dfunc,numSteps,pointsToEvaluate,a_min,a_max)
        startPoints(1,1) = 1.0
        startPoints(2,1)= 5.0
        startVals(1)=splineErr(startPoints(1,:))
        startVals(2)=splineErr(startPoints(2,:))

        CALL amoeba(startPoints,startVals, tol,splinePointer,iter)

        y(1)=startPoints(1,1)
        y(2)=splineErr(startPoints(1,:))
        call splineDestroy()
    END FUNCTION iterateMin

END MODULE hwutil

MODULE  hmwk2
    use nrtype
    use nr
    use hwutil
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
    SUBROUTINE q2a(func1,Dfunc1)
    !-------------------------------------
        PROCEDURE(template_function), POINTER, INTENT(in) :: func1
        PROCEDURE(template_derivative), POINTER, INTENT(in) :: Dfunc1

        INTEGER, PARAMETER :: numSteps=100, numSteps2=200, orderFuncCheb=30, orderFuncCheb2=50
        REAL(DP), DIMENSION(numSteps+1) :: gridPoints, funcAtGrid, splinePoints
        REAL(DP), DIMENSION(numSteps2+1) :: gridPoints2, funcAtGrid2, splinePoints2
        REAL(DP), DIMENSION(orderFuncCheb) :: chebCoeff
        REAL(DP), DIMENSION(orderFuncCheb2) :: chebCoeff2
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
        print *," "

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
        print *," "

        !note: chebev requires as many calls to the utility function, as
        !      the number of grid points. This is the same as linear, if
        !      we are making a lot of interpolations. More, if we are making
        !      just one
        print *,"Chebyshev interpolation"
        chebCoeff = chebft(.01D0,2.05D0,orderFuncCheb,func1)
        chebCoeff2 = chebft(.01D0,2.05D0,orderFuncCheb2,func1)
        print *,"Point                  ",orderFuncCheb,"          ",orderFuncCheb2
        do i=1,floor(2.0D0/.05D0)
            pointToEval = 0.05D0*i
            result1=chebev(.01D0,2.05D0,chebCoeff,pointToEval)
            result11=chebev(.01D0,2.05D0,chebCoeff2,pointToEval)
            pointToEvalTemp(1)=pointToEval
            result2=func1(pointToEvalTemp)
            print *, pointToEval,abs((result2-result1)/result2),abs((result2-result11)/result2)
        end do
        print *," "

    END SUBROUTINE q2a

    !-------------------------------------
    FUNCTION q2b(func1,Dfunc1,numSteps,pointsToEvaluate,a_min,a_max,tol) RESULT (y)
    !-------------------------------------
        PROCEDURE(template_function), POINTER, INTENT(in) :: func1
        PROCEDURE(template_derivative), POINTER, INTENT(in) :: Dfunc1
        INTEGER, INTENT(IN) :: numSteps
        REAL(DP), INTENT(IN) :: a_min, a_max,tol
        REAL(DP), DIMENSION(:), INTENT(IN) :: pointsToEvaluate
        REAL(DP), DIMENSION(2) :: y

        !Note, when solving, need to make sure we don't have two points in the grid
        !that are too close together. Put a fix in grid_gen to avoid this.
        y=iterateMin(func1,Dfunc1,numSteps,pointsToEvaluate,a_min,a_max,tol)

    END FUNCTION q2b

    !-------------------------------------
    SUBROUTINE q2c(func1,Dfunc1)
    !-------------------------------------
        PROCEDURE(template_function), POINTER, INTENT(in) :: func1
        PROCEDURE(template_derivative), POINTER, INTENT(in) :: Dfunc1

    END SUBROUTINE q2c

END MODULE hmwk2

program q2
    use nrtype
    use hmwk2
    implicit none

    PROCEDURE(template_function), POINTER :: func1, func2, func3
    PROCEDURE(template_derivative), POINTER :: Dfunc1, Dfunc2, Dfunc3
    REAL(DP) :: alpha
    INTEGER(I4B), PARAMETER :: numEvalPoints=41
    REAL(DP), DIMENSION(numEvalPoints) :: pointsToEvaluate
    INTEGER(I4B) :: i
    INTEGER(I4B), DIMENSION(10) :: gridPoints = (/ 50,100,150,200,250,500,1000,2500,5000,10000/)

    func1 => u1
    func2 => u2
    func3 => u3

    Dfunc1 => Du1
    Dfunc2 => Du2
    Dfunc3 => Du3

    CALL q2a(func1, Dfunc1)
    CALL q2a(func2, Dfunc2)
    alpha=2.0D0
    print *,"------------------------------"
    print *, "alpha: ",alpha
    print *,"------------------------------"
    CALL q2a(func3, Dfunc3)
    alpha=5.0D0
    print *,"------------------------------"
    print *, "alpha: ",alpha
    print *,"------------------------------"
    CALL q2a(func3, Dfunc3)
    alpha=10.0D0
    print *,"------------------------------"
    print *, "alpha: ",alpha
    print *,"------------------------------"
    CALL q2a(func3, Dfunc3)

    alpha=10.0D0
    do i=1,numEvalPoints
        pointsToEvaluate(i)=0.005*(i)
    end do
    print *,"Grid Points          Gamma                    Max Error"
    do i=1,10
        print *, gridPoints(i),q2b(func3,Dfunc3,gridPoints(i),pointsToEvaluate,.001D0,.05D0*(numEvalPoints+1),.0000000001D0)
    end do
CONTAINS

    !-----------------
    FUNCTION u1(point) RESULT(z)
        !
        ! the function we are trying to minimize
        ! INPUTS: point - the value at which we are evaluating the function
        ! OUTPUTS: z - the value of the function at that point
        REAL(DP), DIMENSION(:), INTENT(IN) :: point
        REAL(DP) :: z
        REAL(DP) :: x
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
        REAL(DP) :: z
        REAL(DP) :: x

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
        REAL(DP) :: x

        x=point(1)
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
        REAL(DP) :: x

        x=point(1)
        z = x**(-alpha)

    END FUNCTION Du3

end program q2
