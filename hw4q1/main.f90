!#define DOPRINTS
module utilFuncs
    use nrtype
    use nrutil
contains
    FUNCTION getSimplexAround(startPoint, distance) RESULT (Y)
        REAL(DP), DIMENSION(2), INTENT(IN) :: startPoint
        REAL(DP), INTENT(IN) :: distance
        REAL(DP), DIMENSION(3,2) :: y
        REAL(DP) :: diffx, diffy, multiplier
        REAL(DP), DIMENSION(3,2) :: startPoint2

        multiplier = distance/2.0D0
        diffy=-multiplier
        diffx=sqrt(3.0D0)*multiplier

        startPoint2(1,1)=startPoint(1)
        startPoint2(1,2)=startPoint(2)+multiplier

        startPoint2(2,1)=startPoint(1)-diffx
        startPoint2(2,2)=startPoint(2)+diffy

        startPoint2(3,1)=startPoint(1)+diffx
        startPoint2(3,2)=startPoint(2)+diffy

        y = startPoint2
    END FUNCTION getSimplexAround

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
        PROCEDURE(template_function3), POINTER, INTENT(in) :: func
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


    ! given arrays x and y of length N containing a tabulated function
    ! given values yp1 and ypn for the first derivative of the interpolating function at points 1 and N
    ! this routine returns an array y2 of length N that contains the second derivatives
    ! of the interpolating function at the tabulated points xi.
    ! If yp1 and/or ypn are equal to 1 × 1030 or larger, the routine is signaled to set the
    ! corresponding boundary condition for a natural spline, with zero second derivative on that
    ! boundary.
    SUBROUTINE spline(x,y,yp1,ypn,y2)
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

    !Given the arrays xa and ya, which tabulate a function (with the xai ’s in increasing or
    !decreasing order), and given the array y2a, which is the output from spline above, and
    !given a value of x, this routine returns a cubic-spline interpolated value. The arrays xa, ya
    !and y2a are all of the same size.
    FUNCTION splint(xa,ya,y2a,x)
        IMPLICIT NONE
        REAL(dp), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
        REAL(dp), INTENT(IN) :: x
        REAL(dp) :: splint
        INTEGER(I4B) :: khi,klo,n
        REAL(dp) :: a,b,h
        n=assert_eq(size(xa),size(ya),size(y2a),'splint')
        klo=max(min(locate(xa,x),n-1),1)
        khi=klo+1
        h=xa(khi)-xa(klo)
        if (h == 0.0) call nrerror('bad xa input in splint')
        a=(xa(khi)-x)/h
        b=(x-xa(klo))/h
        splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_dp
    END FUNCTION splint

    FUNCTION locate(xx,x)
        IMPLICIT NONE
        REAL(dp), DIMENSION(:), INTENT(IN) :: xx
        REAL(dp), INTENT(IN) :: x
        INTEGER(I4B) :: locate
        INTEGER(I4B) :: n,jl,jm,ju
        LOGICAL :: ascnd
        n=size(xx)
        ascnd = (xx(n) >= xx(1))
        jl=0
        ju=n+1
        do
            if (ju-jl <= 1) exit
            jm=(ju+jl)/2
            if (ascnd .eqv. (x >= xx(jm))) then
                jl=jm
            else
                ju=jm
            end if
        end do
        if (x == xx(1)) then
            locate=1
        else if (x == xx(n)) then
            locate=n-1
        else
            locate=jl
        end if
    END FUNCTION locate

    SUBROUTINE tridag(a,b,c,r,u)
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


    SUBROUTINE mov3(a,b,c,d,e,f)
        REAL(dp), INTENT(IN) :: d,e,f
        REAL(dp), INTENT(OUT) :: a,b,c
        a=d
        b=e
        c=f
    END SUBROUTINE mov3

    FUNCTION brent(myfunc, ax,bx,cx,tol,xmin)
        USE nrtype; USE nrutil, ONLY : nrerror
        IMPLICIT NONE
        REAL(dp), INTENT(IN) :: ax,bx,cx,tol
        REAL(dp), INTENT(OUT) :: xmin
        REAL(dp) :: brent
        PROCEDURE(template_function), POINTER :: myfunc
        INTEGER(I4B), PARAMETER :: ITMAX=1000
        REAL(dp), PARAMETER :: CGOLD=0.3819660_dp,ZEPS=1.0e-3_dp*epsilon(ax)
        INTEGER(I4B) :: iter
        REAL(dp) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
        a=min(ax,cx)
        b=max(ax,cx)
        v=bx
        w=v
        x=v
        e=0.0
        fx=myfunc(x)
        fv=fx
        fw=fx
        do iter=1,ITMAX
            xm=0.5_dp*(a+b)
            tol1=tol*abs(x)+ZEPS
            tol2=2.0_dp*tol1
            if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) then
                xmin=x
                brent=fx
                RETURN
            end if
            if (abs(e) > tol1) then
                r=(x-w)*(fx-fv)
                q=(x-v)*(fx-fw)
                p=(x-v)*q-(x-w)*r
                q=2.0_dp*(q-r)
                if (q > 0.0) p=-p
                q=abs(q)
                etemp=e
                e=d
                if (abs(p) >= abs(0.5_dp*q*etemp) .or. &
                    p <= q*(a-x) .or. p >= q*(b-x)) then
                    e=merge(a-x,b-x, x >= xm )
                    d=CGOLD*e
                else
                    d=p/q
                    u=x+d
                    if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
                end if
            else
                e=merge(a-x,b-x, x >= xm )
                d=CGOLD*e
            end if
            u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
            fu=myfunc(u)
            if (fu <= fx) then
                if (u >= x) then
                    a=x
                else
                    b=x
                end if
                call shft(v,w,x,u)
                call shft(fv,fw,fx,fu)
            else
                if (u < x) then
                    a=u
                else
                    b=u
                end if
                if (fu <= fw .or. w == x) then
                    v=w
                    fv=fw
                    w=u
                    fw=fu
                else if (fu <= fv .or. v == x .or. v == w) then
                    v=u
                    fv=fu
                end if
            end if
        end do
        call nrerror('brent: exceed maximum iterations')
    END FUNCTION brent

    FUNCTION linear(func,gridPoints, evalPoint) RESULT(z)
        ! INPUTS: func - a vector of function values
        !         gridPoints - the points at which the function is evaluated
        ! OUTPUTS: evalPoint - the point at which we want to evaluate

        !only makes sense in two dimensions
        REAL(DP),DIMENSION(:),INTENT(in) :: gridPoints
        REAL(DP),DIMENSION(size(gridPoints)), INTENT(in) :: func
        REAL(DP), INTENT(in) :: evalPoint
        REAL(DP) :: z
        REAL :: eps=epsilon(1.0D0)

        !local variables
        REAL(DP) :: slope
        integer, DIMENSION(1) :: indexInGridTemp
        integer :: closest, nextPoint
        LOGICAL :: cond1, cond2

        indexInGridTemp=minloc(abs(evalPoint-gridPoints))
        closest = indexInGridTemp(1)

        !*****
        ! Just return the function value if we are right on the point
        !******
        cond1 = evalPoint .gt. (gridPoints(closest)-eps)
        cond2 = evalPoint .lt. (gridPoints(closest)+eps)
        if(cond1 .and. cond2) then
            z = func(closest)
            return
        end if

        if(evalPoint>gridPoints(closest)) then
            nextPoint = closest + 1
        else
            nextPoint = closest - 1
        end if

        ! if we are past the max grid point, use the left derivative
        if(nextPoint > size(gridPoints))then
            slope=(func(closest)-func(closest-1))/(gridPoints(closest)-gridPoints(closest-1))
        else if(nextPoint < 1)then
            ! if we are before the min grid point, use the right derivative
            slope=(func(2)-func(1))/(gridPoints(2)-gridPoints(1))
        else
            slope=(func(nextPoint)-func(closest))/(gridPoints(nextPoint)-gridPoints(closest))
        endif
        z=func(closest)+slope*(evalPoint-gridPoints(closest))

        if(isNAN(z))then
            print *,"Error. NaN in linear.",evalPoint,func(closest),slope,gridPoints(closest),gridPoints(nextPoint)
            print *,closest, nextPoint
            stop 0
        end if
    END FUNCTION linear

    SUBROUTINE shft(a,b,c,d)
        REAL(dp), INTENT(OUT) :: a
        REAL(dp), INTENT(INOUT) :: b,c
        REAL(dp), INTENT(IN) :: d
        a=b
        b=c
        c=d
    END SUBROUTINE shft

    FUNCTION splineDeriv(gridPoints, fnVal, splinePoints) RESULT(y)
        REAL(DP), DIMENSION(:), INTENT(IN) :: gridPoints
        REAL(DP), DIMENSION(size(gridPoints)), INTENT(IN) :: fnVal
        REAL(DP), DIMENSION(size(GridPoints)), INTENT(IN) :: splinePoints
        REAL(DP), DIMENSION(size(GridPoints)) :: y
        REAL(DP) :: h
        INTEGER :: counter

        y(1) = 0.0D0
        do counter=2,size(gridPoints)
            h= gridPoints(counter)-gridPoints(counter-1)
            y(counter) = (fnVal(counter)-fnVal(counter-1))/h - &
                & (splinePoints(counter)+2.0D0*splinePoints(counter-1))/6.0D0*h
        end do

    end function splineDeriv
end module utilFuncs

module brentWrapper
    USE nrtype
    USE nrutil
    use utilFuncs

    implicit none

    INTEGER :: n_s, currentState, currentCapital
    REAL(DP), allocatable, DIMENSION(:) :: a, s
    REAL(DP), allocatable, DIMENSION(:,:) :: v, y2, transition
    REAL(DP) :: alpha,rho,beta,r,w
    LOGICAL :: doSpline

    PRIVATE n_s,currentState,currentCapital,a,s,v,y2,transition,alpha,beta,rho,r,doSpline,w

contains

    subroutine wrapperCreate(grid,mytransition,mystates,myalpha,myrho,mybeta,myr,myw)
        REAL(DP), DIMENSION(:) ::grid, mystates
        REAL(DP), DIMENSION(:,:) :: mytransition
        REAL(DP), INTENT(IN) :: myalpha,myrho,mybeta,myr,myw

        INTEGER :: tempCounter

        n_s=size(mystates)
        allocate(s(n_s))
        s=mystates

        tempCounter=size(grid)
        allocate(a(tempCounter))
        a=grid

        allocate(transition(n_s,n_s))
        transition=mytransition

        alpha=myalpha
        rho=myrho
        beta=mybeta
        r=myr
        w=myw
    end subroutine wrapperCreate

    subroutine wrapperInit(values,splineParams)
        REAL(DP), DIMENSION(:,:) :: values
        REAL(DP), OPTIONAL, DIMENSION(:,:) :: splineParams

        allocate(v(n_s,size(values,dim=2)))
        v=values

        if(PRESENT(splineParams)) then
            allocate(y2(n_s,size(values,dim=2)))
            y2=splineParams
            doSpline=.TRUE.
        end if
    end subroutine wrapperInit

    subroutine wrapperClean()
        deallocate(v)
        if(allocated(y2)) then
            deallocate(y2)
        end if
    end subroutine wrapperClean

    subroutine wrapperDestroy()
        deallocate(a)
        deallocate(transition)
        deallocate(s)
    end subroutine wrapperDestroy

    function callBrent(state,capital, myfunc, kr1, kr2, kr3, tol, minPoint) RESULT (y)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: state, capital
        REAL(dp), INTENT(IN) :: kr1,kr2,kr3,tol
        REAL(dp), INTENT(OUT) :: minPoint
        PROCEDURE(template_function), POINTER :: myfunc
        REAL(dp) :: y

        currentState = state
        currentCapital = capital

        y=brent(myfunc, kr1, kr2, kr3, tol, minPoint)

    end function callBrent

    function valueFunction(x) RESULT(z)
        IMPLICIT NONE
        REAL(KIND=8), INTENT(IN) :: x
        REAL(KIND=8) :: z
        real*8, dimension(1)        ::  rp
        integer                     ::  jj
        real*8                      :: temp
        real*8, dimension(n_s,1)            ::  interpolated

        do jj=1,n_s
            if (doSpline) then
                interpolated(jj,1)=splint(a,v(jj,:),y2(jj,:),x)
            else
                interpolated(jj,1)=linear(v(jj,:),a,x)
            end if
        end do
        rp=matmul(transition(currentState,:),interpolated**alpha)
        temp = a(currentCapital)*(1+r)+s(currentState)*w-x
        if(temp < 0.0D0) then
            temp=0.0D0
        end if
        z=-(((temp)**rho+beta*rp(1)**(rho/alpha))**(1/rho))
        if(z<-10D50)then
            print *,"value function too low"
            print *,temp,rho,beta,rp(1),alpha
            stop 0
        end if
    end function valueFunction

end module brentWrapper


module aiyagariSolve
    use rowenhorst_module
    USE nrtype
    USE nrutil
    USE brentWrapper
    implicit none

    INCLUDE 'mpif.h'
    real*8, parameter                   ::  phi=.9D0,sigma=.4D0
    integer, parameter                  ::  n_s=9
    integer, parameter                  ::  n_a=301

    !******************
    ! These are here because of screwey splines
    !******************
    integer, parameter                  ::  bottomChop=0   !the number of first points we want to ignore
    integer, parameter                  ::  topChop=0       !the number of end points we want to ignore

    real*8, parameter                   ::  curv=3.0D0,a_max=35.0D0
    real*8, parameter                   ::  beta=0.90
    integer, parameter                  ::  maxit=2000
    real*8,parameter                    ::  toll=1D-8,tol2=1D-8
    LOGICAL                             :: doSpline = .TRUE.

    real*8, dimension(n_s)              ::  s,stationary
    real*8, dimension(n_s,n_s)          ::  transition
    real*8, dimension(n_a)              ::  a
    real*8                              ::  a_min
    real*8                              :: capShare
    real*8                              :: rho,alpha,RRA,EIS,mysigma,wFixed=1.0D0,&
        & rFixed=0.1D0
    real*8, dimension(n_s,n_a,maxit)    ::  v,g
    real*8, dimension(n_s,n_a)          ::  lastStateV
    integer, save                       ::  reportNum, callCount
    PROCEDURE(template_function), POINTER :: funcParam
    PROCEDURE(template_function2), POINTER :: deriv1Func, deriv2Func
    logical                            :: firstCall = .true.

    character(LEN=20)                   :: policyOutput
    character(LEN=20)                   :: distribOutput

    PRIVATE phi, sigma, n_s, n_a, curv, a_max, beta, maxit, toll, doSpline
    PRIVATE s, stationary, a, a_min,rho,alpha,RRA,EIS,mysigma,v,g
    PRIVATE reportNum,funcParam,policyOutput
    PRIVATE aggregateBonds, firstCall

contains
    subroutine setParams(myrra, myeis, func, d1Func, d2Func, file1, file2, every, &
        & capitalShare, doLinear, myR, myW)
        REAL(KIND=8), INTENT(IN) :: myrra, myeis,capitalShare
        PROCEDURE(template_function), POINTER, INTENT(IN) :: func
        PROCEDURE(template_function2), POINTER, INTENT(IN) :: d1Func
        PROCEDURE(template_function2), POINTER, INTENT(IN) :: d2Func
        character(LEN=*),INTENT(IN) :: file1
        character(LEN=*),INTENT(IN) :: file2
        INTEGER, OPTIONAL, INTENT(IN) :: every
        LOGICAL, OPTIONAL, INTENT(IN) :: doLinear
        REAL(DP), OPTIONAL, INTENT(IN) :: myR, myW

        capShare=capitalShare
        RRA=myrra
        EIS=myeis
        funcParam => func
        deriv1Func => d1Func
        deriv2Func => d2Func

        rho=1.0D0-1.0D0/EIS
        alpha=1.0D0-RRA
        reportNum = 50

        if(PRESENT(every)) then
            reportNum=every
        end if

        if(PRESENT(doLinear))then
            doSpline=.not. doLinear
        end if

        if(PRESENT(myW))then
            wFixed=myW
        end if

        if(PRESENT(myR))then
            rFixed=myR
        end if

        policyOutput = file1
        distribOutput = file2
        callCount = 0
        firstCall = .true.

        !**************************************************************************
        ! We use the rowenhorst method to obtain the transition matrix GAMMA,
        ! the stationary probability distribution stationary, and the shocks s.
        !**************************************************************************
        call rouwenhorst(phi,sigma,transition,s,stationary)
        s=1+s/(s(n_s)+0.1D0)
    end subroutine setParams

    function aggregateBondsFixedW(r) RESULT (z)
        ! inputs: r - the interest rate to test
        ! outputs: z - the aggregate level of borrowing
        REAL(KIND=8), INTENT(IN) :: r
        REAL(KIND=8) :: z
        REAL(DP):: totalCapital, temp
        !************
        ! Timing variables
        !************
        real(DP) :: startTime, endTime
        INTEGER :: rank, ierr

        CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        if(rank ==0)then
            call CPU_TIME(startTime)
        end if

        rFixed = r
        temp = (r/capShare)**(1.0D0/(capShare-1))
        wFixed=(1-capShare)*temp**capShare

        totalCapital=aggregateBonds(r,wFixed)

        z=abs(totalCapital-temp)
#ifdef DOPRINTS
        if(rank ==0)then
            call CPU_TIME(endTime)
            print *,"Implied Capital: ",temp, " Calculated Capital: ", totalCapital, "Time to compute: ",endTime-startTime
            flush(6)
        end if
#endif
    end function aggregateBondsFixedW

    function aggregateBondsSetK(k) RESULT (z)
        ! inputs: r - the interest rate to test
        ! outputs: z - the aggregate level of borrowing
        REAL(KIND=8), INTENT(IN) :: k
        REAL(KIND=8) :: z
        REAL(DP):: totalCapital, temp
        !************
        ! Timing variables
        !************
        real(DP) :: startTime, endTime
        INTEGER :: rank, ierr

        CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        if(rank ==0)then
            call CPU_TIME(startTime)
        end if

        rFixed = capShare*k**(capShare-1)
        temp = k
        wFixed=(1-capShare)*temp**capShare
        totalCapital=aggregateBonds(rFixed,wFixed)

        z=abs(totalCapital-temp)
#ifdef DOPRINTS
        if(rank ==0)then
            call CPU_TIME(endTime)
            print *,"Set Capital: ",temp, " Calculated Capital: ", totalCapital, "r: ",rFixed,"Time to compute: ",endTime-startTime
            flush(6)
        end if
#endif
    end function aggregateBondsSetK

    function aggregateBonds(r,w) RESULT (z)
        ! inputs: r - the interest rate to test
        ! outputs: z - the aggregate level of borrowing
        REAL(KIND=8), INTENT(IN) :: r,w
        REAL(KIND=8) :: z

        INTEGER :: iterCount, i,j
        real*8  ::  incr, totalCapital
        real*8, dimension(n_s,n_a-bottomChop-topChop)          ::  steadyStateCapital

        callCount = callCount+1
        a_min = min(-(s(1)*wFixed)/r + 1,0.0D0)

        !**************************************************************************
        ! We set up the grid of asset values based on the curvature, curv
        ! the minimum and maximum values in the grid a_min a_max
        !**************************************************************************
        incr=(a_max-a_min)/(n_a-1)
        a=(/ (    incr*real(i-1,8),i=1,n_a    ) /)
        a=a**curv
        a=a/(a_max-a_min)**(curv-1)+a_min

        !**************************************************************************
        ! we initialize the value function
        ! and we set up an initial guess for it
        !**************************************************************************
        v=0D0
        if(firstCall)then
            do i=1,n_s
                v(i,:,1)=(a-a_min)**2
            !        v(iter,:,1)=0.0
            end do
            firstCall = .false.
        else
            v(:,:,1)=lastStateV
        end if
        g=0D0

        call wrapperCreate(a,transition,s,alpha,rho,beta,r,w)
        iterCount=getPolicyForInterest(r)
        call wrapperDestroy()

        open(unit=1,file=policyOutput)
        write(1,*) a(:)
        do i=1,n_s
            write(1,*) g(i,:,iterCount)
        end do
        close(1)

            !*****************************************************
            ! find the steady state capital
            !*****************************************************
        call findSteadyState(g(:,:,iterCount),steadyStateCapital)

        !**************************************************
        ! now, calculate total capital
        !**************************************************
        totalCapital = dot_product(sum(steadyStateCapital,dim=1),a(bottomChop+1:n_a-topChop))

        z = totalCapital

    end function aggregateBonds

    subroutine findSteadyState(capitalPolicyOrig, statDist)
        ! NOTE: This code is faster single processed than multi-processed
        ! You can confirm this by uncommenting the below define and rebuilding
        !#define MPI_FSS

        !INPUTS: capitalPolicy - the policy function for each state
        !OUTPUTS: statDist - the stationary dist (note: pdf, not cdf)
        INTEGER, parameter :: capitalCount = n_a-bottomChop-topChop
        REAL(DP), dimension(n_s,n_a), INTENT(IN) :: capitalPolicyOrig
        real(DP), dimension(n_s,capitalCount), intent(out) :: statDist
        real(DP), dimension(n_s,capitalCount) ::f_o, f_o_hat, f_n
        real(DP) :: diff, temp
        INTEGER :: i,j, counter
        REAL(DP), dimension(n_s,capitalCount) :: capitalPolicy
        REAL(DP), dimension(capitalCount) :: newCapital


#ifdef MPI_FSS
        !********
        !Variables for MPI
        !**********
        INTEGER :: rank, mysize, source, dest, tag, mpicounter, ierr, it
        REAL(DP), dimension(capitalCount) :: inmsg
        LOGICAL :: allData=.false.
        INTEGER,dimension(MPI_STATUS_SIZE):: stat
#endif
        newCapital = a(bottomChop+1:n_a-topChop)
        capitalPolicy = capitalPolicyOrig(:,bottomChop+1:n_a-topChop)

#ifdef MPI_FSS
        !* check monotonicity of capital policy
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mysize, ierr)

        if(rank == 0) then
#endif
            do i=1,n_s
#ifdef MPI_FSS
                dest = mod(i,mysize);
                if(dest .ne. 0) then
                    tag = i
                    call MPI_SEND(i, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD, ierr)
                else !if (dest .ne. 0)
#endif
                    do j=2,capitalCount
                        if(capitalPolicy(i,j)<capitalPolicy(i,j-1)) then
                            print *,"not monotonic.", i, j
                            print *,capitalPolicy(i,j-1),capitalPolicy(i,j)
                            stop 0
                        end if
                    end do
#ifdef MPI_FSS
                end if
#endif
            end do
#ifdef MPI_FSS
        else !if(rank==0)
            do while(.not. allData)
                source = 0
                call MPI_RECV(i, 1, MPI_INTEGER, source, MPI_ANY_TAG, MPI_COMM_WORLD, stat, ierr)
                do j=2,capitalCount
                    if(capitalPolicy(i,j)<capitalPolicy(i,j-1)) then
                        print *,"not monotonic.", i, j
                        print *,capitalPolicy(i,j-1),capitalPolicy(i,j)
                        stop 0
                    end if
                end do

                !***************
                !check to see if we should expect any more
                !***************
                if(i+mysize>n_s)then
                    allData=.true.
                end if
            end do
        end if

        allData=.false.
#endif
        !setting initial guess to uniform dist across asset grid. First we have it
        !a cdf, and then convert to the appropriate pdf
        do i=1,n_s
            do j=1,capitalCount
                statDist(i,j) = (newCapital(j) - newCapital(1))/(newCapital(capitalCount)-newCapital(1))*dble(1.0_dp/n_s)
            end do
        end do

        f_n=statDist

        ! time to iterate
        diff=100
        counter = 0
        do while((diff>tol2) .and. (counter<maxit))
            counter = counter + 1

            f_o=f_n

            do i=1,n_s
                do j=1,capitalCount
                    f_o_hat(i,j) = linear(f_o(i,:),capitalPolicy(i,:),newCapital(j))
                end do
            end do

            ! need to make sure monotonic and bounded between 0 and 1
#ifdef MPI_FSS
            if(rank == 0) then
#endif
                do i=1,n_s
#ifdef MPI_FSS
                    dest = mod(i,mysize);
                    if(dest .ne. 0) then
                        tag = i
                        call MPI_SEND(i, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD, ierr)
                    else !if (dest .ne. 0)
#endif
                        diff=-1e-5
                        do j=1,capitalCount
                            if(f_o_hat(i,j)>1.0_dp) then
                                f_o_hat(i,j)=1.0_dp
                            else if (f_o_hat(i,j)<0.0_dp) then
                                f_o_hat(i,j)=0.0D0
                            else if (isnan(f_o_hat(i,j))) then
                                print *, "Error: f_o_hat is NAN. Counter:",counter," shock: ",i,"element: ",j
                                print *,"Capital Grid"
                                print *,newCapital(:)
                                print *,"Policy Fn"
                                print *,capitalPolicy(j,:)
                                print *,"f_o"
                                print *,f_o(i,:)
                                print *,"f_o_hat"
                                print *,f_o_hat(i,:)
                                stop 0
                            end if

                            !add a test: if non-monotonic, fail
                            if(diff>f_o_hat(i,j)) then
                                print *, "Error: Non-monotonic cdf function. Counter:",counter," element: ", j, " shock: ",i
                                print *,"value_old: ",diff, " value_new: ",f_o_hat(i,j)
                                print *,"f_o"
                                print *,capitalPolicy(i,:)
                                print *,f_o(i,:)
                                print *,"f_o_hat"
                                print *,newCapital
                                print *,f_o_hat(i,:)
                                stop 0
                            end if
                            diff=f_o_hat(i,j)
                        end do
#ifdef MPI_FSS
                        do mpicounter=i-(mysize-1),i-1
                            source = mod(mpicounter,mysize)
                            call MPI_RECV(f_o_hat(mpicounter,:), capitalCount, MPI_REAL8, source, mpicounter,&
                                &MPI_COMM_WORLD, stat, ierr)
                        end do
                    end if
#endif
                end do
#ifdef MPI_FSS
            else !if (rank == 0)
                do while(.not. allData)
                    source = 0
                    call MPI_RECV(i, 1, MPI_INTEGER, source, MPI_ANY_TAG, MPI_COMM_WORLD, stat, ierr)

                    diff=-1e-5
                    do j=1,capitalCount
                        if(f_o_hat(i,j)>1.0_dp) then
                            f_o_hat(i,j)=1.0_dp
                        else if (f_o_hat(i,j)<0.0_dp) then
                            f_o_hat(i,j)=0.0D0
                        else if (isnan(f_o_hat(i,j))) then
                            print *, "Error: f_o_hat is NAN. Counter:",counter," shock: ",i,"element: ",j
                            print *,"Capital Grid"
                            print *,newCapital(:)
                            print *,"Policy Fn"
                            print *,capitalPolicy(j,:)
                            print *,"f_o"
                            print *,f_o(i,:)
                            print *,"f_o_hat"
                            print *,f_o_hat(i,:)
                            stop 0
                        end if

                        !add a test: if non-monotonic, fail
                        if(diff>f_o_hat(i,j)) then
                            print *, "Error: Non-monotonic cdf function. Counter:",counter," element: ", j, " shock: ",i
                            print *,"value_old: ",diff, " value_new: ",f_o_hat(i,j)
                            print *,"f_o"
                            print *,capitalPolicy(i,:)
                            print *,f_o(i,:)
                            print *,"f_o_hat"
                            print *,newCapital
                            print *,f_o_hat(i,:)
                            stop 0
                        end if
                        diff=f_o_hat(i,j)
                    end do
                    inmsg(1:capitalCount)=f_o_hat(i,:)
                    dest = 0
                    tag = i
                    call MPI_SEND(inmsg, capitalCount, MPI_REAL8, dest, tag, MPI_COMM_WORLD, ierr)

                    !***************
                    !check to see if we should expect any more
                    !***************
                    if(tag+mysize>n_s)then
                        allData=.true.
                    end if
                end do
            end if
#endif

#ifdef MPI_FSS
            if(rank == 0) then
                !*********************
                ! we need to capture things we sent
                !*********************
                mpicounter=mod(n_s,mysize)-1
                do j=n_s-mpicounter,n_s
                    source = mod(j,mysize)
                    call MPI_RECV(f_o_hat(j,:), capitalCount, MPI_REAL8, source, j,&
                        & MPI_COMM_WORLD, stat, ierr)
                end do

                !***************************
                ! Now, we need to communicate the prob function for this iteration across all programs
                !***************************
                do j=1,n_s
                    do it=1,mysize-1
                        call MPI_SEND(f_o_hat(j,:), capitalCount, MPI_REAL8, it, 1000, MPI_COMM_WORLD, ierr)
                    end do
                end do
            else
                source = 0
                do j=1,n_s
                    call MPI_RECV(f_o_hat(j,:), capitalCount, MPI_REAL8, source, MPI_ANY_TAG, MPI_COMM_WORLD, stat, ierr)

                    !Check tag to make sure it is what we have set as "pdf function", i.e. -1
                    if(stat(MPI_TAG) .ne. 1000)then
                        print *,"Error, ", rank, "doesn't know what prob function it received,",stat(MPI_TAG)
                        stop 0
                    end if
                end do
                allData = .false.
            end if
#endif
            f_n = matmul(transpose(transition), f_o_hat)

            !* Fix so that total cdf is 1
            do i=1,n_s
                do j=1,capitalCount
                    temp = f_n(i,j) / f_n(i,capitalCount)

                    if (isnan(temp))then
                        print *,"nan: ",i,j, counter
                        print *,f_n(i,j),f_n(i,capitalCount)
                        print *,f_o_hat(i,:)
                        print *,f_o_hat(1,j),f_o_hat(2,j),f_o_hat(3,j)
                        flush(6)
                        stop 0
                    end if

                    f_n(i,j)=temp
                end do
            end do

            diff = maxval(abs(f_n - f_o))
        end do

        do i=1,n_s
            f_n(i,:)=f_n(i,:)*stationary(i)
        end do

        statDist = f_n
        do i=capitalCount,2,-1
            statDist(:,i) = statDist(:,i) - statDist(:,i-1)
        end do
        do i=1,n_s
            statDist(i,1) = max(0.0D0,stationary(i)-sum(statDist(i,:)))
        end do

        open(unit=1,file=distribOutput)
        do j=1,capitalCount
            write(1,*) a(j+bottomChop),statDist(:,j)
        end do
        close(1)

    end subroutine findSteadyState

    function getPolicyForInterest(r) result(y)
        REAL(KIND=8), INTENT(IN) :: r
        INTEGER :: y
        INTEGER ::  i,it,j,iter
        REAL(DP) :: tempD,tempD2
        real*8, dimension(n_s,n_a)          ::  y2
        real*8                              ::  kr1,kr2,kr3
        !********
        !Variables for MPI
        !**********
        INTEGER :: rank, mysize, source, dest, tag, mpicounter, ierr
        REAL(DP), dimension(n_a*2) :: inmsg
        LOGICAL :: allData=.false.
        INTEGER,dimension(MPI_STATUS_SIZE):: stat


        !**************************************************************************
        ! we begin the iteration
        !**************************************************************************
        do iter=2,maxit
            if(doSpline) then
                do j=1,n_s
                    tempD=(v(j,2,iter-1)-v(j,1,iter-1))/(a(2)-a(1))
                    tempD2=(v(j,n_a,iter-1)-v(j,n_a-1,iter-1))/(a(n_a)-a(n_a-1))
                    call spline(a,v(j,:,iter-1),tempD,tempD2,y2(j,:))
                end do

                call wrapperInit(v(:,:,iter-1),y2)
                CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
                CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mysize, ierr)

                if(rank == 0) then
                    do i=1,n_s
                        dest = mod(i,mysize);
                        if(dest .ne. 0) then
                            tag = i
                            call MPI_SEND(i, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD, ierr)
                        else !if (dest .ne. 0)
                            do it=1,n_a
                                !ensure monotone policy function
                                if(it==1)then
                                    kr1=a(1)
                                else
                                    kr1=g(i,it-1,iter)
                                end if
                                kr3=min(a(it)*(1+r)+s(i)*wFixed,a(n_a))
                                kr2=(kr1+kr3)/2D0
                                v(i,it,iter)=-callBrent(i,it,funcParam,kr1,kr2,kr3,1D-6,g(i,it,iter))
                                if(isNAN(v(i,it,iter)))then
                                    print *,"error, value function is nan"
                                    stop 0
                                end if
                            end do
                            do j=i-(mysize-1),i-1
                                source = mod(j,mysize)
                                call MPI_RECV(inmsg, 2*n_a, MPI_REAL8, source, j, MPI_COMM_WORLD, stat, ierr)
                                v(j,:,iter) = inmsg(1:n_a)
                                g(j,:,iter) = inmsg(n_a+1:2*n_a)
                            end do
                        end if
                    end do
                else !if(rank == 0)
                    do while(.not. allData)
                        source = 0
                        call MPI_RECV(i, 1, MPI_INTEGER, source, MPI_ANY_TAG, MPI_COMM_WORLD, stat, ierr)
                        do it=1,n_a
                             !ensure monotone policy function
                            if(it==1)then
                                kr1=a(1)
                            else
                                kr1=g(i,it-1,iter)
                            end if
                            kr3=min(a(it)*(1+r)+s(i)*wFixed,a(n_a))
                            kr2=(kr1+kr3)/2D0
                            v(i,it,iter)=-callBrent(i,it,funcParam,kr1,kr2,kr3,1D-10,g(i,it,iter))
                            if(isNAN(v(i,it,iter)))then
                                print *,"error, value function is nan"
                                stop 0
                            end if
                        end do
                        inmsg(1:n_a)=v(i,:,iter)
                        inmsg(n_a+1:2*n_a)=g(i,:,iter)
                        dest = 0
                        tag = i
                        call MPI_SEND(inmsg, 2*n_a, MPI_REAL8, dest, tag, MPI_COMM_WORLD, ierr)

                        !***************
                        !check to see if we should expect any more
                        !***************
                        if(tag+mysize>n_s)then
                            allData=.true.
                        end if
                    end do
                end if
                call wrapperClean()
            else
                call wrapperInit(v(:,:,iter-1))
                do i=1,n_s
                    do it=1,n_a
                        kr1=a(1)
                        kr3=min(a(it)*(1+r)+s(i)*wFixed,a(n_a))
                        kr2=(kr1+kr3)/2D0
                        v(i,it,iter)=-callBrent(i,it,funcParam,kr1,kr2,kr3,1D-10,g(i,it,iter))
                        if(isNAN(v(i,it,iter)))then
                            print *,"error, value function is nan"
                            stop 0
                        end if
                    end do
                end do
                call wrapperClean()
            end if
            if(doSpline) then
                if(rank == 0) then
                    !*********************
                    ! we need to capture things we sent
                    !*********************
                    mpicounter=mod(n_s,mysize)-1
                    do j=n_s-mpicounter,n_s
                        source = mod(j,mysize)
                        call MPI_RECV(inmsg, 2*n_a, MPI_REAL8, source, j, MPI_COMM_WORLD, stat, ierr)
                        v(j,:,iter) = inmsg(1:n_a)
                        g(j,:,iter) = inmsg(n_a+1:2*n_a)
                    end do
                    !***************************
                    ! Now, we need to communicate the value function for this iteration across all programs
                    !***************************

                    do j=1,n_s
                        do it=1,mysize-1
                            call MPI_SEND(v(j,:,iter), n_a, MPI_REAL8, it, 1000, MPI_COMM_WORLD, ierr)
                            call MPI_SEND(g(j,:,iter), n_a, MPI_REAL8, it, 1000, MPI_COMM_WORLD, ierr)
                        end do
                    end do
                else
                    source = 0
                    do j=1,n_s
                        call MPI_RECV(v(j,:,iter), n_a, MPI_REAL8, source, MPI_ANY_TAG, MPI_COMM_WORLD, stat, ierr)

                        !Check tag to make sure it is what we have set as "value function", i.e. -1
                        if(stat(MPI_TAG) .ne. 1000)then
                            print *,"Error, don't know what value function we received,",stat(MPI_TAG)
                            stop 0
                        end if

                        call MPI_RECV(g(j,:,iter), n_a, MPI_REAL8, source, MPI_ANY_TAG, MPI_COMM_WORLD, stat, ierr)
                        if(stat(MPI_TAG) .ne. 1000)then
                            print *,"Error, don't know what policy function we received,",stat(MPI_TAG)
                            stop 0
                        end if

                    end do
                    allData = .false.
                end if
            end if
            if( (mod(iter,reportNum)==0) .and. (rank==0))then
                print *,"r:",r,"w:",wFixed,"iter: ",iter,&
                    &"diff: ",maxval(abs(v(:,:,iter)-v(:,:,iter-1)))
                flush(6)
            end if
            if (     maxval(abs(v(:,:,iter)-v(:,:,iter-1)))    .lt.    toll     ) then
#ifdef DOPRINTS
                if(rank == 0) then
                    print*,"done: ",iter, "r: ",r, "w: ",wFixed
                    flush(6)
                end if
#endif
                exit
            end if
        end do

        y=min(iter,maxit)
        lastStateV = v(:,:,y)

    end function         getPolicyForInterest

end module aiyagariSolve

!**************************************************************************
!
!**************************************************************************
program main

    use aiyagariSolve
    use brentWrapper
    REAL(DP) :: RRA, EIS, intDiff,xmin
    REAL(DP) :: capitalShare = 0.33D0
    PROCEDURE(template_function), POINTER :: func
    PROCEDURE(template_function2), POINTER :: d1func, d2func
    PROCEDURE(template_function3), POINTER :: func2
    REAL(DP) :: temp
    REAL(DP), DIMENSION(2) :: startPoint
    REAL(DP), DIMENSION(3,2) :: startPoint2
    REAL(DP), DIMENSION(3) :: startVals
    INTEGER :: temp2,printEvery=500,whichSet=3
    character(LEN=15) :: arg1,arg2,arg3,arg4
    !************
    ! Timing variables
    !************
    real(DP) :: startTime, endTime
    INTEGER :: rank, ierr

    temp2=COMMAND_ARGUMENT_COUNT()
    if(temp2 > 0)then
        call GET_COMMAND_ARGUMENT(1, arg1, whichSet)
        read (arg1,*) whichSet
    end if
    if(temp2 > 1)then
        call GET_COMMAND_ARGUMENT(2, arg1, printEvery)
        read (arg1,*) printEvery
    end if

    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    if(rank ==0)then
        print *,"Question               RRA                    EIS                        r  &
                &                         K                    Error                       Time(s)"        
        flush(6)
    end if

    if(whichSet/=2)then
        !***********************************
        !Question a
        !***********************************
        RRA=2.0D0
        EIS=0.5D0
        func => valueFunction
        d1func => d1prod
        d2func => d2prod

        arg1="policyR2Ep5"
        arg2= "distribR2Ep5"
        call setParams(RRA, EIS, func, d1func, d2func,arg1,arg2 , printEvery, capitalShare,&
            &.FALSE., 0.1D0, 1.0D0)
        func => aggregateBondsFixedW
        if(rank ==0)then
            call CPU_TIME(startTime)
        end if
        intDiff=brent(func,0.01D0,0.09D0,0.12D0,1.0D-8,xmin)
        if(rank ==0)then
            call CPU_TIME(endTime)
            print *,"a             ",RRA,EIS,xmin,impliedCapital(xmin), intDiff, endTime-startTime
            flush(6)
        end if

        !***********************************
        !Question bp1
        !***********************************
        EIS=0.9D0
        do temp2=0,5
            RRA=2.0D0+(20.0D0-2.0D0)/5*temp2
            func => valueFunction
            d1func => d1prod
            d2func => d2prod
            arg1="policyR"
            write(arg2,'(I1)') temp2
            arg3="Ep9"
            arg1=TRIM(arg1)//TRIM(arg2)//TRIM(arg3)
            arg4= "distribR"
            arg2=trim(arg4)//TRIM(arg2)//TRIM(arg3)
            call setParams(RRA, EIS, func, d1func, d2func,arg1,arg2 , printEvery, capitalShare,&
                &.FALSE., 0.1D0, 1.0D0)
            func => aggregateBondsFixedW
            if(rank ==0)then
                call CPU_TIME(startTime)
            end if
            intDiff=brent(func,0.01D0,0.09D0,0.12D0,1.0D-8,xmin)
            if(rank ==0)then
                call CPU_TIME(endTime)
                print *,"b             ",RRA,EIS,xmin,impliedCapital(xmin), intDiff, endTime-startTime
            end if
        end do
    end if

    if(whichSet/=1)then
        !***********************************
        !Question c
        !***********************************
        RRA=2.0D0
        do temp2=0,5
            EIS=0.1D0+(2.0D0-0.1D0)/5*temp2
            func => valueFunction
            d1func => d1prod
            d2func => d2prod
            arg1="policyR2E"
            write(arg2,'(I1)') temp2
            arg1=TRIM(arg1)//TRIM(arg2)
            arg4= "distribR2E"
            arg2=TRIM(arg4)//TRIM(arg2)
            call setParams(RRA, EIS, func, d1func, d2func,arg1,arg2 , printEvery, capitalShare,&
                &.FALSE., 0.1D0, 1.0D0)
            func => aggregateBondsFixedW
            if(rank ==0)then
                call CPU_TIME(startTime)
            end if
            intDiff=brent(func,0.01D0,0.09D0,0.12D0,1.0D-8,xmin)
            if(rank ==0)then
                call CPU_TIME(endTime)
                print *,"c             ",RRA,EIS,xmin,impliedCapital(xmin), intDiff, endTime-startTime
            end if
        end do
    endif
    CALL MPI_FINALIZE(ierr)


contains
    function production(capital,labour) RESULT(y)
        REAL(KIND=8), INTENT(IN) :: capital, labour
        REAL(KIND=8) :: y
        y=capital**capitalShare*labour**(1-capitalShare)
    end function production

    function d1prod(capital,labour) RESULT(y)
        REAL(KIND=8), INTENT(IN) :: capital, labour
        REAL(KIND=8) :: y
        y=capitalShare*capital**(capitalShare-1)*labour**(1-capitalShare)
    end function d1prod

    function d2prod(capital,labour) RESULT(y)
        REAL(KIND=8), INTENT(IN) :: capital, labour
        REAL(KIND=8) :: y
        y=(1-capitalShare)*capital**capitalShare*labour**(-capitalShare)
    end function d2prod

    function impliedCapital(interest) RESULT(y)
        REAL(KIND=8), INTENT(IN) :: interest
        REAL(KIND=8) :: y

        y=(interest/capitalShare)**(1.0D0/(capitalShare-1))
    end function impliedCapital
end program main
