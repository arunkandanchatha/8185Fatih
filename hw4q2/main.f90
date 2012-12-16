!#define DOPRINTS
module utilFuncs
    use nrtype
    use nrutil
contains
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
    use utilFuncs

    implicit none

    INTEGER :: n_s, n_z, currentState, currentAggState, currentCapital
    REAL(DP), allocatable, DIMENSION(:) :: a, s
    REAL(DP), allocatable, DIMENSION(:,:,:,:) :: transition
    REAL(DP), allocatable, DIMENSION(:,:,:) :: v, y2
    REAL(DP) :: beta,r,w,z

    PRIVATE n_s,n_z,currentState,currentCapital,a,s,v,y2,transition,beta,r,w
    PRIVATE currentAggState

contains

    subroutine wrapperCreate(pn_s,pn_z,grid,mytransition,mystates,mybeta,myr,myw)
        INTEGER(I4B), INTENT(IN) :: pn_s, pn_z
        REAL(DP), DIMENSION(:), INTENT(IN) ::grid
        REAL(DP), DIMENSION(:) :: mystates
        REAL(DP), DIMENSION(:,:,:,:) :: mytransition
        REAL(DP), INTENT(IN) :: mybeta,myr,myw

        n_s=pn_s
        n_z=pn_z

        allocate(s(n_s))
        s=mystates

        allocate(a(size(grid)))
        a=grid

        allocate(transition(n_z,n_s,n_z,n_s))
        transition=mytransition

        beta=mybeta
        r=myr
        w=myw
    end subroutine wrapperCreate

    subroutine wrapperInit(values,splineParams)
        REAL(DP), DIMENSION(:,:,:) :: values
        REAL(DP), DIMENSION(:,:,:) :: splineParams

        allocate(v(n_z,n_s,size(values,dim=2)))
        v=values

        allocate(y2(n_z,n_s,size(values,dim=2)))
        y2=splineParams
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

    function callBrent(econState, state,capital, myfunc, kr1, kr2, kr3, tol, minPoint) RESULT (y)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: econState, state, capital
        REAL(dp), INTENT(IN) :: kr1,kr2,kr3,tol
        REAL(dp), INTENT(OUT) :: minPoint
        PROCEDURE(template_function), POINTER :: myfunc
        REAL(dp) :: y

        currentState = state
        currentCapital = capital
        currentAggState = econState
        y=brent(myfunc, kr1, kr2, kr3, tol, minPoint)

    end function callBrent

    function valueFunction(x) RESULT(z)
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: x
        REAL(DP) :: z
        REAL(DP) ::  rp
        integer  ::  i,j
        REAL(DP) :: temp
        REAL(DP), dimension(n_z,n_s)            ::  interpolated

        do i=1,n_z
            do j=1,n_s
                interpolated(i,j)=splint(a,v(i,j,:),y2(i,j,:),x)
            end do
        end do
        rp=1.0D0
        do i=1,n_z
            do j=1,n_s
                rp=rp*transition(currentAggState,currentState,i,j)*interpolated(i,j)
            end do
        end do
        temp = a(currentCapital)*(1+r)+s(currentState)*w-x
        if(temp < 0.0D0) then
            temp=0.0D0
        end if
        z=-(log(temp)+beta*rp)
        if(z<-10D50)then
            print *,"value function too low"
            print *,temp,beta,rp
            stop 0
        end if
    end function valueFunction

end module brentWrapper


module aiyagariSolve
    USE nrtype
    USE brentWrapper
    implicit none

    REAL(DP), parameter                   ::  phi=.9D0,sigma=.4D0
    integer, parameter                  ::  n_a=301

    !******************
    ! These are here because of screwey splines
    !******************
    integer, parameter                  ::  bottomChop=0   !the number of first points we want to ignore
    integer, parameter                  ::  topChop=0       !the number of end points we want to ignore

    REAL(DP), parameter                   ::  curv=3.0D0,a_min=0.0D0, a_max=35.0D0
    REAL(DP), parameter                   ::  k_min=.01D0, k_max = a_max
    REAL(DP), parameter                   ::  beta=0.90
    integer, parameter                  ::  maxit=2000
    REAL(DP),parameter                    ::  toll=1D-8,tol2=1D-8

    integer                  ::  n_z   ! The number of aggregate states
    integer                  ::  n_s   ! The number of employment states
    integer                  ::  n_k   ! The number of capital levels

    REAL(DP), allocatable, dimension(:)              ::  s,stationary
    REAL(DP), allocatable, dimension(:)   ::  zShocks, ssEmployment
    REAL(DP), allocatable, dimension(:,:,:,:)  ::  transition
    REAL(DP), dimension(n_a)              ::  a
    REAL(DP), allocatable, dimension(:)   ::  k
    REAL(DP)                              :: capShare
    REAL(DP)                              :: wFixed=1.0D0,rFixed=0.1D0
    REAL(DP), allocatable, dimension(:,:,:,:,:):: v
    REAL(DP), allocatable, dimension(:,:,:,:)  :: g
    REAL(DP), allocatable, dimension(:,:,:,:)  ::  lastStateV
    integer                       ::  reportNum
    PROCEDURE(template_function), POINTER :: funcParam
    PROCEDURE(template_function2), POINTER :: deriv1Func, deriv2Func
    PROCEDURE(template_function3), POINTER :: capitalCalc
    logical                            :: firstCall = .true.

    character(LEN=20)                   :: policyOutput
    character(LEN=20)                   :: distribOutput

    PRIVATE phi, sigma, n_z, n_a, curv, a_max, beta, maxit, toll
    PRIVATE s, stationary, a, a_min,v,g
    PRIVATE reportNum,funcParam,policyOutput
    PRIVATE aggregateBonds, firstCall

contains
    subroutine setParams(d1Func, d2Func, capitalFunc,file1, file2, every, &
        & capitalShare, myR, myW)
        PROCEDURE(template_function2), POINTER, INTENT(IN) :: d1Func
        PROCEDURE(template_function2), POINTER, INTENT(IN) :: d2Func
        PROCEDURE(template_function3), POINTER, INTENT(IN) :: capitalFunc
        character(LEN=*),INTENT(IN) :: file1
        character(LEN=*),INTENT(IN) :: file2
        INTEGER, OPTIONAL, INTENT(IN) :: every
        REAL(DP), INTENT(IN) :: capitalShare
        REAL(DP), OPTIONAL, INTENT(IN) :: myR, myW

        capShare=capitalShare
        deriv1Func => d1Func
        deriv2Func => d2Func
        capitalCalc => capitalFunc

        reportNum = 50

        if(PRESENT(every)) then
            reportNum=every
        end if

        if(PRESENT(myW))then
            wFixed=myW
        end if

        if(PRESENT(myR))then
            rFixed=myR
        end if

        policyOutput = file1
        distribOutput = file2
        firstCall = .true.
    end subroutine setParams

    subroutine allocateArrays(mn_s,mn_z,mn_k)
        INTEGER(I4B), intent(IN) :: mn_s, mn_z,mn_k

        n_s=mn_s
        n_z=mn_z
        n_k=mn_k
        allocate(s(n_s))
        allocate(k(n_k))
        allocate(stationary(n_s))
        allocate(zShocks(n_z))
        allocate(ssEmployment(n_s))
        allocate(transition(n_z,n_s,n_z,n_s))
        allocate(v(n_z,n_k,n_s,n_a,2))
        allocate(g(n_z,n_k,n_s,n_a))
        allocate(lastStateV(n_z,n_k,n_s,n_a))
    end subroutine allocateArrays

    subroutine deallocateArrays()
        deallocate(s)
        deallocate(k)
        deallocate(stationary)
        deallocate(zShocks)
        deallocate(ssEmployment)
        deallocate(transition)
        deallocate(v)
        deallocate(g)
        deallocate(lastStateV)
    end subroutine deallocateArrays

    subroutine beginKrusellSmith()
        REAL(DP)  ::  incr
        REAL(DP),dimension(2,n_z) :: phi

        INTEGER :: i,j,ii
        REAL(DP) :: temp,xmin
        PROCEDURE(template_function), POINTER :: func

        !**************************************************************************
        ! We set up the grid of asset values based on the curvature, curv
        ! the minimum and maximum values in the grid a_min a_max
        !**************************************************************************
        incr=(a_max-a_min)/(n_a-1)
        a=(/ (    incr*real(i-1,8),i=1,n_a    ) /)
        a=a**curv
        a=a/(a_max-a_min)**(curv-1)+a_min

        !*************************************************************************
        ! First, solve the steady state version of the model with no aggregate shocks.
        !*************************************************************************
        call allocateArrays(1,2,1)
        funcParam => valueFunction
        func => aggregateBondsSetR
        ! we initialize the value function
        ! and we set up an initial guess for it
        v=0D0
        forall(i=1:n_z,j=1:n_k,ii=1:n_s) v(i,j,ii,:,1)=(a-a_min)**2
        g=0D0

        zShocks(1)=1.0D0
        ssEmployment(1)=0.93D0

        !***********************************************************************
        ! This is not the most efficient, but it is clearer to understand what is happening.
        ! The matrix is set up as follows:
        !          A B C D  = A B - current Agg State, current employment state
        !                     C D - next Agg state, next employment state
        !          AggState - 1: only one economy
        !          Employment - 1: unemployed
        !                       2: employed
        transition(1,1,1,1) = 0.07D0
        transition(1,1,1,2) = 0.93D0
        transition(1,2,1,1) = 0.07D0
        transition(1,2,1,2) = 0.93D0

        temp=brent(func,0.01D0,0.09D0,0.12D0,toll,xmin)
        call deallocateArrays()
        stop 0

        call allocateArrays(2,2,50)
        ! Set aggregate capital grid over which we want to evaluate, K
        incr=(k_max-k_min)/(n_k-1)
        k=(/ (    incr*real(i-1,8),i=1,n_k    ) /)
        k=k+k_min

        zShocks=(/0.99D0,1.01D0/)
        s=(0.0D0,1.0D0)
        ssEmployment=(/0.9D0,0.96D0/)

        !***********************************************************************
        ! Now iterate to find stuff.
        !***********************************************************************
        ! This is not the most efficient, but it is clearer to understand what is happening.
        ! These numbers are obtained from Den Haan Judd Juillard (2008) (according to David
        ! Wiczer, at least)
        ! The matrix is set up as follows:
        !          A B C D  = A B - current Agg State, current employment state
        !                     C D - next Agg state, next employment state
        !          AggState - 1: bad economy
        !                     2: good economy
        !          Employment - 1: unemployed
        !              2: employed
        transition(1,1,1,1) = 0.525D0
        transition(1,1,1,2) = 0.35D0
        transition(1,1,2,1) = 0.03125D0
        transition(1,1,2,2) = 0.09375D0
        transition(1,2,1,1) = 0.038889D0
        transition(1,2,1,2) = 0.836111D0
        transition(1,2,2,1) = 0.002083D0
        transition(1,2,2,2) = 0.122917D0
        transition(2,1,1,1) = 0.09375D0
        transition(2,1,1,2) = 0.03125D0
        transition(2,1,2,1) = 0.291667D0
        transition(2,1,2,2) = 0.583333D0
        transition(2,2,1,1) = 0.09375D0
        transition(2,2,1,2) = 0.115885D0
        transition(2,2,2,1) = 0.024306D0
        transition(2,2,2,2) = 0.8506941D0

        !************************************************************************
        ! initial guess for phi0 and phi1
        !************************************************************************
        !        phi(0,:)=log(...)
        phi(1,:)=0.0D0

        ! given K, find K'

        ! interpolate V0(a,K) on K', keeping assets constant.
        ! So now we have V0(a,K') for a grid of points, a

        ! use the bellman equation and interpolate on V0(a,K') (over a) to find a'
        ! and therefore V1(a,K)

        ! Compare V1(a,K) to V0(a,K)

        deallocate(transition)
        deallocate(zShocks)
        deallocate(ssEmployment)
    end subroutine beginKrusellSmith

    function aggregateBondsSetR(r) RESULT (z)
        ! inputs: r - the interest rate to test
        ! outputs: z - the aggregate level of borrowing
        REAL(DP), INTENT(IN) :: r
        REAL(DP) :: z
        REAL(DP):: totalCapital
        INTEGER(I4B) :: states

        do states=1,size(zShocks)
            k(1) = capitalCalc(r,ssEmployment(states),zShocks(states))
            wFixed=deriv2Func(k(1),ssEmployment(states),zShocks(states))
            totalCapital=aggregateBonds(r,wFixed,1)
            z=abs(totalCapital-k(1))
        end do
    end function aggregateBondsSetR

    function aggregateBonds(r,w,currentCap) RESULT (aggK)
        ! inputs: r - the interest rate to test
        !         w - the wages
        !         z - the current aggregate shock level
        !         currentCap - the current aggregate capital level
        ! outputs: aggK - the aggregate level of borrowing
        REAL(DP), INTENT(IN) :: r,w
        INTEGER(I4B), INTENT(IN) :: currentCap
        REAL(DP) :: aggK
        INTEGER :: iterCount, i,j
        REAL(DP)  ::  incr, totalCapital
        REAL(DP), dimension(sizeof(ssEmployment),sizeof(zShocks),n_a-bottomChop-topChop) ::  steadyStateCapital

        rFixed=r
        wFixed=w
        call wrapperCreate(n_s,n_z,a,transition,s,beta,r,w)
        call getPolicyForCapital(currentCap)
        call wrapperDestroy()

        open(unit=1,file=policyOutput)
        write(1,*) a(:)
        do i=1,n_z
            write(1,*) g(i,aggK,:,:)
        end do
        close(1)

            !*****************************************************
            ! find the steady state capital
            !*****************************************************
        call findSteadyState(g(1,:,:),steadyStateCapital)

        !**************************************************
        ! now, calculate total capital
        !**************************************************
        totalCapital = dot_product(sum(steadyStateCapital,dim=1),a(bottomChop+1:n_a-topChop))
        aggK = totalCapital

    end function aggregateBonds


    subroutine getPolicyForCapital(aggK)
        INTEGER, INTENT(IN)                 :: aggK
        INTEGER                             :: i,it,j,iter
        REAL(DP)                            :: tempD,tempD2
        REAL(DP), dimension(n_z,n_s,n_a)    :: y2
        REAL(DP)                            :: kr1,kr2,kr3

        !**************************************************************************
        ! we begin the iteration
        !**************************************************************************
        do iter=2,maxit
            do i=1,n_z
                do j=1,n_s
                    tempD=(v(i,aggK,j,2,1)-v(i,aggK,j,1,1))/(a(2)-a(1))
                    tempD2=(v(i,aggK,j,n_a,1)-v(i,aggK,j,n_a-1,1))/(a(n_a)-a(n_a-1))
                    call spline(a,v(i,aggK,j,:,1),tempD,tempD2,y2(i,j,:))
                end do
            end do
            call wrapperInit(v(:,aggK,:,:,1),y2)
            do i=1,n_z
                do j=1,n_s
                    do it=1,n_a
                        !ensure monotone policy function
                        if(it==1)then
                            kr1=a(1)
                        else
                            kr1=g(i,aggK,j,it-1)
                        end if
                        kr3=min(a(it)*(1+rFixed)+s(j)*wFixed,a(n_a))
                        kr2=(kr1+kr3)/2D0
                        v(i,aggK,j,it,2)=-callBrent(i,j,it,funcParam,kr1,kr2,kr3,toll,g(i,aggK,j,it))
                        if(isNAN(v(i,aggK,j,it,2)))then
                            print *,"error, value function is nan"
                            stop 0
                        end if
                    end do
                end do
            end do
            call wrapperClean()
            if( (mod(iter,reportNum)==0))then
                print *,"r:",rFixed,"w:",wFixed,"iter: ",iter,&
                    &"diff: ",maxval(abs(v(:,aggK,:,:,2)-v(:,aggK,:,:,1)))
                flush(6)
            end if
            if (maxval(abs(v(:,aggK,:,:,2)-v(:,aggK,:,:,1))) .lt. toll) then
#ifdef DOPRINTS
                print*,"done: ",iter, "r: ",r, "w: ",wFixed
                flush(6)
#endif
                exit
            end if

            v(:,aggK,:,:,1)=v(:,aggK,:,:,2)
        end do

        lastStateV(:,aggK,:,:) = v(:,aggK,:,:,2)

    end subroutine  getPolicyForCapital

    subroutine findSteadyState(capitalPolicyOrig, statDist)
        !INPUTS: capitalPolicy - the policy function for each state
        !OUTPUTS: statDist - the stationary dist (note: pdf, not cdf)
        INTEGER, parameter :: capitalCount = n_a-bottomChop-topChop
        REAL(DP), dimension(n_z,n_a), INTENT(IN) :: capitalPolicyOrig
        real(DP), dimension(n_z,capitalCount), intent(out) :: statDist
        real(DP), dimension(n_z,capitalCount) ::f_o, f_o_hat, f_n
        real(DP) :: diff, temp
        INTEGER :: i,j, counter
        REAL(DP), dimension(n_z,capitalCount) :: capitalPolicy
        REAL(DP), dimension(capitalCount) :: newCapital


        newCapital = a(bottomChop+1:n_a-topChop)
        capitalPolicy = capitalPolicyOrig(:,bottomChop+1:n_a-topChop)

        do i=1,n_z
            do j=2,capitalCount
                if(capitalPolicy(i,j)<capitalPolicy(i,j-1)) then
                    print *,"not monotonic.", i, j
                    print *,capitalPolicy(i,j-1),capitalPolicy(i,j)
                    stop 0
                end if
            end do
        end do
        !setting initial guess to uniform dist across asset grid. First we have it
        !a cdf, and then convert to the appropriate pdf
        do i=1,n_z
            do j=1,capitalCount
                statDist(i,j) = (newCapital(j) - newCapital(1))/(newCapital(capitalCount)-newCapital(1))*dble(1.0_dp/n_z)
            end do
        end do

        f_n=statDist

        ! time to iterate
        diff=100
        counter = 0
        do while((diff>tol2) .and. (counter<maxit))
            counter = counter + 1

            f_o=f_n

            do i=1,n_z
                do j=1,capitalCount
                    f_o_hat(i,j) = linear(f_o(i,:),capitalPolicy(i,:),newCapital(j))
                end do
            end do

            ! need to make sure monotonic and bounded between 0 and 1
            do i=1,n_z
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
            end do

! THIS IS INCORRECT            f_n = matmul(transpose(transition), f_o_hat)

            !* Fix so that total cdf is 1
            do i=1,n_z
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

        do i=1,n_z
            f_n(i,:)=f_n(i,:)*stationary(i)
        end do

        statDist = f_n
        do i=capitalCount,2,-1
            statDist(:,i) = statDist(:,i) - statDist(:,i-1)
        end do
        do i=1,n_z
            statDist(i,1) = max(0.0D0,stationary(i)-sum(statDist(i,:)))
        end do

        open(unit=1,file=distribOutput)
        do j=1,capitalCount
            write(1,*) a(j+bottomChop),statDist(:,j)
        end do
        close(1)

    end subroutine findSteadyState

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

    d1func => d1prod
    d2func => d2prod

    arg1="policy"
    arg2= "distrib"
    call setParams(d1func, d2func, impliedCapital, arg1,arg2 , printEvery, capitalShare,&
        &0.1D0, 1.0D0)
    if(rank ==0)then
        call CPU_TIME(startTime)
    end if
    call beginKrusellSmith()
    if(rank ==0)then
        call CPU_TIME(endTime)
        print *,"a             ",RRA,EIS,xmin,impliedCapital(xmin), intDiff, endTime-startTime
        flush(6)
    end if

    CALL MPI_FINALIZE(ierr)


contains
    function production(capital,labour,shock) RESULT(y)
        REAL(DP), INTENT(IN) :: capital, labour,shock
        REAL(DP) :: y
        y=shock*capital**capitalShare*labour**(1-capitalShare)
    end function production

    function d1prod(capital,labour,shock) RESULT(y)
        REAL(DP), INTENT(IN) :: capital, labour,shock
        REAL(DP) :: y
        y=capitalShare*shock*capital**(capitalShare-1)*labour**(1-capitalShare)
    end function d1prod

    function d2prod(capital,labour,shock) RESULT(y)
        REAL(DP), INTENT(IN) :: capital, labour,shock
        REAL(DP) :: y
        y=(1-capitalShare)*shock*capital**capitalShare*labour**(-capitalShare)
    end function d2prod

    function impliedCapital(interest,labour,shock) RESULT(y)
        REAL(DP), INTENT(IN) :: interest, labour, shock
        REAL(DP) :: y
        y = interest/(capitalShare*shock)*&
            1.0D0/labour**(1-capShare)
        y = y**(1.0D0/(capitalShare-1))
    end function impliedCapital

end program main
