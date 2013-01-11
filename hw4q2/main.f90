!#define DOPRINTS
module utilFuncs
    use nrtype
    use nrutil
contains
    ! given arrays x and y of length N containing a tabulated function
    ! given values yp1 and ypn for the first derivative of the interpolating function at points 1 and N
    ! this routine returns an array y2 of length N that contains the second derivatives
    ! of the interpolating function at the tabulated points xi.
    ! If yp1 and/or ypn are equal to 1 � 1030 or larger, the routine is signaled to set the
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

    !Given the arrays xa and ya, which tabulate a function (with the xai �s in increasing or
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

    INTEGER :: n_s, n_z, n_k,n_a,currentState, currentAggState, currentCapital, currentAggCapital
    REAL(DP), allocatable, DIMENSION(:) :: a, s
    REAL(DP), allocatable, DIMENSION(:,:) :: r, w
    REAL(DP), allocatable, DIMENSION(:,:,:,:) :: transition
    REAL(DP), allocatable, DIMENSION(:,:,:,:) :: v, y2
    REAL(DP) :: beta,z,delta

    PRIVATE n_s,n_z,n_k,n_a,currentState,currentCapital,a,s,v,y2,transition,beta,r,w,delta
    PRIVATE currentAggState, currentAggCapital

contains

    subroutine wrapperCreate(pn_s,pn_z,pn_k,grid,mytransition,mystates,mybeta,myr,myw,mydelta)
        INTEGER(I4B), INTENT(IN) :: pn_s, pn_z,pn_k
        REAL(DP), DIMENSION(:), INTENT(IN) ::grid
        REAL(DP), DIMENSION(:), INTENT(IN) :: mystates
        REAL(DP), DIMENSION(:,:,:,:), INTENT(IN) :: mytransition
        REAL(DP), DIMENSION(:,:), INTENT(IN) :: myr,myw
        REAL(DP), INTENT(IN) :: mybeta,mydelta

        n_s=pn_s
        n_z=pn_z
        n_k=pn_k
        n_a=size(grid)
        allocate(s(n_s))
        s=mystates
        allocate(a(n_a))
        a=grid
        allocate(transition(n_z,n_s,n_z,n_s))
        transition=mytransition
        allocate(r(n_z,n_k))
        r=myr
        allocate(w(n_z,n_k))
        w=myw
        beta=mybeta
        delta=mydelta
    end subroutine wrapperCreate

    subroutine wrapperInit(values,splineParams)
        REAL(DP), DIMENSION(:,:,:,:) :: values
        REAL(DP), DIMENSION(:,:,:,:) :: splineParams

        allocate(v(n_z,n_k,n_s,n_a))
        v=values

        allocate(y2(n_z,n_k,n_s,n_a))
        y2=splineParams
    end subroutine wrapperInit

    subroutine wrapperClean()
        deallocate(v)
        deallocate(y2)
    end subroutine wrapperClean

    subroutine wrapperDestroy()
        deallocate(a)
        deallocate(transition)
        deallocate(s)
        deallocate(r)
        deallocate(w)
    end subroutine wrapperDestroy

    function callBrent(econState, aggCap, state,capital, myfunc, kr1, kr2, kr3, tol, minPoint) RESULT (y)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: econState, aggCap, state, capital
        REAL(dp), INTENT(IN) :: kr1,kr2,kr3,tol
        REAL(dp), INTENT(OUT) :: minPoint
        PROCEDURE(template_function), POINTER :: myfunc
        REAL(dp) :: y

        currentState = state
        currentCapital = capital
        currentAggState = econState
        currentAggCapital = aggCap

        y=brent(myfunc, kr1, kr2, kr3, tol, minPoint)

    end function callBrent

    function valueFunction(x) RESULT(z)
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: x
        REAL(DP) :: z
        REAL(DP) ::  rp
        integer  ::  i,j
        REAL(DP) :: temp,eps=epsilon(1.0D0)
        REAL(DP), dimension(n_z,n_s)            ::  interpolated

        do i=1,n_z
            do j=1,n_s
                interpolated(i,j)=splint(a,v(i,currentAggCapital,j,:),y2(i,currentAggCapital,j,:),x)

                if(isNaN(interpolated(i,j)))then
                    print *,"interpolated value in value function is NaN"
                    print *,i,j,x
                    print *,a
                    print *,v(i,currentAggCapital,j,:)
                    print *,y2(i,currentAggCapital,j,:)
                    stop 0
                end if

                if(interpolated(i,j)>10D50)then
                    print *,"interpolated value in value function is too high"
                    print *,i,j,x
                    print *,a
                    print *,v(i,currentAggCapital,j,:)
                    print *,y2(i,currentAggCapital,j,:)
                    stop 0
                end if

                if(interpolated(i,j)<-10D50)then
                    print *,"interpolated value in value function is too low"
                    print *,i,j,x
                    print *,a
                    print *,v(i,currentAggCapital,j,:)
                    print *,y2(i,currentAggCapital,j,:)
                    stop 0
                end if

            end do
        end do
        rp=0.0D0
        do i=1,n_z
            do j=1,n_s
                rp=rp+transition(currentAggState,currentState,i,j)*interpolated(i,j)
            end do
        end do
        temp = a(currentCapital)*(1+r(currentAggState,currentAggCapital)-delta)+&
            s(currentState)*w(currentAggState,currentAggCapital)-x
        if(temp < eps) then
            temp=eps
        end if
        z=-(log(temp)+beta*rp)
        if(z<-10D50)then
            print *,"value function too low"
            print *,temp,beta,rp
            stop 0
        end if
        if(z>10D50)then
            print *,"value function too high"
            print *,temp,beta,rp
            stop 0
        end if
        if(isNaN(z))then
            print *, "value function is NaN"
            print *,temp,beta,rp
            stop 0
        end if
    end function valueFunction

end module brentWrapper

module aiyagariSolve
    USE nrtype
    USE brentWrapper
    implicit none

    INCLUDE 'mpif.h'

    REAL(DP), parameter                   ::  sigma=.4D0
    integer, parameter                  ::  n_a=301
    integer                             :: myseed = 4567
    integer, parameter                  :: periodsForConv = 10000
    integer, parameter                  :: periodsToCut = 1000
    integer, parameter                  :: numHouseholds = 50000

    !******************
    ! These are here because of screwey splines
    !******************
    integer, parameter                  ::  bottomChop=0   !the number of first points we want to ignore
    integer, parameter                  ::  topChop=0       !the number of end points we want to ignore

    REAL(DP), parameter                   ::  curv=3.0D0,a_min=0.01D0/numHouseholds, a_max=150.0D0
    REAL(DP), parameter                   ::  k_min=2.0D0, k_max = 200.0D0
    REAL(DP), parameter                   ::  beta=0.90, delta = 0.025D0
    integer, parameter                  ::  maxit=2000
    REAL(DP),parameter                    ::  toll=1D-8,tol2=1D-5
    REAL(DP),parameter                    ::  lambda = 0.1 ! how much confidence in new value

    integer                  ::  n_z   ! The number of aggregate states
    integer                  ::  n_s   ! The number of employment states
    integer                  ::  n_k   ! The number of capital levels

    REAL(DP), allocatable, dimension(:)              ::  s,stationary,zShocks,k
    REAL(DP), allocatable, dimension(:,:)   ::  ssEmployment, zTransition
    REAL(DP), allocatable, dimension(:,:,:,:)  ::  transition
    REAL(DP), dimension(n_a)              ::  a
    REAL(DP)                              :: capShare
    REAL(DP), allocatable, dimension(:,:):: wFixed,rFixed
    REAL(DP), allocatable, dimension(:,:,:,:,:):: v
    REAL(DP), allocatable, dimension(:,:,:,:)  :: g
    REAL(DP), allocatable, dimension(:,:,:,:)  ::  lastStateV
    TYPE(household), dimension(numHouseholds)  ::  ssDistrib
    REAL(DP), DIMENSION(periodsForConv) :: lastKSeq
    REAL(DP),allocatable,dimension(:,:,:) :: phi


    integer                       ::  reportNum
    PROCEDURE(template_function), POINTER :: funcParam
    PROCEDURE(template_function2), POINTER :: deriv1Func, deriv2Func
    PROCEDURE(template_function3), POINTER :: capitalCalc
    logical                            :: firstCall = .true.

    character(LEN=20)                   :: policyOutput
    character(LEN=20)                   :: distribOutput

    PRIVATE phi, sigma, n_z, n_a, curv, a_max, beta, maxit, toll, delta
    PRIVATE s, stationary, a, a_min,v,g
    PRIVATE reportNum,funcParam,policyOutput
    PRIVATE aggregateBonds, firstCall

contains
    subroutine setParams(d1Func, d2Func, capitalFunc,file1, file2, capitalShare, seedParam, every)
        PROCEDURE(template_function2), POINTER, INTENT(IN) :: d1Func
        PROCEDURE(template_function2), POINTER, INTENT(IN) :: d2Func
        PROCEDURE(template_function3), POINTER, INTENT(IN) :: capitalFunc
        character(LEN=*),INTENT(IN) :: file1
        character(LEN=*),INTENT(IN) :: file2
        REAL(DP), INTENT(IN) :: capitalShare
        INTEGER, OPTIONAL, INTENT(IN) :: seedParam, every

        capShare=capitalShare
        deriv1Func => d1Func
        deriv2Func => d2Func
        capitalCalc => capitalFunc

        reportNum = 500

        if(PRESENT(seedParam)) then
            myseed=seedParam
        end if

        if(PRESENT(every)) then
            reportNum=every
        end if

        policyOutput = file1
        distribOutput = file2
        firstCall = .true.
    end subroutine setParams

    subroutine allocateArrays(mn_z,mn_s,mn_k)
        INTEGER(I4B), intent(IN) :: mn_s, mn_z,mn_k

        n_s=mn_s
        n_z=mn_z
        n_k=mn_k
        allocate(s(n_s))
        allocate(k(n_k))
        allocate(stationary(n_s))
        allocate(zShocks(n_z))
        allocate(ssEmployment(n_z,n_s))
        allocate(transition(n_z,n_s,n_z,n_s))
        allocate(ztransition(n_z,n_z))
        allocate(v(n_z,n_k,n_s,n_a,2))
        allocate(g(n_z,n_k,n_s,n_a))
        allocate(lastStateV(n_z,n_k,n_s,n_a))
        allocate(wFixed(n_z,n_k))
        allocate(rFixed(n_z,n_k))
        allocate(phi(n_z,2,2))
    end subroutine allocateArrays

    subroutine deallocateArrays()
        deallocate(s)
        deallocate(k)
        deallocate(stationary)
        deallocate(zShocks)
        deallocate(ssEmployment)
        deallocate(transition)
        deallocate(ztransition)
        deallocate(v)
        deallocate(g)
        deallocate(lastStateV)
        deallocate(wFixed)
        deallocate(rFixed)
        deallocate(phi)

    end subroutine deallocateArrays

    subroutine beginKrusellSmith(doSS)
        LOGICAL :: doSS
        REAL(DP)  ::  incr, ssErr, ssTot, predicted
        REAL(DP),allocatable,dimension(:,:) :: vals,goodShocks,badShocks
        REAL(DP),allocatable,dimension(:) ::  aggK,aggKp
        REAL(DP), dimension(periodsForConv,2):: aggKHistory
        REAL(DP),dimension(2*periodsForConv*2) :: workArray
        LOGICAL :: iterComplete
        INTEGER :: i,j,ii,iter, whichState
        REAL(DP) :: temp,xmin, avgK
        PROCEDURE(template_function), POINTER :: func
        !************
        ! MPI vars
        !************
        INTEGER rank, ierr, mysize
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mysize, ierr)

        funcParam => valueFunction

        !**************************************************************************
        ! We set up the grid of asset values based on the curvature, curv
        ! the minimum and maximum values in the grid a_min a_max
        !**************************************************************************
        incr=(a_max-a_min)/(n_a-1)
        a=(/ (    incr*real(i-1,8),i=1,n_a    ) /)
        a=a**curv
        a=a/(a_max-a_min)**(curv-1)+a_min

        if(doSS)then
            !*************************************************************************
            ! First, solve the steady state version of the model with no aggregate shocks.
            !*************************************************************************
            call allocateArrays(1,2,1)
            s(1)=0.0D0
            s(2)=1.0D0
            func => aggregateBondsSetR

            ! we initialize the value function
            ! and we set up an initial guess for it
            allocate(vals(n_z,2))
            allocate(aggK(n_z))
            allocate(aggKp(n_z))
            v=0D0
            forall(i=1:n_z,j=1:n_k,ii=1:n_s) v(i,j,ii,:,1)=(a-a_min)**2
            g=0D0

            zShocks(1)=1.0D0
            ssEmployment(1,1)=0.07D0
            ssEmployment(1,2)=0.93D0

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
            ztransition(1,1) = 1.0D0

            temp=brent(func,0.05D0,0.13D0,0.20D0,tol2,xmin)

            !*************
            !Let's print out the ss distrib
            !*************
            open(unit=1,file=distribOutput)
            do i=1,numHouseholds
                write(1,*) ssDistrib(i)%employmentState
                write(1,*) ssDistrib(i)%capital
            end do
            close(1)

            call deallocateArrays()
            deallocate(vals)
            deallocate(aggK)
            deallocate(aggKp)
        else
            open(unit=1,file=distribOutput)
            do i=1,numHouseholds
                read(1,*) ssDistrib(i)%employmentState
                read(1,*) ssDistrib(i)%capital
            end do
        end if

        !**********************************************************************
        ! Now do proper Krusell-Smith algorithm
        !********************************************************************

        call allocateArrays(2,2,100)
        allocate(vals(n_z,2))
        allocate(aggK(n_z))
        allocate(aggKp(n_z))

        ! Set aggregate capital grid over which we want to evaluate, K
        ! curving this allows convergence. Without, it doesn't. Damn I hate splines!
        incr=(k_max-k_min)/(n_k-1)
        k=(/ (    incr*real(i-1,8),i=1,n_k    ) /)
        k=k**curv
        k=k/(k_max-k_min)**(curv-1)+k_min

        zShocks=(/0.99D0,1.01D0/)
        s(1)=0.0D0
        s(2)=1.0D0
        ssEmployment=reshape((/0.1D0,0.04D0,0.90D0,0.96D0/),(/2,2/))

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

        zTransition(1,1) = transition(1,1,1,1)+transition(1,1,1,2)
        zTransition(1,2) = 1-ztransition(1,1)
        zTransition(2,2) = transition(2,2,2,1)+transition(2,2,2,2)
        zTransition(2,1) = 1-ztransition(2,2)

        !************************************************************************
        ! initial guess for phi0 and phi1
        !************************************************************************
        aggK(:)=sum(ssDistrib(:)%capital)/numHouseholds
        !cheating, using pre-calculated values from Maliar
        phi(1,1,1)=0.15107110616449188D0
        phi(1,2,1)=0.93731513060912675D0
        phi(2,1,1)=0.10818954905023634D0
        phi(2,2,1)=0.93747819452287962D0
        !phi(:,1,1)=log(aggK(1))
        !phi(:,2,1)=0.1D0
        vals(1,:)=1.0D0
        vals(2,:)=log(aggK(1))

        ! given K, find K'
        do i=1,n_z
            aggKp(i)=exp(dot_product(phi(i,:,1),vals(:,i)))
        end do

        do i=1,n_z
            do j=1,n_k
                wFixed(i,j)=deriv2Func(k(j),ssEmployment(i,2),zShocks(i))
                rFixed(i,j)=deriv1Func(k(j),ssEmployment(i,2),zShocks(i))
            end do
        end do

        iter = 0
        iterComplete=.false.
        do while ( (iter<maxit) .and. (.not. iterComplete))
            iter = iter+1

            ! initialize the value function
            ! and we set up an initial guess for it
            ! do this each time, since previous value functions are probably
            ! garbage
            v=0D0
            forall(i=1:n_z,j=1:n_k,ii=1:n_s) v(i,j,ii,:,1)=(a-a_min)**2
            g=0D0

            call wrapperCreate(n_s,n_z,n_k,a,transition,s,beta,rFixed,wFixed, delta)
            call getAllPolicy()
            call wrapperDestroy()

            !*****************************************************
            ! find the steady state capital
            !*****************************************************
            aggKHistory=findSteadyStateCapital(.false.)

            !****************************************************
            ! Split the data into two sets, one for good shocks and
            ! one for bad shocks
            !****************************************************
            j=0
            do i=1,periodsForConv
                if(aggKHistory(i,1)>1.0D0)then
                    j=j+1
                end if
            end do

            allocate(goodShocks(j,3))
            allocate(badShocks(periodsForConv-j,3))

            j=0
            ii=0
            do i=1,periodsForConv
                if(aggKHistory(i,1)>1.0D0)then
                    j=j+1
                    goodShocks(j,1)=1.0D0
                    if(i==1)then
                        goodShocks(j,2)=log(aggK(1))
                    else
                        goodShocks(j,2)=log(aggKHistory(i-1,2))
                    end if
                    goodShocks(j,3)=log(aggKHistory(i,2))
                else
                    ii=ii+1
                    badShocks(ii,1)=1.0D0
                    if(i==1)then
                        badShocks(ii,2)=log(aggK(1))
                    else
                        badShocks(ii,2)=log(aggKHistory(i-1,2))
                    end if
                    badShocks(ii,3)=log(aggKHistory(i,2))
                end if
            end do

            !**************************************************
            ! Find phi(0) and phi(1) using OLS
            !**************************************************

            !First for good shocks
            j=size(goodShocks,dim=1)
            call dgels('N', j, 2, 1, goodShocks(:,1:2), j, goodShocks(:,3),j, workArray,size(workArray),i)

            phi(1,1,2)=goodShocks(1,3)
            phi(1,2,2)=goodShocks(2,3)
            j=size(badShocks,dim=1)
            call dgels('N', j, 2, 1, badShocks(:,1:2), j, badShocks(:,3),j, workArray,size(workArray),i)
            phi(2,1,2)=badShocks(1,3)
            phi(2,2,2)=badShocks(2,3)

            if( (mod(iter,1)==0) .and. (rank==0))then
                print *,"KS:" ,iter,maxval(abs(phi(:,:,2)-phi(:,:,1)))
                print *,"OLD:"
                print *,"1:",phi(1,:,1)
                print *,"2:",phi(2,:,1)
                print *,"NEW:"
                print *,"1:",phi(1,:,2)
                print *,"2:",phi(2,:,2)
                flush(6)
            end if

            if(maxval(abs(phi(:,:,2)-phi(:,:,1)))<1D0-6)then
                iterComplete = .true.
            end if

            phi(:,:,1)=lambda*phi(:,:,2)+(1.0D0-lambda)*phi(:,:,1)

            deallocate(goodShocks)
            deallocate(badShocks)

            ssErr=0.0D0
            ssTot=0.0D0
            avgK = sum(aggKHistory(:,2))/periodsForConv
            do i=1,periodsForConv
                ssTot=ssTot+(aggKHistory(i,2)-avgK)**2
                whichState=floor(aggKHistory(i,1))
                if(i==1)then
                    predicted = exp(phi(whichState,1,1)+phi(whichState,2,1)*log(aggK(1)))
                else
                    predicted = exp(phi(whichState,1,1)+phi(whichState,2,1)*log(aggKHistory(i-1,2)))
                end if
                ssErr=ssErr+(aggKHistory(i,2)-predicted)**2
            end do

            if(rank == 0)then
                print *,"R-squared: ",1.0D0-ssErr/ssTot
            end if

        end do

        ssErr=0.0D0
        ssTot=0.0D0
        avgK = sum(aggKHistory(:,2))/periodsForConv
        do i=1,periodsForConv
            ssTot=ssTot+(aggKHistory(i,2)-avgK)**2
            whichState=floor(aggKHistory(i,1))
            if(i==1)then
                predicted = exp(phi(whichState,1,1)+phi(whichState,2,1)*log(aggK(1)))
            else
                predicted = exp(phi(whichState,1,1)+phi(whichState,2,1)*log(aggKHistory(i-1,2)))
            end if
            ssErr=ssErr+(aggKHistory(i,2)-predicted)**2
        end do

        if(rank == 0)then
            print *,"R-squared: ",1.0D0-ssErr/ssTot
        end if
        deallocate(vals)
        deallocate(aggK)
        deallocate(aggKp)
        deallocate(goodShocks)
        deallocate(badShocks)

    end subroutine beginKrusellSmith

    function aggregateBondsSetR(r) RESULT (z)
        !*************************************
        !Only call this when we have no aggregate shocks
        !*************************************
        ! inputs: r - the interest rate to test
        ! outputs: z - the aggregate level of borrowing
        REAL(DP), INTENT(IN) :: r
        REAL(DP) :: z
        REAL(DP):: totalCapital
        REAL(DP),dimension(periodsForConv) :: aggK
        !************
        ! MPI vars
        !************
        INTEGER rank, ierr, mysize
!        INTEGER,dimension(MPI_STATUS_SIZE):: stat
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mysize, ierr)

        k(1) = capitalCalc(r,ssEmployment(1,2),zShocks(1))
        aggK=aggregateBonds(.true.)
        totalCapital=aggK(periodsForConv)
        z=abs(totalCapital-k(1))

        if(rank == 0)then
            print *,"Implied Capital: ",k(1),"Actual Capital: ",totalCapital
            flush(6)
        end if
    end function aggregateBondsSetR

    function aggregateBonds(noAggShocks) RESULT (aggK)
        ! inputs: r - the interest rate to test
        !         w - the wages
        !         z - the current aggregate shock level
        !         currentCap - the current aggregate capital level
        ! outputs: aggK - the aggregate level of borrowing
        LOGICAL, INTENT(IN) :: noAggShocks
        REAL(DP),dimension(periodsForConv) :: aggK
        REAL(DP),dimension(periodsForConv,2) :: aggKtemp
        INTEGER :: i,j
        !************
        ! MPI vars
        !************
        INTEGER rank, ierr, mysize
!        INTEGER,dimension(MPI_STATUS_SIZE):: stat
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mysize, ierr)

        do i=1,n_z
            do j=1,n_k
                wFixed(i,j)=deriv2Func(k(j),ssEmployment(i,2),zShocks(i))
                rFixed(i,j)=deriv1Func(k(j),ssEmployment(i,2),zShocks(i))
            end do
        end do

        call wrapperCreate(n_s,n_z,n_k,a,transition,s,beta,rFixed,wFixed, delta)
        call getAllPolicy(noAggShocks)
        call wrapperDestroy()

        if(rank == 0)then
            open(unit=1,file=policyOutput)
            write(1,*) a(:)
            do j=1,n_s
                do i=1,n_z
                    write(1,*) g(i,1,j,:)
                end do
            end do
            close(1)
        end if

        !*****************************************************
        ! find the steady state capital
        !*****************************************************
        aggKTemp=findSteadyStateCapital(noAggShocks)
        aggK=aggKTemp(:,2)
    end function aggregateBonds

    subroutine getAllPolicy(inSS)
        LOGICAL, OPTIONAL, INTENT(IN)       :: inSS
        LOGICAL                             :: inSS2

        INTEGER                             :: i,it,j,iter,aggK,ii
        REAL(DP)                            :: tempD,tempD2
        REAL(DP), dimension(n_z,n_k,n_s,n_a)    :: y2
        REAL(DP)                            :: kr1,kr2,kr3
        real(DP),dimension(n_z,n_k,n_s,n_a) :: tempV,tempY,tempG, tempV2
        REAL(DP),dimension(n_z,n_k) :: kprime
        LOGICAL :: exitLoop

        !************
        ! MPI vars
        !************
        INTEGER rank, ierr, mysize
!        INTEGER,dimension(MPI_STATUS_SIZE):: stat

        !************
        ! Timing variables
        !************
        real(DP) :: startTime, endTime

        CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mysize, ierr)

        if(present(inSS))then
            inSS2=inSS
        else
            inSS2=.false.
        end if

        if(rank==0)then
            call CPU_TIME(startTime)
            if(reportNum<maxit)then
                if(inSS2)then
                    print *," "
                    print *,"Value Function iteration",rFixed(1,1)
                else
                    print *," "
                    print *,"Value Function iteration"
                end if
                print *,"--------------------------------------"
                flush(6)
            end if
        end if

        if(inSS2)then
            do i=1,n_z
                kprime(i,:)=k
            end do
        else
            do j=1,n_k
                do i=1,n_z
                    kprime(i,j)=max(k(1),exp(phi(i,1,1)+phi(i,2,1)*log(k(j))))
                    kprime(i,j)=min(k(n_k),kprime(i,j))
                end do
            end do
        end if

        !**************************************************************************
        ! we begin the iteration
        !**************************************************************************
        exitLoop = .false.
        do iter=1,maxit
            if(exitLoop)then
                exit
            end if

            if (.not. inSS2)then
                !Get splines for each a, across all k
                do ii = 1,n_a
                    do j=1,n_s
                        do i=1,n_z
                            tempD=(v(i,2,j,ii,1)-v(i,1,j,ii,1))/(k(2)-k(1))
                            tempD2=(v(i,n_k,j,ii,1)-v(i,n_k-1,j,ii,1))/(k(n_k)-k(n_k-1))
                            call spline(k,v(i,:,j,ii,1),tempD,tempD2,y2(i,:,j,ii))
                        end do
                    end do
                end do
            end if

            !interpolate values at K'. This is what we need when we evaluate policy functions
            do it=1,n_a
                do j=1,n_s
                    do aggK=1,n_k
                        do i=1,n_z
                            if(.not. inSS2)then
                                tempV(i,aggK,j,it)=splint(k,v(i,:,j,it,1),y2(i,:,j,it),kprime(i,aggK))
                            else
                                tempV(i,aggK,j,it)=v(i,aggK,j,it,1)
                            end if
                        end do
                    end do
                end do
            end do

            !now, find a  policy function that is valid across each k'
            do aggK=1,n_k
                do j=1,n_s
                    do i=1,n_z
                        tempD = (tempV(i,aggK,j,2)-tempV(i,aggK,j,1))/(a(2)-a(1))
                        tempD2 = (tempV(i,aggK,j,n_a)-tempV(i,aggK,j,n_a-1))/(a(n_a)-a(n_a-1))
                        call spline(a,tempV(i,aggK,j,:),tempD,tempD2,tempY(i,aggK,j,:))
                    end do
                end do
            end do

            !use this single policy function for each k' (actually one for each future agg state and emloyment
            !level (so really four in basic K-S) for evaluating next period value function
            call wrapperInit(tempV,tempY)

            tempG=0
            tempV2=0
            !find next period's policy function
            do it=rank+1,n_a,mysize
                do j=1,n_s
                    do aggK = 1,n_k
                        do i=1,n_z
                            kr1=a(1)
                            kr3=min(a(it)*(1+rFixed(i,aggK)-delta)+s(j)*wFixed(i,aggK),a(n_a))
                            kr2=(kr1+kr3)/2D0
                            tempV2(i,aggK,j,it)=-callBrent(i,aggK,j,it,funcParam,kr1,kr2,kr3,toll,tempG(i,aggK,j,it))
                            if(isNAN(tempV2(i,aggK,j,it)))then
                                print *,"error, value function is nan"
                                print *,i,j,it
                                print *,kr1,kr2,kr3
                                flush(6)
                                stop 0
                            end if
                            if(tempV2(i,aggK,j,it)<-10D50)then
                                print *,"error, value function is too low"
                                print *,i,j,it
                                print *,kr1,kr2,kr3
                                print *,v(i,aggK,j,it,1),tempV2(i,aggK,j,it)
                                print *,s(j)
                                flush(6)
                                stop 0
                            end if
                        end do
                    end do
                end do
            end do

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            call MPI_ALLREDUCE(tempV2,v(:,:,:,:,2),n_s*n_k*n_z*n_a,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(tempG,g,n_s*n_k*n_z*n_a,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

            call wrapperClean()
            if(rank == 0)then
                if( (mod(iter,reportNum)==0))then
                    call CPU_TIME(endTime)
                    print *,iter,maxval(abs(v(:,:,:,:,2)-v(:,:,:,:,1))),endTime-startTime
                    flush(6)
                end if
            end if

            if (maxval(abs(v(:,:,:,:,2)-v(:,:,:,:,1))) .lt. toll) then
#ifdef DOPRINTS
                print*,"done: ",iter
                flush(6)
#endif
                exitLoop = .true.
            end if
            v(:,:,:,:,1)=v(:,:,:,:,2)
        end do

        if(rank == 0)then
            if(reportNum<maxit)then
                print *," "
            end if
        end if

        lastStateV(:,:,:,:) = v(:,:,:,:,1)
    end subroutine  getAllPolicy

    FUNCTION findSteadyStateCapital(noShocks) RESULT(y)
        LOGICAL, INTENT(IN) :: noShocks
        REAL(DP), dimension(periodsForConv,2) :: y

        INTEGER(I4B) :: i,j,ii,jj
        REAL(DP) :: tempD, tempD2
        REAL(DP), DIMENSION(n_s,n_a) :: y2
        INTEGER(I4B), DIMENSION(periodsForConv+periodsToCut) :: z
        REAL(DP), DIMENSION(periodsForConv+periodsToCut) :: aggK
        TYPE(household), DIMENSION(numHouseholds) :: hhs
        REAL(DP), DIMENSION(n_s,n_a) :: gint,tempInt
        REAL(DP),DIMENSION(numHouseholds) :: tempArray,tempArray2,num
        INTEGER, DIMENSION(numHouseholds) :: tempEmp,tempEmp2

        !************
        ! MPI vars
        !************
        INTEGER rank, ierr, mysize
!        INTEGER,dimension(MPI_STATUS_SIZE):: stat

        !************
        ! Timing variables
        !************
        real(DP) :: startTime, endTime

        CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mysize, ierr)

        call CPU_TIME(startTime)
        !first, draw T states for the economy
        num= rand(myseed)
        z(1)=1

        do i=2,periodsForConv+periodsToCut
            num(1)=rand(0)
            flush(6)
            z(i)=1
            do j=1,n_z-1
                if(num(1)>sum(zTransition(z(i-1),1:j)))then
                    z(i)=j+1
                else
                    exit
                end if
            end do
        end do

        !we now have our series of economic shocks. Now let's set initial distribution
        if(noshocks)then
            do i=1,numHouseholds,2
                hhs(i)%employmentState=1
                hhs(i)%capital=k(1)
                hhs(i+1)%employmentState=2
                hhs(i+1)%capital=k(1)
            end do
        else
            hhs(:)=ssDistrib
        end if

        !calculate aggregate capital
        aggK(1) = sum(hhs(:)%capital)/numHouseholds

        !now let's update
        if((reportNum<maxit).and.(rank == 0))then
            print *, "SS"
            print *,"------------------------"
            print *,"      iter      Capital                  Time"
            flush(6)
        end if

        do i=2,periodsForConv+periodsToCut
            !interpolate all the policy functions to the current aggregate capital level
            if(noshocks)then
                gint = g(1,1,:,:)
            else
                gint=0
                tempInt=0
                do ii=rank+1,n_a,mysize
                    do j=1,n_s
                        tempInt(j,ii)=linear(g(z(i-1),:,j,ii),k,aggK(i-1))
                    end do
                end do
                call MPI_BARRIER(MPI_COMM_WORLD,ierr)
                call MPI_ALLREDUCE(tempInt,gInt,n_s*n_a,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
            end if

            !set this period's employment and capital levels
            do j=1,n_s
                tempD=(gint(j,2)-gint(j,1))/(a(2)-a(1))
                tempD2=(gint(j,n_a)-gint(j,n_a-1))/(a(n_a)-a(n_a-1))
                call spline(a,gint(j,:),tempD,tempD2,y2(j,:))
            end do

            tempArray=0
            tempEmp=0
            tempArray2=0
            tempEmp2=0
            do j=1,numHouseholds
                num(j)=rand(0)
            end do

            do j=rank+1,numHouseholds,mysize
                ii = hhs(j)%employmentState
                tempArray(j)=splint(a,gint(ii,:),y2(ii,:),hhs(j)%capital)
                tempEmp(j)=1
                do jj=1,n_s-1
                    if(num(j)>sum(transition(z(i-1),ii,z(i),1:jj)))then
                        tempEmp(j)=jj+1
                    else
                        exit
                    end if
                end do
            end do
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            call MPI_ALLREDUCE(tempEmp,tempEmp2,numHouseholds,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(tempArray,tempArray2,numHouseholds,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

            do j=1,numHouseholds
                hhs(j)%capital=min(max(tempArray2(j),a_min),a_max)
                hhs(j)%employmentState=tempEmp2(j)
            end do

            aggK(i)=sum(hhs(:)%capital)/numHouseholds

            if( (mod(i,reportNum)==0) .and. (rank ==0))then
                call CPU_TIME(endTime)
                print *,i,aggK(i),endTime-startTime
                flush(6)
            end if
        end do
        if(noShocks)then
            ssDistrib=hhs
        end if
        y(:,1)=z(periodsToCut+1:periodsForConv+periodsToCut)
        y(:,2)=aggK(periodsToCut+1:periodsForConv+periodsToCut)
    end function findSteadyStateCapital

end module aiyagariSolve

!**************************************************************************
!
!**************************************************************************
program main

    use aiyagariSolve
    use brentWrapper
    REAL(DP) :: intDiff,xmin
    REAL(DP) :: capitalShare = 0.36D0
    PROCEDURE(template_function), POINTER :: func
    PROCEDURE(template_function2), POINTER :: d1func, d2func
    PROCEDURE(template_function3), POINTER :: func2
    INTEGER :: temp,temp2,printEvery=500,whichSet=3,ierr
    character(LEN=15) :: arg1,arg2
    logical :: readFromFile

    !************
    ! Timing variables
    !************
    real(DP) :: startTime, endTime

    readFromFile=.false.
    temp2=COMMAND_ARGUMENT_COUNT()
    if(temp2 > 0)then
        call GET_COMMAND_ARGUMENT(1, arg1, whichSet)
        read (arg1,*) whichSet
    end if
    if(temp2 > 1)then
        call GET_COMMAND_ARGUMENT(2, arg1, printEvery)
        read (arg1,*) printEvery
    end if
    if(temp2 > 2)then
        call GET_COMMAND_ARGUMENT(3, arg1, temp)
        read (arg1,*) readFromFile
    end if

    d1func => d1prod
    d2func => d2prod
    func2 => impliedCapital

    arg1="policy"
    arg2= "distrib"
    call setParams(d1func, d2func, func2, arg1,arg2 , capitalShare, whichSet, printEvery)
    call CPU_TIME(startTime)
    CALL MPI_INIT(ierr)
    call beginKrusellSmith(.not. readFromFile)
    CALL MPI_FINALIZE(ierr)
    call CPU_TIME(endTime)
    print *,"a             ",xmin, intDiff, endTime-startTime
    flush(6)

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
