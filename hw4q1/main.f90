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
        INTEGER(I4B), PARAMETER :: ITMAX=100
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

        !local variables
        REAL(DP) :: baseGrid, nextGrid, fBase, fNext
        REAL(DP) :: slope
        integer, DIMENSION(1) :: indexInGridTemp
        integer :: indexInGrid, nextPoint

        indexInGridTemp=minloc(abs(gridPoints-evalPoint))
        indexInGrid = indexInGridTemp(1)
        if(evalPoint>gridPoints(indexInGrid)) then
            nextPoint = indexInGrid + 1
        else
            nextPoint = indexInGrid
            indexInGrid = indexInGrid - 1
        end if

        ! if we are past the max grid point, use the left derivative
        if( (indexInGrid==size(gridPoints)) .and. (evalPoint>gridPoints(indexInGrid)) )then
            nextPoint=indexInGrid
            indexInGrid=indexInGrid-1
        endif

        ! if we are before the min grid point, use the right derivative
        if( (indexInGrid==1) .and. (evalPoint<gridPoints(indexInGrid)) )then
            nextPoint=indexInGrid
            indexInGrid=indexInGrid+1
        endif

        baseGrid=gridPoints(indexInGrid)
        nextGrid=gridPoints(nextPoint)

        fBase = func(indexInGrid)
        fNext = func(nextPoint)
        slope = (fNext-fBase)/(nextGrid-baseGrid)

        z=fBase+slope*(evalPoint-baseGrid)
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
    REAL(DP) :: alpha,rho,beta,r
    LOGICAL :: doSpline

    PRIVATE n_s,currentState,currentCapital,a,s,v,y2,transition,alpha,beta,rho,r,doSpline

contains

    subroutine wrapperCreate(grid,mytransition,mystates,myalpha,myrho,mybeta,myr)
        REAL(DP), DIMENSION(:) ::grid, mystates
        REAL(DP), DIMENSION(:,:) :: mytransition
        REAL(DP), INTENT(IN) :: myalpha,myrho,mybeta,myr

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
        temp = a(currentCapital)*(1+r)+s(currentState)-x
        if(temp < 0.0D0) then
            temp=0.0D0
        end if
        z=-(((temp)**rho+beta*rp(1)**(rho/alpha))**(1/rho))
    end function valueFunction

end module brentWrapper


module aiyagariSolve
    use rowenhorst_module
    USE nrtype
    USE nrutil
    USE brentWrapper
    implicit none

    real*8, parameter                   ::  phi=.9D0,sigma=.4D0
    integer, parameter                  ::  n_s=3
    integer, parameter                  ::  n_a=101
    real*8, parameter                   ::  curv=3.0D0,a_max=45.0D0
    real*8, parameter                   :: beta=0.90
    integer, parameter                  ::  maxit=2000
    real*8,parameter                    ::  toll=1D-8,tol2=1D-8
    LOGICAL                             :: doSpline = .TRUE.

    real*8, dimension(n_s)              ::  s,stationary
    real*8, dimension(n_s,n_s)          ::  transition
    real*8, dimension(n_a)              ::  a
    real*8                              ::  a_min


    real*8                              :: rho,alpha,RRA,EIS,mysigma
    real*8, dimension(n_s,n_a,maxit)    ::  v,g
    integer                             ::  reportNum
    PROCEDURE(template_function), POINTER :: funcParam

    character(LEN=20)                   :: policyOutput

    PRIVATE phi, sigma, n_s, n_a, curv, a_max, beta, maxit, toll, doSpline
    PRIVATE s, stationary, a, a_min,rho,alpha,RRA,EIS,mysigma,v,g
    PRIVATE reportNum,funcParam,policyOutput

contains
    subroutine setParams(myrra, myeis, func, file1, every, doLinear)
        REAL(KIND=8), INTENT(IN) :: myrra, myeis
        PROCEDURE(template_function), POINTER, INTENT(IN) :: func
        character(LEN=*),INTENT(IN) :: file1
        INTEGER, OPTIONAL, INTENT(IN) :: every
        LOGICAL, OPTIONAL, INTENT(IN) :: doLinear

        RRA=myrra
        EIS=myeis
        funcParam => func

        rho=1-1.0D0/EIS
        alpha=1-RRA

        reportNum = 50
        if(PRESENT(every)) then
            reportNum=every
        end if

        if(PRESENT(doLinear))then
            doSpline=.not. doLinear
        end if
        policyOutput = file1

        !**************************************************************************
        ! We use the rowenhorst method to obtain the transition matrix GAMMA,
        ! the stationary probability distribution stationary, and the shocks s.
        !**************************************************************************
        call rouwenhorst(phi,sigma,transition,s,stationary)
        s=1+s

    end subroutine setParams

    function aggregateBonds(r) RESULT (z)
        ! inputs: r - the interest rate to test
        ! outputs: z - the aggregate level of borrowing
        REAL(KIND=8), INTENT(IN) :: r
        REAL(KIND=8) :: z
        INTEGER :: iterCount, i
        real*8                              ::  incr
        real*8, dimension(n_s,n_a)          ::  steadyStateCapital

        a_min = max(-s(1)/r+1.0D0,-3.0D0)

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
        do iterCount=1,n_s
            v(iterCount,:,1)=(a-a_min)**2
        end do
        g=0D0

        call wrapperCreate(a,transition,s,alpha,rho,beta,r)
        iterCount=getPolicyForInterest(r)
        call wrapperDestroy()

        !*****************************************************
        ! find the steady state capital
        !*****************************************************
        call findSteadyState(g(:,:,iterCount),steadyStateCapital)

        !**************************************************
        ! now, calculate total borrowing using policy functions
        !**************************************************
        !TO DO

        !***************************************************
        ! return the absolute value of total borrowing
        !***************************************************
        !TO DO
        z=0.0D0
    end function aggregateBonds

    subroutine findSteadyState(capitalPolicy, statDist)
        !INPUTS: capitalPolicy - the policy function for each state
        !OUTPUTS: statDist - the stationary dist (note: pdf, not cdf)
        REAL(DP), dimension(n_s,n_a), INTENT(IN) :: capitalPolicy
        real(DP), dimension(n_s,n_a), intent(out) :: statDist
        real(DP), dimension(n_s,n_a) ::f_o, f_o_hat, f_n
        real(DP) :: diff
        INTEGER :: i,j, counter

        !setting initial guess to uniform dist across asset grid. First we have it
        !a cdf, and then convert to the appropriate pdf
        do i=1,n_s
            do j=1,n_a
                statDist(i,j) = (a(j) - a(1))/(a(n_a)-a(1))
            end do
        end do

        f_n=statDist

        !normalize so that we can compare
        do i=1,n_s
            f_n(i,:)=f_n(i,:)*stationary(i)/f_n(i,n_a)
        end do

        ! time to iterate
        diff=100
        counter = 0
        do while((diff>tol2) .and. (counter<100))
            counter = counter + 1
            f_o=f_n

            do i=1,n_s
                do j=1,n_a
                    f_o_hat(i,j) = linear(f_o(i,:),capitalPolicy(i,:),a(j))
                end do
            end do

            ! need to make sure monotonic and bounded between 0 and 1
            do i=1,n_s
                diff=-1e-5
                do j=1,n_a
                    if(f_o_hat(i,j)>1.0_dp) then
                        f_o_hat(i,j)=1.0_dp
                    else if (f_o_hat(i,j)<0.0_dp) then
                        if(j==1)then
                            f_o_hat(i,j)=0.0D0
                        else
                            f_o_hat(i,j)=f_o_hat(i,j-1)+epsilon(1.0D0)
                        end if
                    else if (isnan(f_o_hat(i,j))) then
                        print *, "Error: f_o_hat is NAN. Counter:",counter," shock: ",i,"element: ",j
                        print *,"Capital Grid"
                        print *,a(:)
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
                        print *,a
                        print *,f_o_hat(i,:)
                        stop 0
                    end if
                    diff=f_o_hat(i,j)
                end do
            end do

            f_n = matmul(transition, f_o_hat)

            !normalize to account for rounding errors
            do i=1,n_s
                f_n(i,:)=f_n(i,:)*stationary(i)/f_n(i,n_a)
            end do

            diff = maxval(abs(f_n - f_o))
            if (mod(counter,50)==0) then
                print*,"findSteadyState Iteration: ",counter, " diff: ",diff
            end if
        end do

        statDist = f_n
        print *,"done: ",counter
        do i=n_a,2,-1
            statDist(:,i) = statDist(:,i) - statDist(:,i-1)
        end do
        statDist(:,1) = stationary-sum(statDist(:,2:n_a),dim=2)

    end subroutine findSteadyState

    function getPolicyForInterest(r) result(y)
        REAL(KIND=8), INTENT(IN) :: r
        INTEGER :: y
        INTEGER ::  i,it,j,iter
        REAL(DP) :: tempD,tempD2
        real*8, dimension(n_s,n_a)          ::  y2
        real*8                              ::  kr1,kr2,kr3

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
                do i=1,n_s
                    do it=1,n_a
                        kr1=a(1)
                        kr3=min(a(it)*(1+r)+s(i),a(n_a))
                        kr2=(kr1+kr3)/2D0
                        !                        splineD=splineDeriv(a, v(i,:,iter-1), y2(i,:))
                        !                        v(i,it,iter)=-dbrent(funcParam,kr1,kr2,kr3,1D-10,g(i,it,iter))
                        v(i,it,iter)=-callBrent(i,it,funcParam,kr1,kr2,kr3,1D-10,g(i,it,iter))
                    end do
                end do
                call wrapperClean()
            else
                call wrapperInit(v(:,:,iter-1))
                do i=1,n_s
                    do it=1,n_a
                        kr1=a(1)
                        kr3=min(a(it)*(1+r)+s(i),a(n_a))
                        kr2=(kr1+kr3)/2D0
                        v(i,it,iter)=-callBrent(i,it,funcParam,kr1,kr2,kr3,1D-10,g(i,it,iter))
                    end do
                end do
                call wrapperClean()
            end if
            if(mod(iter,reportNum)==0)then
                print *,"r:",r,"iter: ",iter,"diff: ",maxval(maxval(abs(v(:,:,iter)-v(:,:,iter-1)),1),2)
                flush(6)
            end if
            if (     maxval(maxval(abs(v(:,:,iter)-v(:,:,iter-1)),1),2)    .lt.    toll     ) then
                print*,"done: ",iter
                flush(6)
                exit
            end if
        end do

        y=min(iter,maxit)
    end function         getPolicyForInterest

end module aiyagariSolve

!**************************************************************************
!
!**************************************************************************
program main

    use aiyagariSolve
    use brentWrapper

    REAL(DP) :: RRA, EIS
    PROCEDURE(template_function), POINTER :: func

    RRA=2.0D0
    EIS=2.0D0
    func => valueFunction
    call setParams(RRA, EIS, func, "policyR2E2", 50, .FALSE.)
    print *,aggregateBonds(0.05D0)
end program main
