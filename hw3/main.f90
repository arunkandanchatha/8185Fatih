!**************************************************************************
!
!**************************************************************************
program main
    use rowenhorst_module
    USE nrtype
    USE nrutil
    implicit none
    !**************************************************************************
    ! Parameters and variables needed for the rowenhorst method
    !**************************************************************************
    real*8, parameter                   ::  phi=.9D0,sigma=.4D0
    integer, parameter                  ::  n_s=3
    real*8, dimension(n_s)              ::  s,stationary
    real*8, dimension(n_s,n_s)          ::  GAMMA
    !**************************************************************************
    ! Preparing the grid for the asset a
    !**************************************************************************
    integer, parameter                  ::  n_a=101
    real*8, dimension(n_a)              ::  a
    real*8, parameter                   ::  curv=3.0D0,a_max=50.0D0
     real*8                              ::  incr, a_min
    !**************************************************************************
    ! convergence parameter for supnorm
    !**************************************************************************
    real*8                              ::  toll=1D-8
    !**************************************************************************
    ! variables necessary for the spline interpolation
    !**************************************************************************
    real*8, dimension(n_s,n_a)          ::  y2
    real*8, dimension(n_s,1)            ::  interpolated
    !**************************************************************************
    ! housekeeping inside the loop that does the iteration
    !**************************************************************************
    real*8                              ::  kr1,kr2,kr3
    !**************************************************************************
    ! other parameters of the model
    !**************************************************************************
    real*8, parameter                   ::  r=0.1D0,beta=0.90
    real*8                              :: rho,alpha,RRA,EIS
    !**************************************************************************
    ! Preparing the value function and the policy function arrays
    !**************************************************************************
    integer, parameter                  ::  maxit=2000
    real*8, dimension(n_s,n_a,maxit)    ::  v,g
    !**************************************************************************
    ! counters
    !**************************************************************************
    integer                             ::  i,iter,it,j,jj
    PROCEDURE(template_function), POINTER :: epsZin,crraFn
    LOGICAL                             :: doSpline = .TRUE.

    epsZin => func1
    crraFn => func2

    RRA=2.0D0
    EIS=2.0D0
    rho=1-1.0D0/EIS
    alpha=1-RRA
#if 1
    call question(epsZin, "policy1a", "value1a")
    rho=RRA
    call question(crraFn, "policy1b", "value1b")
#endif

    EIS=0.1D0
    rho=1-1/EIS
    alpha=1-RRA
#if 1
    call question(epsZin, "policy2a", "value2a")
    rho=RRA
    call question(crraFn, "policy2b", "value2b")
#endif

    RRA=10.0D0
    EIS=2.0D0
    rho=1-1/EIS
    alpha=1-RRA
#if 1
    call question(epsZin, "policy3a", "value3a",5)
    rho=RRA
    call question(crraFn, "policy3b", "value3b",5)
#endif
    EIS=0.1D0
    rho=1-1/EIS
    alpha=1-RRA
#if 1
    call question(epsZin, "policy4a", "value4a")
    rho=RRA
    call question(crraFn, "policy4b", "value4b")
#endif



contains

    SUBROUTINE question(funcParam, file1, file2, every)
        IMPLICIT NONE
        PROCEDURE(template_function), POINTER, INTENT(IN) :: funcParam
        character(LEN=*),INTENT(IN) :: file1
        character(LEN=*),INTENT(IN) :: file2
        INTEGER, OPTIONAL, INTENT(IN) :: every

        real(DP) tempD2,tempD
        integer :: reportNum=50

        if(PRESENT(every)) then
            reportNum=every
        end if

    !**************************************************************************
    ! We use the rowenhorst method to obtain the transition matrix GAMMA,
    ! the stationary probability distribution stationary, and the shocks s.
    !**************************************************************************
    call rouwenhorst(phi,sigma,GAMMA,s,stationary)
    s=1+s
    a_min = -s(1)/r+1.0D0

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
    do iter=1,n_s
        v(iter,:,1)=(a-a_min)**2
!        v(iter,:,1)=0.0
    end do
    g=0D0

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
        do i=1,n_s
            do it=1,n_a
                kr1=a(1)
                kr3=min(a(it)*(1+r)+s(i),a(n_a))
                kr2=(kr1+kr3)/2D0
                v(i,it,iter)=-brent(funcParam,kr1,kr2,kr3,1D-10,g(i,it,iter))
            end do
        end do
else
        do i=1,n_s
            do it=1,n_a
                kr1=a(1)
                kr3=min(a(it)*(1+r)+s(i),a(n_a))
                kr2=(kr1+kr3)/2D0
                v(i,it,iter)=-brent(funcParam,kr1,kr2,kr3,1D-10,g(i,it,iter))
            end do
        end do
end if
        if(mod(iter,reportNum)==0)then
            print *,"count: ",iter,file1,maxval(maxval(abs(v(:,:,iter)-v(:,:,iter-1)),1),2)
            flush(6)
        end if
        if (     maxval(maxval(abs(v(:,:,iter)-v(:,:,iter-1)),1),2)    .lt.    toll     ) then
            print*,"done: ",iter
            flush(6)
            exit
        end if
    end do

    !**************************************************************************
    ! we print the results of the iterations in a file
    !**************************************************************************
    open(unit=1,file=file1)
    open(unit=2,file=file2)
    open(unit=3,file="asset")
    open(unit=4,file="shock")
    open(unit=5,file="gamma")

    do it=1,n_s
        write(1,*) g(it,:,iter)
        write(2,*) v(it,:,iter)
    end do
    do jj=1,n_a
        write(3,*) a(jj)
    end do
    do jj=1,n_s
        write(4,*) s(jj)
    end do
    do jj=1,n_s
        write(5,*) GAMMA(jj,:)
    end do

    close(1)
    close(2)
    close(3)
    close(4)
    close(5)

    END SUBROUTINE question
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

    !    FUNCTION dbrent(ax,bx,cx,tol,xmin)
    !        IMPLICIT NONE
    !        REAL(dp), INTENT(IN) :: ax,bx,cx,tol
    !        REAL(dp), INTENT(OUT) :: xmin
    !        REAL(dp) :: dbrent
    !        INTEGER(I4B), PARAMETER :: ITMAX=100
    !        REAL(dp), PARAMETER :: ZEPS=1.0e-3_dp*epsilon(ax)
    !        INTEGER(I4B) :: iter
    !        REAL(dp) :: a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2,&
    !            u,u1,u2,v,w,x,xm
    !        LOGICAL :: ok1,ok2
    !        a=min(ax,cx)
    !        b=max(ax,cx)
    !        v=bx
    !        w=v
    !        x=v
    !        e=0.0
    !        fx=func(x)
    !        fv=fx
    !        fw=fx
    !        dx=dfunc(x)
    !        dv=dx
    !        dw=dx
    !        do iter=1,ITMAX
    !            xm=0.5_dp*(a+b)
    !            tol1=tol*abs(x)+ZEPS
    !            tol2=2.0_dp*tol1
    !            if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) exit
    !            if (abs(e) > tol1) then
    !                d1=2.0_dp*(b-a)
    !                d2=d1
    !                if (dw /= dx) d1=(w-x)*dx/(dx-dw)
    !                if (dv /= dx) d2=(v-x)*dx/(dx-dv)
    !                u1=x+d1
    !                u2=x+d2
    !                ok1=((a-u1)*(u1-b) > 0.0) .and. (dx*d1 <= 0.0)
    !                ok2=((a-u2)*(u2-b) > 0.0) .and. (dx*d2 <= 0.0)
    !                olde=e
    !                e=d
    !                if (ok1 .or. ok2) then
    !                    if (ok1 .and. ok2) then
    !                        d=merge(d1,d2, abs(d1) < abs(d2))
    !                    else
    !                        d=merge(d1,d2,ok1)
    !                    end if
    !                    if (abs(d) <= abs(0.5_dp*olde)) then
    !                        u=x+d
    !                        if (u-a < tol2 .or. b-u < tol2) &
    !                            d=sign(tol1,xm-x)
    !                    else
    !                        e=merge(a,b, dx >= 0.0)-x
    !                        d=0.5_dp*e
    !                    end if
    !                else
    !                    e=merge(a,b, dx >= 0.0)-x
    !                    d=0.5_dp*e
    !                end if
    !            else
    !                e=merge(a,b, dx >= 0.0)-x
    !                d=0.5_dp*e
    !            end if
    !            if (abs(d) >= tol1) then
    !                u=x+d
    !                fu=func(u)
    !            else
    !                u=x+sign(tol1,d)
    !                fu=func(u)
    !                if (fu > fx) exit
    !            end if
    !            du=dfunc(u)
    !            if (fu <= fx) then
    !                if (u >= x) then
    !                    a=x
    !                else
    !                    b=x
    !                end if
    !                call mov3(v,fv,dv,w,fw,dw)
    !                call mov3(w,fw,dw,x,fx,dx)
    !                call mov3(x,fx,dx,u,fu,du)
    !            else
    !                if (u < x) then
    !                    a=u
    !                else
    !                    b=u
    !                end if
    !                if (fu <= fw .or. w == x) then
    !                    call mov3(v,fv,dv,w,fw,dw)
    !                    call mov3(w,fw,dw,u,fu,du)
    !                else if (fu <= fv .or. v == x .or. v == w) then
    !                    call mov3(v,fv,dv,u,fu,du)
    !                end if
    !            end if
    !        end do
    !        if (iter > ITMAX) call nrerror('dbrent: exceeded maximum iterations')
    !        xmin=x
    !        dbrent=fx
    !    END FUNCTION dbrent

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

    FUNCTION linear(func,evalPoint,gridPoints) RESULT(z)
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
            nextPoint = indexInGrid -1
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

    function func1(x) RESULT(z)
        IMPLICIT NONE
            REAL(KIND=8), INTENT(IN) :: x
            REAL(KIND=8) :: z
        real*8, dimension(1)        ::  rp
        integer                     ::  jj
        real*8                      :: temp
        do jj=1,n_s
if (doSpline) then
            interpolated(jj,1)=splint(a,v(jj,:,iter-1),y2(jj,:),x)
else
            interpolated(jj,1)=linear(v(jj,:,iter-1),x,a)
end if
        end do
        rp=matmul(GAMMA(i,:),interpolated**alpha)
        temp = a(it)*(1+r)+s(i)-x
        if(temp < 0.0D0) then
            temp=0.0D0
        end if
            z=-(((temp)**rho+beta*rp(1)**(rho/alpha))**(1/rho))
    end function func1

    function func2(x) RESULT(z)
        IMPLICIT NONE
            REAL(KIND=8), INTENT(IN) :: x
            REAL(KIND=8) :: z
        real*8, dimension(1)        ::  rp
        integer                     ::  jj
        real*8                      :: temp
        do jj=1,n_s
if (doSpline) then
            interpolated(jj,1)=splint(a,v(jj,:,iter-1),y2(jj,:),x)
else
            interpolated(jj,1)=linear(v(jj,:,iter-1),x,a)
end if
        end do

        rp=matmul(GAMMA(i,:),interpolated)
        temp = a(it)*(1+r)+s(i)-x
        if(temp < 0.0D0) then
            print *,a(it),temp,x,s(i)
            stop 0
            temp=epsilon(1.0D0)
        end if
        z=-(temp**(1-rho)/(1-rho)+beta*rp(1))

    end function func2


end program main
