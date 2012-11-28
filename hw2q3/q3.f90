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

    FUNCTION brent(ax,bx,cx,func,tol,xmin)
        USE nrtype; USE nrutil, ONLY : nrerror
        IMPLICIT NONE
        REAL(dp), INTENT(IN) :: ax,bx,cx,tol
        REAL(dp), INTENT(OUT) :: xmin
        REAL(dp) :: brent
        PROCEDURE(template_function), POINTER :: func
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
        fx=func(x)
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
            fu=func(u)
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
    CONTAINS
        !BL
        SUBROUTINE shft(a,b,c,d)
            REAL(dp), INTENT(OUT) :: a
            REAL(dp), INTENT(INOUT) :: b,c
            REAL(dp), INTENT(IN) :: d
            a=b
            b=c
            c=d
        END SUBROUTINE shft
    END FUNCTION brent

end module nr

program q3
use nr
use nrtype
    implicit none

        REAL(DP), parameter :: A=1.0D0
        REAL(DP), parameter :: alf=0.4D0
        REAL(DP), parameter :: bet=0.9D0
        REAL(DP), parameter :: cbar=0.0D0
        integer, parameter :: NUM_GRID_POINTS = 1001
        integer, parameter :: NUM_ITER = 500
        integer, parameter :: NUM_POINTS=NUM_GRID_POINTS*NUM_ITER

        REAL(DP), dimension(NUM_GRID_POINTS) :: k, splinePoints, RHS, PPF
        REAL(DP) :: mydel
        REAL(DP), dimension(NUM_GRID_POINTS,NUM_ITER) ::v_iter, k_tilde
        INTEGER(I4B) :: i,j
        PROCEDURE(template_function), POINTER :: func

        data v_iter /NUM_POINTS * 0.0D0/

        func => splineMinVal

        mydel=1.0D0
        call policy(mydel,"q3del1p0.csv")
        mydel=0.1D0
        call policy(mydel,"q3del0p1.csv")

contains
    subroutine policy(del,outputFile)
        real(KIND=8),INTENT(IN) :: del
        character(LEN=*),INTENT(IN) :: outputFile

        real(KIND=8) :: kss, kd, ku,eps, err, slope1, slope2, tol=10e-3
        real(KIND=8), dimension(NUM_GRID_POINTS,1) :: one
        real(KIND=8), dimension(NUM_GRID_POINTS,NUM_ITER) ::c_policy
        REAL(DP), dimension(NUM_GRID_POINTS) :: newRHS

        real(KIND=8), dimension(NUM_GRID_POINTS,1) :: colPtr1, colPtr2
        integer  :: tempInt1

        !some basic initialization
        data c_policy /NUM_POINTS * 0.0D0/
        data one /NUM_GRID_POINTS * 1.0D0 /
        print *,"Starting calculation."
        flush(6)

        !calculating the kss
        kss=(1.0D0/(bet*alf*A)-(1.0D0-del)/(alf*A))**(1.0D0/(alf-1.0D0))
        !lowest possible capital expenditure
        kd  = .01*kss
        ku  = 1.5*kss

        call sub_grid_generation(k,kd,ku,3.0D0, 1)

        !An array where cell i,j indicates how much we would produce if we had
        !ki capital level this period and chose kj capital level for next
        !period (note: kj is independent of output, so all elements in a given
        !row are the same)
        PPF = (1.0D0-del)*k+A*k**alf

        eps=epsilon(1.0D0)
        RHS = log(PPF)+100

        ! Numerical value function iterations
        i=1
        err = 100.0D0
        do while ((err > 1.0e-6) .and. (i<NUM_ITER))
            i=i+1
            print *,i
            flush(6)
            !RHS gives us the value function for capital on the current
            !grid. But this is not neccessarily the maximum value. So,
            !we will use cubic spline to estimate the value function.
            !Then, we will find the maximum using the brent algorithm.
            !This will give us k_tilde for the next period at each
            !starting k, with which we can then get new estimates of
            !the value function.
            slope1=(RHS(2)-RHS(1))/(k(2)-k(1))
            slope2=(RHS(NUM_GRID_POINTS)-RHS(NUM_GRID_POINTS-1))/(k(NUM_GRID_POINTS)-k(NUM_GRID_POINTS-1))
            call spline(k,RHS,slope1,slope2,splinePoints)

            do j=1,NUM_GRID_POINTS

                newRHS = -(log(max(PPF(j)-k,eps))+bet*RHS)

                tempInt1=minloc(newRHS,1)
                if (tempInt1 .eq. 1) tempInt1=2
                if (tempInt1 .eq. NUM_GRID_POINTS) tempInt1=NUM_GRID_POINTS-1
                v_iter(j,i)= -brent(k(tempInt1-1), k(tempInt1), k(tempInt1+1), &
                                        &func, tol,k_tilde(j,i))
            end do
            c_policy(:,i) = max(PPF - k_tilde(:,i),eps)

            ! need to pretend each of these columns is an array
            colPtr1 = reshape(v_iter(:,i-1), (/NUM_GRID_POINTS , 1/))
            colPtr2 = reshape(v_iter(:,i), (/NUM_GRID_POINTS , 1/))
            ! checking the norm of the difference between successive iterations:
            ! this is convergence criterion
            err = maxval(abs(colPtr2 - colPtr1))
            RHS = v_iter(:,i)
         end do

        print *,"Finished calculation. ",i," iterations"

        open(1, file=outputFile)
        do j=1,NUM_GRID_POINTS
            write(1,*) k(j),",",k_tilde(j,i),",",k(j)
        end do
        close(1)

    end subroutine policy

    !-------------------------------------
    FUNCTION splineMinVal(x) RESULT (z)
        !-------------------------------------
        REAL(DP),INTENT(IN) :: x
        REAL(DP) :: z
        z=-splineMaxVal(x)
    END FUNCTION splineMinVal

    !-------------------------------------
    FUNCTION splineMaxVal(x) RESULT (z)
        !-------------------------------------
        REAL(DP), INTENT(IN) :: x
        REAL(DP) :: z
        REAL(DP),parameter :: eps=epsilon(1.0D0)
        REAL(DP) :: fixedC
        fixedC=max(PPF(j)-x,eps)
        z=log(fixedC)+bet*splint(k,RHS,splinePoints,x)
    END FUNCTION splineMaxVal

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

end program q3
