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
        INTEGER(I4B), PARAMETER :: ITMAX=50!00
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
                            if (i /= ilo) then
                             y(i)=func(p(i,:))
                            end if
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

end module nr

program q3
use nr
use nrtype
    implicit none

        REAL(DP), parameter :: A=1.0D0
        REAL(DP), parameter :: alf=0.4D0
        REAL(DP), parameter :: bet=0.9D0
        REAL(DP), parameter :: cbar=0.0D0
        integer, parameter :: NUM_GRID_POINTS = 201
        integer, parameter :: NUM_ITER = 100
        integer, parameter :: NUM_POINTS=NUM_GRID_POINTS*NUM_ITER

        REAL(DP), dimension(NUM_GRID_POINTS) :: k, k_tilde,splinePoints
        REAL(DP) :: mydel
        REAL(DP), dimension(NUM_GRID_POINTS,NUM_ITER) ::v_iter
        REAL(DP), dimension(NUM_GRID_POINTS,NUM_GRID_POINTS) :: RHS
        INTEGER(I4B) :: i,j
        PROCEDURE(template_function), POINTER :: func
        REAL(DP), dimension(2,1) :: pointToEval
        REAL(DP), dimension(2) :: valAtPoint

        data v_iter /NUM_POINTS * 0.0D0/

        func => splineMinVal

        mydel=1.0D0
    call policy(mydel,"q7del1p0.csv")
!        mydel=0.1D0
 !   call policy(mydel,"q7del0p1.csv")

contains
    subroutine policy(del,outputFile)
        real(KIND=8),INTENT(IN) :: del
        character(LEN=*),INTENT(IN) :: outputFile

        real(KIND=8) :: kss, kd, ku,eps, err, slope1, slope2, temp2, tol=10e-3
        real(KIND=8), dimension(NUM_GRID_POINTS,1) :: one
        integer, dimension(NUM_GRID_POINTS) :: ind_k_tilde
        real(KIND=8), dimension(NUM_GRID_POINTS,NUM_ITER) ::c_policy
        real(KIND=8), dimension(NUM_GRID_POINTS,NUM_GRID_POINTS) ::Ut, PPF, c

        real(KIND=8), dimension(NUM_GRID_POINTS,1) :: colPtr1, colPtr2
        integer  :: iter

        !some basic initialization
        data c_policy /NUM_POINTS * 0.0D0/
        data one /NUM_GRID_POINTS * 1.0D0 /
        print *,"Starting calculation."
        flush(6)

        eps=epsilon(1.0D0)

        !calculating the kss
        kss=(1.0D0/(bet*alf*A)-(1.0D0-del)/(alf*A))**(1.0D0/(alf-1.0D0))
        !lowest possible capital expenditure
        kd  = .01*kss
        ku  = 1.5*kss

        k(1)=epsilon(1.0D0)
        do i=2,NUM_GRID_POINTS-1
            k(i) = kd+(i-2)*(ku-kd)/(NUM_GRID_POINTS-1)
        end do
        k(NUM_GRID_POINTS)=ku

        !An array where cell i,j indicates how much we would produce if we had
        !ki capital level this period and chose kj capital level for next
        !period (note: kj is independent of output, so all elements in a given
        !row are the same)
        colPtr1 = reshape(k, (/NUM_GRID_POINTS , 1/))
        PPF = transpose(matmul(one,transpose((1.0D0-del)*colPtr1+A*colPtr1**alf)))

        !consumption in period. Cell i,j is consumption if we had used ki to
        !produce and kj is next period's capital choice
        c = PPF-matmul(one,transpose(colPtr1))

        ! we have a minimum level of consumption (to avoid crazy values)
        c = max(c,-cbar+eps)

        !so, we can now find utility for each consumption/capital pair
        !where Ut(i,j) is the utility we obtain by consumption where we produced
        !with ki capital and saved kj capital for next period.
        Ut = log(c+cbar)

        !v_iter contains our initial guess at the value function. So
        !RHS contains our initial guess at  the value of the bellman equation
        colPtr2 = reshape(v_iter(:,1), (/NUM_GRID_POINTS , 1/))
        RHS = Ut+bet*matmul(one,transpose(colPtr2))

        ! Numerical value function iterations
        i=1
        err = 100.0D0
        do while ((err > 1.0e-6) .and. (i<NUM_ITER))

            i=i+1
            !RHS gives us the value function for capital on the current
            !grid. But this is not neccessarily the maximum value. So,
            !we will use cubic spline to estimate the value function.
            !Then, we will find the maximum using the amoeba algorithm.
            !This will give us k_tilde for the next period at each
            !starting k, with which we can then get new estimates of
            !the value function.

            k_tilde(1)=k(1)
            do j=2,NUM_GRID_POINTS
                pointToEval(1,1) = k(2)
                valAtPoint(1) = splineMinVal(pointToEval(1,:))
                pointToEval(2,1) = k(NUM_GRID_POINTS)
                valAtPoint(2) = splineMinVal(pointToEval(2,:))

                !okay, let's estimate the value function using cubic spline
                slope1 = (RHS(j,3)-RHS(j,2))/(k(3)-k(2))
                slope1 = (RHS(j,NUM_GRID_POINTS)-RHS(j,NUM_GRID_POINTS-1))/(k(NUM_GRID_POINTS)-k(NUM_GRID_POINTS-1))

                call spline(k,RHS(j,:),slope1,slope2,splinePoints)

                CALL amoeba(pointToEval,valAtPoint, tol,func,iter)

                v_iter(j,i) = -valAtPoint(1)
                k_tilde(j)  = pointToEval(1,1)

                temp2 = k(j)**alf + (1.0D0-del)*k(j) - k_tilde(j)
                if (temp2 < 0.0D0)then
                    print *,i,j
                    print *, "Dumb result: ",k(j),":",temp2,":",k_tilde(j)
                    stop 1
                end if
            end do
            c_policy(:,i) = k(:)**alf +(1.0D0-del)*k(:) - k_tilde(:)

            print *,"k: ",k
            print *,"c: ",c_policy(:,i)
            print *,"k':",k_tilde
            print *,"v':",v_iter(:,i)
            print *," "
            ! need to pretend each of these columns is an array
            colPtr1 = reshape(v_iter(:,i-1), (/NUM_GRID_POINTS , 1/))
            colPtr2 = reshape(v_iter(:,i), (/NUM_GRID_POINTS , 1/))
            ! checking the norm of the difference between successive iterations:
            ! this is convergence criterion
            err = maxval(abs(colPtr2 - colPtr1))
            RHS = Ut + bet*matmul(one,transpose(colPtr2))

         end do

        print *,"Finished calculation. ",i," iterations"

        do i=1,NUM_GRID_POINTS
!            write(1,*) k(i),",",k_tilde(i),",",k(i)
        end do

    end subroutine policy

    !-------------------------------------
    FUNCTION splineMinVal(x) RESULT (z)
        !-------------------------------------
        REAL(DP), DIMENSION(:), INTENT(IN) :: x
        REAL(DP) :: z
        z=-splineMaxVal(x)
    END FUNCTION splineMinVal

    !-------------------------------------
    FUNCTION splineMaxVal(x) RESULT (z)
        !-------------------------------------
        REAL(DP), DIMENSION(:), INTENT(IN) :: x
        REAL(DP) :: z
        REAL(DP) :: result1, c, eps=epsilon(1.0D0), prod, diff, tempc
        INTEGER :: which
        prod = k(j)**alf +(1.0D0-mydel)*k(j)
        tempc = prod - x(1)
        if(tempc<eps)then
            which = 1
            c=eps
            z=-10.0e20
        else if (tempc>prod) then
            which = 2
            c=prod-eps
            z=-10.0e20
        else if(x(1)>k(NUM_GRID_POINTS)) then
            which = 3
            c=prod-k(NUM_GRID_POINTS)
            z=-10.0e20
        else
            which = 3
            c=tempc
        diff=c-tempc
        result1=splint(k,RHS(j,:),splinePoints,x(1)-diff)
        z=log(c)+bet*result1
        if(z>0.0D0)then
            print *,"func:",k
            print *,"RHS: ",RHS(j,:)
            print *,k(j),x(1),result1,z
            flush(6)
            stop 0
        end if
        end if

!        print *,"func: ","k:",k(j),"p:",prod,"tx:",x(1),"x:",x(1)-diff,"c:",c,z
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
