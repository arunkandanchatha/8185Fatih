
MODULE nr
    use nrtype
    implicit none

contains
    SUBROUTINE lnsrch(xold,fold,g,p,x,f,stpmax,check,func)
        USE nrtype; USE nrutil, ONLY : assert_eq,nrerror,vabs
        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(IN) :: xold,g
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: p
        REAL(DP), INTENT(IN) :: fold,stpmax
        REAL(DP), DIMENSION(:), INTENT(OUT) :: x
        REAL(DP), INTENT(OUT) :: f
        LOGICAL(LGT), INTENT(OUT) :: check
        PROCEDURE(template_function), POINTER, INTENT(in) :: func
        REAL(DP), PARAMETER :: ALF=1.0e-4_dp,TOLX=epsilon(x)
        INTEGER(I4B) :: ndum
        REAL(DP) :: a,alam,alam2,alamin,b,disc,f2,pabs,rhs1,rhs2,slope,tmplam
        ndum=assert_eq(size(g),size(p),size(x),size(xold),'lnsrch')
        check=.false.
        pabs=vabs(p(:))
        if (pabs > stpmax) p(:)=p(:)*stpmax/pabs
        slope=dot_product(g,p)
        if (slope >= 0.0) call nrerror('roundoff problem in lnsrch')
        alamin=TOLX/maxval(abs(p(:))/max(abs(xold(:)),1.0_dp))
        alam=1.0
        do
            x(:)=xold(:)+alam*p(:)
            f=func(x)
            if (alam < alamin) then
                x(:)=xold(:)
                check=.true.
                RETURN
            else if (f <= fold+ALF*alam*slope) then
                RETURN
            else
                if (alam == 1.0) then
                    tmplam=-slope/(2.0_dp*(f-fold-slope))
                else
                    rhs1=f-fold-alam*slope
                    rhs2=f2-fold-alam2*slope
                    a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
                    b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
                        (alam-alam2)
                    if (a == 0.0) then
                        tmplam=-slope/(2.0_dp*b)
                    else
                        disc=b*b-3.0_dp*a*slope
                        if (disc < 0.0) then
                            tmplam=0.5_dp*alam
                        else if (b <= 0.0) then
                            tmplam=(-b+sqrt(disc))/(3.0_dp*a)
                        else
                            tmplam=-slope/(b+sqrt(disc))
                        end if
                    end if
                    if (tmplam > 0.5_dp*alam) tmplam=0.5_dp*alam
                end if
            end if
            alam2=alam
            f2=f
            alam=max(tmplam,0.1_dp*alam)
        end do
    END SUBROUTINE lnsrch

    SUBROUTINE dfpmin(p,gtol,iter,fret,func,dfunc)
        USE nrtype; USE nrutil, ONLY : nrerror,outerprod,unit_matrix,vabs
        IMPLICIT NONE
        INTEGER(I4B), INTENT(OUT) :: iter
        REAL(DP), INTENT(IN) :: gtol
        REAL(DP), INTENT(OUT) :: fret
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: p
        PROCEDURE(template_function), POINTER, INTENT(in) :: func
        PROCEDURE(template_derivative), POINTER, INTENT(in) :: dfunc
        INTEGER(I4B), PARAMETER :: ITMAX=200
        REAL(DP), PARAMETER :: STPMX=100.0_dp,EPS=epsilon(p),TOLX=4.0_dp*EPS
        INTEGER(I4B) :: its
        LOGICAL :: check
        REAL(DP) :: den,fac,fad,fae,fp,stpmax,sumdg,sumxi
        REAL(DP), DIMENSION(size(p)) :: dg,g,hdg,pnew,xi
        REAL(DP), DIMENSION(size(p),size(p)) :: hessin

        fp=func(p)
        g=dfunc(p)
        call unit_matrix(hessin)
        xi=-g
        stpmax=STPMX*max(vabs(p),real(size(p),dp))
        do its=1,ITMAX
            iter=its
            call lnsrch(p,fp,g,xi,pnew,fret,stpmax,check,func)
            fp=fret
            xi=pnew-p
            p=pnew
            if (maxval(abs(xi)/max(abs(p),1.0_sp)) < TOLX) RETURN
            dg=g
            g=dfunc(p)
            den=max(fret,1.0_sp)
            if (maxval(abs(g)*max(abs(p),1.0_sp)/den) < gtol) RETURN
            dg=g-dg
            hdg=matmul(hessin,dg)
            fac=dot_product(dg,xi)
            fae=dot_product(dg,hdg)
            sumdg=dot_product(dg,dg)
            sumxi=dot_product(xi,xi)
            if (fac > sqrt(EPS*sumdg*sumxi)) then
                fac=1.0_sp/fac
                fad=1.0_sp/fae
                dg=fac*xi-fad*hdg
                hessin=hessin+fac*outerprod(xi,xi)-&
                    fad*outerprod(hdg,hdg)+fae*outerprod(dg,dg)
            end if
            xi=-matmul(hessin,g)
        end do
        call nrerror('dfpmin: too many iterations')
    END SUBROUTINE dfpmin

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

MODULE  hmwk2
    use nrtype
    use nr
    implicit none
    INTEGER(I8B), parameter :: NUM_ITER=10000000
    INTEGER(KIND=8), save:: thetaCalls=0
    REAL(DP), allocatable, dimension(:,:) :: intervals
    LOGICAL, save :: firstTime = .true.

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
    SUBROUTINE q1b(func,dfunc,startPoint,fret)
        !-------------------------------------

        !use a quasi-Newton based algorithm, NR’s dfpmin.f90 code
        !which implements BFGS.
        !
        ! INPUTS:
        !   1) func - the function we are trying to find the min value of
        !   2) startPoint - the co-ordinates of the starting point
        PROCEDURE(template_function), POINTER, INTENT(in) :: func
        PROCEDURE(template_derivative), POINTER, INTENT(in) :: dfunc
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: startPoint
        REAL(DP), INTENT(OUT) :: fret
        INTEGER(I4B) :: iter
        REAL(DP):: tol=1.0e-12


        if(maxval(abs(dfunc(startPoint))) > 0.0D0) then
            CALL dfpmin(startPoint,tol,iter,fret,func,dfunc)
        end if
    END SUBROUTINE q1b

    !-------------------------------------
    SUBROUTINE q1c(func,startPoint,startVals)
        !-------------------------------------

        !use nelder-meade
        !
        ! INPUTS:
        !   1) func - the function we are trying to find the min value of
        !   2) startPoint - the co-ordinates of the starting point simplex
        !   3) startVals - the value of the function at the starting points
        PROCEDURE(template_function), POINTER, INTENT(in) :: func
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: startPoint
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: startVals
        INTEGER(I4B) :: iter
        REAL(DP):: tol=1.0e-12

        CALL amoeba(startPoint,startVals, tol,func,iter)

    END SUBROUTINE q1c

    !-------------------------------------
    SUBROUTINE q1d(func,startPoints,startVals,distance,dfunc)
        !-------------------------------------

        !use nelder-meade
        !
        ! INPUTS:
        !   1) func - the function we are trying to find the min value of
        !   2) startPoint - the co-ordinates of the starting point simplex
        !   3) startVals - the value of the function at the starting points
        PROCEDURE(template_function), POINTER, INTENT(in) :: func
        PROCEDURE(template_derivative), OPTIONAL, POINTER, INTENT(in) :: dfunc
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: startPoints
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: startVals
        REAL(DP), INTENT(IN) :: distance
        INTEGER(I4B) :: iter, counter=0, i, bestCount
        REAL(DP):: tol=1.0e-2,tol2=1.0e-10,diff
        REAL(DP), DIMENSION(size(startPoints,dim=2)) :: currentBest, currentStart, currentVal
        REAL(DP) :: currentBestVal,fret

        diff=1.0D0
        currentBest = startPoints(1,:)
        currentBestVal = startVals(1)
        bestCount = 0
        DO WHILE ((counter<NUM_ITER).and.(diff>tol2))
            counter = counter + 1

            if(mod(counter,1000000)==1)then
                print *, "q1d:",counter
                flush(6)
            end if

            CALL amoeba(startPoints,startVals, tol,func,iter)

            !option to do a quasi-newton method when we are close
            if(present(dfunc))then
                currentVal = startPoints(1,:)
                fret=maxval(abs(dfunc(currentVal)))
                if(fret > 0.0D0) then
                    CALL dfpmin(currentVal,tol2,iter,startVals(1),func,dfunc)
                end if
            end if

            if (startVals(1)<currentBestVal) then
                diff = abs(currentBestVal-startVals(1))
                currentBest = startPoints(1,:)
                currentBestVal = startVals(1)
                bestCount = counter
            end if

            currentStart = getNextPoint(currentBest)
            startPoints = getSimplexAround(currentStart,distance)

            DO i=1,size(startVals)
                startVals(i) = func(startPoints(i,:))
            END DO
        END DO

        print *,"Final count: ",bestCount," Point: ",currentBest," Val: ",currentBestVal
    END SUBROUTINE q1d

    FUNCTION getNextPoint(currentBest) RESULT (y)
        REAL(DP), DIMENSION(:), INTENT(IN) :: currentBest
        REAL(DP), DIMENSION(size(currentBest)) :: y,nextVal
        REAL(DP) :: mytheta

        mytheta = nextTheta()
        if(firstTime)then
            call sobseq(nextVal,1)
            firstTime = .false.
        end if
        call sobseq(nextVal)
        nextVal = (intervals(:,2)-intervals(:,1)) * nextVal+intervals(:,1)
        y=mytheta*currentBest+(1-mytheta)*nextVal
    END FUNCTION

    FUNCTION nextTheta() RESULT (y)
        REAL(DP) :: y
        thetaCalls = thetaCalls+1
        y=(1.0D0*min(thetaCalls,NUM_ITER))/NUM_ITER
    END FUNCTION nextTheta

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
END MODULE hmwk2

PROGRAM main
    use nrtype
    use hmwk2
    implicit none
    REAL(DP), DIMENSION(2,2) :: startPoints
    REAL(DP), DIMENSION(2) :: startPoint
    REAL(DP), DIMENSION(3,2) :: startPoint2
    REAL(DP), DIMENSION(3) :: startVals
    REAL(DP) :: fret
    PROCEDURE(template_function), POINTER :: func
    PROCEDURE(template_derivative), POINTER :: deriv
    INTEGER :: i,j
    REAL(DP), DIMENSION(3) :: distance

    allocate(intervals(2,2))
    intervals(1,1)=-50.0D0
    intervals(1,2)=50.0D0
    intervals(2,1)=-20.0D0
    intervals(2,2)=20.0D0

    func => myfunction
    deriv => myderivative

    distance(1)=5.0D0
    distance(2)=10.0D0
    distance(3)=15.0D0

    startpoints(1,1)=-25.0D0
    startpoints(1,2)=25.0D0

    startpoints(2,1)=15.0D0
    startpoints(2,2)=15.0D0

    do i=1,2
        startPoint(1)=startpoints(i,1)
        startPoint(2)=startpoints(i,2)

        print *, "Iteration ",i

        print *, "dfpmin"
        print *, "Initial Point: ",startPoint, " Value: ",func(startPoint)
        CALL q1b(func, deriv, startPoint,fret)
        print *,"Min: ",startPoint," Value: ",fret

        do j=1,3

            startPoint2 = getSimplexAround(startPoint, distance(j))

            startVals(1)=func(startPoint2(1,:))
            startVals(2)=func(startPoint2(2,:))
            startVals(3)=func(startPoint2(3,:))

            print *, "Amoeba"
            print *, "Initial Simplex, Distance-",distance(j)
            print *, "1: ",startPoint2(1,:), " Value: ",startVals(1)
            print *, "2: ",startPoint2(2,:), " Value: ",startVals(2)
            print *, "3: ",startPoint2(3,:), " Value: ",startVals(3)

            CALL q1c(func, startPoint2, startVals)

            print *, "Final"
            print *, "1: ",startPoint2(1,:), " Value: ",startVals(1)
            print *, "2: ",startPoint2(2,:), " Value: ",startVals(2)
            print *, "3: ",startPoint2(3,:), " Value: ",startVals(3)
        end do
    end do

    startPoint2 = getSimplexAround(startpoints(1,:), 5.0D0)
    startVals(1)=func(startPoint2(1,:))
    startVals(2)=func(startPoint2(2,:))
    startVals(3)=func(startPoint2(3,:))

    !why this distance? It seems that the larger the size of the original simplex,
    !the better the algorithm does.
    CALL q1d(func, startPoint2, startVals,distance(3),deriv)

    deallocate(intervals)
CONTAINS

    !-----------------
    FUNCTION myfunction(point) RESULT(z)
        !
        ! the function we are trying to minimize
        ! INPUTS: point - the value at which we are evaluating the function
        ! OUTPUTS: y - the value of the function at that point
        REAL(DP), DIMENSION(:), INTENT(IN) :: point
        REAL(DP) :: x,y,z

        IF(size(point)/=2)THEN
            CALL sub_mystop("myfunction: incorrect dimensions in initial point.&
                & Should be two."            )
        END IF
        x=point(1)
        if(abs(x)<epsilon(1.0D0))then
            z=1000
        else
            y=point(2)
            z=5.0D0*sin(x)/x * max(20.0D0-abs(y),0.0D0) ** 1.2D0
        end if
    END FUNCTION myfunction

    FUNCTION myderivative(point) RESULT(z)
        !
        ! the derivative of the function we are trying to minimize
        ! INPUTS: point - the value at which we are evaluating the function
        ! OUTPUTS: y - the value of the partial derivatives of the function
        !              at that point
        REAL(DP), DIMENSION(:), INTENT(IN) :: point
        REAL(DP), DIMENSION(size(point,dim=1)) :: z
        REAL(DP):: x,y

        IF(size(point)/=2)THEN
            print *,size(point),point
            CALL sub_mystop("myderivative: incorrect dimensions in initial point.&
                & Should be two."            )
        END IF
        x=point(1)
        y=point(2)

        IF (abs(y)<20.0) THEN
            z(1)=5.0D0/x*(cos(x)-sin(x)/x)*(20-abs(y))**1.2D0
            z(2)=-(5.0D0*sin(x)/x)*(1.2D0*(20.0D0-abs(y))**0.2D0)
        ELSE
            z(1)=0.0D0
            z(2)=0.0D0
        END IF
    END FUNCTION myderivative

end program main
