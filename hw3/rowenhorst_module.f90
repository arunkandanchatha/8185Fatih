module rowenhorst_module
    implicit none
contains

    SUBROUTINE rouwenhorst(rho,sigma,p,y,s)
        ! rowenhorst method to approximate univariate AR(1) process by Markov chain
        !      y(t) = rho y(t-1)+ sigma sqrt(1-rho^2) e(t),   e(t)~N(0,1)
        ! INPUTS: rho - serial correlation coefficient,
        !         sigma - coefficient of variation
        ! OUTPUT: P is an n-by-n matrix of Markov transition probabilities
        !         y is an n-by-1 vector of symmetric and evenly-spaced Markov state space
        !         s is an n-by-1 vector of stationary distribution (binomial)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: rho, sigma
        REAL*8, DIMENSION(:,:), INTENT(OUT) :: p
        REAL*8, DIMENSION(:), INTENT(OUT) :: y,s
        REAL*8 :: ybar, q
        INTEGER :: n
        n=size(y)
        IF (size(p,dim=1)/=n .or. size(p,dim=2)/=n) THEN
            PRINT '(a,i3,a,i3)', 'rouwenhorst: p must be a square matrix of size ',n,' x ',n
            STOP 'program terminated by rouwenhorst'
        END IF
        IF (size(s)/=n) THEN
            PRINT '(a,i3)', 'rouwenhorst: y and s must be vectors of the same size ',n
            STOP 'program terminated by rouwenhorst'
        END IF
        ybar=sigma*sqrt(real(n-1,8))
        q=(1+rho)/2
        CALL rhmat(p)
        CALL grid(y,-ybar,ybar,1.0D0)
        CALL binom(s)
    CONTAINS
        RECURSIVE SUBROUTINE rhmat(p)
            IMPLICIT NONE
            REAL*8, DIMENSION(:,:), INTENT(OUT) :: p
            REAL*8, DIMENSION(size(p,dim=1)-1,size(p,dim=2)-1) :: p1
            INTEGER :: h
            h=size(p,dim=1)
            IF (size(p,dim=2)/=h) STOP 'P must be a square matrix'
            IF (h<2) STOP 'P must be at least 2-by-2 matrix'
            IF (h==2) THEN
                p=reshape((/q,1-q,1-q,q/),(/2,2/))
            ELSE
                CALL rhmat(p1)
                p=0
                p(1:h-1,1:h-1)=q*p1
                p(1:h-1,2:h)=(1-q)*p1+p(1:h-1,2:h)
                p(2:h,1:h-1)=(1-q)*p1+p(2:h,1:h-1)
                p(2:h,2:h)=q*p1+p(2:h,2:h)
                p(2:h-1,:)=p(2:h-1,:)/2
            END IF
        END SUBROUTINE rhmat
    END SUBROUTINE rouwenhorst

    SUBROUTINE binom(f)
        ! Binomial probability mass function with p=1/2
        IMPLICIT NONE
        REAL*8, DIMENSION(:), INTENT(OUT) :: f
        INTEGER :: n,k
        n=size(f)
        f(1)=2.0D0**(1-n)
        DO k=1,n-1
            f(k+1)=f(k)*(n-k)/k
        END DO
    END SUBROUTINE binom

    SUBROUTINE grid(x,xmin,xmax,s)
        ! Purpose: Generate grid x on [xmin,xmax] using spacing parameter s set as follows:
        ! s=1       linear spacing
        ! s>1       left skewed grid spacing with power s
        ! 0<s<1     right skewed grid spacing with power s
        ! s<0       geometric spacing with distances changing by a factor -s^(1/(n-1)), (>1 grow, <1 shrink)
        ! s=-1      logarithmic spacing with distances changing by a factor (xmax-xmin+1)^(1/(n-1))
        ! s=0       logarithmic spacing with distances changing by a factor (xmax/xmin)^(1/(n-1)), only if xmax,xmin>0
        IMPLICIT NONE
        REAL*8, DIMENSION(:), INTENT(OUT) :: x
        REAL*8, INTENT(IN) :: xmin,xmax,s
        REAL*8 :: c ! growth rate of grid subintervals for logarithmic spacing
        INTEGER :: n,i
        n=size(x)
        FORALL(i=1:n) x(i)=(i-1)/real(n-1,8)
        IF (s>0.0D0) THEN
            x=x**s*(xmax-xmin)+xmin
            IF (s==1.0D0) THEN
            !           PRINT '(a,i8,a,f6.3,a,f6.3,a)', 'Using ',n,' equally spaced grid points over domain [',xmin,',',xmax,']'
            ELSE
            !           PRINT '(a,i8,a,f6.3,a,f6.3,a,f6.3,a)', 'Using ',n,' skewed spaced grid points with power ',s,' over domain [',xmin,',',xmax,']'
            END IF
        ELSE
            IF (s==-1.0D0) THEN
                c=xmax-xmin+1
            !       ELSEIF (s==0.0D0) THEN
            !           IF (xmin>0.0D0) THEN
            !               c=xmax/xmin
            !           ELSE
            !               STOP 'grid: can not use logarithmic spacing for nonpositive values'
            !           END IF
            ELSE
                c=-s
            END IF
            !       PRINT '(a,i8,a,f6.3,a,f6.3,a,f6.3,a)', 'Using ',n,' logarithmically spaced grid points with growth rate ',c,' over domain [',xmin,',',xmax,']'
            x=((xmax-xmin)/(c-1))*(c**x)-((xmax-c*xmin)/(c-1))
        END IF
    END SUBROUTINE grid

end module rowenhorst_module
