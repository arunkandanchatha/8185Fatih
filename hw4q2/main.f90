module aiyagariSolve
    USE nrtype
    USE QSHEP2D

    implicit none

    REAL(DP), parameter                   ::  sigma=.4D0
    integer, parameter                  ::  n_a=200, n_z = 2, n_s = 2, n_k = 10
    integer, parameter                  :: myseed = 45678
    integer, parameter                  :: periodsForConv = 10001
    integer, parameter                  :: periodsToCut = 1000
    integer, parameter                  :: numHouseholds = 10000
    integer, parameter                  :: maxit = 500

    REAL(DP),parameter                   :: alpha=0.36
    REAL(DP), parameter                   ::  curv=2.0D0,a_min=0.001D0/numHouseholds, a_max=1000.0D0
    REAL(DP), parameter                   ::  k_min=30.0D0, k_max = 50.0D0
    REAL(DP), parameter                   ::  beta=0.90, delta = 0.025D0,c_gamma=1.0D0
    REAL(DP),parameter                    ::  criter_a=1D-9,criter_B=1D-9
    REAL(DP),parameter                  :: update_a = 0.7, update_B = 0.3
    REAL(DP),parameter                    ::  mu = 0.15 !unemp benefits

    REAL(DP), dimension(n_z)              ::  zShocks
    REAL(DP), dimension(n_k)              ::  k
    REAL(DP), dimension(n_z,n_s)   ::  ssEmployment
    REAL(DP), dimension(n_z,n_s,n_z,n_s)  ::  transition
    REAL(DP), dimension(n_a)              ::  a
    REAL(DP),dimension(n_z,2,2) :: phi
    TYPE(household), dimension(:,:), allocatable  ::  ssDistrib
    integer, dimension(periodsForConv) :: z

    PROCEDURE(template_function), POINTER :: funcParam
    PROCEDURE(template_function2), POINTER :: deriv1Func, deriv2Func
    PROCEDURE(template_function3), POINTER :: capitalCalc

    character(LEN=20)                   :: policyOutput
    character(LEN=20)                   :: distribOutput

contains
    subroutine setParams(d1Func, d2Func, capitalFunc,file1, file2)
        PROCEDURE(template_function2), POINTER, INTENT(IN) :: d1Func
        PROCEDURE(template_function2), POINTER, INTENT(IN) :: d2Func
        PROCEDURE(template_function3), POINTER, INTENT(IN) :: capitalFunc
        character(LEN=*),INTENT(IN) :: file1
        character(LEN=*),INTENT(IN) :: file2

        deriv1Func => d1Func
        deriv2Func => d2Func
        if(.not.associated(capitalFunc))then
            print *,"capital func in setParams is null."
            stop 0
        end if
        capitalCalc => capitalFunc
        if(.not.associated(capitalCalc,target=capitalFunc))then
            print *,"capitalCalc in setParams is null."
            stop 0
        end if

        policyOutput = file1
        distribOutput = file2
    end subroutine setParams

    subroutine myshocks()
        integer i,j,ii
        REAL(DP), dimension(n_z,n_z) :: ztransition
        REAL, dimension(periodsForConv, numHouseholds) :: harvest
        REAL, dimension(periodsForConv) :: aggShocks

        zTransition(1,1) = transition(1,1,1,1)+transition(1,1,1,2)
        zTransition(1,2) = 1-ztransition(1,1)
        zTransition(2,2) = transition(2,2,2,1)+transition(2,2,2,2)
        zTransition(2,1) = 1-ztransition(2,2)

        call srand(myseed)
        do j=1,numHouseholds
            do i=1,periodsForConv
                if(j==1)then
                    aggShocks(i) = rand(0)
                end if
                harvest(i,j) = rand(0)
            end do
        end do

        do i=1,numHouseholds
            if(harvest(1,i) < ssEmployment(1,1)) then
                ssDistrib(1,i)%employmentState = 1
            else
                ssDistrib(1,i)%employmentState = 2
            end if
        end do

        !z(i) is the shock received in period i
        z(1)=1

        !set states for all future periods
        do i=2,periodsForConv
            z(i)=1
            do j=1,n_z-1
                if(aggShocks(i-1)>sum(zTransition(z(i-1),1:j)))then
                    z(i)=j+1
                else
                    exit
                end if
            end do
        end do

        ! set employment for each household in each period
        do ii=1,numHouseholds
            print *,ii
            do j=2,periodsForConv
                ssDistrib(j,ii)%employmentState=1
                do i=1,n_s-1
                    if(harvest(j,ii)>sum(transition(z(j-1),ssDistrib(j-1,ii)%employmentState,z(j),1:i)))then
                        ssDistrib(j,ii)%employmentState=i+1
                    else
                        exit
                    end if
                end do
            end do
        end do
    end subroutine myshocks

    subroutine beginKrusellSmith()
        LOGICAL :: iterComplete
        INTEGER :: i,j,ii,jj,iter,iter2
        REAL(DP) :: diff_a, eps=epsilon(1.0_dp), kSS, incr, minc=0.0_dp
        PROCEDURE(template_function), POINTER :: func
        REAL(DP),dimension(n_a,n_k,n_z,n_s) :: aprime, prob_bu, prob_be, prob_gu, prob_ge
        REAL(DP),dimension(n_a,n_k,n_z,n_s) :: aaux,kaux,aglabor,agshock_aux,idshock_aux
        REAL(DP),dimension(n_a,n_k,n_z,n_s) :: irateaux,wageaux,wealth,kprime
        REAL(DP),dimension(n_a,n_k,n_z,n_s) :: irate_b, irate_g, wagerate_b, wagerate_g
        REAL(DP),dimension(n_a,n_k,n_z,n_s) :: a2prime_bu,cprime_bu,muprime_bu
        REAL(DP),dimension(n_a,n_k,n_z,n_s) :: a2prime_be,cprime_be,muprime_be
        REAL(DP),dimension(n_a,n_k,n_z,n_s) :: a2prime_gu,cprime_gu,muprime_gu
        REAL(DP),dimension(n_a,n_k,n_z,n_s) :: a2prime_ge,cprime_ge,muprime_ge
        REAL(DP),dimension(n_a,1,1,1) :: a2prime2_xu, a2prime2_xe, temp1, temp2
        REAL(DP),dimension(n_a,n_s) :: aprimet
        REAL(DP),dimension(n_a,n_k,n_z,n_s) :: expec, cn, aprimen, c
        REAL(DP),dimension(numHouseholds,1,1,1) :: temp1a
        REAL(DP),dimension(n_s,1,1,1) :: temp2a
        REAL(DP),dimension(numHouseholds) :: across
        REAL(DP),dimension(numHouseholds,1,1,1) :: acrossn
        REAL(DP),dimension(n_s) :: s
        REAL(DP),dimension(periodsForConv) :: kmts

        allocate(ssDistrib(periodsForConv,numHouseholds))

        !**************************************************************************
        ! We set up the grid of asset values based on the curvature, curv
        ! the minimum and maximum values in the grid a_min a_max
        !**************************************************************************
        incr=(a_max-a_min)/(n_a-1)
        a=(/ (    incr*real(i-1,8),i=1,n_a    ) /)
        a=a**curv
        a=a/(a_max-a_min)**(curv-1)+a_min

        ! Set aggregate capital grid over which we want to evaluate, K
        ! curving this allows convergence. Without, it doesn't. Damn I hate splines!
        incr=(k_max-k_min)/(n_k-1)
        k=(/ (    incr*real(i-1,8),i=1,n_k    ) /)
        k=k**curv
        k=k/(k_max-k_min)**(curv-1)+k_min

        ! Initial capital function
        forall(i=1:n_z,j=1:n_k,ii=1:n_s) aprime(:,j,i,ii)=0.9*a

        kSS=((1/beta-(1-delta))/alpha)**(1/(alpha-1))
        across(:)=kSS

        phi(1,1,:)= 0.0_dp
        phi(1,2,:)= 1.0_dp
        phi(2,1,:)= 0.0_dp
        phi(2,2,:)= 1.0_dp

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
        transition(2,2,1,1) = 0.009115D0
        transition(2,2,1,2) = 0.115885D0
        transition(2,2,2,1) = 0.024306D0
        transition(2,2,2,2) = 0.8506941D0

        ssEmployment=reshape((/0.1D0,0.04D0,0.90D0,0.96D0/),(/2,2/))

        call myshocks()

        zShocks=(/0.99D0,1.01D0/)
        s=(/0.0D0,1.0D0/)

        !next state probabilities
        prob_bu(:,:,1,1)=transition(1,1,1,1)
        prob_bu(:,:,1,2)=transition(1,2,1,1)
        prob_bu(:,:,2,1)=transition(2,1,1,1)
        prob_bu(:,:,2,2)=transition(2,2,1,1)

        prob_be(:,:,1,1)=transition(1,1,1,2)
        prob_be(:,:,1,2)=transition(1,2,1,2)
        prob_be(:,:,2,1)=transition(2,1,1,2)
        prob_be(:,:,2,2)=transition(2,2,1,2)

        prob_gu(:,:,1,1)=transition(1,1,2,1)
        prob_gu(:,:,1,2)=transition(1,2,2,1)
        prob_gu(:,:,2,1)=transition(2,1,2,1)
        prob_gu(:,:,2,2)=transition(2,2,2,1)

        prob_ge(:,:,1,1)=transition(1,1,2,2)
        prob_ge(:,:,1,2)=transition(1,2,2,2)
        prob_ge(:,:,2,1)=transition(2,1,2,2)
        prob_ge(:,:,2,2)=transition(2,2,2,2)

        forall(i=1:n_z,j=1:n_k,ii=1:n_s) aaux(:,j,i,ii)=a
        forall(i=1:n_z,j=1:n_a,ii=1:n_s) kaux(j,:,i,ii)=k
        forall(i=1:n_a,j=1:n_k,ii=1:n_s) aglabor(i,j,1,ii)=ssEmployment(1,2)
        forall(i=1:n_a,j=1:n_k,ii=1:n_s) aglabor(i,j,2,ii)=ssEmployment(2,2)
        forall(i=1:n_a,j=1:n_k,ii=1:n_s) agshock_aux(i,j,1,ii)=zShocks(1)
        forall(i=1:n_a,j=1:n_k,ii=1:n_s) agshock_aux(i,j,2,ii)=zShocks(2)
        forall(i=1:n_a,j=1:n_k,ii=1:n_s) idshock_aux(i,j,1,ii)=s(1)
        forall(i=1:n_a,j=1:n_k,ii=1:n_s) idshock_aux(i,j,2,ii)=s(2)

        forall(i=1:n_a,j=1:n_k,ii=1:n_z,jj=1:n_s)
            irateaux(i,j,ii,jj)=deriv1Func(k(j),ssEmployment(jj,2)/0.9_dp,zShocks(ii))
        end forall
        forall(i=1:n_a,j=1:n_k,ii=1:n_z,jj=1:n_s)
            wageaux(i,j,ii,jj)=deriv2Func(k(j),ssEmployment(jj,2)/0.9_dp,zShocks(ii))
        end forall
        forall(i=1:n_a,j=1:n_k,ii=1:n_z,jj=1:n_s)
            wealth(i,j,ii,jj)=irateaux(i,j,ii,jj)*kaux(i,j,ii,jj) + &
                wageaux(i,j,ii,jj)*idshock_aux(i,j,ii,jj)*0.9 + &
                mu*wageaux(i,j,ii,jj)*(1_dp-idshock_aux(i,j,ii,jj))
        end forall


        iter = 0
        iterComplete = .false.
        do while ( (iter<maxit) .and. (.not. iterComplete))
            iter = iter+1

            !individual stuff
            forall(i=1:n_a,j=1:n_k,ii=1:n_z,jj=1:n_s)
                kprime(i,j,ii,jj) = exp(phi(ii,1,1)+phi(ii,2,1)*log(kaux(i,j,ii,jj)))
            end forall

            kprime=min(max(k_min,kprime),k_max)

            forall(i=1:n_a,j=1:n_k,ii=1:n_z,jj=1:n_s)
                irate_b(i,j,ii,jj)=deriv1Func(kprime(i,j,ii,jj),&
                    ssEmployment(1,2)/0.9_dp,zShocks(1))
                irate_g(i,j,ii,jj)=deriv1Func(kprime(i,j,ii,jj),&
                    ssEmployment(2,2)/0.9_dp,zShocks(2))
                wagerate_b(i,j,ii,jj)=deriv2Func(kprime(i,j,ii,jj),&
                    ssEmployment(1,2)/0.9_dp,zShocks(1))
                wagerate_g(i,j,ii,jj)=deriv2Func(kprime(i,j,ii,jj),&
                    ssEmployment(2,2)/0.9_dp,zShocks(2))
            end forall

            diff_a=1
            iter2=1
            do while ((diff_a>criter_a) .and. (iter2<maxit))
                iter2 = iter2+1
                a2prime_bu=interpolate(n_a,a,n_k,k,aprime(:,:,1,1),n_a,n_k,n_z,n_s,aprime,kprime)
                a2prime_be=interpolate(n_a,a,n_k,k,aprime(:,:,1,2),n_a,n_k,n_z,n_s,aprime,kprime)
                a2prime_gu=interpolate(n_a,a,n_k,k,aprime(:,:,2,1),n_a,n_k,n_z,n_s,aprime,kprime)
                a2prime_ge=interpolate(n_a,a,n_k,k,aprime(:,:,2,2),n_a,n_k,n_z,n_s,aprime,kprime)
                forall(i=1:n_a,j=1:n_k,ii=1:n_z,jj=1:n_s)
                    cprime_bu(i,j,ii,jj) = max(irate_b(i,j,ii,jj)*aprime(i,j,ii,jj)+&
                        mu*wagerate_b(i,j,ii,jj)+(1-delta)*aprime(i,j,ii,jj)-&
                        a2prime_bu(i,j,ii,jj),minc+eps)
                    cprime_be(i,j,ii,jj) = max(irate_b(i,j,ii,jj)*aprime(i,j,ii,jj)+&
                        wagerate_b(i,j,ii,jj)/0.9+(1-delta)*aprime(i,j,ii,jj)-&
                        mu*wagerate_b(i,j,ii,jj)*(ssEmployment(1,1)/(1-ssEmployment(1,1)))&
                        -a2prime_be(i,j,ii,jj),minc+eps)
                    cprime_gu(i,j,ii,jj) = max(irate_g(i,j,ii,jj)*aprime(i,j,ii,jj)+&
                        mu*wagerate_g(i,j,ii,jj)+(1-delta)*aprime(i,j,ii,jj)-&
                        a2prime_gu(i,j,ii,jj),minc+eps)
                    cprime_ge(i,j,ii,jj) = max(irate_g(i,j,ii,jj)*aprime(i,j,ii,jj)+&
                        wagerate_g(i,j,ii,jj)/0.9+(1-delta)*aprime(i,j,ii,jj)-&
                        mu*wagerate_b(i,j,ii,jj)*(ssEmployment(1,1)/(1-ssEmployment(1,1)))&
                        -a2prime_ge(i,j,ii,jj),minc+eps)
                end forall
                forall(i=1:n_a,j=1:n_k,ii=1:n_z,jj=1:n_s)
                    muprime_bu(i,j,ii,jj) = (cprime_bu(i,j,ii,jj) - minc)**(-c_gamma)
                    muprime_be(i,j,ii,jj) = (cprime_be(i,j,ii,jj) - minc)**(-c_gamma)
                    muprime_gu(i,j,ii,jj) = (cprime_gu(i,j,ii,jj) - minc)**(-c_gamma)
                    muprime_ge(i,j,ii,jj) = (cprime_ge(i,j,ii,jj) - minc)**(-c_gamma)
                end forall

                forall(i=1:n_a,j=1:n_k,ii=1:n_z,jj=1:n_s)
                    expec(i,j,ii,jj) = &
                        muprime_bu(i,j,ii,jj)*(1-delta+irate_b(i,j,ii,jj))*prob_bu(i,j,ii,jj) +&
                        muprime_be(i,j,ii,jj)*(1-delta+irate_b(i,j,ii,jj))*prob_be(i,j,ii,jj) +&
                        muprime_gu(i,j,ii,jj)*(1-delta+irate_g(i,j,ii,jj))*prob_gu(i,j,ii,jj) +&
                        muprime_ge(i,j,ii,jj)*(1-delta+irate_g(i,j,ii,jj))*prob_ge(i,j,ii,jj)
                end forall

                forall(i=1:n_a,j=1:n_k,ii=1:n_z,jj=1:n_s)
                    cn(i,j,ii,jj) = (beta*expec(i,j,ii,jj))**(-1/c_gamma)+minc
                end forall

                aprimen = min(max(wealth-cn,k_min),k_max)

                diff_a = maxval(abs(aprimen-aprime))

                print *,diff_a

                aprime = update_a*aprimen + (1-update_a)*aprime
            end do
            c=wealth-kprime

            !aggregate stuff
            do iter2=1,periodsForConv
                kmts(iter2) = sum(across)/numHouseholds
                kmts(iter2) = min(max(kmts(iter2),k_min),k_max)

                temp1(:,1,1,1) = a
                temp2(:,1,1,1) = kmts(iter2)

                a2prime2_xu=interpolate(n_a,a,n_k,k,aprime(:,:,z(iter2),1),n_a,1,1,1,temp1,temp2)
                a2prime2_xe=interpolate(n_a,a,n_k,k,aprime(:,:,z(iter2),2),n_a,1,1,1,temp1,temp2)

                aprimet(:,1) = a2prime2_xu(:,1,1,1)
                aprimet(:,2) = a2prime2_xe(:,1,1,1)

                temp1a(:,1,1,1) = across
                temp2a(:,1,1,1) = (/1.0D0,2.0D0/)
                acrossn=interpolate(n_a,a,n_s,temp2a,aprimet,numHouseholds,1,1,1,temp1a,&
                    DBLE(ssDistrib(iter2,:)%employmentState))

                acrossn=min(max(acrossn,a_min),a_max)
                across=acrossn(:,1,1,1)

                if(mod(iter2,10)==0)then
                    print *,iter2
                end if
            end do

            if(maxval(abs(phi(:,:,2)-phi(:,:,1)))<criter_B)then
                iterComplete = .true.
            end if
        end do

        deallocate(ssDistrib)

    end subroutine beginKrusellSmith

    function interpolate(l1, dim1, l2, dim2, fn, s1, s2, s3, s4, xpoints, ypoints) RESULT (interps)
        integer, intent(in) :: l1, l2, s1, s2, s3, s4
        REAL(DP), dimension(l1), intent(IN) :: dim1
        REAL(DP), dimension(l2), intent(IN) :: dim2
        REAL(DP), dimension(l1,l2), intent(IN) :: fn
        REAL(DP), dimension(s1,s2,s3,s4), intent(IN) :: xpoints,ypoints
        REAL(DP), dimension(s1,s2,s3,s4) :: interps

        integer :: n
        REAL(DP), dimension(:), allocatable :: x, y, f, rsq
        integer, dimension(:), allocatable :: lnext
        REAL(DP) :: xmin, ymin, dx, dy, rmax
        integer, dimension(:,:), allocatable :: lcell
        REAL(DP), dimension(:,:), allocatable :: a

        integer :: i,j, ii, jj, ier

        n=l1*l2
        allocate(x(n),y(n),f(n),rsq(n),lcell(n,n),a(5,n))

        do i = 1,l1
            do j = 1,l2
                x((i-1)*n_k+j)=dim1(i)
                y((i-1)*n_k+j)=dim2(j)
                f((i-1)*n_k+j)=fn(i,j)
            end do
        end do

        call qshep2 ( n, x, y, f, 13, 19, ceiling(sqrt(DBLE(n)/3_dp)), lcell, lnext, xmin, ymin, &
            dx, dy, rmax, rsq, a, ier )

        if ( ier /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'interpolate - Error!'
            write ( *, '(a,i6)' ) '  Error in interpoloate, IER = ', ier
            stop
        end if

        do jj=1,s4
            do ii=1,s3
                do j=1,s2
                    do i=1,s1
                        interps(i,j,ii,jj) = qs2val ( xpoints(i,j,ii,jj), ypoints(i,j,ii,jj), n, x, y, f, 13, lcell, lnext, xmin, &
                            ymin, dx, dy, rmax, rsq, a )
                    end do
                end do
            end do
        end do

        deallocate(x,y,f,rsq,lcell,a)

    end function interpolate

end module aiyagariSolve

!**************************************************************************
!
!**************************************************************************
program main

    use aiyagariSolve

    PROCEDURE(template_function), POINTER,save :: func
    PROCEDURE(template_function2), POINTER,save :: d1func, d2func
    PROCEDURE(template_function3), POINTER,save :: func2a
    character(LEN=15) :: arg1,arg2

    !************
    ! Timing variables
    !************
    real(DP) :: startTime, endTime

    flush(6)
    d1func => d1prod
    d2func => d2prod
    func2a => impliedCapital
    arg1="policy"
    arg2= "distrib"
    call setParams(d1func, d2func, func2a, arg1,arg2)
    call CPU_TIME(startTime)
    call beginKrusellSmith()
    call CPU_TIME(endTime)
    print *,"Time: ",endTime-startTime
    flush(6)
contains
    function production(capital,labour,shock) RESULT(y)
        use nrtype
        implicit none
        REAL(DP), INTENT(IN) :: capital, labour,shock
        REAL(DP) :: y
        y=shock*capital**alpha*labour**(1-alpha)
    end function production

    pure function d1prod(capital,labour,shock) RESULT(y)
        use nrtype
        implicit none
        REAL(DP), INTENT(IN) :: capital, labour,shock
        REAL(DP) :: y
        y=alpha*shock*capital**(alpha-1)*labour**(1-alpha)
    end function d1prod

    pure function d2prod(capital,labour,shock) RESULT(y)
        use nrtype
        implicit none
        REAL(DP), INTENT(IN) :: capital, labour,shock
        REAL(DP) :: y
        y=(1-alpha)*shock*capital**alpha*labour**(-alpha)
    end function d2prod

    function impliedCapital(interest,labour,shock) RESULT(y)
        use nrtype
        implicit none
        REAL(DP), INTENT(IN) :: interest, labour, shock
        REAL(DP) :: y
        y = interest/(alpha*shock)*&
            1.0D0/labour**(1-alpha)
        y = y**(1.0D0/(alpha-1))
    end function impliedCapital

end program main
