module aiyagariSolve
    USE nrtype
    USE QSHEP2D
    USE lib_array, only : interp2d, locate
    implicit none

    integer, parameter                  ::  n_a=200, n_z = 2, n_s = 2, n_k = 4
    integer, parameter                  :: myseed = 4567890
    integer, parameter                  :: periodsForConv = 6000
    integer, parameter                  :: periodsToCut = 1000
    integer, parameter                  :: numHouseholds = 10000
    integer, parameter                  :: maxit = 4000

    REAL(DP),parameter                   :: alpha=0.36D0
    REAL(DP), parameter                   ::  curv=7.0D0,a_min=0.0D0, a_max=1000.0D0
    REAL(DP), parameter                   ::  km_min=30.0D0, km_max = 50.0D0
    REAL(DP), parameter                   ::  beta=0.99_dp, delta = 0.025D0,c_gamma=1.0D0
    REAL(DP),parameter                    ::  criter_a=1D-8,criter_B=1D-8
    REAL(DP),parameter                  :: update_a = 0.7_dp, lambda = 0.3_dp
    REAL(DP),parameter                    ::  mu = 0.0_dp !unemp benefits

    REAL(DP), dimension(n_z)              ::  zShocks
    REAL(DP), dimension(n_k)              ::  k
    REAL(DP), dimension(n_z,n_s)   ::  ssEmployment
    REAL(DP), dimension(n_z,n_s,n_z,n_s)  ::  transition
    REAL(DP), dimension(n_a)              ::  a
    REAL(DP),dimension(n_z,2,2) :: phi
    TYPE(household), dimension(:,:), allocatable  ::  ssDistrib
    integer, dimension(periodsForConv) :: z
    REAL(DP)                   :: l_bar

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
        integer i,j,ii, n
        REAL(DP), dimension(n_z,n_z) :: ztransition
        REAL, dimension(periodsForConv, numHouseholds) :: harvest
        REAL, dimension(periodsForConv) :: aggShocks
        INTEGER, dimension(:), allocatable :: seed

        zTransition(1,1) = transition(1,1,1,1)+transition(1,1,1,2)
        zTransition(1,2) = 1-ztransition(1,1)
        zTransition(2,2) = transition(2,2,2,1)+transition(2,2,2,2)
        zTransition(2,1) = 1-ztransition(2,2)

        call random_seed(size = n)
        allocate(seed(n))
        seed(1)=myseed
        call random_seed(put = seed)
        call random_number(harvest)
        call random_number(aggShocks)
        deallocate(seed)

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
        INTEGER :: i,j,ii,jj,iter,iter2, interChoice
        REAL(DP) :: diff_a, eps=10D-10, kSS, incr, minc=0.1_dp, diff_b
        REAL(DP) :: ssErrG, ssTotG, ssErrB, ssTotB
        PROCEDURE(template_function), POINTER :: func
        REAL(DP),dimension(n_a,n_k,n_z,n_s) :: aprime, prob_bu, prob_be, prob_gu, prob_ge
        REAL(DP),dimension(n_a,n_k,n_z,n_s) :: aaux,kmaux,aglabor,agshock_aux,idshock_aux
        REAL(DP),dimension(n_a,n_k,n_z,n_s) :: irateaux,wageaux,wealth,kmprime
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
        REAL(DP),dimension(numHouseholds) :: across, acrosstemp
        REAL(DP),dimension(numHouseholds,1,1,1) :: acrossn
        REAL(DP),dimension(n_s) :: s
        REAL(DP),dimension(periodsForConv) :: kmts, gini
        REAL(DP),dimension(:,:),allocatable :: goodShocks,badShocks
        REAL(DP),dimension(2*periodsForConv*2) :: workArray

        INTEGER, dimension(4) :: maxError

        allocate(ssDistrib(periodsForConv,numHouseholds))

        !**************************************************************************
        ! We set up the grid of asset values based on the curvature, curv
        ! the minimum and maximum values in the grid a_min a_max
        !**************************************************************************
        incr=DBLE(0.5_dp/(n_a-1))
        a=(/ (    incr*real(i-1,8),i=1,n_a    ) /)
        a = a**curv
        a=a/maxval(a)
        a=a_min+(a_max-a_min)*a

        ! Set aggregate capital grid over which we want to evaluate, K
        ! curving this allows convergence. Without, it doesn't. Damn I hate splines!
        incr=(km_max-km_min)/(n_k-1)
        k=(/ (    incr*real(i-1,8),i=1,n_k    ) /)
        !        k=k**2
        !        k=k/(km_max-km_min)**(1)+km_min
        k=k+km_min
        ! Initial capital function
        forall(i=1:n_z,j=1:n_k,ii=1:n_s) aprime(:,j,i,ii)=0.9_dp*a

        kSS=((1.0_dp/beta-(1.0_dp-delta))/alpha)**(1.0_dp/(alpha-1.0_dp))
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
        transition(2,2,2,2) = 0.850694D0

        ssEmployment=reshape((/0.1D0,0.04D0,0.90D0,0.96D0/),(/2,2/))
        l_bar = 1.0_dp/ssEmployment(1,2)

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
        forall(i=1:n_z,j=1:n_a,ii=1:n_s) kmaux(j,:,i,ii)=k
        forall(i=1:n_a,j=1:n_k,ii=1:n_s) aglabor(i,j,1,ii)=ssEmployment(1,2)
        forall(i=1:n_a,j=1:n_k,ii=1:n_s) aglabor(i,j,2,ii)=ssEmployment(2,2)
        forall(i=1:n_a,j=1:n_k,ii=1:n_s) agshock_aux(i,j,1,ii)=zShocks(1)
        forall(i=1:n_a,j=1:n_k,ii=1:n_s) agshock_aux(i,j,2,ii)=zShocks(2)
        forall(i=1:n_a,j=1:n_k,ii=1:n_z) idshock_aux(i,j,ii,1)=s(1)
        forall(i=1:n_a,j=1:n_k,ii=1:n_z) idshock_aux(i,j,ii,2)=s(2)

        forall(i=1:n_a,j=1:n_k,ii=1:n_z,jj=1:n_s)
            irateaux(i,j,ii,jj)=deriv1Func(kmaux(i,j,ii,jj),ssEmployment(ii,2)*l_bar,zShocks(ii))
        end forall
        forall(i=1:n_a,j=1:n_k,ii=1:n_z,jj=1:n_s)
            wageaux(i,j,ii,jj)=deriv2Func(kmaux(i,j,ii,jj),ssEmployment(ii,2)*l_bar,zShocks(ii))
        end forall
        forall(i=1:n_a,j=1:n_k,ii=1:n_z,jj=1:n_s)
            wealth(i,j,ii,jj)=irateaux(i,j,ii,jj)*aaux(i,j,ii,jj) + &
                wageaux(i,j,ii,jj)*idshock_aux(i,j,ii,jj)*l_bar + &
                mu*wageaux(i,j,ii,jj)*(1.0_dp-idshock_aux(i,j,ii,jj)) + &
                (1.0_dp-delta)*aaux(i,j,ii,jj)-&
                mu*wageaux(i,j,ii,jj)*(1.0_dp-aglabor(i,j,ii,jj))/aglabor(i,j,ii,jj)*idshock_aux(i,j,ii,jj)
        end forall

        iter = 0
        iterComplete = .false.
        print *,"        iter               diff"
        do while ( (iter<maxit) .and. (.not. iterComplete))
            iter = iter+1

            !individual stuff
            forall(i=1:n_a,j=1:n_k,ii=1:n_z,jj=1:n_s)
                kmprime(i,j,ii,jj) = exp(phi(ii,1,1)+phi(ii,2,1)*log(kmaux(i,j,ii,jj)))
            end forall

            kmprime=min(max(km_min,kmprime),km_max)

            forall(i=1:n_a,j=1:n_k,ii=1:n_z,jj=1:n_s)
                irate_b(i,j,ii,jj)=deriv1Func(kmprime(i,j,ii,jj), ssEmployment(1,2)*l_bar,zShocks(1))
                irate_g(i,j,ii,jj)=deriv1Func(kmprime(i,j,ii,jj),&
                    ssEmployment(2,2)*l_bar,zShocks(2))
                wagerate_b(i,j,ii,jj)=deriv2Func(kmprime(i,j,ii,jj),&
                    ssEmployment(1,2)*l_bar,zShocks(1))
                wagerate_g(i,j,ii,jj)=deriv2Func(kmprime(i,j,ii,jj),&
                    ssEmployment(2,2)*l_bar,zShocks(2))
            end forall

            diff_a=1
            iter2=0
            do while ((diff_a>criter_a) .and. (iter2<maxit))
                iter2 = iter2+1

                interChoice=5
                a2prime_bu=interpolate(interChoice,n_a,a,n_k,k,aprime(:,:,1,1),n_a,n_k,n_z,n_s,aprime,kmprime)
                a2prime_be=interpolate(interChoice,n_a,a,n_k,k,aprime(:,:,1,2),n_a,n_k,n_z,n_s,aprime,kmprime)
                a2prime_gu=interpolate(interChoice,n_a,a,n_k,k,aprime(:,:,2,1),n_a,n_k,n_z,n_s,aprime,kmprime)
                a2prime_ge=interpolate(interChoice,n_a,a,n_k,k,aprime(:,:,2,2),n_a,n_k,n_z,n_s,aprime,kmprime)

                forall(i=1:n_a,j=1:n_k,ii=1:n_z,jj=1:n_s)
                    cprime_bu(i,j,ii,jj) = max(irate_b(i,j,ii,jj)*aprime(i,j,ii,jj)+&
                        mu*wagerate_b(i,j,ii,jj)+(1.0_dp-delta)*aprime(i,j,ii,jj)-&
                        a2prime_bu(i,j,ii,jj),minc+eps)
                    cprime_be(i,j,ii,jj) = max(irate_b(i,j,ii,jj)*aprime(i,j,ii,jj)+&
                        wagerate_b(i,j,ii,jj)*l_bar+(1.0_dp-delta)*aprime(i,j,ii,jj)-&
                        mu*wagerate_b(i,j,ii,jj)*(ssEmployment(1,1)/(1.0_dp-ssEmployment(1,1)))&
                        -a2prime_be(i,j,ii,jj),minc+eps)
                    cprime_gu(i,j,ii,jj) = max(irate_g(i,j,ii,jj)*aprime(i,j,ii,jj)+&
                        mu*wagerate_g(i,j,ii,jj)+(1.0_dp-delta)*aprime(i,j,ii,jj)-&
                        a2prime_gu(i,j,ii,jj),minc+eps)
                    cprime_ge(i,j,ii,jj) = max(irate_g(i,j,ii,jj)*aprime(i,j,ii,jj)+&
                        wagerate_g(i,j,ii,jj)*l_bar+(1.0_dp-delta)*aprime(i,j,ii,jj)-&
                        mu*wagerate_g(i,j,ii,jj)*(ssEmployment(2,1)/(1.0_dp-ssEmployment(2,1)))&
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
                        muprime_bu(i,j,ii,jj)*(1.0_dp-delta+irate_b(i,j,ii,jj))*prob_bu(i,j,ii,jj) +&
                        muprime_be(i,j,ii,jj)*(1.0_dp-delta+irate_b(i,j,ii,jj))*prob_be(i,j,ii,jj) +&
                        muprime_gu(i,j,ii,jj)*(1.0_dp-delta+irate_g(i,j,ii,jj))*prob_gu(i,j,ii,jj) +&
                        muprime_ge(i,j,ii,jj)*(1.0_dp-delta+irate_g(i,j,ii,jj))*prob_ge(i,j,ii,jj)
                end forall

                forall(i=1:n_a,j=1:n_k,ii=1:n_z,jj=1:n_s)
                    cn(i,j,ii,jj) = (beta*expec(i,j,ii,jj))**(-1.0_dp/c_gamma)+minc
                end forall

                aprimen = min(max(wealth-cn,a_min),a_max)

                diff_a = maxval(abs(aprimen-aprime))
                maxError = maxloc(abs(aprimen-aprime))

                aprime = update_a*aprimen + (1.0_dp-update_a)*aprime

            end do

            c=wealth-aprime

            !aggregate stuff
            acrosstemp=across
            do iter2=1,periodsForConv
                kmts(iter2) = sum(across)/numHouseholds
                kmts(iter2) = min(max(kmts(iter2),km_min),km_max)

                temp1(:,1,1,1) = a
                temp2(:,1,1,1) = kmts(iter2)

                interChoice=5

                a2prime2_xu=interpolate(interChoice,n_a,a,n_k,k,aprime(:,:,z(iter2),1),n_a,1,1,1,temp1,temp2)
                a2prime2_xe=interpolate(interChoice,n_a,a,n_k,k,aprime(:,:,z(iter2),2),n_a,1,1,1,temp1,temp2)

                aprimet(:,1) = a2prime2_xu(:,1,1,1)
                aprimet(:,2) = a2prime2_xe(:,1,1,1)

                temp1a(:,1,1,1) = across
                temp2a(:,1,1,1) = (/1.0D0,2.0D0/)
                acrossn=interpolate(interChoice,n_a,a,n_s,temp2a,aprimet,numHouseholds,1,1,1,temp1a,&
                    DBLE(ssDistrib(iter2,:)%employmentState))

                acrossn=min(max(acrossn,a_min),a_max)
                across=acrossn(:,1,1,1)
                ssDistrib(iter2,:)%capital=across

            end do

            !****************************************************
            ! Split the data into two sets, one for good shocks and
            ! one for bad shocks
            !****************************************************
            j=0
            do i=periodsToCut+1,periodsForConv-1
                if(z(i)>1)then
                    j=j+1
                end if
            end do

            allocate(goodShocks(j,3))
            allocate(badShocks(periodsForConv-j-periodsToCut-1,3))

            j=0
            ii=0
            do i=periodsToCut+1,periodsForConv-1
                if(z(i)>1)then
                    j=j+1
                    goodShocks(j,1)=1.0D0
                    goodShocks(j,2)=log(kmts(i))
                    goodShocks(j,3)=log(kmts(i+1))
                else
                    ii=ii+1
                    badShocks(ii,1)=1.0D0
                    badShocks(ii,2)=log(kmts(i))
                    badShocks(ii,3)=log(kmts(i+1))
                end if
            end do

            !**************************************************
            ! Find phi(0) and phi(1) using OLS
            !**************************************************

            !First for good shocks, which are state 2
            j=size(goodShocks,dim=1)
            incr=sum(goodShocks(:,3))/max(1,size(goodShocks(:,3)))
            ssTotG = sum((goodShocks(:,3)-incr)**2)
            call dgels('N', j, 2, 1, goodShocks(:,1:2), j, goodShocks(:,3),j, workArray,size(workArray),i)
            if(i/=0)then
                print *,"error regressing good shocks.",i
                stop 0
            end if
            ssErrG = sum(goodShocks(3:,3)**2)
            phi(2,1,2)=goodShocks(1,3)
            phi(2,2,2)=goodShocks(2,3)


            !And for bad shocks, which are state 1
            j=size(badShocks,dim=1)
            incr=sum(badShocks(:,3))/max(1,size(badShocks(:,3)))
            ssTotB = sum((badShocks(:,3)-incr)**2)
            call dgels('N', j, 2, 1, badShocks(:,1:2), j, badShocks(:,3),j, workArray,size(workArray),i)
            if(i/=0)then
                print *,"error regressing bad shocks.",i
                stop 0
            end if
            ssErrB = sum(badShocks(3:,3)**2)
            phi(1,1,2)=badShocks(1,3)
            phi(1,2,2)=badShocks(2,3)

            diff_B=maxval(abs(phi(:,:,2)-phi(:,:,1)))

            ! we use the terminal distribution of the current iteration as initial
            ! distribution for a subsequent iteration. When the solution is sufficiently
            ! accurate, dif_B<(criter_B*100), we stop such an updating and hold the
            ! distribution "kcross" fixed for the rest of iterations. ·
            if (diff_B<(criter_B*100)) then
                across=acrosstemp;             ! don't update
            end if

            if (mod(iter,1)==0)then
                print *,iter,maxval(abs(phi(:,:,2)-phi(:,:,1)))
                flush(6)
            end if
            phi(:,:,1)=lambda*phi(:,:,2)+(1.0D0-lambda)*phi(:,:,1)

            if(diff_B<criter_B)then
                iterComplete = .true.
            end if

            deallocate(goodShocks, badShocks)
        end do

        print *,"G:",phi(2,:,1),1.0_dp-ssErrG/ssTotG
        print *,"B:",phi(1,:,1),1.0_dp-ssErrB/ssTotB

        call calculateGini(ssDistrib,gini)
        print *,gini,z
        deallocate(ssDistrib)

    end subroutine beginKrusellSmith

    subroutine calculateGini(distrib, ginis)
        TYPE(household), dimension(periodsForConv,numHouseholds), INTENT(IN) ::  distrib
        REAL(DP), dimension(periodsForConv), INTENT(OUT) :: ginis
        INTEGER, dimension(periodsForConv,n_a) :: periodDistrib
        REAL(DP), dimension(periodsForConv, n_a) :: statDist, cdf, percWealth
        REAL(DP), dimension(periodsForConv) :: sn
        INTEGER :: i,j ! iterators
        INTEGER :: ind

        periodDistrib = 0
        do j=1,numHouseholds
            do i=1,periodsForConv
                ind=locate(a,distrib(i,j)%capital)
                periodDistrib(i,ind) = periodDistrib(i,ind) + 1
            end do
        end do

        statDist = periodDistrib/DBLE(numHouseholds)
        forall(i=1:periodsForConv) sn(i) = sum(statDist(i,:)*a)
        forall(i=1:periodsForConv,j=1:n_a) cdf(i,j) = sum(statDist(i,1:j))
        forall(i=1:periodsForConv,j=1:n_a) percWealth(i,j) = sum(statDist(i,1:j)*a(1:i))/sn(i)
        forall(i=1:periodsForConv) ginis(i) = sum(cdf(i,:)-percWealth(i,:))/sum(cdf(i,:))

#if 1
        print *, statDist(1,:)
        print *,"=============================="
        print *, sum(statDist(1,:))
        print *,"=============================="
        print *, sn(1)
        print *,"=============================="
        print *, cdf(1,:)
        print *,"=============================="
        print *, percWealth(1,:)
        print *,"=============================="
        print *, ginis(1)
        stop 0
#endif
    end subroutine calculateGini

    function interpolate(which, l1, dim1, l2, dim2, fn, s1, s2, s3, s4, xpoints, ypoints) RESULT (interps)
        integer, intent(in) :: which, l1, l2, s1, s2, s3, s4
        REAL(DP), dimension(l1), intent(IN) :: dim1
        REAL(DP), dimension(l2), intent(IN) :: dim2
        REAL(DP), dimension(l1,l2), intent(IN) :: fn
        REAL(DP), dimension(s1,s2,s3,s4), intent(IN) :: xpoints,ypoints
        REAL(DP), dimension(s1,s2,s3,s4) :: interps

        integer :: n, n2
        REAL(DP), dimension(:), allocatable :: x, y, f, rsq, x2, y2, f2
        integer, dimension(:), allocatable :: lnext
        REAL(DP) :: xmin, ymin, dx, dy, rmax
        integer, dimension(:,:), allocatable :: lcell
        REAL(DP), dimension(:,:), allocatable :: a
        REAL(DP), dimension(:,:,:), allocatable :: work

        integer :: i,j, ii, jj, ier

        if(which == 1) then
            n=l1*l2
            allocate(x(n),y(n),f(n))

            do i = 1,l1
                do j = 1,l2
                    x((i-1)*l2+j)=dim1(i)
                    y((i-1)*l2+j)=dim2(j)
                    f((i-1)*l2+j)=fn(i,j)
                end do
            end do


            allocate(rsq(n),lnext(n),lcell(n,n),a(5,n))
            call qshep2 ( n, x, y, f, 13, 19, ceiling(sqrt(DBLE(n)/3.0_dp)), lcell, lnext, xmin, ymin, &
                dx, dy, rmax, rsq, a, ier )

            if ( ier /= 0 ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'interpolate - Error!'
                write ( *, '(a,i6)' ) '  Error in qshep2, IER = ', ier
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

            deallocate(rsq,lnext,lcell,a)
            deallocate(x,y,f)

        else if(which == 2) then

            n2 = s1*s2*s3*s4
            allocate(x2(n2),y2(n2),f2(n2),work(3,l1,l2))

            ier=1
            do jj = 1,s4
                do ii = 1,s3
                    do j = 1,s2
                        do i = 1,s1
                            x2(ier)=xpoints(i,j,ii,jj)
                            y2(ier)=ypoints(i,j,ii,jj)
                            ier = ier+1
                        end do
                    end do
                end do
            end do

            CALL RGBI3P(1,l1,l2,dim1,dim2,fn,n2,x2,y2,f2,ier,work)

            if ( ier /= 0 ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'interpolate - Error!'
                write ( *, '(a,i6)' ) '  Error in SDBI3P, IER = ', ier
                stop
            end if

            ier=1
            do jj = 1,s4
                do ii = 1,s3
                    do j = 1,s2
                        do i = 1,s1
                            interps(i,j,ii,jj)=f2(ier)
                            ier = ier+1
                        end do
                    end do
                end do
            end do

            deallocate(x2,y2,f2,work)
        else if (which == 3) then

            allocate(work(3,l1,l2))
            do jj = 1,s4
                do ii = 1,s3
                    do j = 1,s2
                        CALL RGSF3P(1,l1,l2,dim1,dim2,fn,s1,xpoints(:,j,ii,jj),1,ypoints(1,j,ii,jj),&
                            interps(:,j,ii,jj), ier, work)
                    end do
                end do
            end do
            if ( ier /= 0 ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'interpolate - Error!'
                write ( *, '(a,i6)' ) '  Error in SDSF3P, IER = ', ier
                stop
            end if
            deallocate(work)
        else if (which == 4) then
            n2 = s1*s2*s3*s4
            allocate(x2(n2),y2(n2),f2(n2))

            ier=1
            do jj = 1,s4
                do ii = 1,s3
                    do j = 1,s2
                        do i = 1,s1
                            x2(ier)=xpoints(i,j,ii,jj)
                            y2(ier)=ypoints(i,j,ii,jj)
                            ier = ier+1
                        end do
                    end do
                end do
            end do

            call pwl_interp_2d ( l1, l2, dim1, dim2, fn, n2, x2, y2, f2 )
            ier=1
            do jj = 1,s4
                do ii = 1,s3
                    do j = 1,s2
                        do i = 1,s1
                            interps(i,j,ii,jj)=f2(ier)
                            ier = ier+1
                        end do
                    end do
                end do
            end do

        else if(which == 5) then

            ier=1
            do jj = 1,s4
                do ii = 1,s3
                    do j = 1,s2
                        do i = 1,s1
                            interps(i,j,ii,jj)=interp2d(dim1,dim2,fn,xpoints(i,j,ii,jj),ypoints(i,j,ii,jj),bounds_error=.true.)
                            ier = ier+1
                        end do
                    end do
                end do
            end do

        else
            print *, "error. Unknown interpolation type: ",which
            stop 0
        end if
    end function interpolate

    subroutine print4DArray(a)
        REAL(DP), dimension(:,:,:,:), INTENT(IN) :: a
        INTEGER :: s1, s2, s3, s4
        INTEGER :: s,t,u,v

        s1 = size(a,1)
        s2 = size(a,2)
        s3 = size(a,3)
        s4 = size(a,4)

        do v = 1,s4
            do u = 1,s3
                do t = 1,s2
                    do s = 1,s1
                        print *,s,t,u,v,a(s,t,u,v)
                    end do
                end do
            end do
        end do

        print *,"================================"
    end subroutine print4DArray

end module aiyagariSolve

!**************************************************************************
! This code is essentially a copy of the Maliar and Maliar (2010) code but
! written in fortran.
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
