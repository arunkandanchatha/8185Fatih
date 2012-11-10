program main
    implicit none

    call policy(0.1D0,"q7del0p1.csv")
    call policy(0.0D0,"q7del0p0.csv")

contains
    subroutine policy(del,outputFile)

        real(KIND=8), parameter :: A=1.0D0
        real(KIND=8), parameter :: alf=0.4D0
        real(KIND=8), parameter :: bet=0.9D0
        real(KIND=8), parameter :: cbar=0.0D0

        integer, parameter :: NUM_GRID_POINTS = 2001
        integer, parameter :: NUM_ITER = 1000
        integer, parameter :: NUM_POINTS=2*NUM_GRID_POINTS*NUM_ITER

        real(KIND=8),INTENT(IN) :: del
        character(LEN=*),INTENT(IN) :: outputFile

        real(KIND=8) :: kss, kd, ku,eps, err
        real(KIND=8), dimension(NUM_GRID_POINTS,1) :: k, k1, one
        real(KIND=8), dimension(NUM_GRID_POINTS,1) :: k_tilde
        integer, dimension(NUM_GRID_POINTS) :: ind_k_tilde
        real(KIND=8), dimension(NUM_GRID_POINTS,NUM_ITER) ::c_policy, v_iter
        real(KIND=8), dimension(NUM_GRID_POINTS,NUM_GRID_POINTS) ::Ut, RHS, PPF, c

        real(KIND=8), dimension(NUM_GRID_POINTS,1) :: colPtr1, colPtr2
        integer i,j

        !some basic initialization
        data c_policy,v_iter /NUM_POINTS * 0.0D0/
        data one /NUM_GRID_POINTS * 1.0D0 /
        print *,"Starting calculation."
        flush(6)

        eps=epsilon(1.0D0)

        !calculating the kss
        kss=(1.0D0/(bet*alf*A)-(1.0D0-del)/alf*A)**(1.0D0/(alf-1.0D0))
        !lowest possible capital expenditure
        kd  = .01*kss
        ku  = 1.5*kss

        do i=1,NUM_GRID_POINTS-1
            k(i,1) = kd+(i-1)*(ku-kd)/(NUM_GRID_POINTS-1)
        end do
        k(NUM_GRID_POINTS,1)=ku
        k1=k

        !An array where cell i,j indicates how much we would produce if we had
        !ki capital level this period and chose kj capital level for next
        !period (note: kj is independent of output, so all elements in a given
        !row are the same)
        PPF = transpose(matmul(one,transpose((1.0D0-del)*k+A*k**alf)))

        !consumption in period. Cell i,j is consumption if we had used ki to
        !produce and kj is next period's capital choice
        c = PPF-matmul(one,transpose(k))

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

            !we only care about the maximum utility
            !gives the maximum in each row of RHS,
            !saves values in v1 and index of k_tilde in ind_k_tilde
            ind_k_tilde=maxloc(RHS,dim=2)

            do j=1,NUM_GRID_POINTS
                v_iter(j,i) = RHS(j,ind_k_tilde(j))

                ! capital next period, they are column vectors
                k_tilde(j,1)  = k(ind_k_tilde(j),1)
            end do

            c_policy(:,i) = k(:,1)**alf +(1.0D0-del)*k(:,1) - k_tilde(:,1)

            ! need to pretend each of these columns is an array
            colPtr1 = reshape(v_iter(:,i-1), (/NUM_GRID_POINTS , 1/))
            colPtr2 = reshape(v_iter(:,i), (/NUM_GRID_POINTS , 1/))
            ! checking the norm of the difference between successive iterations:
            ! this is convergence criterion
            err = maxval(abs(colPtr2 - colPtr1))
            RHS = Ut + bet*matmul(one,transpose(colPtr2))
        end do

        print *,"Finished calculation. ",i," iterations"

        open(1, file=outputFile)
        do i=1,NUM_GRID_POINTS
            write(1,*) k(i,1),",",k_tilde(i,1),",",k(i,1)
        end do
        close(1)

    end subroutine policy

end program main
