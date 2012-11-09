program q3
    implicit none
    integer,dimension(10) :: n=(/-1,-2,-3,-4,-5,-6,-7,-8,-9,-10/)
    real(KIND=8),dimension(10) :: estimate
    real(KIND=8) :: derivative,eps,p=1.5d0
    integer,dimension(1) :: j
    integer :: i

    do i=1,10
        eps=10.0D0**n(i)
        estimate(i)=(0.5D0*((p+eps)**(-0.5D0)-p**(-0.5D0))+0.5D0*((p+eps)**(-0.2D0)-p**(-0.2D0)))/eps
        print *,"eps: ",eps," estimate: ",estimate(i)
    end do
    derivative=-0.25D0*p**(-1.5d0)-0.1D0*p**(-1.2D0)
    print *, "Real: ",derivative

    j=minloc(abs(estimate-derivative))
    print *,"Most accurate: ",j(1)

end program q3
