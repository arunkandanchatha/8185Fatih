program q2
    implicit none

    real(KIND=8), dimension(2) :: phi = (/1.0D0,0.61803398D0/)
    real(KIND=8) :: temp,realPhi = (sqrt(5.D0)-1.0D0)/2.0D0
    integer :: i

    do i=2,20
        temp = phi(1)-phi(2)
        print *,"i: ",i,"  calculated: ",temp, "actual: ", realPhi**i
        phi(1)=phi(2)
        phi(2)=temp
    end do

end program q2
