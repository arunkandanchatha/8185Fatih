program q1
    implicit none
    integer,dimension(8) :: n=(/1,2,3,4,5,6,7,8/)
    real(KIND=8),dimension(8,2) :: roots1,roots2
    real(KIND=8) :: a,b,c,q
    integer :: i

    a=1.D0
    b=100000.D0
    do i=1,8
        c=10.D0**(-(n(i)))
        roots1(i,1)=(-b+sqrt(b**2-4*a*c))/2.D0
        roots1(i,2)=(-b+sqrt(b**2-4*a*c))/2.D0
        if(b<0.0) then
            q=-(b-sqrt(b**2-4*a*c))/2.D0
        else
            q=-(b+sqrt(b**2-4*a*c))/2.D0
        end if
        print *,q
        roots2(i,1)=q/a
        roots2(i,2)=c/q
    end do

    print *,"     Method 1                          Method2"
    do i=1,8
        print *,max(roots1(i,1),roots1(i,2)),"    ",max(roots2(i,1),roots2(i,2))
    end do
end program q1
