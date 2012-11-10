program q5
    implicit none

    real(KIND=8) :: a=5.0D0,b=1.0D0,diff=100.0D0,x,valA,valB,valX
    integer :: counter=0

    valA=0.5D0*a**(-0.5D0)+0.5D0*a**(-0.2D0)-0.75D0
    valB=0.5D0*b**(-0.5D0)+0.5D0*b**(-0.2D0)-0.75D0
    if(valA*valB>0.0D0) then
        print *,"a: ",a," b: ", b, "f(a): ",valA,"f(b): ",valB
        stop 0
    end if

    do while ((diff>10.0e-6).and.(counter<10000))
        counter = counter +1
        x = (a+b)/2.0D0
        valX=0.5D0*x**(-0.5D0)+0.5D0*x**(-0.2D0)-0.75D0

        if(valX*valA>0.0D0) then
            a=x
        else
            b=x
        end if

        diff=abs(valX)
    end do

    print *,"Bijection Method. root: ",x, " iterations: ",counter

    !reinitialize for secant method
    a=5.0D0
    b=1.0D0
    valA=0.5D0*a**(-0.5D0)+0.5D0*a**(-0.2D0)-0.75D0
    valB=0.5D0*b**(-0.5D0)+0.5D0*b**(-0.2D0)-0.75D0
    counter = 0
    diff = 100.0D0

    do while ((diff>10.0e-6).and.(counter<10000))
        counter = counter +1
        x=b-valB*(b-a)/(valB-valA)
        a=b
        b=x
        valA=valB
        valB=0.5D0*b**(-0.5D0)+0.5D0*b**(-0.2D0)-0.75D0

        diff=abs(valB)
    end do

    print *,"Secant Method. root: ",b, " iterations: ",counter

    !reinitialize for newton method

    !a represents x(n)
    a=5.0D0
    !use valA to represent the value of function at x(n)
    valA=0.5D0*a**(-0.5D0)+0.5D0*a**(-0.2D0)-0.75D0
    diff = abs(valA)
    counter = 0

    do while ((diff>10.0e-6).and.(counter<10000))
        counter = counter +1

        !use valB to represent the derivative at x(n)
        valB=-0.25D0*a**(-1.5d0)-0.1D0*a**(-1.2D0)

        !find x(n+1)
        a=a-valA/valB

        !find x(n+1)
        valA=0.5D0*a**(-0.5D0)+0.5D0*a**(-0.2D0)-0.75D0

        !how far are we from 0?
        diff=abs(valA)
    end do

    print *,"Newton Method. root: ",a, " iterations: ",counter

end program q5
