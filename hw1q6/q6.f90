program q6
    implicit none

    real(KIND=8) :: a=5.0D0,b=1.0D0,valA,valB,tmp,c,valC,d=0.0D0
    real(KIND=8) :: s, valS, toler=10.0e-6
    integer :: counter=0
    logical mflag, cond1, cond2, cond3, cond4, cond5

    valA=0.5D0*a**(-0.5D0)+0.5D0*a**(-0.2D0)-0.75D0
    valB=0.5D0*b**(-0.5D0)+0.5D0*b**(-0.2D0)-0.75D0

    if(valA*valB>0.0D0) then
        print *,"a: ",a," b: ", b, "f(a): ",valA,"f(b): ",valB
        stop 0
    end if

    if(abs(valA)<abs(valB)) then
        tmp=b
        b=a
        a=tmp

        tmp=valB
        valB=valA
        valA=tmp
    end if

    c=a
    valC=valA
    mflag = .true.

    do while ((abs(valB)>toler) .and. (abs(b-a)>toler) .and. (counter<10000))

        counter = counter+1

        if((valA/=valC).and.(valB/=valC)) then
            s=(a*valA*valC)/((valA-valB)*(valA-valC))+(b*valA*valC)/((valB-valA)*(valB-valC))+&
                &(c*valA*valB)/((valC-valA)*(valC-valB))
        else
            s=b-valB*(b-a)/(valB-valA)
        end if

        cond1 = ((s<((3.0D0*a+b)/4.0D0)) .and. (s>b)) .or. ((s>((3.0D0*a+b)/4.0D0)) .and. (s<b))
        cond2 = mflag .and. (abs(s-b) .ge. (abs(b-c)/2.0D0))
        cond3 = (.not. mflag) .and. (abs(s-b) .ge. (abs(c-d)/2.0D0))
        cond4 = mflag .and. (abs(b-c) < abs(toler))
        cond5 = (.not. mflag) .and. (abs(c-d)<abs(toler))

        if(cond1 .or. cond2 .or. cond3 .or. cond4 .or. cond5) then
            s=(a+b)/2
            mflag = .true.
        else
            mflag = .false.
        end if

        valS = 0.5D0*s**(-0.5D0)+0.5D0*s**(-0.2D0)-0.75D0
        d = c
        c = b
        valC = valB

        if (valA*valS < 0.0D0) then
            b = s
            valB = valS
        else
            a = s
            valA = valS
        end if

        if (abs(a)<abs(b)) then
            tmp=b
            b=a
            a=tmp

            tmp=valB
            valB=valA
            valA=tmp
        end if

    end do

    print *,"Brent Method. root: ",b, " iterations: ",counter


end program q6
