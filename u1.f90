program u13d
    implicit none

    ! defining our types
    integer, parameter :: r8=selected_real_kind(15)
    integer, parameter :: i4=selected_int_kind(8)

    ! defining our parameters
    complex(kind=r8), parameter :: cr=(1.0_r8,0.0_r8)
    complex(kind=r8), parameter :: ci=(0.0_r8,1.0_r8)
    complex(kind=r8), parameter :: cz=(0.0_r8,0.0_r8)
    real(kind=r8), parameter :: pi=dacos(-1.0_r8)
    integer(kind=i4), parameter :: d=4 ! dimensions of our euclidian space-time

    ! defining our matrices
    complex(kind=r8), allocatable, dimension(:,:,:,:,:) :: u
    real(kind=r8), allocatable, dimension(:,:) :: w

    ! defining our variables
    real(kind=r8) :: beta, delta, ts, te
    integer(kind=i4) :: imc, nr, nmc, nterm, i, j, a, b

    call cpu_time(ts)
    ! reading the input value
    open(unit=1,file='u1-in.dat')    
    read(1,*) nr ! number of sites
    read(1,*) beta ! beta
    read(1,*) nmc ! number of monte carlo sweeps
    read(1,*) nterm ! number of termalization updates
    read(1,*) delta ! delta for the metropolis
    read(1,*) a, b  ! size of the wilson loops
    close(1)
    
    ! allocate the memory
    allocate(u(nr,nr,nr,nr,d),w(a,b))
    u=cr
    ! thermalize
    do i=1,nterm
    call metropolis(u,beta,delta)
    enddo

    ! make the measuments
    do i=1,nmc
        call metropolis(u,beta,delta)
        call wilsonloop(u,w,a,b)
        !call measurewilson(u,a,b,w)
        do j=1,a
            write(100+j,*) i, w(j,:)
        enddo
    enddo

    ! deallocate the used memory
    deallocate(u)
    call cpu_time(te)
    print*,(te-ts)/(60.0_r8),"minutes"
    contains
    ! subroutine that return us the six neighbors
    subroutine neighbors(n,npmi,nmmi,npni,nmni,npmimni,mi,ni)
        integer(kind=i4), dimension(d) :: n, npmi, nmmi, npni, nmni, npmimni
        integer(kind=i4) :: mi, ni
        npmi=n
        nmmi=n
        npni=n
        nmni=n

        ! calculating the vectors with the neighbors positionr
        npmi(mi)=npmi(mi)+1
        nmmi(mi)=nmmi(mi)-1
        npni(ni)=npni(ni)+1
        nmni(ni)=nmni(ni)-1

        npmimni=npmi
        npmimni(ni)=npmimni(ni)-1

        ! applyng the periodic bounduary conditionr
        if((npmi(mi).gt.nr))then
            npmi(mi)=1
            npmimni(mi)=1
        endif
        if((npni(ni).gt.nr))then
            npni(ni)=1
        endif
        if(nmmi(mi).lt.1)then
            nmmi(mi)=nr
        endif
        if(nmni(ni).lt.1)then
            nmni(ni)=nr
            npmimni(ni)=nr
        endif
    endsubroutine neighbors

    ! subroutine that compute the sum of the stamples
    function sumstamples(u,n,mi)
        complex(kind=r8) :: sumstamples, temp, u(nr,nr,nr,nr,d)
        integer(kind=i4), dimension(d) :: n,npmi,nmmi,npni,nmni,npmimni
        integer(kind=i4) :: mi, ni
        temp=cz
        do ni=1,d
            if(ni.ne.mi)then
                call neighbors(n,npmi,nmmi,npni,nmni,npmimni,mi,ni)

                temp=temp+u(npmi(1),npmi(2),npmi(3),npmi(4),ni)*&
                &conjg(u(npni(1),npni(2),npni(3),npni(4),mi))*conjg(u(n(1),n(2),n(3),n(4),ni))

                temp=temp+conjg(u(npmimni(1),npmimni(2),npmimni(3),npmimni(4),ni))*&
                &conjg(u(nmni(1),nmni(2),nmni(3),nmni(4),mi))*u(nmni(1),nmni(2),nmni(3),nmni(4),ni)
            endif
        enddo
        sumstamples=temp
    endfunction sumstamples

    ! subroutine that thermalize the lattice
    subroutine metropolis(u,beta,delta)
        complex(kind=r8), intent(inout) :: u(nr,nr,nr,nr,d)
        complex(kind=r8) :: uprime, x
        real(kind=r8) :: alpha, beta, delta, zeta, ds
        integer(kind=i4) :: e1, e2, e3, e4, mi

        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=1,nr
            do mi=1,d
                ! we make a metropolis step on this link
                ! first, we make a x random number closy to the unity
                call random_number(alpha)
                alpha=1.0_r8-2.0_r8*alpha   ! put this number in the interval [-1,1)
                alpha=delta*alpha           ! put this numbr in the interval [-delta,delta)
                x=(cr+alpha*ci)/dsqrt(1.0_r8-alpha*alpha)

                ! we crete a candidate for a new link from x
                uprime=x*u(e1,e2,e3,e4,mi)

                ! and we compute the change in the action due to this uprime
                ds=-beta*dreal((uprime-u(e1,e2,e3,e4,mi))*sumstamples(u,(/e1,e2,e3,e4/),mi))
                !print*,ds
                ! we accept it if the action is 
                if(ds.lt.0.0_r8)then
                    u(e1,e2,e3,e4,mi)=uprime/abs(uprime)
                else
                    ! or we accept it following the boltzman distribution
                    call random_number(zeta)
                    if(zeta.lt.dexp(-ds))then
                        u(e1,e2,e3,e4,mi)=uprime/abs(uprime)
                    endif
                endif
            enddo
        enddo
        enddo
        enddo
        enddo
    endsubroutine metropolis

    ! we compute the observables
    subroutine wilsonloop(u,w,a,b)
        complex(kind=r8), intent(in) :: u(nr,nr,nr,nr,d)
        real(kind=r8) :: w(a,b)
        integer(kind=i4), intent(in) :: a, b
        complex(kind=r8) :: umini
        real(kind=r8) :: sum
        integer(kind=i4) :: e1, e2, e3, e4, mi, ni
        integer(kind=i4) :: l, lmi, lni, n(4)
        do lmi=1,a
        do lni=1,b
            ! compute the loop of size (lmi, lni)
            sum=0.0_r8
            ! visit each lattice site
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nr
                n=(/e1,e2,e3,e4/)
                ! visit all planes
                do mi=1,3
                do ni=mi+1,4
                    umini=cr
                    
                    do l=1,lmi
                    umini=umini*u(n(1),n(2),n(3),n(4),mi)
                    n(mi)=n(mi)+1
                    if(n(mi).gt.nr) n(mi)=1
                    enddo

                    do l=1,lni
                    umini=umini*u(n(1),n(2),n(3),n(4),ni)
                    n(ni)=n(ni)+1
                    if(n(ni).gt.nr) n(ni)=1
                    enddo

                    do l=1,lmi
                    n(mi)=n(mi)-1
                    if(n(mi).lt.1) n(mi)=nr
                    umini=umini*conjg(u(n(1),n(2),n(3),n(4),mi))
                    enddo

                    do l=1,lni
                    n(ni)=n(ni)-1
                    if(n(ni).lt.1) n(ni)=nr
                    umini=umini*conjg(u(n(1),n(2),n(3),n(4),ni))
                    enddo
                    sum=sum+dreal(umini)
                enddo
                enddo
            enddo
            enddo
            enddo
            enddo
            w(lmi,lni)=sum/(6.0_r8*nr**4.0_r8)
        enddo
        enddo
    endsubroutine wilsonloop
endprogram u13d