module qlsum_deps

contains

    subroutine zpow_dlg ( n, z, iharm, zout )

        implicit none
        integer, intent(in) :: n, iharm
        complex, intent(in) :: z(n)
        complex, intent(out) :: zout(n)

        integer i, nharm 
        complex zin
        complex zk(n)

        zout = 1d0
        
        if ( iharm .eq. 0 ) return

        zk = z

        nharm = abs(iharm)
        do while (nharm .gt. 0)

            if ( mod(nharm,2).eq.1 ) then
                    
                zout    = zout * zk

            endif
            zk = zk**2
            nharm = nharm / 2
        enddo


        iharm_lt0: &
        if (iharm.lt.0) then

                do i=1,n
                    zin = zout(i)
                    call zdiv( zin, zout(i) )
                enddo

        endif iharm_lt0

        return

    end subroutine zpow_dlg
 
    subroutine zpow(n, z, iharm, zout )
    implicit none
    integer n, iharm
    complex z(n), zout(n)

    integer i
    complex one, zero
        complex zin

        integer nharm
        logical isodd
        intrinsic mod

        logical use_zdiv
        parameter(use_zdiv=.true.)

        integer nb
        parameter(nb=1024*4*4)
        complex zk(nb)
        integer istart,iend,isize

        one = 1.0d0
        zero = 0.0d0

        if (iharm.eq.0) then
             do i=1,n
                zout(i) = one
             enddo
             return
        endif


        do istart=1,n,nb

           iend = min(n,istart+nb-1)
           isize = iend-istart+1

      do i=1,isize
        zout(istart-1+i) = one
          enddo


         do i=1,isize
           zk(i) = z(istart-1+i)
         enddo

        nharm = abs(iharm)
        do while (nharm .gt. 0)
           isodd = (mod(nharm,2).eq.1)
           if (isodd) then
               do i=1,isize
                 zout(istart-1+i) = zout(istart-1+i) * zk(i)
               enddo
           endif
           do i=1,isize
              zk(i) = zk(i) * zk(i)
           enddo
           nharm = int( nharm/2 )
        enddo



    if (iharm.lt.0) then
           if (use_zdiv) then
         do i=1,isize
                zin = zout(istart-1+i)
                call zdiv( zin, zout(istart-1+i) )
            enddo
           else
             do i=1,isize
                zin = zout(istart-1+i)
                zout(istart-1+i) = one/zin
             enddo
           endif
        endif

        enddo

    return
    end subroutine 
    

    subroutine zdiv( zin, zout )
        implicit none
        complex zin, zout
        real a, b
        real d

        real one
        parameter(one=1.0d0)
        real rd, a_over_b, b_over_a


        a = real(zin)
        b = aimag(zin)

!       z = (a + i * b)
!       1/z =  a/(a^2 + b^2) - i * b/(a^2 + b^2)
!
!       or    1/(a + (b/a)*b) - i * (b/a) / (a + (b/a)*b)
!       or    (a/b)/( (a/b)*a + b ) - i * 1/( (a/b)*a + b )
!        
        if (abs(a).gt.abs(b)) then
            b_over_a = b/a
            d = a + (b_over_a)*b
            rd = one/d
            zout = cmplx( rd, -(b_over_a)*rd )

            rd  = a**2 + b**2
            zout    = cmplx ( a / rd, -b / rd )
        else
            a_over_b = a/b
            d = (a_over_b)*a + b
            rd = one/d
            zout = cmplx( (a_over_b)*rd, -rd )
        endif

        return
    end subroutine
    
    subroutine zdiv_dlg( zin, zout )
        implicit none
        complex zin, zout
        real a, b
        real d

        real one
        parameter(one=1.0d0)
        real rd, a_over_b, b_over_a


        a = real(zin)
        b = aimag(zin)

        if (abs(a).gt.abs(b)) then
            rd  = a**2 + b**2
            zout    = cmplx ( a / rd, -b / rd )
        else
            a_over_b = a/b
            d = (a_over_b)*a + b
            rd = one/d
            zout = cmplx( (a_over_b)*rd, -rd )
        endif

        return
    end subroutine zdiv_dlg

end module qlsum_deps
     

