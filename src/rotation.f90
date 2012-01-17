modUle rotation

implicit none

real :: sqx

contains

    subroutine init_rotation ( g )

        use grid

        implicit none

        type(gridBlock), intent(inout) :: g

        integer :: i,j 
        real, allocatable :: sqr(:,:)
        real, allocatable :: det(:,:)

        allocate ( sqr ( g%nR, g%nZ ) )

        ! Create a cylindrical coordinate version of the 
        ! rotation matrix so we have some consistency in the
        ! naming conventions of coordinates and hence preserve
        ! what is left of my sanity

        allocate ( &
            g%Urr( g%nR, g%nZ ), & 
            g%Urt( g%nR, g%nZ ), &
            g%Urz( g%nR, g%nZ ), &
            g%Utr( g%nR, g%nZ ), & 
            g%Utt( g%nR, g%nZ ), &
            g%Utz( g%nR, g%nZ ), &
            g%Uzr( g%nR, g%nZ ), & 
            g%Uzt( g%nR, g%nZ ), &
            g%Uzz( g%nR, g%nZ ) )

        sqr     = sqrt ( 1d0 - g%bR_unit**2 )

        g%Urr     = sqr
        g%Urt    = -g%bR_unit * g%bT_unit / sqr
        g%Urz     = -g%bR_unit * g%bZ_unit / sqr

        g%Utr    = 0
        g%Utt   = g%bZ_unit / sqr
        g%Utz    = -g%bT_unit / sqr

        g%Uzr     = g%bR_unit
        g%Uzt    = g%bT_unit
        g%Uzz     = g%bZ_unit


        ! Check the determinant, should = 1
        ! ---------------------------------
    
        allocate(det(g%nR,g%nZ))

        det = g%Urr * g%Utt * g%Uzz &
            + g%Urt * g%Utz * g%Uzr &
            + g%Urz * g%Utr * g%Uzt &
            - g%Urr * g%Utz * g%Uzt &
            - g%Urt * g%Utr * g%Uzz &
            - g%Urz * g%Utt * g%Uzr

        if ( any(abs(1-det)>5e-1) ) then

            write(*,*) 'rotation.f90: ERROR det != 1, it is = ', 1-det
            stop

        endif

        deallocate(det)

        allocate ( g%U_RTZ_to_ABb(g%nR,g%nZ,3,3) )

        g%U_RTZ_to_ABb(:,:,1,1)  = g%Urr
        g%U_RTZ_to_ABb(:,:,1,2)  = g%Urt
        g%U_RTZ_to_ABb(:,:,1,3)  = g%Urz
        
        g%U_RTZ_to_ABb(:,:,2,1)  = g%Utr
        g%U_RTZ_to_ABb(:,:,2,2)  = g%Utt
        g%U_RTZ_to_ABb(:,:,2,3)  = g%Utz

        g%U_RTZ_to_ABb(:,:,3,1)  = g%Uzr
        g%U_RTZ_to_ABb(:,:,3,2)  = g%Uzt
        g%U_RTZ_to_ABb(:,:,3,3)  = g%Uzz

    end sUbroUtine init_rotation


    subroutine deriv_rotation ( g )

        use aorsa2din_mod, &
        only: nZFUn, r0
        use grid
        use derivatives 

        implicit none

        type(gridBlock), intent(inout) :: g

        integer :: i, j
        real :: distance


        allocate ( g%gradPrlB (g%nR,g%nZ), g%sinTh(g%nR,g%nZ) )
        allocate ( g%bPol(g%nR,g%nZ) )

        ! NOTE: CHECK THIS, i.e., unit or absolute?
        g%bPol = sqrt ( g%bR_unit**2 + g%bZ_unit**2 ) * g%bMag

        g%gradPrlB = 0

        do i = 1, g%nR
            do j = 1, g%nZ

                ! Brambillas approximation to dBdPar
                ! ----------------------------------

                distance = sqrt ( (g%R(i)-r0)**2 + g%Z(j)**2 )
                if ( distance > 0 ) then

                    g%sinTh(i,j) =  g%Z(j) / distance
                    g%gradPrlB(i,j) = g%bMag(i,j) / g%R(i) * abs ( g%bPol(i,j) * g%sinTh(i,j) )

                else ! on magnetic axis (r=r0, y=y0=0)

                    g%gradPrlB(i,j) = 0

                endif

                !if (nzfUn == 0) gradPrlB(i,j) = 1.0e-10

           enddo
        enddo

        ! R,Th,z version
        ! --------------

        allocate ( &
            g%drUrr(g%nR,g%nZ), g%drrUrr(g%nR,g%nZ), &
            g%drUrt(g%nR,g%nZ), g%drrUrt(g%nR,g%nZ), &
            g%drUrz(g%nR,g%nZ), g%drrUrz(g%nR,g%nZ), &
            g%drUtr(g%nR,g%nZ), g%drrUtr(g%nR,g%nZ), &
            g%drUtt(g%nR,g%nZ), g%drrUtt(g%nR,g%nZ), &
            g%drUtz(g%nR,g%nZ), g%drrUtz(g%nR,g%nZ), &
            g%drUzr(g%nR,g%nZ), g%drrUzr(g%nR,g%nZ), &
            g%drUzt(g%nR,g%nZ), g%drrUzt(g%nR,g%nZ), &
            g%drUzz(g%nR,g%nZ), g%drrUzz(g%nR,g%nZ) )

        allocate ( &
            g%dzUrr(g%nR,g%nZ), g%dzzUrr(g%nR,g%nZ), &
            g%dzUrt(g%nR,g%nZ), g%dzzUrt(g%nR,g%nZ), &
            g%dzUrz(g%nR,g%nZ), g%dzzUrz(g%nR,g%nZ), &
            g%dzUtr(g%nR,g%nZ), g%dzzUtr(g%nR,g%nZ), &
            g%dzUtt(g%nR,g%nZ), g%dzzUtt(g%nR,g%nZ), &
            g%dzUtz(g%nR,g%nZ), g%dzzUtz(g%nR,g%nZ), &
            g%dzUzr(g%nR,g%nZ), g%dzzUzr(g%nR,g%nZ), &
            g%dzUzt(g%nR,g%nZ), g%dzzUzt(g%nR,g%nZ), &
            g%dzUzz(g%nR,g%nZ), g%dzzUzz(g%nR,g%nZ) )

        allocate ( &
            g%drzUrr(g%nR,g%nZ), g%drzUrt(g%nR,g%nZ), g%drzUrz(g%nR,g%nZ), &
            g%drzUtr(g%nR,g%nZ), g%drzUtt(g%nR,g%nZ), g%drzUtz(g%nR,g%nZ), &
            g%drzUzr(g%nR,g%nZ), g%drzUzt(g%nR,g%nZ), g%drzUzz(g%nR,g%nZ) ) 

        do i = 1, g%nR
            do j = 1, g%nZ

                call deriv_x(g%R, g%Urr, i, j, dfdx = g%drUrr(i,j), d2fdx2 = g%drrUrr(i,j))
                call deriv_x(g%R, g%Urt, i, j, dfdx = g%drUrt(i,j), d2fdx2 = g%drrUrt(i,j))
                call deriv_x(g%R, g%Urz, i, j, dfdx = g%drUrz(i,j), d2fdx2 = g%drrUrz(i,j))

                call deriv_x(g%R, g%Utr, i, j, dfdx = g%drUtr(i,j), d2fdx2 = g%drrUtr(i,j))
                call deriv_x(g%R, g%Utt, i, j, dfdx = g%drUtt(i,j), d2fdx2 = g%drrUtt(i,j))
                call deriv_x(g%R, g%Utz, i, j, dfdx = g%drUtz(i,j), d2fdx2 = g%drrUtz(i,j))

                call deriv_x(g%R, g%Uzr, i, j, dfdx = g%drUzr(i,j), d2fdx2 = g%drrUzr(i,j))
                call deriv_x(g%R, g%Uzt, i, j, dfdx = g%drUzt(i,j), d2fdx2 = g%drrUzt(i,j))
                call deriv_x(g%R, g%Uzz, i, j, dfdx = g%drUzz(i,j), d2fdx2 = g%drrUzz(i,j))

                call deriv_y(g%z, g%Urr, i, j, dfdy = g%drUrr(i,j), d2fdy2 = g%drrUrr(i,j))
                call deriv_y(g%z, g%Urt, i, j, dfdy = g%drUrt(i,j), d2fdy2 = g%drrUrt(i,j))
                call deriv_y(g%z, g%Urz, i, j, dfdy = g%drUrz(i,j), d2fdy2 = g%drrUrz(i,j))

                call deriv_y(g%z, g%Utr, i, j, dfdy = g%drUtr(i,j), d2fdy2 = g%drrUtr(i,j))
                call deriv_y(g%z, g%Utt, i, j, dfdy = g%drUtt(i,j), d2fdy2 = g%drrUtt(i,j))
                call deriv_y(g%z, g%Utz, i, j, dfdy = g%drUtz(i,j), d2fdy2 = g%drrUtz(i,j))

                call deriv_y(g%z, g%Uzr, i, j, dfdy = g%drUzr(i,j), d2fdy2 = g%drrUzr(i,j))
                call deriv_y(g%z, g%Uzt, i, j, dfdy = g%drUzt(i,j), d2fdy2 = g%drrUzt(i,j))
                call deriv_y(g%z, g%Uzz, i, j, dfdy = g%drUzz(i,j), d2fdy2 = g%drrUzz(i,j))

                call deriv_xy(g%R, g%z, g%Urr, i, j, d2fdxy = g%drzUrr(i,j))
                call deriv_xy(g%R, g%z, g%Urt, i, j, d2fdxy = g%drzUrt(i,j))
                call deriv_xy(g%R, g%z, g%Urz, i, j, d2fdxy = g%drzUrz(i,j))

                call deriv_xy(g%R, g%z, g%Utr, i, j, d2fdxy = g%drzUtr(i,j))
                call deriv_xy(g%R, g%z, g%Utt, i, j, d2fdxy = g%drzUtt(i,j))
                call deriv_xy(g%R, g%z, g%Utz, i, j, d2fdxy = g%drzUtz(i,j))

                call deriv_xy(g%R, g%z, g%Uzr, i, j, d2fdxy = g%drzUzr(i,j))
                call deriv_xy(g%R, g%z, g%Uzt, i, j, d2fdxy = g%drzUzt(i,j))
                call deriv_xy(g%R, g%z, g%Uzz, i, j, d2fdxy = g%drzUzz(i,j))

           enddo
        enddo

        !drUrr =0 
        !drUrt =0 
        !drUrz =0 

        !drUtr =0 
        !drUtt =0 
        !drUtz =0 

        !drUzr =0 
        !drUzt =0 
        !drUzz =0 

        !drUrr =0 
        !drUrt =0 
        !drUrz =0 

        !drUtr =0 
        !drUtt =0 
        !drUtz =0 

        !drUzr =0 
        !drUzt =0 
        !drUzz =0 

        !drrUrr=0 
        !drrUrt=0
        !drrUrz=0

        !drrUtr=0
        !drrUtt=0
        !drrUtz=0

        !drrUzr=0
        !drrUzt=0
        !drrUzz=0

        !drrUrr=0
        !drrUrt=0
        !drrUrz=0

        !drrUtr=0
        !drrUtt=0
        !drrUtz=0

        !drrUzr=0
        !drrUzt=0
        !drrUzz=0

        !drzUrr=0 
        !drzUrt=0
        !drzUrz=0

        !drzUtr=0
        !drzUtt=0
        !drzUtz=0

        !drzUzr=0
        !drzUzt=0
        !drzUzz=0

    end subroutine deriv_rotation

end module rotation
