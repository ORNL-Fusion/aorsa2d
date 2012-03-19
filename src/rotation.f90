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

        ! Alternate approach
        real :: zAxis(0:2),perp(0:2),PerpUnit(0:2),bUnitCar(0:2)
        real :: PerpDotB,PerpMag,bUnitMag,CheckTh,Theta
        real :: q0,q1,q2,q3,RotQ(3,3),InvRotQ(3,3)
 
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

        allocate ( g%U_RTZ_to_ABb(g%nR,g%nZ,3,3) )

        !sqr     = sqrt ( 1d0 - g%bR_unit**2 )

        !g%Urr    = sqr
        !g%Urt    = -g%bR_unit * g%bT_unit / sqr
        !g%Urz    = -g%bR_unit * g%bZ_unit / sqr

        !g%Utr   = 0
        !g%Utt   = g%bZ_unit / sqr
        !g%Utz   = -g%bT_unit / sqr

        !g%Uzr    = g%bR_unit
        !g%Uzt    = g%bT_unit
        !g%Uzz    = g%bZ_unit


        !! Check the determinant, should = 1
        !! ---------------------------------
    
        !allocate(det(g%nR,g%nZ))

        !det = g%Urr * g%Utt * g%Uzz &
        !    + g%Urt * g%Utz * g%Uzr &
        !    + g%Urz * g%Utr * g%Uzt &
        !    - g%Urr * g%Utz * g%Uzt &
        !    - g%Urt * g%Utr * g%Uzz &
        !    - g%Urz * g%Utt * g%Uzr

        !if ( any(abs(1-det)>5e-1) ) then

        !    write(*,*) 'rotation.f90: ERROR det != 1, it is = ', 1-det
        !    stop

        !endif

        !deallocate(det)


        !g%U_RTZ_to_ABb(:,:,1,1)  = g%Urr
        !g%U_RTZ_to_ABb(:,:,1,2)  = g%Urt
        !g%U_RTZ_to_ABb(:,:,1,3)  = g%Urz
        !
        !g%U_RTZ_to_ABb(:,:,2,1)  = g%Utr
        !g%U_RTZ_to_ABb(:,:,2,2)  = g%Utt
        !g%U_RTZ_to_ABb(:,:,2,3)  = g%Utz

        !g%U_RTZ_to_ABb(:,:,3,1)  = g%Uzr
        !g%U_RTZ_to_ABb(:,:,3,2)  = g%Uzt
        !g%U_RTZ_to_ABb(:,:,3,3)  = g%Uzz


        ! OK, build the rotation matrix a different way.
        ! Remember, this is the rotation to and from 
        ! alp,bet,prl and r,t,z
        !
        ! Directions are defined as follows?
        !
        ! prl - in the direction of B
        ! alp - perp to both z and B
        ! bet - alp x prl

        do i=1,g%nR
            do j=1,g%nZ

                ! Get vector perp to both z axis and b

                zAxis = 0.0
                zAxis(2) = 1.0

                bUnitCar = (/ g%bR_unit(i,j),g%bT_unit(i,j),g%bZ_unit(i,j) /)

                Perp(0) = +(zAxis(1)*bUnitCar(2)-zAxis(2)*bUnitCar(1))
                Perp(1) = -(zAxis(0)*bUnitCar(2)-zAxis(2)*bUnitCar(0))
                Perp(2) = +(zAxis(0)*bUnitCar(1)-zAxis(1)*bUnitCar(0)) 

                ! Check angle between Perp and b

                PerpDotB    = Perp(0)*bUnitCar(0)+Perp(1)*bUnitCar(1)+Perp(2)*bUnitCar(2)
                PerpMag = sqrt ( sum( Perp**2 ) )
                bUnitMag    = sqrt ( sum( bUnitCar**2 ) )
                CheckTh = aCos ( PerpDotB / ( PerpMag * bUnitMag ) ) * 180.0/pi

                PerpUnit    = Perp / PerpMag

                ! Get angle between z axis and b

                theta   = aCos ( bUnitCar(2) )

                ! Calculate the quaternions

                q0  = cos ( theta / 2.0 )
                q1  = sin ( theta / 2.0 ) * PerpUnit(0)
                q2  = sin ( theta / 2.0 ) * PerpUnit(1) 
                q3  = sin ( theta / 2.0 ) * PerpUnit(2)

                ! Construct the rotation matrix

                RotQ(1,1:3) = (/ q0**2+q1**2-q2**2-q3**2, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2) /)
                RotQ(2,1:3) = (/ 2*(q2*q1+q0*q3), q0**2-q1**2+q2**2-q3**2, 2*(q2*q3-q0*q1) /)
                RotQ(3,1:3) = (/ 2*(q3*q1-q0*q2), 2*(q3*q2+q0*q1), q0**2-q1**2-q2**2+q3**2 /)

                InvRotQ    = transpose ( rotQ )

                g%U_RTZ_to_ABb(i,j,1:3,1:3)  = RotQ

                g%Urr = g%U_RTZ_to_ABb(i,j,1,1) 
                g%Urt = g%U_RTZ_to_ABb(i,j,1,2) 
                g%Urz = g%U_RTZ_to_ABb(i,j,1,3) 

                g%Utr = g%U_RTZ_to_ABb(i,j,2,1) 
                g%Utt = g%U_RTZ_to_ABb(i,j,2,2) 
                g%Utz = g%U_RTZ_to_ABb(i,j,2,3) 

                g%Uzr = g%U_RTZ_to_ABb(i,j,3,1) 
                g%Uzt = g%U_RTZ_to_ABb(i,j,3,2) 
                g%Uzz = g%U_RTZ_to_ABb(i,j,3,3) 

        enddo
    enddo

    end sUbroUtine init_rotation


    subroutine deriv_rotation ( g )

        use aorsaNamelist, &
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

                call deriv_y(g%z, g%Urr, i, j, dfdy = g%dzUrr(i,j), d2fdy2 = g%dzzUrr(i,j))
                call deriv_y(g%z, g%Urt, i, j, dfdy = g%dzUrt(i,j), d2fdy2 = g%dzzUrt(i,j))
                call deriv_y(g%z, g%Urz, i, j, dfdy = g%dzUrz(i,j), d2fdy2 = g%dzzUrz(i,j))

                call deriv_y(g%z, g%Utr, i, j, dfdy = g%dzUtr(i,j), d2fdy2 = g%dzzUtr(i,j))
                call deriv_y(g%z, g%Utt, i, j, dfdy = g%dzUtt(i,j), d2fdy2 = g%dzzUtt(i,j))
                call deriv_y(g%z, g%Utz, i, j, dfdy = g%dzUtz(i,j), d2fdy2 = g%dzzUtz(i,j))

                call deriv_y(g%z, g%Uzr, i, j, dfdy = g%dzUzr(i,j), d2fdy2 = g%dzzUzr(i,j))
                call deriv_y(g%z, g%Uzt, i, j, dfdy = g%dzUzt(i,j), d2fdy2 = g%dzzUzt(i,j))
                call deriv_y(g%z, g%Uzz, i, j, dfdy = g%dzUzz(i,j), d2fdy2 = g%dzzUzz(i,j))

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
