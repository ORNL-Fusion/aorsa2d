modUle rotation

implicit none

real :: sqx

contains

    function RotMatHere (bR_unit,bT_unit,bZ_unit)

        use constants

        implicit none

        real, intent(in) :: bR_unit, bT_unit, bZ_unit

        real :: RotMatHere(3,3)

        ! Alternate approach
        real :: zAxis(0:2),perp(0:2),PerpUnit(0:2),bUnitCar(0:2)
        real :: PerpDotB,PerpMag,bUnitMag,CheckTh,Theta
        real :: q0,q1,q2,q3,InvRotQ(3,3)

        ! Get vector perp to both z axis and b

        zAxis = 0.0
        zAxis(2) = 1.0

        !bUnitCar = (/ g%bR_unit(w),g%bT_unit(w),g%bZ_unit(w) /)
        bUnitCar = (/ bR_unit,bT_unit,bZ_unit /)

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

        RotMatHere(1,1:3) = (/ q0**2+q1**2-q2**2-q3**2, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2) /)
        RotMatHere(2,1:3) = (/ 2*(q2*q1+q0*q3), q0**2-q1**2+q2**2-q3**2, 2*(q2*q3-q0*q1) /)
        RotMatHere(3,1:3) = (/ 2*(q3*q1-q0*q2), 2*(q3*q2+q0*q1), q0**2-q1**2-q2**2+q3**2 /)

        InvRotQ    = transpose ( RotMatHere )

    end function RotMatHere

    subroutine init_rotation ( g )

        use grid

        implicit none

        type(gridBlock), intent(inout) :: g

        integer :: i,j,w 
        !real, allocatable :: sqr(:)
        real, allocatable :: det(:,:)
        real :: RotQ(3,3)


        !allocate ( sqr ( g%nR, g%nZ ) )
        !allocate ( sqr ( size(g%wl) ) )

        ! Create a cylindrical coordinate version of the 
        ! rotation matrix so we have some consistency in the
        ! naming conventions of coordinates and hence preserve
        ! what is left of my sanity

        !allocate ( &
        !    g%Urr( g%nR, g%nZ ), & 
        !    g%Urt( g%nR, g%nZ ), &
        !    g%Urz( g%nR, g%nZ ), &
        !    g%Utr( g%nR, g%nZ ), & 
        !    g%Utt( g%nR, g%nZ ), &
        !    g%Utz( g%nR, g%nZ ), &
        !    g%Uzr( g%nR, g%nZ ), & 
        !    g%Uzt( g%nR, g%nZ ), &
        !    g%Uzz( g%nR, g%nZ ) )
        allocate ( &
            g%Urr( size(g%pt) ), & 
            g%Urt( size(g%pt) ), &
            g%Urz( size(g%pt) ), &
            g%Utr( size(g%pt) ), & 
            g%Utt( size(g%pt) ), &
            g%Utz( size(g%pt) ), &
            g%Uzr( size(g%pt) ), & 
            g%Uzt( size(g%pt) ), &
            g%Uzz( size(g%pt) ) )

        !allocate ( g%U_RTZ_to_ABb(g%nR,g%nZ,3,3) )
        allocate ( g%U_RTZ_to_ABb(size(g%pt),3,3) )

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

        do w=1,size(g%pt)

            !i=g%wl(w)%i
            !j=g%wl(w)%j

        !do i=1,g%nR
        !    do j=1,g%nZ

                RotQ = RotMatHere (g%bR_unit(w),g%bT_unit(w),g%bZ_unit(w))

                g%U_RTZ_to_ABb(w,1:3,1:3)  = RotQ

                g%Urr = g%U_RTZ_to_ABb(w,1,1) 
                g%Urt = g%U_RTZ_to_ABb(w,1,2) 
                g%Urz = g%U_RTZ_to_ABb(w,1,3) 

                g%Utr = g%U_RTZ_to_ABb(w,2,1) 
                g%Utt = g%U_RTZ_to_ABb(w,2,2) 
                g%Utz = g%U_RTZ_to_ABb(w,2,3) 

                g%Uzr = g%U_RTZ_to_ABb(w,3,1) 
                g%Uzt = g%U_RTZ_to_ABb(w,3,2) 
                g%Uzz = g%U_RTZ_to_ABb(w,3,3) 

    !    enddo
    !enddo
    enddo

    end sUbroUtine init_rotation


    subroutine deriv_rotation ( g )

        use aorsaNamelist, &
        only: nZFUn, r0
        use grid
        use derivatives 
        use interp, only: dlg_interpB 

        implicit none

        type(gridBlock), intent(inout) :: g

        integer :: i, j, w
        real :: distance

        real :: bTmp(3),bMagTmp
        real :: bRu,bTu,bZu,Rot(3,3),RotLf(3,3),RotRt(3,3),&
            RotDn(3,3),RotUp(3,3),RotRt_(3,3),RotLf_(3,3),&
            RotDn_(3,3),RotUp_(3,3)
        real :: dx, dy
        real :: dfdx(3,3), d2fdx2(3,3), dfdy(3,3), d2fdy2(3,3),d2fdxy(3,3)
 
        !allocate ( g%gradPrlB (g%nR,g%nZ), g%sinTh(g%nR,g%nZ) )
        !allocate ( g%bPol(g%nR,g%nZ) )
        allocate ( g%gradPrlB(size(g%pt)), g%sinTh(size(g%pt)) )
        allocate ( g%bPol(size(g%pt)) )

        ! NOTE: CHECK THIS, i.e., unit or absolute?
        g%bPol = sqrt ( g%bR_unit**2 + g%bZ_unit**2 ) * g%bMag

        g%gradPrlB = 0

        do w=1,size(g%pt)
            i=g%pt(w)%i
            j=g%pt(w)%j

        !do i = 1, g%nR
        !    do j = 1, g%nZ

                ! Brambillas approximation to dBdPar
                ! ----------------------------------

                distance = sqrt ( (g%R(i)-r0)**2 + g%Z(j)**2 )
                if ( distance > 0 ) then

                    g%sinTh(w) =  g%Z(j) / distance
                    g%gradPrlB(w) = g%bMag(w) / g%R(i) * abs ( g%bPol(w) * g%sinTh(w) )

                else ! on magnetic axis (r=r0, y=y0=0)

                    g%gradPrlB(w) = 0

                endif

                !if (nzfUn == 0) gradPrlB(i,j) = 1.0e-10

        !   enddo
        !enddo
        enddo

        ! R,Th,z version
        ! --------------

        !allocate ( &
        !    g%drUrr(g%nR,g%nZ), g%drrUrr(g%nR,g%nZ), &
        !    g%drUrt(g%nR,g%nZ), g%drrUrt(g%nR,g%nZ), &
        !    g%drUrz(g%nR,g%nZ), g%drrUrz(g%nR,g%nZ), &
        !    g%drUtr(g%nR,g%nZ), g%drrUtr(g%nR,g%nZ), &
        !    g%drUtt(g%nR,g%nZ), g%drrUtt(g%nR,g%nZ), &
        !    g%drUtz(g%nR,g%nZ), g%drrUtz(g%nR,g%nZ), &
        !    g%drUzr(g%nR,g%nZ), g%drrUzr(g%nR,g%nZ), &
        !    g%drUzt(g%nR,g%nZ), g%drrUzt(g%nR,g%nZ), &
        !    g%drUzz(g%nR,g%nZ), g%drrUzz(g%nR,g%nZ) )

        !allocate ( &
        !    g%dzUrr(g%nR,g%nZ), g%dzzUrr(g%nR,g%nZ), &
        !    g%dzUrt(g%nR,g%nZ), g%dzzUrt(g%nR,g%nZ), &
        !    g%dzUrz(g%nR,g%nZ), g%dzzUrz(g%nR,g%nZ), &
        !    g%dzUtr(g%nR,g%nZ), g%dzzUtr(g%nR,g%nZ), &
        !    g%dzUtt(g%nR,g%nZ), g%dzzUtt(g%nR,g%nZ), &
        !    g%dzUtz(g%nR,g%nZ), g%dzzUtz(g%nR,g%nZ), &
        !    g%dzUzr(g%nR,g%nZ), g%dzzUzr(g%nR,g%nZ), &
        !    g%dzUzt(g%nR,g%nZ), g%dzzUzt(g%nR,g%nZ), &
        !    g%dzUzz(g%nR,g%nZ), g%dzzUzz(g%nR,g%nZ) )

        !allocate ( &
        !    g%drzUrr(g%nR,g%nZ), g%drzUrt(g%nR,g%nZ), g%drzUrz(g%nR,g%nZ), &
        !    g%drzUtr(g%nR,g%nZ), g%drzUtt(g%nR,g%nZ), g%drzUtz(g%nR,g%nZ), &
        !    g%drzUzr(g%nR,g%nZ), g%drzUzt(g%nR,g%nZ), g%drzUzz(g%nR,g%nZ) ) 

        allocate ( &
            g%drUrr(size(g%pt)), g%drrUrr(size(g%pt)), &
            g%drUrt(size(g%pt)), g%drrUrt(size(g%pt)), &
            g%drUrz(size(g%pt)), g%drrUrz(size(g%pt)), &
            g%drUtr(size(g%pt)), g%drrUtr(size(g%pt)), &
            g%drUtt(size(g%pt)), g%drrUtt(size(g%pt)), &
            g%drUtz(size(g%pt)), g%drrUtz(size(g%pt)), &
            g%drUzr(size(g%pt)), g%drrUzr(size(g%pt)), &
            g%drUzt(size(g%pt)), g%drrUzt(size(g%pt)), &
            g%drUzz(size(g%pt)), g%drrUzz(size(g%pt)) )

        allocate ( &
            g%dzUrr(size(g%pt)), g%dzzUrr(size(g%pt)), &
            g%dzUrt(size(g%pt)), g%dzzUrt(size(g%pt)), &
            g%dzUrz(size(g%pt)), g%dzzUrz(size(g%pt)), &
            g%dzUtr(size(g%pt)), g%dzzUtr(size(g%pt)), &
            g%dzUtt(size(g%pt)), g%dzzUtt(size(g%pt)), &
            g%dzUtz(size(g%pt)), g%dzzUtz(size(g%pt)), &
            g%dzUzr(size(g%pt)), g%dzzUzr(size(g%pt)), &
            g%dzUzt(size(g%pt)), g%dzzUzt(size(g%pt)), &
            g%dzUzz(size(g%pt)), g%dzzUzz(size(g%pt)) )

        allocate ( &
            g%drzUrr(size(g%pt)), g%drzUrt(size(g%pt)), g%drzUrz(size(g%pt)), &
            g%drzUtr(size(g%pt)), g%drzUtt(size(g%pt)), g%drzUtz(size(g%pt)), &
            g%drzUzr(size(g%pt)), g%drzUzt(size(g%pt)), g%drzUzz(size(g%pt)) ) 


        do w=1,size(g%pt)
            i=g%pt(w)%i
            j=g%pt(w)%j

        !do i = 1, g%nR
        !    do j = 1, g%nZ

                bTmp = dlg_interpB ( (/g%R(i),0.0,g%z(j)/), bMagHere = bMagTmp )  
                bRu = bTmp(1)/bMagTmp
                bTu = bTmp(2)/bMagTmp
                bZu = bTmp(3)/bMagTmp
                Rot = RotMatHere(bRu,bTu,bZu)
      
                ! R direction derivatives
                if(i.ne.1.and.i.ne.g%nR)then
                    bTmp = dlg_interpB ( (/g%R(i-1),0.0,g%z(j)/), bMagHere = bMagTmp )  
                    bRu = bTmp(1)/bMagTmp
                    bTu = bTmp(2)/bMagTmp
                    bZu = bTmp(3)/bMagTmp
                    RotLf = RotMatHere(bRu,bTu,bZu)

                    bTmp = dlg_interpB ( (/g%R(i+1),0.0,g%z(j)/), bMagHere = bMagTmp )  
                    bRu = bTmp(1)/bMagTmp
                    bTu = bTmp(2)/bMagTmp
                    bZu = bTmp(3)/bMagTmp
                    RotRt = RotMatHere(bRu,bTu,bZu)

                    dx = ( (g%R(i)-g%R(i-1)) + (g%R(i+1)-g%R(i)) ) / 2
                    dfdx = (RotRt - RotLf) / (2.0*dx)
                    d2fdx2 = (RotRt - 2*Rot + RotLf) / dx**2
                endif
                if(i.eq.1)then

                    bTmp = dlg_interpB ( (/g%R(i+1),0.0,g%z(j)/), bMagHere = bMagTmp )  
                    bRu = bTmp(1)/bMagTmp
                    bTu = bTmp(2)/bMagTmp
                    bZu = bTmp(3)/bMagTmp
                    RotRt = RotMatHere(bRu,bTu,bZu)

                    bTmp = dlg_interpB ( (/g%R(i+2),0.0,g%z(j)/), bMagHere = bMagTmp )  
                    bRu = bTmp(1)/bMagTmp
                    bTu = bTmp(2)/bMagTmp
                    bZu = bTmp(3)/bMagTmp
                    RotRt_ = RotMatHere(bRu,bTu,bZu)

                    dx = ( (g%R(i+2)-g%R(i+1)) + (g%R(i+1)-g%R(i)) ) / 2
                    dfdx = (RotRt - Rot) / dx
                    d2fdx2 = (RotRt_ - 2.0 * RotRt + Rot) / dx**2
                endif
                if(i.eq.g%nR)then

                    bTmp = dlg_interpB ( (/g%R(i-1),0.0,g%z(j)/), bMagHere = bMagTmp )  
                    bRu = bTmp(1)/bMagTmp
                    bTu = bTmp(2)/bMagTmp
                    bZu = bTmp(3)/bMagTmp
                    RotLf = RotMatHere(bRu,bTu,bZu)

                    bTmp = dlg_interpB ( (/g%R(i-2),0.0,g%z(j)/), bMagHere = bMagTmp )  
                    bRu = bTmp(1)/bMagTmp
                    bTu = bTmp(2)/bMagTmp
                    bZu = bTmp(3)/bMagTmp
                    RotLf_ = RotMatHere(bRu,bTu,bZu)

                    dx = ( (g%R(i)-g%R(i-1)) + (g%R(i-1)-g%R(i-2)) ) / 2
                    dfdx = (Rot - RotLf) / dx
                    d2fdx2 = (Rot - 2.0 * RotLf + RotLf_) / dx**2

                endif

                ! Z direction derivatives 
                if (j .ne. 1 .and. j .ne. g%nZ)then
                    bTmp = dlg_interpB ( (/g%R(i),0.0,g%z(j-1)/), bMagHere = bMagTmp )  
                    bRu = bTmp(1)/bMagTmp
                    bTu = bTmp(2)/bMagTmp
                    bZu = bTmp(3)/bMagTmp
                    RotDn = RotMatHere(bRu,bTu,bZu)

                    bTmp = dlg_interpB ( (/g%R(i),0.0,g%z(j+1)/), bMagHere = bMagTmp )  
                    bRu = bTmp(1)/bMagTmp
                    bTu = bTmp(2)/bMagTmp
                    bZu = bTmp(3)/bMagTmp
                    RotUp = RotMatHere(bRu,bTu,bZu)

                    dy = ( (g%z(j)-g%z(j-1)) + (g%z(j+1)-g%z(j)) ) / 2
                    dfdy = (RotUp - RotDn) / (2.0*dy)
                    d2fdy2 = (RotUp - 2*Rot + RotDn) / dy**2
                endif
                if(j.eq.1)then
                    bTmp = dlg_interpB ( (/g%R(i),0.0,g%z(j+1)/), bMagHere = bMagTmp )  
                    bRu = bTmp(1)/bMagTmp
                    bTu = bTmp(2)/bMagTmp
                    bZu = bTmp(3)/bMagTmp
                    RotUp = RotMatHere(bRu,bTu,bZu)

                    bTmp = dlg_interpB ( (/g%R(i),0.0,g%z(j+2)/), bMagHere = bMagTmp )  
                    bRu = bTmp(1)/bMagTmp
                    bTu = bTmp(2)/bMagTmp
                    bZu = bTmp(3)/bMagTmp
                    RotUp_ = RotMatHere(bRu,bTu,bZu)

                    dy = ( (g%z(j+2)-g%z(j+1)) + (g%z(j+1)-g%z(j)) ) / 2
                    dfdy = (RotUp - Rot) / dy
                    d2fdy2 = (RotUp_ - 2.0 * RotUp + Rot) / dy**2
                endif
                if(j.eq.g%nZ)then
                    bTmp = dlg_interpB ( (/g%R(i),0.0,g%z(j-1)/), bMagHere = bMagTmp )  
                    bRu = bTmp(1)/bMagTmp
                    bTu = bTmp(2)/bMagTmp
                    bZu = bTmp(3)/bMagTmp
                    RotDn = RotMatHere(bRu,bTu,bZu)

                    bTmp = dlg_interpB ( (/g%R(i),0.0,g%z(j-2)/), bMagHere = bMagTmp )  
                    bRu = bTmp(1)/bMagTmp
                    bTu = bTmp(2)/bMagTmp
                    bZu = bTmp(3)/bMagTmp
                    RotDn_ = RotMatHere(bRu,bTu,bZu)

                    dy = ( (g%z(j)-g%z(j-1)) + (g%z(j-1)-g%z(j-2)) ) / 2
                    dfdy = (Rot - RotDn) / dy
                    d2fdy2 = (Rot - 2.0 * RotDn + RotDn_) / dy**2
                endif
                

                ! R/Z derivatives
                d2fdxy = 0
                if (i /= 1 .and. i /= g%nR .and. j /= 1 .and. j /= g%nZ) then

                    dx = ( (g%R(i)-g%R(i-1)) + (g%R(i+1)-g%R(i)) ) / 2
                    dy = ( (g%z(j)-g%z(j-1)) + (g%z(j+1)-g%z(j)) ) / 2

                    bTmp = dlg_interpB ( (/g%R(i+1),0.0,g%z(j+1)/), bMagHere = bMagTmp )  
                    bRu = bTmp(1)/bMagTmp
                    bTu = bTmp(2)/bMagTmp
                    bZu = bTmp(3)/bMagTmp
                    RotRt_ = RotMatHere(bRu,bTu,bZu)

                    bTmp = dlg_interpB ( (/g%R(i-1),0.0,g%z(j+1)/), bMagHere = bMagTmp )  
                    bRu = bTmp(1)/bMagTmp
                    bTu = bTmp(2)/bMagTmp
                    bZu = bTmp(3)/bMagTmp
                    RotLf_ = RotMatHere(bRu,bTu,bZu)

                    bTmp = dlg_interpB ( (/g%R(i+1),0.0,g%z(j-1)/), bMagHere = bMagTmp )  
                    bRu = bTmp(1)/bMagTmp
                    bTu = bTmp(2)/bMagTmp
                    bZu = bTmp(3)/bMagTmp
                    RotDn_ = RotMatHere(bRu,bTu,bZu)

                    bTmp = dlg_interpB ( (/g%R(i-1),0.0,g%z(j-1)/), bMagHere = bMagTmp )  
                    bRu = bTmp(1)/bMagTmp
                    bTu = bTmp(2)/bMagTmp
                    bZu = bTmp(3)/bMagTmp
                    RotUp_ = RotMatHere(bRu,bTu,bZu)

                    d2fdxy = (RotRt_ - RotLf_ - RotDn_ + RotUp_ ) &
                      / (4.0 * dx * dy)

                endif


                g%drUrr(w) = dfdx(1,1)
                g%drUrt(w) = dfdx(1,2)
                g%drUrz(w) = dfdx(1,3)

                g%drUtr(w) = dfdx(2,1)
                g%drUtt(w) = dfdx(2,2)
                g%drUtz(w) = dfdx(2,3)

                g%drUzr(w) = dfdx(3,1)
                g%drUzt(w) = dfdx(3,2)
                g%drUzz(w) = dfdx(3,3)

                g%drrUrr(w) = d2fdx2(1,1)
                g%drrUrt(w) = d2fdx2(1,2)
                g%drrUrz(w) = d2fdx2(1,3)

                g%drrUtr(w) = d2fdx2(2,1)
                g%drrUtt(w) = d2fdx2(2,2)
                g%drrUtz(w) = d2fdx2(2,3)

                g%drrUzr(w) = d2fdx2(3,1)
                g%drrUzt(w) = d2fdx2(3,2)
                g%drrUzz(w) = d2fdx2(3,3)

                g%dzUrr(w) = dfdy(1,1)
                g%dzUrt(w) = dfdy(1,2)
                g%dzUrz(w) = dfdy(1,3)

                g%dzUtr(w) = dfdy(2,1)
                g%dzUtt(w) = dfdy(2,2)
                g%dzUtz(w) = dfdy(2,3)

                g%dzUzr(w) = dfdy(3,1)
                g%dzUzt(w) = dfdy(3,2)
                g%dzUzz(w) = dfdy(3,3)

                g%dzzUrr(w) = d2fdy2(1,1)
                g%dzzUrt(w) = d2fdy2(1,2)
                g%dzzUrz(w) = d2fdy2(1,3)

                g%dzzUtr(w) = d2fdy2(2,1)
                g%dzzUtt(w) = d2fdy2(2,2)
                g%dzzUtz(w) = d2fdy2(2,3)

                g%dzzUzr(w) = d2fdy2(3,1)
                g%dzzUzt(w) = d2fdy2(3,2)
                g%dzzUzz(w) = d2fdy2(3,3)

                g%drzUrr(w) = d2fdxy(1,1)
                g%drzUrt(w) = d2fdxy(1,2)
                g%drzUrz(w) = d2fdxy(1,3)

                g%drzUtr(w) = d2fdxy(2,1)
                g%drzUtt(w) = d2fdxy(2,2)
                g%drzUtz(w) = d2fdxy(2,3)

                g%drzUzr(w) = d2fdxy(3,1)
                g%drzUzt(w) = d2fdxy(3,2)
                g%drzUzz(w) = d2fdxy(3,3)

                !call deriv_x(g%R, g%Urr, i, j, dfdx = g%drUrr(i,j), d2fdx2 = g%drrUrr(i,j))
                !call deriv_x(g%R, g%Urt, i, j, dfdx = g%drUrt(i,j), d2fdx2 = g%drrUrt(i,j))
                !call deriv_x(g%R, g%Urz, i, j, dfdx = g%drUrz(i,j), d2fdx2 = g%drrUrz(i,j))

                !call deriv_x(g%R, g%Utr, i, j, dfdx = g%drUtr(i,j), d2fdx2 = g%drrUtr(i,j))
                !call deriv_x(g%R, g%Utt, i, j, dfdx = g%drUtt(i,j), d2fdx2 = g%drrUtt(i,j))
                !call deriv_x(g%R, g%Utz, i, j, dfdx = g%drUtz(i,j), d2fdx2 = g%drrUtz(i,j))

                !call deriv_x(g%R, g%Uzr, i, j, dfdx = g%drUzr(i,j), d2fdx2 = g%drrUzr(i,j))
                !call deriv_x(g%R, g%Uzt, i, j, dfdx = g%drUzt(i,j), d2fdx2 = g%drrUzt(i,j))
                !call deriv_x(g%R, g%Uzz, i, j, dfdx = g%drUzz(i,j), d2fdx2 = g%drrUzz(i,j))

                !call deriv_y(g%z, g%Urr, i, j, dfdy = g%dzUrr(i,j), d2fdy2 = g%dzzUrr(i,j))
                !call deriv_y(g%z, g%Urt, i, j, dfdy = g%dzUrt(i,j), d2fdy2 = g%dzzUrt(i,j))
                !call deriv_y(g%z, g%Urz, i, j, dfdy = g%dzUrz(i,j), d2fdy2 = g%dzzUrz(i,j))

                !call deriv_y(g%z, g%Utr, i, j, dfdy = g%dzUtr(i,j), d2fdy2 = g%dzzUtr(i,j))
                !call deriv_y(g%z, g%Utt, i, j, dfdy = g%dzUtt(i,j), d2fdy2 = g%dzzUtt(i,j))
                !call deriv_y(g%z, g%Utz, i, j, dfdy = g%dzUtz(i,j), d2fdy2 = g%dzzUtz(i,j))

                !call deriv_y(g%z, g%Uzr, i, j, dfdy = g%dzUzr(i,j), d2fdy2 = g%dzzUzr(i,j))
                !call deriv_y(g%z, g%Uzt, i, j, dfdy = g%dzUzt(i,j), d2fdy2 = g%dzzUzt(i,j))
                !call deriv_y(g%z, g%Uzz, i, j, dfdy = g%dzUzz(i,j), d2fdy2 = g%dzzUzz(i,j))

                !call deriv_xy(g%R, g%z, g%Urr, i, j, d2fdxy = g%drzUrr(i,j))
                !call deriv_xy(g%R, g%z, g%Urt, i, j, d2fdxy = g%drzUrt(i,j))
                !call deriv_xy(g%R, g%z, g%Urz, i, j, d2fdxy = g%drzUrz(i,j))

                !call deriv_xy(g%R, g%z, g%Utr, i, j, d2fdxy = g%drzUtr(i,j))
                !call deriv_xy(g%R, g%z, g%Utt, i, j, d2fdxy = g%drzUtt(i,j))
                !call deriv_xy(g%R, g%z, g%Utz, i, j, d2fdxy = g%drzUtz(i,j))

                !call deriv_xy(g%R, g%z, g%Uzr, i, j, d2fdxy = g%drzUzr(i,j))
                !call deriv_xy(g%R, g%z, g%Uzt, i, j, d2fdxy = g%drzUzt(i,j))
                !call deriv_xy(g%R, g%z, g%Uzz, i, j, d2fdxy = g%drzUzz(i,j))

        !   enddo
        !enddo
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
