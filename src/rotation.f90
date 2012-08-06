modUle rotation

implicit none

real :: sqx

contains

    function Cross (v1,v2)

        real, intent(in) :: v1(:), v2(:)
        real :: Cross(0:2)

        Cross(0) = +(v1(1)*v2(2)-v1(2)*v2(1))
        Cross(1) = -(v1(0)*v2(2)-v1(2)*v2(0))
        Cross(2) = +(v1(0)*v2(1)-v1(1)*v2(0)) 

    end function Cross
    
    function Dot (v1,v2)

        real, intent(in) :: v1(:), v2(:)
        real :: Dot

        Dot = v1(0)*v2(0)+v1(1)*v2(1)+v1(2)*v2(2)

    end function Dot

    function RotMatHere (bR_unit,bT_unit,bZ_unit)

        use constants

        implicit none

        real, intent(in) :: bR_unit, bT_unit, bZ_unit

        real :: RotMatHere(3,3)

        ! Alternate approach
        real :: theta
        real :: q0,q1,q2,q3,InvRotQ(3,3)

        real :: zu_rtz(0:2), au_rtz(0:2), bu_rtz(0:2), pu_rtz(0:2)

        ! Vectors are:
        !
        ! a,b,p: alpha, beta, parallel
        ! r,t,z: cylindical, right handed
        !
        ! z cross b = a 
        !
        ! u are unit vectors

        ! Get vector perp to both z axis and b

        zu_rtz = 0
        zu_rtz(2) = 1

        pu_rtz = (/ bR_unit, bT_unit, bZ_unit /)

        au_rtz = Cross(zu_rtz,pu_rtz)
        bu_rtz = Cross(pu_rtz,au_rtz)

        ! Get angle between z axis and b

        theta = aCos ( Dot(zu_rtz,pu_rtz) )

        ! Calculate the quaternions

        q0  = cos ( theta / 2.0 )
        q1  = sin ( theta / 2.0 ) * bu_rtz(0)
        q2  = sin ( theta / 2.0 ) * bu_rtz(1) 
        q3  = sin ( theta / 2.0 ) * bu_rtz(2)

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
        real :: RotQ(3,3)
        real :: det

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

        allocate ( g%U_RTZ_to_ABb(size(g%pt),3,3) )

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

            RotQ = RotMatHere (g%bR_unit(w),g%bT_unit(w),g%bZ_unit(w))

            !det = RotQ(1,1) * RotQ(2,2) * RotQ(3,3) &
            !    + RotQ(1,2) * RotQ(2,3) * RotQ(3,1) &
            !    + RotQ(1,3) * RotQ(2,1) * RotQ(3,2) &
            !    - RotQ(1,1) * RotQ(2,3) * RotQ(3,2) &
            !    - RotQ(1,2) * RotQ(2,1) * RotQ(3,3) &
            !    - RotQ(1,3) * RotQ(2,2) * RotQ(3,1)
            !write(*,*) det

            g%U_RTZ_to_ABb(w,1:3,1:3)  = RotQ

        enddo

        g%Urr = g%U_RTZ_to_ABb(:,1,1) 
        g%Urt = g%U_RTZ_to_ABb(:,1,2) 
        g%Urz = g%U_RTZ_to_ABb(:,1,3) 

        g%Utr = g%U_RTZ_to_ABb(:,2,1) 
        g%Utt = g%U_RTZ_to_ABb(:,2,2) 
        g%Utz = g%U_RTZ_to_ABb(:,2,3) 

        g%Uzr = g%U_RTZ_to_ABb(:,3,1) 
        g%Uzt = g%U_RTZ_to_ABb(:,3,2) 
        g%Uzz = g%U_RTZ_to_ABb(:,3,3)

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
 
        allocate ( g%gradPrlB(size(g%pt)), g%sinTh(size(g%pt)) )
        allocate ( g%bPol(size(g%pt)) )

        ! NOTE: CHECK THIS, i.e., unit or absolute?
        g%bPol = sqrt ( g%bR_unit**2 + g%bZ_unit**2 ) * g%bMag

        g%gradPrlB = 0

        do w=1,size(g%pt)
            i=g%pt(w)%i
            j=g%pt(w)%j

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

        enddo

        ! R,Th,z version
        ! --------------

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

            bTmp = dlg_interpB ( (/g%R(i),0.0,g%z(j)/), bMagHere = bMagTmp )  
            bRu = bTmp(1)/bMagTmp
            bTu = bTmp(2)/bMagTmp
            bZu = bTmp(3)/bMagTmp
            Rot = RotMatHere(bRu,bTu,bZu)
      
            ! R direction derivatives
            dfdx = 0
            d2fdx2 = 0
            if(i.ne.1.and.i.ne.g%nR.and.g%nR>1)then
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
            if(i.eq.1.and.g%nR>1)then

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
            if(i.eq.g%nR.and.g%nR>1)then

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
            dfdy = 0
            d2fdy2 = 0
            if (j .ne. 1 .and. j .ne. g%nZ .and. g%nZ>1)then
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
            if(j==1.and.g%nZ>1)then
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
            if(j==g%nZ.and.g%nZ>1)then
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
            if (i /= 1 .and. i /= g%nR .and. j /= 1 .and. j /= g%nZ .and.g%nR>1.and.g%nZ>1) then

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
             
        enddo

    end subroutine deriv_rotation

end module rotation
