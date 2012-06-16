module solve

use constants
#ifdef par
use scalapack_mod
#endif
#ifdef USE_PGESVR
use pgesvr_mod
#endif

use timer_mod
implicit none


contains

    subroutine solve_lsq ()

        use aorsaNamelist, &
        only: square, magma
        use mat_fill, &
        only: aMat!, aMat_
        use antenna, &
        only: brhs, brhs_global

        implicit none

        integer :: info, nRow, nCol
        integer, allocatable, dimension(:) :: ipiv

        !   lapack lsq solve variables
        
        integer :: M_, N_, NRHS, LDA, LDB, MN_, LWORK, RANK
        real :: RCOND
#ifndef dblprec   
        real, allocatable :: rWork(:)
        complex, allocatable :: work(:)
#else
        real(kind=dbl), allocatable :: rWork(:)
        complex(kind=dbl), allocatable :: work(:)
#endif
        integer, allocatable, dimension(:) :: jpvt

        ! Cuda
        integer(kind=long) :: k1, k2, nb
        integer, external :: magma_get_zgetrf_nb
        integer :: magma_solve, magStat
        integer :: dA_dim

        nRow    = size ( aMat, 2 ) !nPtsX * nPtsY * 3
        nCol    = size ( aMat, 1 ) !nModesX * nModesY * 3

        if(square) then

            allocate ( ipiv(nRow) )
            ipiv = 0
#ifndef dblprec
            call cgesv ( nRow, 1, aMat, nRow, ipiv, brhs, nRow, info )
#else            
            if(magma)then

                !! MAGMA Cuda solve
                !! ----------------

                !nb  = magma_get_zgetrf_nb (nRow)
                !k1  = 32 - mod ( maxVal ( (/nRow,nCol/) ), 32 )
                !k2  = 32 - mod ( nCol, 32 ) 
                !if(k1==32)k1=0
                !if(k2==32)k2=0

                !lWork = nRow * nb
                !dA_dim = ( nRow + k1 )**2 &
                !    + (nCol + k2 ) * nb + 2 * nb**2

                !magStat = magma_solve ( dA_dim, lWork, aMat, ipiv, nRow )

                !call zgetrs ( 'N', nRow, 1, aMat, nRow, ipiv, brhs, nRow, info )

            else 

                ! CPU LAPACK
                ! ----------

                call zgesv ( nRow, 1, aMat, nRow, ipiv, brhs, nRow, info )

            endif

#endif

        else

            M_  = nRow
            N_  = nCol
            NRHS   = 1
            LDA = maxVal ( (/ 1, M_ /) )
            LDB = maxVal ( (/ 1, M_, N_ /) )
            MN_ = minVal ( (/ M_, N_ /) )
            LWORK   = MN_ + maxVal ( (/ 2 * MN_, N_ + 1, MN_ + NRHS /) )
            RCOND   = 1E-12

            allocate ( WORK ( maxVal ( (/ 1, LWORK /) ) ) )
            allocate ( RWORK ( 2 * N_ ) )
            allocate ( JPVT ( N_ ) ) 

            JPVT    = 0
            WORK    = 0
            RWORK   = 0

            write(*,*) shape ( amat ), M_*N_, size(amat)
#ifndef dblprec
            call cgelsy ( M_, N_, NRHS, aMat, LDA, brhs, LDB, JPVT, RCOND, RANK, &
                                WORK, LWORK, RWORK, info )
#else
            call zgelsy ( M_, N_, NRHS, aMat, LDA, brhs, LDB, JPVT, RCOND, RANK, &
                                WORK, LWORK, RWORK, info )
#endif

        endif

        if  (info/=0) then

            write(*,*) '    solve.F90: ERROR, solve did not complete successfully'
            write(*,*) '    info: ', info
            stop

        else
            write(*,*) '    LAPACK status: ', info
        endif 
        
        deallocate ( aMat )

        brhs_global = brhs

    end subroutine solve_lsq

#ifdef par

    subroutine solve_lsq_parallel ()

        use aorsaNamelist, &
        only: npRow, npCol, square
        use mat_fill, &
        only: aMat
        use antenna, &
        only: brhs, brhs_global
        use parallel

        implicit none

        character :: trans
        integer :: m, n, nrhs
        integer :: ia, ja
        complex :: b
        integer :: ib, jb
#ifndef dblprec
        complex, allocatable :: work(:)
#else
        complex(kind=dbl), allocatable :: work(:)
#endif
        integer :: lWork, info
        integer :: ltau, lwf, lws, MpA0, NqA0, NRHSqB0, MpB0
        integer :: IROFFA, IAROW, IACOL, IROFFB, IBROW, IBCOL
        integer :: mb_a, nb_a, csrc_a, csrc_b, icoffa, mb_b
        integer :: icoffb, nb_b, Npb0, rsrc_a, rsrc_b

        integer, external :: ILCM, NUMROC, INDXG2P
        integer, allocatable :: ipiv(:)

        integer :: ii, gg
        type(timer) :: tSolve0,tSolve0_onlyfactor

#ifdef USE_GPU
        integer :: memsize,ivalue,istatus
        character(len=127) :: memsize_str
#endif
#ifdef USE_ROW_SCALING
        logical, parameter :: use_row_scaling = .true.
        real*8 :: rnorm,rnorm_inv
        integer :: iia, i, incx
        type(timer) :: tRowScale

        logical, parameter :: use_column_scaling = .false.
        integer :: j,jja,iib
        real*8 :: cnorm, cnorm_inv
        real*8, dimension(:), allocatable :: cnorm_inv_array
#endif
#ifdef USE_RCOND
        character, parameter :: norm = '1'
        integer :: lzwork,lrwork
        real*8, allocatable,dimension(:) :: rwork
        complex*16, allocatable, dimension(:) :: zwork
        complex*16 :: zwork1(1024*1024)
        real*8 :: rwork1(1024*1024)
        real*8 :: anorm, rcond
        type(timer) :: tRCond
        real :: TimeRCond

        interface
          double precision function pzlange(norm,m,n,A,ia,ja,descA,work)
          character norm
          integer ia,ja,m,n
          integer descA(*)
          double precision work(*)
          complex*16 A(*)
          end function pzlange
        end interface
#endif

        trans   = 'N'

        mb_a    = rowBlockSize
        nb_a    = colBlockSize

        mb_b    = rowBlockSize
        nb_b    = 1

        rsrc_a  = rowStartProc
        csrc_a  = colStartProc
       
        rsrc_b  = rowStartProc 
        csrc_b  = colStartProc

        m   = nRow
        n   = nCol
        nrhs    = 1
        ia  = 1!myRow * rowBlockSize + 1
        ja  = 1!myCol * colBlockSize + 1 
        ib  = 1!myRow * rowBlockSize + 1
        jb  = 1

        iroffa  = mod( ia-1, mb_a )
        icoffa  = mod( ja-1, nb_a )
        iarow   = indxg2p( ia, mb_a, myRow, rsrc_a, npRow )
        iacol   = indxg2p( ja, nb_a, myCol, csrc_a, npCol )
        MpA0    = numroc( m+iroffa, mb_a, myrow, iarow, nprow )
        NqA0    = numroc( n+icoffa, nb_a, mycol, iacol, npcol )

        iroffb  = mod( ib-1, mb_b )
        icoffb  = mod( jb-1, nb_b )
        ibrow   = indxg2p( ib, mb_b, myRow, rsrc_b, npRow )
        ibcol   = indxg2p( jb, nb_b, myCol, csrc_b, npCol )
        Mpb0    = numroc( m+iroffb, mb_b, myRow, ibrow, npRow )
        Npb0    = numroc( n+iroffb, mb_b, myRow, ibrow, npRow )
        nrhsqb0 = numroc( nrhs+icoffb, nb_b, myCol, ibcol, npCol )

        ltau    = numroc( ja+min(m,n)-1, nb_a, myCol, csrc_a, npCol )
        lwf     = nb_a * ( Mpa0 + Nqa0 + nb_a )
        lws     = max( (nb_a*(nb_a-1))/2, (nrhsqb0 + Mpb0)*nb_a ) + &
                    nb_a * nb_a

        lWork   = ltau + max ( lwf, lws )


        ! Least squares or LU decomp for square array
        ! -------------------------------------------

        if(square) then 

            if(iAm==0) write(*,*) &
            "   using LU decomposition for square system"

            allocate ( ipiv ( numroc ( m, mb_a, myRow, rsrc_a, npRow ) + mb_a ) )
            ipiv = 0

#ifdef USE_ROW_SCALING
             if (use_row_scaling) then
               if(iAm==0)write(*,*) '		Using row scaling'
               call start_timer( tRowScale )

               do i=1,n
                 iia = (ia-1) + i
                 incx = descriptor_aMat(M_)
                 call pdznrm2( n,rnorm,aMat,iia,ja,descriptor_aMat, incx )

                 rnorm_inv = 1.0
                 if (abs(rnorm) .gt. epsilon(rnorm)) then
                     rnorm_inv = 1.0/rnorm
                 endif
                 call pzdscal(n, rnorm_inv, aMat,iia,ja,descriptor_aMat,incx)

                 incx = descriptor_brhs(M_)
                 iib = (ib-1) + i
                 call pzdscal(nrhs,rnorm_inv, brhs,iib,jb,descriptor_brhs,incx)
                enddo

                if (iAm == 0) then
                  write(*,*) 'Row Scaling took ', end_timer( tRowScale )
                endif
               endif

               if (use_column_scaling) then
                 if(iAm==0)write(*,*) '		Using col scaling'
                 allocate( cnorm_inv_array(n) )
                 do j=1,n
                   jja = (ja-1) + j
                   incx = 1
                   call pdznrm2(n,cnorm,aMat,ia,jja,descriptor_aMat,incx)

                   cnorm_inv = 1.0
                   if (abs(cnorm) .gt. epsilon(cnorm)) then
                       cnorm_inv = 1.0/cnorm
                   endif
                   cnorm_inv_array(j) = cnorm_inv

                   incx = 1
                   call pzdscal(n,cnorm_inv, aMat,ia,jja,descriptor_aMat,incx)
                  enddo
                 endif
#endif

#ifdef USE_RCOND
! -------------------------------------------
! compute norm(A) needed later in estimating
! condition number
! -------------------------------------------
             TimeRCond = 0
             call start_timer(tRCond)
             lrwork = 2*n
             allocate(rwork(lrwork))
             anorm = pzlange(norm,n,n,aMat,ia,ja,descriptor_aMat,rwork)
             deallocate(rwork)
             TimeRCond = end_timer(tRCond)
#endif

#ifndef USE_GPU

#ifndef dblprec
            call pcgesv ( n, nrhs, aMat, ia, ja, descriptor_aMat, ipiv, &
                    brhs, ib, jb, descriptor_brhs, info ) 
#else
#ifdef USE_PGESVR
            if(iAm==0)write(*,*) '		Using iterative refinement PZGESVR'
            call start_timer(tSolve0)
            call pzgesvr ( n, nrhs, aMat, ia, ja, descriptor_aMat, ipiv, &
                    brhs, ib, jb, descriptor_brhs, info ) 
            if(iAm==0)write(*,*) '    Time to solve (0): ', end_timer(tSolve0)
#else
            if(iAm==0)write(*,*) '		Using regular PZGESV complex*16'
            call start_timer(tSolve0)
            call pzgesv ( n, nrhs, aMat, ia, ja, descriptor_aMat, ipiv, &
                    brhs, ib, jb, descriptor_brhs, info ) 
            if(iAm==0)write(*,*) '    Time to solve (0): ', end_timer(tSolve0)
#endif
#endif

#else
! use GPU version
! 4 GBytes on  1 MPI task, memsize = 256*1024*1024
! 4 GBytes on  4 MPI tasks, memsize = 64*1024*1024
! 4 GBytes on  8 MPI tasks, memsize = 32*1024*1024
! 4 GBytes on 16 MPI tasks, memsize = 16*1024*1024
!

             memsize = 16*1024*1024
             call getenv("MEMSIZE",memsize_str)
             if (len(trim(memsize_str)).ge.1) then
               ivalue = 32
               read(memsize_str,*,iostat=istatus) ivalue
               if ((istatus.eq.0).and.(ivalue.ge.1)) then
                 memsize = ivalue*1024*1024
               endif
             endif
             if (iAm == 0) then
               write(*,*) 'memsize = ', memsize
             endif


#ifdef USE_PGESVR
            if(iAm==0)write(*,*) '		Using iterative refinement PZGESVR on GPU'
            call start_timer(tSolve0)
            call pzgesvr ( n, nrhs, aMat, ia, ja, descriptor_aMat, ipiv, &
                    brhs, ib, jb, descriptor_brhs, info ) 
            if(iAm==0)write(*,*) '    Time to solve (0): ', end_timer(tSolve0)
#else

!debug       call pzgetrf(n,n,aMat, ia,ja,descriptor_aMat,ipiv,  info)
            if(iAm==0)write(*,*) '		Using OOC PZGESV complex*16 on GPU'
            call start_timer(tSolve0)
            call start_timer(tSolve0_onlyfactor)
            call pzgetrf_ooc2(n,n,aMat, ia,ja,descriptor_aMat,ipiv,  &
                     memsize, info )
            if(iAm==0)write(*,*) '    Time to solve (0 - onlyfactor): ', &
                end_timer(tSolve0_onlyfactor)


             if ((info.ne.0) .and. (iAm == 0)) then
               write(*,*) 'pzgetrf_ooc status: ',info
             endif

            call pzgetrs( 'N',n,nrhs,aMat,ia,ja,descriptor_aMat,ipiv, &
                    brhs, ib, jb, descriptor_brhs, info)
            if(iAm==0)write(*,*) '    Time to solve (0): ', end_timer(tSolve0)

             if ((info.ne.0) .and. (iAm == 0)) then
               write(*,*) 'pzgetrs status: ',info
             endif
#endif



#endif

#ifdef USE_ROW_SCALING
             if (use_column_scaling) then
               do i=1,n
                  iib = (ib-1) + i
                  cnorm_inv = cnorm_inv_array(i)
                  incx = descriptor_brhs(M_)
                  call pzdscal( nrhs, cnorm_inv, &
                         brhs,iib,jb,descriptor_brhs,incx)
                enddo
                deallocate( cnorm_inv_array )
              endif
#endif

#ifdef USE_RCOND
             if(iAm==0)write(*,*) '        Estimating rCond'
             call start_timer(tRCond)
             lzwork = -1
             lrwork = -1
             call pzgecon( norm,n, aMat,ia,ja,descriptor_aMat,  &
                  anorm,rcond, zwork1,lzwork,rwork1,lrwork,info)

             lzwork = int(abs(zwork1(1))) + 1
             lrwork = int( abs(rwork1(1)) ) + 1
             allocate( zwork(lzwork), rwork(lrwork) )
             call pzgecon( '1',n, aMat,ia,ja,descriptor_aMat,  &
                  anorm,rcond, zwork,lzwork,rwork,lrwork,info)
             deallocate( zwork, rwork )

             if ((info.ne.0) .and. (iAm == 0)) then
               write(*,*) 'pzgecon status ',info
             endif
             if (iAm == 0) then
               write(*,*) 'estimated rcond = ',rcond
             endif
             TimeRCond = TimeRCond + end_timer(tRCond)
             if(iAm==0)write(*,*) '    Time to estimate RCond: ', TimeRCond

#endif


            if (iAm==0) then 
                write(*,*) '    pcgesv status: ', info
            endif

            deallocate ( ipiv )


          else ! if (square)
            if(iAm==0) write(*,*) &
            "   using least squares for n x m system"

            allocate ( work ( lWork ) )
            work = 0

#ifndef dblprec
            call pcgels ( trans, m, n, nrhs, aMat, ia, ja, descriptor_aMat, &
                        brhs, ib, jb, descriptor_brhs, work, lWork, info ) 
#else
            call pzgels ( trans, m, n, nrhs, aMat, ia, ja, descriptor_aMat, &
                        brhs, ib, jb, descriptor_brhs, work, lWork, info ) 
#endif
            if (iAm==0) then 
                write(*,*) '    pcgels status: ', info
                write(*,*) '    actual/optimal lwork: ', lwork, real(work(1))
            endif

            deallocate ( work )

          endif  ! if (square)
     
        if  (info/=0) then

            write(*,*) 'solve.F90: ERROR, solve did not complete successfully'
            stop

        endif 
        
        call gather_coeffs ()

    end subroutine solve_lsq_parallel


    subroutine gather_coeffs ()

        use aorsaNamelist, &
        only: npRow, npCol
        use antenna, &
        only: brhs, brhs_global
        use parallel

        implicit none

        integer :: ii, gg

        !   Gather the solution vector from all processors by
        !   creating a global solution vector of the full, not
        !   local size, and have each proc fill in its piece.
        !   Then do a global sum on all procs such that all procs
        !   will have a complete copy of the solution.
        !   NOT SURE IF THIS WILL WORK FOR A NON 2x2 grid yet.

        if (myCol==0) then 
            do ii=1,nRowLocal

                gg   =  myRow*rowBlockSize+1 &
                        +mod(ii-1,rowBlockSize) &
                        +(ii-1)/rowBlockSize * npRow*rowBlockSize
                brhs_global(gg) = brhs(ii)

                !write(*,*) iAm, brhs(ii), brhs_global(gg)

            enddo
        endif

        call cgSum2D ( iContext, 'All', ' ', nRow, 1, &
                brhs_global, nRow, -1, -1 )

    end subroutine gather_coeffs

#endif

    subroutine extract_coeffs ( g )

        use grid
        use antenna, &
        only: brhs_global
 
        implicit none

        type(gridBlock), intent(inout) :: g

        integer :: m, n, iCol

        !   Extract the k coefficents from the solution
        !   for each field component
        !   -------------------------------------------

        allocate ( &
            g%ealphak(g%nMin:g%nMax,g%mMin:g%mMax), &
            g%ebetak(g%nMin:g%nMax,g%mMin:g%mMax), &
            g%eBk(g%nMin:g%nMax,g%mMin:g%mMax) )

        do n=g%nMin,g%nMax
            do m=g%mMin,g%mMax

                iCol = (m-g%mMin) * 3 + (n-g%nMin) * g%nModesZ * 3 + 1
                iCol = iCol + ( g%StartCol-1 )
       
                g%ealphak(n,m)    = brhs_global(iCol+0)
                g%ebetak(n,m)     = brhs_global(iCol+1)
                g%eBk(n,m)        = brhs_global(iCol+2)

            enddo
        enddo

    end subroutine extract_coeffs

end module solve
