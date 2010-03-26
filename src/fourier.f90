module fourier_mod

contains

        subroutine sft_inv_t_right( n,nseq, A, lda, &
                nkx1,nkx2, nnodex, xx_save )
        implicit none

        integer n,nseq,lda
        complex A(lda,*)
        integer nkx1,nkx2,nnodex
        complex xx_save(1:(nkx2-nkx1+1))
!
        integer i,j, mshift
        complex cvec(n)
!       -------------------------------------
!       compute inv_xx_t(i,n) * A(1:n,1:nseq)
!
!       inv_xx_t is (1/n) * inv(xx)
!       inv(xx) = inv( D * inv(S) * B * S )
!               = inv(S) * inv(B) * S * inv(D)
!               = inv(S) * F * S * inv(D)
!
!       where S is cshift, F is forward fourier transform
!       D is twiddle factors
!       -------------------------------------

        mshift = -nkx1
        do j=1,nseq
          do i=1,n
            cvec(i) = A(i,j) * xx_save(i)
          enddo
          cvec = cshift( cvec, mshift )
          do i=1,n
            A(i,j) = cvec(i)
          enddo
        enddo

        call vzfftf( n,nseq, A, lda )

        do j=1,nseq
           do i=1,n
            cvec(i) = A(i,j)
           enddo
           cvec = cshift(cvec, -mshift )
           do i=1,n
            A(i,j) = cvec(i)
           enddo
        enddo

        return
        end subroutine sft_inv_t_right

        subroutine sft_inv_left( n, nseq, A, lda, &
                nkx1,nkx2, nnodex, xx_save )
        implicit none

        integer n,nseq,lda
        integer nkx1,nkx2,nnodex
        complex A(lda,*)
        complex xx_save(1:(nkx2-nkx1+1))
!
        integer i,j
        complex A_t( n, nseq )
!       -------------------------------------
!       compute Z(1:nseq,1:n) = A(1:nseq, 1:n) * inv_xx(i,n)
!
!       this is similar to computing
!            A_t(1:n,1:nseq) = transpose(A(1:nseq,1:n))
!
!            Z_t = inv_xx_t(i,n) * A_t(1:n,1:nseq)
!
!       ----------------------------------------



!       ------------------
!       transpose of input
!       ------------------
        do j=1,nseq
        do i=1,n
          A_t(i,j) = A(j,i)
        enddo
        enddo

        call sft_inv_t_right( n,nseq, A_t, size(A_t,1), &
                nkx1,nkx2, nnodex, xx_save )

!       -------------------
!       transpose of output
!       -------------------
        do i=1,n
        do j=1,nseq
          A(j,i) = A_t(i,j)
        enddo
        enddo


        return
        end subroutine sft_inv_left


        subroutine sft_setup( xx, nkdim1,nkdim2, nxdim, &
                nkx1,nkx2, nnodex, xx_save )
        implicit none

        integer nkdim1,nkdim2,nxdim, nnodex
        integer nkx1,nkx2
        complex xx(nkdim1 : nkdim2, nxdim)
        complex xx_save( 1:(nkx2-nkx1+1) )
!
!       --------------------------
!       precompute twiddle factors
!
!       xx * v = (D * inv(S) * B * S ) * v
!       where S is a shift, B is backward fourier transform
!       and  D contain the twiddle factors
!
!       --------------------------
        real ar,ac
        integer n,nseq,lda
        complex, dimension(:,:), allocatable :: A, XA

        integer i,j, mshift
        complex tmp
        complex inv_xx(nkx1:nkx2,1:nnodex)
        complex inv_xx_t(1:nnodex,nkx1:nkx2)

        integer, parameter :: idebug = 0

        complex v(1:nnodex), xxv(1:(nkx2-nkx1+1))
!      complex v(nnodex), xxv(nkx2-nkx1+1)

        complex zi
        real err, tol
        logical isok

        tol = 1.0e-6
        zi = cmplx(0.0,1.0)

        isok = nnodex .eq. (nkx2-nkx1+1)
        if (.not.isok) then
           write(*,*) 'sft_setup: error with nkx1,nkx2,nnodex '
           write(*,*) 'nkx1,nkx2,nnodex ', &
                  nkx1,nkx2,nnodex
           stop '** error ** '
        endif


!       ---------------------
!       some arbitrary values
!       ---------------------
        n = nnodex
        do i=1,nnodex
          v(i) = real(i)/real(n) + zi * real(n-i+1)/real(n)
        enddo


!       ----------------------------------------------
!efd    internal compiler error with intel 8 compiler
!       rewrite matmul as do loop
!       ----------------------------------------------
!       xxv(1:(nkx2-nkx1+1))  =
!     &        matmul( xx(nkx1:nkx2,1:nnodex), v(1:nnodex) )


        do i=nkx1,nkx2
           tmp = 0.0
           do j=1,nnodex
             tmp = tmp + xx(i,j) * v(j)
           enddo
           xxv(i-nkx1+1) = tmp
        enddo

!       --------------------------------------------
!       note cshift in f90 is the opposite direction
!       of cshift in matlab
!       --------------------------------------------
        mshift = -nkx1
        v = cshift( v, mshift )
        call vzfftb( nnodex, 1, v, nnodex )
        v = cshift( v, -mshift )
        do i=1,(nkx2-nkx1+1)
          xx_save(i) = xxv(i) / v(i)
        enddo

!       ---------------------------------------
!       store the reciprical of twiddle factors
!       to avoid very costly complex divisions
!       ---------------------------------------
        do i=1,(nkx2-nkx1+1)
           xx_save(i) = 1.0/xx_save(i)
        enddo

        if (idebug.ge.1) then

!       ------------
!       double check
!       ------------
        do i=1,n
          v(i) = exp( zi*real(i)/real(n) )
        enddo
        xxv(1:(nkx2-nkx1+1)) = &
              matmul( xx(nkx1:nkx2,1:nnodex), v(1:nnodex))
        v = cshift( v, mshift )
        call vzfftb( nnodex, 1, v, nnodex )
        v = cshift( v, -mshift )
        do i=1,(nkx2-nkx1+1)
           v(i) = v(i) / xx_save(i)
        enddo

        err = maxval(abs(xxv(1:(nkx2-nkx1+1)) - v(1:(nkx2-nkx1+1))))
        if (idebug.ge.1) then
          write(*,*) 'err for xxv and v is ', err
        endif

        if (err.gt. tol) then
           write(*,9100)  err
 9100      format('sft_setup: err = ', 1pe14.4)

           do i=1,(nkx2-nkx1+1)
            write(*,*) ' i, xxv(i), v(i) ', i,xxv(i),v(i)
           enddo

           do i=1,(nkx2-nkx1+1)
            write(*,*) 'i,xx_save(i) ', i, xx_save(i)
           enddo
        endif



!       --------------
!       further checks with sft_inv_t_right
!       --------------
        do j=1,nnodex
        do i=nkx1,nkx2
           inv_xx(i,j) = 1.0/xx(i,j)
        enddo
        enddo

        do j=1,nnodex
        do i=nkx1,nkx2
           inv_xx_t(j,i) = inv_xx(i,j)
        enddo
        enddo

        n = nkx2-nkx1+1
        nseq = n
        lda = n
        allocate( A(lda, nseq ) )
        allocate( XA(lda, nseq ) )

        do j=1,nseq
        do i=1,n
           ar = real(i+(j-1)*n)/real(n*nseq)
           ac = real(j + (i-1)*nseq)/real(n*nseq)
          A(i,j) = exp(ar + zi*ac)
        enddo
        enddo

        XA(1:nnodex,1:nseq) = &
             matmul( inv_xx_t(1:nnodex,nkx1:nkx2), A(1:n,1:nseq) )


        call sft_inv_t_right( n,nseq, A, lda, &
                nkx1,nkx2, nnodex, xx_save )

        err = maxval( abs(XA(1:nnodex,1:nseq) - A(1:nnodex,1:nseq)) )
        if (idebug.ge.1) then
          write(*,*) 'err with sft_inv_xx_t_right ', err
        endif

        if (err .gt. tol) then
           write(*,*) 'sft_setup: err with sft_inv_t_right is ', err

           do j=1,nseq
           do i=1,nnodex
            err = abs( XA(i,j) - A(i,j) )

            write(*, 9010) i,j, &
              real(XA(i,j)),aimag(XA(i,j)), &
              real(A(i,j)), aimag(A(i,j)),err
 9010       format(' i,j ',2(1x,i3), &
              ' XA ',2(1x,1pe11.3), &
              ' A ',2(1x,1pe11.3),' err ',1pe11.3)
           enddo
           enddo

           do i=nkx1,nkx2
           do j=1,nnodex
             write(*,*) 'j,i ', j,i, 'inv_xx_t(j,i) ', inv_xx_t(j,i)
           enddo
           enddo

           stop '** stop for error ** '
        endif


        deallocate( A )
        deallocate( XA )

!       --------------------------------
!       further checks with sft_inv_left
!       --------------------------------
        nseq = nkx2-nkx1+1
        lda = nseq
        allocate( A(nseq, n), XA(nseq,n) )

        do j=1,n
        do i=1,nseq
           ar = real(i+(j-1)*n)/real(n*nseq)
           ac = real(j + (i-1)*nseq)/real(n*nseq)
          A(i,j) = ar + zi*ac
        enddo
        enddo

        XA = matmul( A(1:nseq,1:n), inv_xx(nkx1:nkx2, 1:nnodex) )

        call sft_inv_left( n,nseq, A, lda, &
                   nkx1,nkx2, nnodex, xx_save )

        err = maxval( abs( XA(1:nseq,1:n) - A(1:nseq,1:n) ) )
        if (idebug.ge.1) then
          write(*,*) 'err with sft_inv_left ', err
        endif


        if (err.gt.tol) then
           write(*,*) 'sft_setup: err with sft_inv_left is ', err
           stop '** stop for error ** '
        endif

        deallocate( A )
        deallocate( XA )

        endif

        return
        end subroutine sft_setup



      subroutine vzfftf( n,nseq, C, ldc )
!       -------------------------------------------
!       forward transform, using zfftf in dfftpack
!       -------------------------------------------
      integer  n,nseq,ldc
      complex*8 C(ldc,*)

      real*8 wsave(4*(4*n+15))
      complex*16 cvec(n)
      integer  i,j

        call zffti( n, wsave )
      do j=1,nseq
        do i=1,n
          cvec(i) = C(i,j)
        enddo
        call zfftf( n, cvec, wsave )
        do i=1,n
          C(i,j) = cvec(i)
        enddo
      enddo

      return
      end subroutine vzfftf

      subroutine vzfftb( n,nseq, C, ldc )
!       -------------------------------------------
!       backward transform, using zfftb in dfftpack
!       -------------------------------------------
      integer  n,nseq,ldc
      complex*8 C(ldc,*)

      real*8 wsave(2*(4*n+15))
      complex*16 cvec(n)
      integer  i,j

        call zffti( n, wsave )
      do j=1,nseq
        do i=1,n
          cvec(i) = C(i,j)
        enddo
        call zfftb( n, cvec, wsave )
        do i=1,n
          C(i,j) = cvec(i)
        enddo
      enddo

      return
      end subroutine vzfftb



      subroutine convert2d_row(row_in, row_out, &
         xlen, ylen, nnodex, nnodey, &
         nxdim, nydim, nkdim1, nkdim2, mkdim1, mkdim2, &
         nkx1, nkx2, nky1, nky2, xx, yy, dx, dy, ndfmax, &
         nmodesx, nmodesy,  use_fft)

      implicit none

      integer nnodex, nnodey, nxdim, nydim, nkdim1, nkdim2, &
         mkdim1, mkdim2, i, j, n, m, nkx1, nkx2, nky1, nky2, &
         ndfmax, ncol, nmodesx, nmodesy

      real xlen, ylen, dx, dy


      complex xx(nkdim1 : nkdim2, nxdim), yy(mkdim1 : mkdim2, nydim)
      complex row_in(ndfmax), row_out(ndfmax)

!       ---------------
!       local variables
!       ---------------
        logical::  use_fft
        complex, dimension(:), allocatable :: xx_save, yy_save

        logical, parameter :: use_gemm = .true.
        integer :: mm,nn,kk, ld1,ld2,ld3
        complex :: alpha,beta

        integer :: irnc
        integer :: j_lo,j_hi,j_size
        integer :: i_lo,i_hi,i_size
        integer :: m_lo,m_hi,m_size
        integer :: n_lo,n_hi,n_size

        real :: dx_over_xlen, dy_over_ylen

        complex, dimension(:,:), allocatable :: inv_xx, inv_xx_t
        complex, dimension(:,:), allocatable :: inv_yy, inv_yy_t
        complex, dimension(:,:), allocatable :: ealpha, ealphak
        complex, dimension(:,:), allocatable :: ebeta, ebetak
        complex, dimension(:,:), allocatable :: eb, ebk
        complex, dimension(:,:), allocatable :: tmp_i_m

!       -------------------------------------------------------
!       Forward transformation used in sft2d_mat from
!       (i,j) configuration space to (n,m) fourier space is
!
!       F( (n,m); (i,j) )  = kron( inv_yy(m,j),  inv_xx(n,i) )
!
!        ealphak(n,m) = F( (n,m);(i,j) ) * ealpha(i,j)
!
!       Now we need the action of transpose(F) * vector
!
!       Let F_t( (i,j); (n,m) ) = transpose( F( (n,m); (i,j) )
!
!       Let inv_xx_t(i,n) = transpose( inv_xx(n,i) )
!
!       Let inv_yy_t(j,m) = transpose( inv_yy(m,j) )
!
!       Then
!
!       F_t(  (i,j);(n,m) )  = kron(inv_yy_t(j,m), inv_xx_t(i,n))
!
!       ealpha(i,j) = F_t( (i,j);(n,m) ) * ealphak(n,m)
!       -------------------------------------------------------

        i_lo = 1
        i_hi = nnodex
        i_size = i_hi-i_lo+1

        j_lo = 1
        j_hi = nnodey
        j_size = j_hi-j_lo+1


        n_lo = nkx1
        n_hi = nkx2
        n_size = n_hi-n_lo+1

        m_lo = nky1
        m_hi = nky2
        m_size = m_hi-m_lo+1

!       -----------------
!       initialize arrays
!       -----------------


      allocate( tmp_i_m(i_lo:i_hi,m_lo:m_hi) )

      allocate( ealphak(n_lo:n_hi,m_lo:m_hi) )
      allocate( ebetak(n_lo:n_hi,m_lo:m_hi) )
      allocate( ebk(n_lo:n_hi,m_lo:m_hi) )

      allocate( ealpha(i_lo:i_hi,j_lo:j_hi) )
      allocate( ebeta(i_lo:i_hi,j_lo:j_hi) )
      allocate( eb(i_lo:i_hi,j_lo:j_hi) )


      dx_over_xlen = dx/xlen
      dy_over_ylen = dy/ylen

      if (use_fft) then

!       ----------------------
!       initialization for FFT
!       ----------------------
          allocate( xx_save( n_size + i_size ) )
          allocate( yy_save( m_size + j_size ) )

          xx_save(:) = 0.0
          yy_save(:) = 0.0

          call sft_setup( xx, nkdim1,nkdim2, nxdim, &
                nkx1,nkx2, nnodex, xx_save )

          call sft_setup( yy, mkdim1,mkdim2, nydim, &
                nky1,nky2, nnodey, yy_save )
      else

      allocate( inv_xx(n_lo:n_hi,i_lo:i_hi) )
      allocate( inv_xx_t(i_lo:i_hi,n_lo:n_hi) )

      allocate( inv_yy(m_lo:m_hi,j_lo:j_hi) )
      allocate( inv_yy_t(j_lo:j_hi,m_lo:m_hi) )



      do i=i_lo,i_hi
      do n=n_lo,n_hi
        inv_xx(n,i) = (dx_over_xlen) / xx(n,i)
        inv_xx_t(i,n) = inv_xx(n,i)
      enddo
      enddo


      do j=j_lo,j_hi
      do m=m_lo,m_hi
         inv_yy(m , j) = (dy_over_ylen) / yy(m,j)
         inv_yy_t(j,m) = inv_yy(m,j)
      enddo
      enddo

      endif


!       ------------------------------
!       split row_in into  components
!       ealphak(n,m), ebetak(n,m), ebk(n,m)
!       ------------------------------

      do m = m_lo,m_hi
      do n = n_lo,n_hi

            irnc = (m - m_lo) * 3 + (n - n_lo) * (3*m_size) + 1

            ealphak(n,m) = row_in(irnc)
            ebetak(n,m)  = row_in(irnc + 1)
            ebk(n,m)     = row_in(irnc + 2)

      end do
      end do


!       ----------------------------------------------------------
!       F_t(  (i,j);(n,m) )  = kron(inv_yy_t(j,m), inv_xx_t(i,n))
!
!       ealpha(i,j) = F_t( (i,j);(n,m) ) * ealphak(n,m)
!
!                   = inv_xx_t(i,n) * ealphak(n,m) * inv_yy(m,j)
!       -----------------------------------------------------------






!       --------------------
!       compute with ealphak
!       --------------------



!       ---------------------------------
!       compte tmp_i_m = matmul( inv_xx * ealphak )
!       ---------------------------------

        if (use_fft) then


           tmp_i_m(i_lo:i_size,m_lo:m_hi) = ealphak(n_lo:n_hi,m_lo:m_hi)

           call sft_inv_t_right( i_size,m_size, &
                  tmp_i_m, size(tmp_i_m,1), &
                  nkx1,nkx2, nnodex, xx_save )

           do m=m_lo,m_hi
           do i=i_lo,i_hi
             tmp_i_m(i,m) = dx_over_xlen * tmp_i_m(i,m)
           enddo
           enddo

        elseif  (use_gemm) then
           mm = i_size
           nn = m_size
           kk = n_size
           ld1 = size(inv_xx_t,1)
           ld2 = size(ealphak,1)
           ld3 = size(tmp_i_m,1)
           alpha = 1.0
           beta = 0.0

           call zgemm( 'N','N', mm,nn,kk, &
                       alpha, inv_xx_t,ld1, ealphak,ld2, &
                       beta,  tmp_i_m, ld3 )

        else

        tmp_i_m(i_lo:i_hi,m_lo:m_hi) = &
          matmul( inv_xx_t(i_lo:i_hi,n_lo:n_hi), &
                  ealphak(n_lo:n_hi,m_lo:m_hi) )

        endif


!debug-begin
!       do m=m_lo,m_hi
!       do i=i_lo,i_hi
!          if (use_fft) then
!          write(*,9101) i,m,real(tmp_i_m(i,m)),aimag(tmp_i_m(i,m))
!         else
!          write(*,9102) i,m,real(tmp_i_m(i,m)),aimag(tmp_i_m(i,m))
!         endif
!       enddo
!       enddo
!
! 9101  format('fft:i,m ',2i4,' tmp_i_m ',2(1x,1pe14.4))
! 9102  format('mat:i,m ',2i4,' tmp_i_m ',2(1x,1pe14.4))
!debug-end


!       ---------------------------------
!       compute ealpha = matmul( tmp_i_m, inv_yy )
!       ---------------------------------
        if (use_fft) then

          call sft_inv_left( m_size, i_size, tmp_i_m, size(tmp_i_m,1), &
                nky1,nky2, nnodey, yy_save )

          ealpha(i_lo:i_hi,j_lo:j_hi) = &
               dy_over_ylen * tmp_i_m(i_lo:i_hi,m_lo:m_hi)

        elseif (use_gemm) then

          mm = i_size
          nn = j_size
          kk = m_size
          ld1 = size(tmp_i_m,1)
          ld2 = size(inv_yy,1)
          ld3 = size(ealpha,1)
          alpha = 1.0
          beta = 0.0

          call zgemm( 'N', 'N', mm,nn,kk, &
                alpha, tmp_i_m, ld1, inv_yy, ld2, &
                beta,  ealpha, ld3 )

        else

        ealpha(i_lo:i_hi,j_lo:j_hi) = &
          matmul( tmp_i_m(i_lo:i_hi,m_lo:m_hi), &
                  inv_yy(m_lo:m_hi,j_lo:j_hi)  )

        endif

!       -------------------
!       compute with ebetak
!       -------------------

        if (use_fft) then


           tmp_i_m(i_lo:i_size,m_lo:m_hi) = ebetak(n_lo:n_hi,m_lo:m_hi)

           call sft_inv_t_right( i_size,m_size, &
                tmp_i_m, size(tmp_i_m,1), &
                nkx1,nkx2, nnodex, xx_save )

           do m=m_lo,m_hi
           do i=i_lo,i_hi
             tmp_i_m(i,m) = dx_over_xlen * tmp_i_m(i,m)
           enddo
           enddo


        elseif (use_gemm) then

           mm = i_size
           nn = m_size
           kk = n_size
           ld1 = size(inv_xx_t,1)
           ld2 = size(ebetak,1)
           ld3 = size(tmp_i_m,1)
           alpha = 1.0
           beta = 0.0

           call zgemm( 'N','N', mm,nn,kk, &
                       alpha, inv_xx_t,ld1, ebetak,ld2, &
                       beta,  tmp_i_m, ld3 )


        else

        tmp_i_m(i_lo:i_hi,m_lo:m_hi) = &
          matmul( inv_xx_t(i_lo:i_hi,n_lo:n_hi), &
                  ebetak(n_lo:n_hi,m_lo:m_hi) )
        endif



!       ---------------------------------
!       compute ebeta = matmul( tmp_i_m, inv_yy )
!       ---------------------------------
        if (use_fft) then

          call sft_inv_left( m_size, i_size, tmp_i_m, size(tmp_i_m,1), &
                nky1,nky2, nnodey, yy_save )

          ebeta(i_lo:i_hi,j_lo:j_hi) = &
              dy_over_ylen * tmp_i_m(i_lo:i_hi,m_lo:m_hi)

        elseif (use_gemm) then


          mm = i_size
          nn = j_size
          kk = m_size
          ld1 = size(tmp_i_m,1)
          ld2 = size(inv_yy,1)
          ld3 = size(ebeta,1)
          alpha = 1.0
          beta = 0.0

          call zgemm( 'N', 'N', mm,nn,kk, &
                alpha, tmp_i_m, ld1, inv_yy, ld2, &
                beta,  ebeta, ld3 )


        else

        ebeta(i_lo:i_hi,j_lo:j_hi) = &
          matmul( tmp_i_m(i_lo:i_hi,m_lo:m_hi), &
                  inv_yy(m_lo:m_hi,j_lo:j_hi)  )

        endif



!       ----------------
!       compute with ebk
!       ----------------


        if (use_fft) then


           tmp_i_m(i_lo:i_size,m_lo:m_hi) = ebk(n_lo:n_hi,m_lo:m_hi)

           call sft_inv_t_right( i_size,m_size, &
                tmp_i_m, size(tmp_i_m,1), &
                nkx1,nkx2, nnodex, xx_save )

           do m=m_lo,m_hi
           do i=i_lo,i_hi
            tmp_i_m(i,m) = dx_over_xlen * tmp_i_m(i,m)
           enddo
           enddo

        elseif (use_gemm) then


           mm = i_size
           nn = m_size
           kk = n_size
           ld1 = size(inv_xx_t,1)
           ld2 = size(ebk,1)
           ld3 = size(tmp_i_m,1)
           alpha = 1.0
           beta = 0.0

           call zgemm( 'N','N', mm,nn,kk, &
                       alpha, inv_xx_t,ld1, ebk,ld2, &
                       beta,  tmp_i_m, ld3 )


        else

        tmp_i_m(i_lo:i_hi,m_lo:m_hi) = &
          matmul( inv_xx_t(i_lo:i_hi,n_lo:n_hi), &
                  ebk(n_lo:n_hi,m_lo:m_hi) )

        endif



!       ---------------------------------
!       compute eb = matmul( tmp_i_m, inv_yy )
!       ---------------------------------
        if (use_fft) then

          call sft_inv_left( m_size, i_size, tmp_i_m, size(tmp_i_m,1), &
                nky1,nky2, nnodey, yy_save )

          eb(i_lo:i_hi,j_lo:j_hi) = &
               dy_over_ylen * tmp_i_m(i_lo:i_hi,m_lo:m_hi)

        elseif (use_gemm) then


          mm = i_size
          nn = j_size
          kk = m_size
          ld1 = size(tmp_i_m,1)
          ld2 = size(inv_yy,1)
          ld3 = size(eb,1)
          alpha = 1.0
          beta = 0.0

          call zgemm( 'N', 'N', mm,nn,kk, &
                alpha, tmp_i_m, ld1, inv_yy, ld2, &
                beta,  eb, ld3 )


        else

        eb(i_lo:i_hi,j_lo:j_hi) = &
          matmul( tmp_i_m(i_lo:i_hi,m_lo:m_hi), &
                  inv_yy(m_lo:m_hi,j_lo:j_hi)  )

        endif
!       --------------------------------------
!       copy transformed components to row_out
!       --------------------------------------
         do j = j_lo,j_hi
         do i = i_lo,i_hi

            irnc = (i - i_lo) * (3*j_size) + (j - j_lo) * 3 + 1

            row_out(irnc)   = ealpha(i, j)
            row_out(irnc+1) = ebeta(i,j)
            row_out(irnc+2) = eb(i, j)

         end do
         end do

!       ------------------
!       deallocate storage
!       ------------------


        deallocate( ealphak )
        deallocate( ebetak )
        deallocate( ebk )

        deallocate( ealpha )
        deallocate( ebeta )
        deallocate( eb )

        deallocate( tmp_i_m )

        if (use_fft) then
          deallocate( xx_save )
          deallocate( yy_save )
        else

        deallocate( inv_xx )
        deallocate( inv_xx_t )
        deallocate( inv_yy )
        deallocate( inv_yy_t )

        endif

        return
        end subroutine convert2d_row
!
!*******************************************************************
!

      subroutine sft_row(row, rowk, xlen, ylen, nnodex, nnodey, &
         nxdim, nydim, nkdim1, nkdim2, mkdim1, mkdim2, &
         nkx1, nkx2, nky1, nky2, xx, yy, dx, dy, ndfmax, &
         nmodesx, nmodesy)

!     ----------------------------------------------
!     performs a Fourier transform on one row of the
!     collocation matrix, p_amat
!        row(ncol)  = input  vector of dimension ndfmax
!        rowk(ncol) = output vector of dimension ndfmax
!     ----------------------------------------------------

      implicit none

      integer nnodex, nnodey, nxdim, nydim, nkdim1, nkdim2, &
         mkdim1, mkdim2, i, j, n, m, nkx1, nkx2, nky1, nky2, &
         ndfmax, ncol, nmodesx, nmodesy

      real xlen, ylen, dx, dy

      complex ealphak(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              ebetak(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              ebk(nkdim1 : nkdim2, mkdim1 : mkdim2)

      complex ealpha(nxdim, nydim), &
              ebeta(nxdim, nydim), &
              eb(nxdim, nydim)

      complex xx(nkdim1 : nkdim2, nxdim), yy(mkdim1 : mkdim2, nydim)
      complex row(ndfmax), rowk(ndfmax)

      do i = 1, nnodex
         do j = 1, nnodey
            ncol = (i - 1) * 3 * nnodey + (j - 1) * 3 + 1
            ealpha(i,j) = row(ncol)
            ebeta(i,j)  = row(ncol + 1)
            eb(i,j)     = row(ncol + 2)
         end do
      end do

      call sft2d_mat(ealpha, xlen, ylen, nnodex, nnodey, &
         nxdim, nydim, nkdim1, nkdim2, mkdim1, mkdim2, &
         nkx1, nkx2, nky1, nky2, xx, yy, ealphak, dx, dy)

      call sft2d_mat(ebeta, xlen, ylen, nnodex, nnodey, &
         nxdim, nydim, nkdim1, nkdim2, mkdim1, mkdim2, &
         nkx1, nkx2, nky1, nky2, xx, yy, ebetak, dx, dy)

      call sft2d_mat(eb, xlen, ylen, nnodex, nnodey, &
         nxdim, nydim, nkdim1, nkdim2, mkdim1, mkdim2, &
         nkx1, nkx2, nky1, nky2, xx, yy, ebk, dx, dy)




      do n = nkx1, nkx2
         do m = nky1, nky2

            ncol = (n - nkx1) * 3 * nmodesy + (m - nky1) * 3 + 1

            rowk(ncol)     = ealphak(n,m)
            rowk(ncol + 1) = ebetak(n,m)
            rowk(ncol + 2) = ebk(n,m)


         end do
      end do

      return
      end subroutine sft_row

!
!*******************************************************************
!

      subroutine sft2d_mat(f, xlen, ylen, nnodex, nnodey, &
         nxdim, nydim, nkdim1, nkdim2, mkdim1, mkdim2, &
         nkx1, nkx2, nky1, nky2, xx, yy, fk, dx, dy)

      implicit none

      integer nnodex, nnodey, nxdim, nydim, nkdim1, nkdim2, &
          mkdim1, mkdim2, i, j, n, m, nkx1, nkx2, nky1, nky2

      real xlen, ylen, dx, dy

      complex fk(nkdim1 : nkdim2, mkdim1 : mkdim2), xint, cexpkxky
      complex f(nxdim, nydim)
      complex xx(nkdim1 : nkdim2, nxdim), yy(mkdim1 : mkdim2, nydim)

!       ---------------
!       local variables
!       ---------------
        logical, parameter :: use_gemm = .true.
        integer :: n_size, m_size, i_size, j_size
        integer :: mm,nn,kk,  ld1,ld2,ld3
        integer :: n_lo,n_hi, m_lo,m_hi
        complex :: alpha,beta

        complex, dimension(:,:), allocatable :: fk_tmp
        complex, dimension(:,:), allocatable :: inv_xx
        complex, dimension(:,:), allocatable :: inv_yy
        complex, dimension(:,:), allocatable :: inv_xx_f

        real :: dx_over_xlen, dy_over_ylen

!       interface
!       subroutine zgemm( transa, transb, m,n,k,
!     &         alpha, A, lda, B,ldb,
!     &         beta, C, ldc )
!
!       character transa, transb
!       integer m,n,k,  lda,ldb,ldc
!       complex*16 alpha,beta, A(*),B(*),C(*)
!       end subroutine zgemm
!       end interface



      if (use_gemm) then

      n_lo = nkx1
      n_hi = nkx2
      m_lo = nky1
      m_hi = nky2

      n_size = n_hi - n_lo + 1
      m_size = m_hi - m_lo + 1
      i_size = nnodex
      j_size = nnodey

      allocate( fk_tmp(n_size, m_size) )
      allocate( inv_xx(n_size, i_size) )
      allocate( inv_yy(m_size, j_size) )
      allocate( inv_xx_f(n_size,j_size) )

      dx_over_xlen = dx/xlen

      do i=1,nnodex
      do n=nkx1,nkx2
        inv_xx(  n-nkx1+1, i) = (dx_over_xlen) / xx(n,i)
      enddo
      enddo

      dy_over_ylen = dy/ylen

      do j=1,nnodey
      do m=nky1,nky2
         inv_yy( m-nky1+1, j) = (dy_over_ylen) / yy(m,j)
      enddo
      enddo

!       -------------------------------------------------
!       fk_tmp(n,m) =  inv_xx(n,i)* f(i,j) * inv_yy(m,j)
!       -------------------------------------------------


!       -----------------------------------------------------
!       compute inv_xx_f(n,j) = matmul( inv_xx(n,i), f(i,j) )
!       -----------------------------------------------------


        mm = n_size
        nn = j_size
        kk = i_size

        ld1 = size(inv_xx,1)
        ld2 = size(f,1)
        ld3 = size(inv_xx_f,1)

        alpha = 1.0
        beta = 0.0

        call zgemm( 'N','N', mm,nn,kk, &
                alpha, inv_xx, ld1, f, ld2, &
                beta, inv_xx_f, ld3 )


!       --------------------------------------------------
!       compute fk_tmp(n,m) = inv_xx_f(n,j) * inv_yy(m,j)
!       --------------------------------------------------
        mm = n_size
        nn = m_size
        kk = j_size

        ld1 = size(inv_xx_f,1)
        ld2 = size(inv_yy,1)
        ld3 = size(fk_tmp,1)

        alpha = 1.0
        beta = 0.0

        call zgemm( 'N', 'T', mm,nn,kk, &
            alpha,  inv_xx_f, ld1,  inv_yy, ld2, &
            beta,   fk_tmp, ld3 )


!       -----------------------
!       copy results to fk(n,m)
!       -----------------------
        do m=nky1,nky2
        do n=nkx1,nkx2
          fk(n,m) =  fk_tmp( n-nkx1+1, m-nky1+1)
        enddo
        enddo

!       -------------------
!       deallocate storage
!       -------------------
        deallocate( fk_tmp )
        deallocate( inv_xx_f )
        deallocate( inv_yy )
        deallocate( inv_xx )


      else

      do n = nkx1, nkx2
         do m = nky1, nky2

            fk(n, m) = 0.0

            do i = 1, nnodex
               do j = 1, nnodey

!                 ----------------------------------------------------
!                 cexpkxky = exp(zi * (xkx(n) * x(i) + xky(m) * y(j)))
!                 ----------------------------------------------------
                  cexpkxky = xx(n, i) * yy(m, j)

                  xint = f(i, j) / cexpkxky

                  fk(n, m) = fk(n, m) + xint * dx / xlen * dy / ylen


               end do
            end do

         end do
      end do


      endif

      return
      end subroutine sft2d_mat


!
!*******************************************************************
!

    subroutine sftinv2d( a, f, fx, fy )

        use aorsa2din_mod, &
        only: nPtsX, nPtsY, nModesX, nModesY
        use grid, &
        only: kxL, kxR, kyL, kyR, xkxsav, xkysav, xx, yy
 
        implicit none
        
        complex, intent(in) :: a(:,:)
        complex, intent(inout), optional, allocatable :: &
            f(:,:), fx(:,:), fy(:,:)

        complex :: cexpkxky
        integer :: i, j, n, m

        if (.not. allocated ( f ) ) allocate ( f(nPtsX,nPtsY) )
        if ( present ( fx ) ) then
            if (.not. allocated ( fx ) ) allocate ( fx(nPtsX,nPtsY) )
        endif
        if ( present ( fy ) ) then 
            if (.not. allocated ( fy ) ) allocate ( fy(nPtsX,nPtsY) )
        endif

        f = 0
        if (present(fx)) fx = 0
        if (present(fy)) fy = 0

        do i = 1, nPtsX
            do j = 1, nPtsY
        
                do n = 1, nModesX 
                    do m = 1, nModesY 

                      !----------------------------------------------------
                      !cexpkxky = exp(zi * (xkx(n) * x(i) + xky(m) * y(j)))
                      !----------------------------------------------------

                      cexpkxky = xx(kxL+n-1, i) * yy(kyL+m-1, j)

                      f(i,j) = f(i,j) + a(n,m) * cexpkxky
                      if (present(fx)) fx(i,j) = f(i,j) + xkxsav(kxL+n-1) * a(n,m) * cexpkxky
                      if (present(fy)) fy(i,j) = f(i,j) + xkysav(kyL+m-1) * a(n,m) * cexpkxky

                    enddo
                enddo

            enddo
        enddo
    
    
    end subroutine sftinv2d
!
!*******************************************************************
!

      subroutine deriv_fourier_r(x, y, xkx, xky, &
         a, nxdim, nydim, &
         nkdim1, nkdim2, mkdim1, mkdim2, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, xx, yy, &
         f, fx, fxx, fy, fyy, fxy)

      implicit none

      integer nxdim, nydim, n, m, i, j, &
         nkdim1, nkdim2, mkdim1, mkdim2, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2

      real x(nxdim), y(nydim)
      real xkx(nkdim1 : nkdim2), xky(mkdim1 : mkdim2)
      complex zi, &
         a(nkdim1 : nkdim2, mkdim1 : mkdim2)
      real f(nxdim, nydim), fx(nxdim, nydim), fxx(nxdim, nydim), &
        fy(nxdim, nydim), fyy(nxdim, nydim), fxy(nxdim, nydim)

      complex cexpkxky

      complex xx(nkdim1 : nkdim2, nxdim), &
              yy(mkdim1 : mkdim2, nydim)



      zi = cmplx(0.0, 1.0)


      do i = 1, nnodex
         do j = 1, nnodey

            f(i,j)   = 0.0
            fx(i,j)  = 0.0
            fxx(i,j) = 0.0
            fy(i,j)  = 0.0
            fyy(i,j) = 0.0
            fxy(i,j) = 0.0

            do n = nkx1, nkx2
               do m = nky1, nky2

!                 ----------------------------------------------------
!                 cexpkxky = exp(zi * (xkx(n) * x(i) + xky(m) * y(j)))
!                 ----------------------------------------------------
                  cexpkxky = xx(n, i) * yy(m, j)

                  f(i,j) = f(i,j) + real(a(n,m) * cexpkxky)

                  fx(i,j) = fx(i,j) &
                              + real(zi * xkx(n) * a(n,m) * cexpkxky)
                  fxx(i,j) = fxx(i,j) &
                                - real(xkx(n)**2 * a(n,m) * cexpkxky)

                  fy(i,j) = fy(i,j) &
                              + real(zi * xky(m) * a(n,m) * cexpkxky)
                  fyy(i,j) = fyy(i,j) &
                                - real(xky(m)**2 * a(n,m) * cexpkxky)

                  fxy(i,j) = fxy(i,j) &
                            - real(xkx(n) * xky(m) * a(n,m) * cexpkxky)


               end do
            end do

         end do
      end do


      return
      end subroutine deriv_fourier_r
!
!*******************************************************************
!

      subroutine deriv_fourier(x, y, xkx, xky, &
         a, nxdim, nydim, &
         nkdim1, nkdim2, mkdim1, mkdim2, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, xx, yy, &
         f, fx, fxx, fy, fyy, fxy)

      implicit none

      integer nxdim, nydim, n, m, i, j, &
         nkdim1, nkdim2, mkdim1, mkdim2, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2

      real x(nxdim), y(nydim)
      real xkx(nkdim1 : nkdim2), xky(mkdim1 : mkdim2)
      complex zi, &
         a(nkdim1 : nkdim2, mkdim1 : mkdim2)
      complex f(nxdim, nydim), fx(nxdim, nydim), fxx(nxdim, nydim), &
        fy(nxdim, nydim), fyy(nxdim, nydim), fxy(nxdim, nydim)

      complex cexpkxky

      complex xx(nkdim1 : nkdim2, nxdim), &
              yy(mkdim1 : mkdim2, nydim)



      zi = cmplx(0.0, 1.0)


      do i = 1, nnodex
         do j = 1, nnodey

            f(i,j)   = 0.0
            fx(i,j)  = 0.0
            fxx(i,j) = 0.0
            fy(i,j)  = 0.0
            fyy(i,j) = 0.0
            fxy(i,j) = 0.0

            do n = nkx1, nkx2
               do m = nky1, nky2

!                 ----------------------------------------------------
!                 cexpkxky = exp(zi * (xkx(n) * x(i) + xky(m) * y(j)))
!                 ----------------------------------------------------
                  cexpkxky = xx(n, i) * yy(m, j)

                  f(i,j) = f(i,j) + a(n,m) * cexpkxky

                  fx(i,j) = fx(i,j) + zi * xkx(n) * a(n,m) * cexpkxky
                  fxx(i,j) = fxx(i,j) - xkx(n)**2 * a(n,m) * cexpkxky

                  fy(i,j) = fy(i,j) + zi * xky(m) * a(n,m) * cexpkxky
                  fyy(i,j) = fyy(i,j) - xky(m)**2 * a(n,m) * cexpkxky

                  fxy(i,j) = fxy(i,j) &
                                - xkx(n) * xky(m) * a(n,m) * cexpkxky


               end do
            end do

         end do
      end do


      return
      end subroutine deriv_fourier
!
!*******************************************************************
!

      subroutine sft2d(f, xlen, ylen, nnodex, nnodey, &
         nxdim, nydim, nkdim1, nkdim2, mkdim1, mkdim2, &
         nkx1, nkx2, nky1, nky2, xx_inv, yy_inv, fk, dx, dy, &
         myid, nproc, icontxt)

      implicit none

      integer nnodex, nnodey, nxdim, nydim, nkdim1, nkdim2, &
          mkdim1, mkdim2, i, j, n, m, nkx1, nkx2, nky1, nky2

      integer  myid, nproc, id, ngrid, icontxt

      real xlen, ylen, dx, dy, fact

      complex fk(nkdim1 : nkdim2, mkdim1 : mkdim2), xint, cexpkxky
      complex f(nxdim, nydim)
      complex xx_inv(nkdim1 : nkdim2, nxdim), &
              yy_inv(mkdim1 : mkdim2, nydim)

      fact = dx / xlen * dy / ylen


      do n = nkx1, nkx2
         do m = nky1, nky2

            fk(n, m) = 0.0

            do i = 1, nnodex
               do j = 1, nnodey

!                 ----------------------------------------------------------
!                 cexpkxky = 1.0 / exp(zi * (xkx(n) * x(i) + xky(m) * y(j)))
!                 ----------------------------------------------------------
                  cexpkxky = xx_inv(n, i) * yy_inv(m, j)

                  xint = f(i, j) * cexpkxky

                  fk(n, m) = fk(n, m) + xint * fact


               end do
            end do

         end do
      end do

      return
      end subroutine sft2d

!
!!*******************************************************************
!!
!
!      subroutine sft2d_parallel(f, xlen, ylen, nnodex, nnodey, &
!         nxdim, nydim, nkdim1, nkdim2, mkdim1, mkdim2, &
!         nkx1, nkx2, nky1, nky2, xx_inv, yy_inv, fk, dx, dy, &
!         myid, nproc, icontxt)
!
!      implicit none
!
!      integer nnodex, nnodey, nxdim, nydim, nkdim1, nkdim2, &
!          mkdim1, mkdim2, i, j, n, m, nkx1, nkx2, nky1, nky2
!
!      integer  myid, nproc, id, ngrid, icontxt
!
!      real xlen, ylen, dx, dy, fact
!
!      complex fk(nkdim1 : nkdim2, mkdim1 : mkdim2), xint, cexpkxky
!      complex f(nxdim, nydim), fksav(nxdim, nydim)
!      complex xx_inv(nkdim1 : nkdim2, nxdim), &
!              yy_inv(mkdim1 : mkdim2, nydim)
!
!      fksav(:,:) = (0.0,0.0)
!
!      fact = dx / xlen * dy / ylen
!
!!      print*,myid,"fourier: nnodex  nnodey",nnodex,nnodey
!
!      call blacs_barrier(icontxt, 'All')
!
!
!      do n = nkx1, nkx2
!         do m = nky1, nky2
!
!            ngrid = (n - nkx1) * nnodey + (m - nky1 + 1)
!            id = mod(ngrid, nproc)
!            if(id .eq. myid)then
!
!               fksav(n - nkx1 + 1, m - nky1 + 1) = 0.0
!
!               do i = 1, nnodex
!                  do j = 1, nnodey
!
!!                    ----------------------------------------------------------
!!                    cexpkxky = 1.0 / exp(zi * (xkx(n) * x(i) + xky(m) * y(j)))
!!                    ----------------------------------------------------------
!!                     cexpkxky = xx_inv(n, i) * yy_inv(m, j)
!!
!!                     xint = f(i, j) * cexpkxky
!!
!!                     fksav(n - nkx1 + 1, m - nky1 + 1)
!!     .                 = fksav(n - nkx1 + 1, m - nky1 + 1) + xint * fact
!                     fksav(n - nkx1 + 1, m - nky1 + 1) &
!                       = fksav(n - nkx1 + 1, m - nky1 + 1) + &
!                         f(i,j) * xx_inv(n,i) * yy_inv(m,j) *fact
!
!                  end do
!               end do
!
!            end if
!
!         end do
!      end do
!
!      call blacs_barrier(icontxt, 'All')
!
!
!      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, fksav, &
!         nxdim, -1, -1)
!
!
!      do n = nkx1, nkx2
!         do m = nky1, nky2
!            fk(n, m) = fksav(n - nkx1 + 1, m - nky1 + 1)
!         end do
!      end do
!
!
!
!
!      return
!      end subroutine sft2d_parallel
!
!!
!!*******************************************************************
!




      subroutine sft2d_old(x, y, f, xlen, ylen, nnodex, nnodey, &
         nxdim, nydim, nkdim1, nkdim2, mkdim1, mkdim2, &
         nkx1, nkx2, nky1, nky2, xx, yy, fk, dx, dy)

      implicit none

      integer nnodex, nnodey, nxdim, nydim, nkdim1, nkdim2, &
          mkdim1, mkdim2, &
         i, j, n, m, nkx1, nkx2, nky1, nky2
      real x(nxdim), y(nydim), xlen, ylen
      real f(nxdim, nydim), dx, dy
      complex fk(nkdim1 : nkdim2, mkdim1 : mkdim2), xint

      complex xx(nkdim1 : nkdim2, nxdim), yy(mkdim1 : mkdim2, nydim)


      do n = nkx1, nkx2
         do m = nky1, nky2

            fk(n, m) = 0.0

            do i = 1, nnodex - 1
               do j = 1, nnodey -1
                  xint = (f(i, j)     / (xx(n, i)   * yy(m, j)) &
                        + f(i+1, j)   / (xx(n, i+1) * yy(m, j)) &
                        + f(i, j+1)   / (xx(n, i)   * yy(m, j+1)) &
                        + f(i+1, j+1) / (xx(n, i+1) * yy(m, j+1)) &
                                )/ 4.0
                  fk(n, m) = fk(n, m) + xint * (x(i+1) - x(i)) / xlen &
                                             * (y(j+1) - y(j)) / ylen


               end do
            end do

         end do
      end do

      return
      end subroutine sft2d_old
!
!*******************************************************************
!


      subroutine sft(x, fx, xlen, nx, nxdim, nkdim1, nkdim2, &
         nk1, nk2, an)

      implicit none

      integer nx, nxdim, nkdim1, nkdim2, n, nk, nk1, nk2
      real x(nxdim), xlen, xkn, pi
      complex fx(nxdim), an(nkdim1:nkdim2), zi, xint

      pi = 3.14159
      zi = cmplx(0.0, 1.0)

      do nk = nk1, nk2
         xkn = 2.0 * pi/xlen * nk
         an(nk) = 0.0
         do n = 1, nx-1
            xint = (fx(n)   * exp(-zi * xkn * x(n)) &
                  + fx(n+1) * exp(-zi * xkn * x(n+1)) )/ 2.0
            an(nk) = an(nk) + xint * (x(n+1) - x(n)) / xlen
         end do
      end do

      return
      end subroutine sft
!
!*******************************************************************
!
      subroutine sftinv(x, fx, xlen, nx, nxdim, nkdim1, nkdim2, &
          nk1, nk2, an)

      implicit none

      integer nx, nxdim, nkdim1, nkdim2, n, nk, nk1, nk2
      real x(nxdim), xlen, xkn, pi
      complex zi, fx(nxdim), an(nkdim1 : nkdim2)

      pi = 3.14159
      zi = cmplx(0.0, 1.0)

      do n = 1, nx
         fx(n) = 0.0
         do nk = nk1, nk2
            xkn = 2.0 * pi * nk / xlen
            fx(n) = fx(n) + an(nk) * exp(zi * xkn * x(n))
         end do
      end do

      return
      end subroutine sftinv
!
!*******************************************************************
!
      subroutine fft(nx, work, nxdim, nkdim1, nkdim2, fx, an)

      implicit none

      integer nhalf, nx, nxdim, nkdim1, nkdim2, n, nk
      complex fx(nxdim), work(nxdim), an(nkdim1 : nkdim2)
      real x(nxdim)

      nhalf = nx / 2
      do n = 1, nx
         work(n) = fx(n)
      end do

      call four1(real(work), nx, -1)

      do n = 1, nx
         nk = n - 1
         if(nk .gt. nhalf) nk = nk - nx
         an(nk) = work(n) / nx
      end do
      return
      end subroutine fft
!
!***************************************************************************
!
      subroutine fftinv(nx, work, nxdim, nkdim1, nkdim2, fx, an)

      implicit none

      integer nhalf, nx, nxdim, nkdim1, nkdim2, n, nk
      complex fx(nxdim), work(nxdim), an(nkdim1 : nkdim2)
      real x(nxdim)

      nhalf = nx / 2

      do n = 1, nx
         nk = n - 1
         if(nk .gt. nhalf) nk = nk - nx
         work(n) = an(nk)
      end do

      call four1(real(work), nx, 1)

      do n = 1, nx
         fx(n) = work(n)
      end do

      return
      end subroutine fftinv
!
!***************************************************************************
!
      subroutine fftinvr(nx, work, nxdim, nkdim1, nkdim2, fx, an)

      implicit none

      integer nhalf, nx, nxdim, nkdim1, nkdim2, n, nk
      real fx(nxdim), an(nkdim1 : nkdim2)
      real x(nxdim)
      complex work(nxdim)

      nhalf = nx / 2

      do n = 1, nx
         nk = n - 1
         if(nk .gt. nhalf) nk = nk - nx
         work(n) = cmplx(an(nk), 0.0)
!         write(6, 100) n, work(n)
      end do

      call four1(real(work), nx, 1)

      do n = 1, nx
         fx(n) = real(work(n))
      end do
  100 format (i10, 1p8e12.4)
      return
      end subroutine fftinvr
!
!*******************************************************************
!
      subroutine four1(data, nn, isign)

      implicit none

      integer isign, nn
      real data(2*nn)
!         Replaces data(1:2*nn) by its discrete Fourier transform, if isign
!         is input as 1; or replaces
!         data(1:2*nn) by nn times its inverse discrete Fourier transform,
!         if isign is input as -1.
!         data is a complex array of length nn or, equivatlently a real
!         array of length 2*nn.  nn MUST be an integer power of 2
!         (this is not checked for!).
      integer i, istep, j, m, mmax, n
      real tempi, tempr
!          Double precision for the trigonometic recurrences
      double precision theta, wi, wpi, wpr, wr, wtemp

      n = 2 * nn
      j = 1
!          This is bit reversal section of the routine
      do i = 1, n, 2
         if (j .gt. i) then
            tempr = data(j)
            tempi = data(j+1)
            data(j) = data(i)
            data(j+1) = data(i+1)
            data(i) = tempr
            data(i+1) = tempi
         endif
         m = n/2

    1    if (m.ge.2 .and. j.gt.m) then
            j = j - m
            m = m/2
            go to 1
         endif
         j = j + m
      end do
      mmax = 2
    2 if (n .gt. mmax) then
         istep = 2*mmax
         theta = 6.28318530717959D0/(isign*mmax)
         wpr=-2.d0*dsin(0.5d0*theta)**2

         wpi=dsin(theta)

         wr=1.d0

         wi=0.d0

         do m = 1, mmax, 2
            do i = m, n, istep
               j = i + mmax
               tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
               tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)

               data(j) = data(i) - tempr
               data(j+1) = data(i+1) - tempi
               data(i) = data(i) + tempr
               data(i+1) = data(i+1) + tempi
            end do
            wtemp = wr
            wr = wr*wpr - wi*wpi + wr
            wi = wi*wpr + wtemp*wpi + wi
         end do
         mmax = istep
         go to 2
      endif
      return
      end subroutine four1

      end module fourier_mod


