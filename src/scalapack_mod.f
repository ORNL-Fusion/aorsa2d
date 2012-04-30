        module scalapack_mod

        use assert_mod
        use t4_mod_tprof

        implicit none

        public

!     any way to combine this interface and the next?
        interface pcopy2d
        module procedure pcopy2d_z, pcopy2d_d, pcopy2d_c, pcopy2d_s,     &
     &       pcopy2d_z2, pcopy2d_d2, pcopy2d_c2, pcopy2d_s2,                 &
     &       pcopy2d_gl_c,pcopy2d_gl_z,pcopy2d_gl_d,pcopy2d_gl_s,            &
     &       pcopy2d_gl_c2,pcopy2d_gl_z2,pcopy2d_gl_d2,pcopy2d_gl_s2
        end interface
        
        interface pcopy2d_lg
        module procedure pcopy2d_lg_c,pcopy2d_lg_z,pcopy2d_lg_d,         &
     &       pcopy2d_lg_s
        end interface

        interface pcopy
        module procedure pcopy_z,pcopy_d,pcopy_c,pcopy_s

        end interface

!       ---------------------
!       contains constants and
!       interface information
!       ---------------------
        interface pgemv
        module procedure pgemv_c111,pgemv_d111,pgemv_s111, pgemv_z111,      &
     &                   pgemv_c211,pgemv_d211,pgemv_s211, pgemv_z211,      &
     &                   pgemv_c221,pgemv_d221,pgemv_s221, pgemv_z221
        end interface

        interface pgeadd
        module procedure pgeadd_c, pgeadd_d, pgeadd_s, pgeadd_z,         &
     &    pgeadd_c21, pgeadd_d21, pgeadd_s21, pgeadd_z21,                &
     &    pgeadd_c12, pgeadd_d12, pgeadd_s12, pgeadd_z12,                &
     &    pgeadd_c22, pgeadd_d22, pgeadd_s22, pgeadd_z22
        end interface

        interface pgemr2d
        module procedure pgemr2d_s,pgemr2d_d,pgemr2d_c,pgemr2d_z,        &
     &       pgemr2d_s1,pgemr2d_d1,pgemr2d_c1,pgemr2d_z1
        end interface

        interface gsum2d
        module procedure gsum2d_i,gsum2d_s,gsum2d_d,gsum2d_c,gsum2d_z,    &
     &       gsum2d_i1,gsum2d_s1,gsum2d_d1,gsum2d_c1,gsum2d_z1,                &
     &       gsum2d_i0,gsum2d_s0,gsum2d_d0,gsum2d_c0,gsum2d_z0
        end interface

        interface pelget
        module procedure pelget_c, pelget_d, pelget_s, pelget_z,          &
     &       pelget_c2, pelget_d2, pelget_s2, pelget_z2
        end interface

        interface pelset
        module procedure pelset_c, pelset_d, pelset_s, pelset_z,         &
     &       pelset_c2, pelset_d2, pelset_s2, pelset_z2
        end interface

!     gebX2d_Y1 are for A of rank 1:
        interface gebs2d
        module procedure gebs2d_c, gebs2d_d, gebs2d_s, gebs2d_z,          &
     &       gebs2d_c1, gebs2d_d1, gebs2d_s1, gebs2d_z1,                  &
     &       gebs2d_c0, gebs2d_d0, gebs2d_s0, gebs2d_z0 
        end interface

        interface gebr2d
        module procedure gebr2d_c, gebr2d_d, gebr2d_s, gebr2d_z,          &
     &       gebr2d_c1, gebr2d_d1, gebr2d_s1, gebr2d_z1,                  &
     &       gebr2d_c0, gebr2d_d0, gebr2d_s0, gebr2d_z0
        end interface

        interface gesd2d
        module procedure gesd2d_c, gesd2d_d, gesd2d_s, gesd2d_z,          &
     &       gesd2d_c1, gesd2d_d1, gesd2d_s1, gesd2d_z1,                  &
     &       gesd2d_c0, gesd2d_d0, gesd2d_s0, gesd2d_z0
        end interface

        interface gerv2d
        module procedure gerv2d_c, gerv2d_d, gerv2d_s, gerv2d_z,          &
     &       gerv2d_c1, gerv2d_d1, gerv2d_s1, gerv2d_z1,                  &
     &       gerv2d_c0, gerv2d_d0, gerv2d_s0, gerv2d_z0
        end interface

        interface pgemm
        module procedure pgemm_c, pgemm_z, pgemm_s, pgemm_d,             &
     &       pgemm_c2, pgemm_z2, pgemm_s2, pgemm_d2
        end interface

        interface pger
        module procedure pger_c, pger_z, pger_s, pger_d,                  &
     &       pger_c2, pger_z2, pger_s2, pger_d2
        end interface

        interface plaprnt
        module procedure plaprnt_c,plaprnt_z,plaprnt_s,plaprnt_d,         &
     &       plaprnt_c2, plaprnt_z2, plaprnt_s2, plaprnt_d2
        end interface


        interface plange
        module procedure plange_c, plange_z, plange_s, plange_d,          &
     &             plange_c2, plange_z2, plange_s2, plange_d2
        end interface

        interface placp3
        module procedure placp3_c, placp3_z, placp3_s, placp3_d
        end interface

        interface ptran
        module procedure ptran_s, ptran_d
        end interface

        interface pdot
        module procedure pdot_s,pdot_d, pdot_s22, pdot_d22,               &
     &       pdotc_c,pdotc_z
        end interface

        interface pdotu
        module procedure pdotu_c, pdotu_z
        end interface



        integer, parameter :: block_cyclic_2d = 1
        integer, parameter :: dlen_ = 9 
        integer, parameter :: dtype_ = 1 
        integer, parameter :: ctxt_ = 2 
        integer, parameter :: m_ = 3 
        integer, parameter :: n_ = 4 
        integer, parameter :: mb_ = 5 
        integer, parameter :: nb_ = 6 
        integer, parameter :: rsrc_ = 7 
        integer, parameter :: csrc_ = 8 
        integer, parameter :: lld_ = 9 

!       ---------------------
!       Parallel BLAS (PBLAS) library
!       ---------------------

        interface 

        subroutine pddot(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        implicit none
        integer n,ix,jx,descx(*),incx,iy,jy,descy(*),incy
        real*8 dot,x(*),y(*)
        end subroutine pddot

        subroutine psdot(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        implicit none
        integer n,ix,jx,descx(*),incx,iy,jy,descy(*),incy
        real*4 dot,x(*),y(*)
        end subroutine psdot

        subroutine pcdotc(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        implicit none
        integer n,ix,jx,descx(*),incx,iy,jy,descy(*),incy
        complex*8 dot,x(*),y(*)
        end subroutine pcdotc

        subroutine pcdotu(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        integer n,ix,jx,descx(*),incx,iy,jy,descy(*),incy
        complex*8 dot,x(*),y(*)
        end subroutine pcdotu

        subroutine pzdotc(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        implicit none
        integer n,ix,jx,descx(*),incx,iy,jy,descy(*),incy
        complex*16 dot,x(*),y(*)
        end subroutine pzdotc

        subroutine pzdotu(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        implicit none
        integer n,ix,jx,descx(*),incx,iy,jy,descy(*),incy
        complex*16 dot,x(*),y(*)
        end subroutine pzdotu

        end interface


        interface
        subroutine pzgeru( m, n,                                             &
     &       alpha,   X,ix,jx,descX,incX,                                    &
     &                Y,iy,jy,descY,incY,                                    &
     &                A,ia,ja,descA )
        integer m,n
        integer ix,jx,descX(*),incX
        integer iy,jy,descY(*),incY
        integer ia,ja,descA(*)
        complex*16 alpha, X(*),Y(*),A(*)
        end subroutine pzgeru



        subroutine pcgeru( m, n,                                             &
     &       alpha,   X,ix,jx,descX,incX,                                    &
     &                Y,iy,jy,descY,incY,                                    &
     &                A,ia,ja,descA )
        integer m,n
        integer ix,jx,descX(*),incX
        integer iy,jy,descY(*),incY
        integer ia,ja,descA(*)
        complex*8 alpha, X(*),Y(*),A(*)
        end subroutine pcgeru


        subroutine pzgerc( m, n,                                             &
     &       alpha,   X,ix,jx,descX,incX,                                    &
     &                Y,iy,jy,descY,incY,                                    &
     &                A,ia,ja,descA )
        integer m,n
        integer ix,jx,descX(*),incX
        integer iy,jy,descY(*),incY
        integer ia,ja,descA(*)
        complex*16 alpha, X(*),Y(*),A(*)
        end subroutine pzgerc



        subroutine pcgerc( m, n,                                             &
     &       alpha,   X,ix,jx,descX,incX,                                    &
     &                Y,iy,jy,descY,incY,                                    &
     &                A,ia,ja,descA )
        integer m,n
        integer ix,jx,descX(*),incX
        integer iy,jy,descY(*),incY
        integer ia,ja,descA(*)
        complex*8 alpha, X(*),Y(*),A(*)
        end subroutine pcgerc

        subroutine pdger( m, n,                                              &
     &       alpha,   X,ix,jx,descX,incX,                                    &
     &                Y,iy,jy,descY,incY,                                    &
     &                A,ia,ja,descA )
        integer m,n
        integer ix,jx,descX(*),incX
        integer iy,jy,descY(*),incY
        integer ia,ja,descA(*)
        real*8 alpha, X(*),Y(*),A(*)
        end subroutine pdger

        subroutine psger( m, n,                                              &
     &       alpha,   X,ix,jx,descX,incX,                                    &
     &                Y,iy,jy,descY,incY,                                    &
     &                A,ia,ja,descA )
        integer m,n
        integer ix,jx,descX(*),incX
        integer iy,jy,descY(*),incY
        integer ia,ja,descA(*)
        real*4 alpha, X(*),Y(*),A(*)
        end subroutine psger

        end interface


        interface 

        subroutine pztrsm(side,uplo,trans,diag,                               &
     &       m,n,alpha,   A,ia,ja,descA, B,ib,jb,descB )
        implicit none
        character side, uplo, trans, diag
        integer m,n,  ia,ja, ib,jb
        integer descA(*), descB(*)
        complex*16 alpha
        complex*16 A(*), B(*)
        end subroutine pztrsm

        subroutine pzscal(n,alpha,X,ix,jx,descX,incx)
        implicit none
        integer n,ix,jx,descX(*),incx
        complex*16 alpha, X(*)
        end subroutine pzscal

        subroutine pdznrm2(n,norm2,X,ix,jx,descX,incx)
        implicit none
        integer n,ix,jx,descX(*),incx
        complex*16 X(*)
        real*8 norm2
        end subroutine pdznrm2

        subroutine pxerbla(icontext,srname,info)
        character*(*) srname
        integer icontext,info
        end subroutine pxerbla







        subroutine pzcopy( n,X,ix,jx,descX,incx,Y,iy,jy,descY,incy)
        implicit none
        integer n,ix,jx,incx,  iy,jy,incy
        integer descX(*),descY(*)
        complex*16 X(*),Y(*)
        end subroutine pzcopy

        end interface


!       -----------------
!       generic interface
!       -----------------
        interface


        subroutine pzgemv(trans, m,n, alpha, A,ia,ja,descA,                  &
     &                                       X,ix,jx,descX,incx,             &
     &                                beta,  Y,iy,jy,descY,incy)
        character trans
        integer m,n,  ia,ja,descA(*) 
        integer ix,jx,descX(*),incx
        integer iy,jy,descY(*),incy
        complex*16 alpha,beta
        complex*16 A(*),X(*),Y(*)
        end subroutine pzgemv



        subroutine pcgemv(trans, m,n, alpha, A,ia,ja,descA,                  &
     &                                       X,ix,jx,descX,incx,             &
     &                                beta,  Y,iy,jy,descY,incy)
        character trans
        integer m,n,  ia,ja,descA(*) 
        integer ix,jx,descX(*),incx
        integer iy,jy,descY(*),incy
        complex*8 alpha,beta
        complex*8 A(*),X(*),Y(*)
        end subroutine pcgemv




        subroutine pdgemv(trans, m,n, alpha, A,ia,ja,descA,                   &
     &                                       X,ix,jx,descX,incx,              &
     &                                beta,  Y,iy,jy,descY,incy)
        character trans
        integer m,n,  ia,ja,descA(*) 
        integer ix,jx,descX(*),incx
        integer iy,jy,descY(*),incy
        real*8 alpha,beta
        real*8 A(*),X(*),Y(*)
        end subroutine pdgemv




        subroutine psgemv(trans, m,n, alpha, A,ia,ja,descA,                   &
     &                                       X,ix,jx,descX,incx,              &
     &                                beta,  Y,iy,jy,descY,incy)
        character trans
        integer m,n,  ia,ja,descA(*) 
        integer ix,jx,descX(*),incx
        integer iy,jy,descY(*),incy
        real*4 alpha,beta
        real*4 A(*),X(*),Y(*)
        end subroutine psgemv

        end interface



        interface 

        subroutine pzelset(A,ia,ja,descA,alpha)
        integer ia,ja,descA(*)
        complex*16 alpha
        complex*16 A(*)
        end subroutine pzelset

        subroutine pcelset(A,ia,ja,descA,alpha)
        integer ia,ja,descA(*)
        complex*8 alpha
        complex*8 A(*)
        end subroutine pcelset

        subroutine pdelset(A,ia,ja,descA,alpha)
        integer ia,ja,descA(*)
        real*8 alpha
        real*8 A(*)
        end subroutine pdelset

        subroutine pselset(A,ia,ja,descA,alpha)
        integer ia,ja,descA(*)
        real*4 alpha
        real*4 A(*)
        end subroutine pselset

        end interface


        interface

        subroutine pzlaprnt(m,n,A,ia,ja,descA,irprnt,icprnt,                  &
     &       cmatnm,nout,work)
        implicit none
        integer m,n,ia,ja,descA(*),irprnt,icprnt,nout
        complex*16 A(*),work(*)
        character*(*) cmatnm
        end subroutine pzlaprnt

        subroutine pclaprnt(m,n,A,ia,ja,descA,irprnt,icprnt,                  &
     &       cmatnm,nout,work)
        implicit none
        integer m,n,ia,ja,descA(*),irprnt,icprnt,nout
        complex*8 A(*),work(*)
        character*(*) cmatnm
        end subroutine pclaprnt

        subroutine pdlaprnt(m,n,A,ia,ja,descA,irprnt,icprnt,                  &
     &       cmatnm,nout,work)
        implicit none
        integer m,n,ia,ja,descA(*),irprnt,icprnt,nout
        real*8 A(*),work(*)
        character*(*) cmatnm
        end subroutine pdlaprnt

        subroutine pslaprnt(m,n,A,ia,ja,descA,irprnt,icprnt,                  &
     &       cmatnm,nout,work)
        implicit none
        integer m,n,ia,ja,descA(*),irprnt,icprnt,nout
        real*4 A(*),work(*)
        character*(*) cmatnm
        end subroutine pslaprnt

        end interface

        interface

        real*4 function pclange( norm,m,n,A,ia,ja,descA,work)
        implicit none
        character norm
        integer m,n,ia,ja,descA(*)
        complex*8 A(*)
        real*4 work(*)
        end function pclange



        real*8 function pzlange( norm,m,n,A,ia,ja,descA,work)
        implicit none
        character norm
        integer m,n,ia,ja,descA(*)
        complex*16 A(*)
        real*8 work(*)
        end function pzlange


        real*4 function pslange( norm,m,n,A,ia,ja,descA,work)
        implicit none
        character norm
        integer m,n,ia,ja,descA(*)
        real*4 A(*)
        real*4 work(*)
        end function pslange


        real*8 function pdlange( norm,m,n,A,ia,ja,descA,work)
        implicit none
        character norm
        integer m,n,ia,ja,descA(*)
        real*8 A(*)
        real*8 work(*)
        end function pdlange

        end interface


        interface


        subroutine psgemm( transA, transB,                                   &
     &       m,n,k,                                                          & 
     &       alpha,   A,ia,ja,descA,                                         & 
     &                B,ib,jb,descB,                                         & 
     &       beta,    C,ic,jc,descC )

        implicit none
        character transA,transB
        integer m,n,k
        integer ia,ja,  ib,jb,  ic,jc
        integer descA(*),descB(*),descC(*)
        real*4 alpha,beta
        real*4  A(*),B(*),C(*)
        end subroutine psgemm

        subroutine pdgemm( transA, transB,                                   &
     &       m,n,k,                                                          &
     &       alpha,   A,ia,ja,descA,                                         &
     &                B,ib,jb,descB,                                         &
     &       beta,    C,ic,jc,descC )

        implicit none
        character transA,transB
        integer m,n,k
        integer ia,ja,  ib,jb,  ic,jc
        integer descA(*),descB(*),descC(*)
        real*8 alpha,beta
        real*8  A(*),B(*),C(*)
        end subroutine pdgemm

        subroutine pcgemm( transA, transB,                                    &
     &       m,n,k,                                                           &
     &       alpha,   A,ia,ja,descA,                                          &
     &                B,ib,jb,descB,                                          & 
     &       beta,    C,ic,jc,descC )

        implicit none
        character transA,transB
        integer m,n,k
        integer ia,ja,  ib,jb,  ic,jc
        integer descA(*),descB(*),descC(*)
        complex*8 alpha,beta
        complex*8  A(*),B(*),C(*)
        end subroutine pcgemm



      subroutine pzgemm( transA, transB,                                   &
     &       m,n,k,                                                          &
     &       alpha,   A,ia,ja,descA,                                         & 
     &                B,ib,jb,descB,                                         &
     &       beta,    C,ic,jc,descC )

        implicit none
        character transA,transB
        integer m,n,k
        integer ia,ja,  ib,jb,  ic,jc
        integer descA(*),descB(*),descC(*)
        complex*16 alpha,beta
        complex*16  A(*),B(*),C(*)
        end subroutine pzgemm





        end interface


        interface

        subroutine pzelget( scope, top, alpha, A, ia, ja, descA )
        character scope
        character top
        integer ia,ja,descA(*)
        complex*16 alpha
        complex*16 A(*)
        end subroutine pzelget

        subroutine pcelget( scope, top, alpha, A, ia, ja, descA )
        character scope
        character top
        integer ia,ja,descA(*)
        complex*8 alpha
        complex*8 A(*)
        end subroutine pcelget

        subroutine pdelget( scope, top, alpha, A, ia, ja, descA )
        character scope
        character top
        integer ia,ja,descA(*)
        real*8 alpha
        real*8 A(*)
        end subroutine pdelget

        subroutine pselget( scope, top, alpha, A, ia, ja, descA )
        character scope
        character top
        integer ia,ja,descA(*)
        real*4 alpha
        real*4 A(*)
        end subroutine pselget





      subroutine pzgeadd(trans,m,n,alpha, A,ia,ja,descA,                   & 
     &    beta, C,ic,jc,descC)
        character trans
        integer m,n,  ia,ja,descA(*),   ic,jc,descC(*)
        complex*16 alpha,beta
        complex*16 A(*), C(*)
        end subroutine pzgeadd

      subroutine pcgeadd(trans,m,n,alpha, A,ia,ja,descA,                    &
     &    beta, C,ic,jc,descC)
        character trans
        integer m,n,  ia,ja,descA(*),   ic,jc,descC(*)
        complex*8 alpha,beta
        complex*8 A(*), C(*)
        end subroutine pcgeadd

      subroutine pdgeadd(trans,m,n,alpha, A,ia,ja,descA,                    &
     &    beta, C,ic,jc,descC)
        character trans
        integer m,n,  ia,ja,descA(*),   ic,jc,descC(*)
        real*8 alpha,beta
        real*8 A(*), C(*)
        end subroutine pdgeadd

      subroutine psgeadd(trans,m,n,alpha, A,ia,ja,descA,                    &
     &    beta, C,ic,jc,descC)
        character trans
        integer m,n,  ia,ja,descA(*),   ic,jc,descC(*)
        real*4 alpha,beta
        real*4 A(*), C(*)
        end subroutine psgeadd

        end interface 



!       ---------------
!       Scalapack TOOLS
!       ---------------
        interface

        integer function numroc(n,nb,iproc,isrcproc,nprocs)
        implicit none
        integer n,nb,iproc,isrcproc,nprocs
        end function numroc

        subroutine descinit(desc,m,n,mb,nb,                                   &
     &             irsrc,icsrc,ictxt,lld,info)
        implicit none
        integer m,n,mb,nb,irsrc,icsrc,ictxt,lld,info
        integer desc(*)
        end subroutine descinit


!        EFD: iceil already defined on cray x1
!
!        integer function iceil(inum,idenom)
!        implicit none
!        integer inum,idenom
!        end function iceil

        subroutine infog2l(grindx,gcindx,desc,nprow,npcol,                    &
     &           myrow,mycol,lrindx,lcindx,   rsrc,csrc )
        implicit none
        integer grindx,gcindx,nprow,npcol
        integer desc(*)
        integer myrow,mycol,lrindx,lcindx,   rsrc,csrc
        end subroutine infog2l

        subroutine infog1l( gindx,nb,nprocs,myroc,                            &
     &           isrcproc,lindx,rocsrc) 
        implicit none
        integer gindx,nb,nprocs,myroc
        integer isrcproc,lindx,rocsrc
        end subroutine infog1l

        integer function indxg2l(indxglob,nb,iproc,isrcproc,nprocs)
        implicit none
        integer indxglob,nb,iproc,isrcproc,nprocs
        end function indxg2l

        integer function indxg2p(indxglob,nb,iproc,isrcproc,nprocs)
        implicit none
        integer indxglob,nb,iproc,isrcproc,nprocs
        end function indxg2p

        integer function indxl2g(indxloc,nb,iproc,isrcproc,nprocs)
        implicit none
        integer indxloc,nb,iproc,isrcproc,nprocs
        end function indxl2g

        end interface


!       complex*16 interface

        interface

        subroutine pzmatadd(m,n,alpha,  A,ia,ja,descA,                        &
     &     beta, C,ic,jc,descC )
        implicit none
        integer m,n,   ia,ja,  ic,jc
        integer descA(*),descC(*)
        complex*16 alpha,beta
        complex*16 A(*),C(*)
        end subroutine pzmatadd




        subroutine pzelset2( alpha, A, ia, ja, descA, beta )
        integer ia,ja,descA(*)
        complex*16 alpha,beta
        complex*16 A(*)
        end subroutine pzelset2





        subroutine pzamax( n, amax, indx, X, ix, jx, descX, incX ) 
        implicit none
        integer n, indx, ix, jx, descX(*), incX
        complex*16 amax
        complex*16 X(*)
        end subroutine pzamax


        subroutine pzdscal( n, dalpha, X, ix, jx, descX, incX )
        implicit none
        integer n,  ix, jx, descX(*), incX
        real*8 dalpha
        complex*16 X(*)
        end subroutine pzdscal



        end interface


!       ----------------
!       Scalapack solver
!       ----------------
        interface

        subroutine pzgeqrf(m,n, A,ia,ja,descA,tau,work,lwork,info)
        implicit none
        integer m,n,   ia,ja, lwork, info
        integer descA(*)
        complex*16 A(*), tau(*), work(*)
        end subroutine pzgeqrf

        subroutine pzunmqr(side,trans,m,n,k,  A,ia,ja,descA, tau,             &
     &           B,ib,jb,descB,   work,lwork,info)
        implicit none
        character side,trans
        integer m,n,k,   ia,ja, ib,jb, lwork,info
        integer descA(*), descB(*)
        complex*16 A(*), tau(*), B(*), work(*)
        end subroutine pzunmqr

        subroutine pzlaset(uplo,m,n,alpha,beta,A,ia,ja,descA)
        implicit none
        character uplo
        integer m,n,ia,ja
        complex*16 alpha,beta
        integer descA(*)
        complex*16 A(*)
        end subroutine pzlaset

        

        subroutine pzgetrf(m,n,A,ia,ja,descA,ipiv,info)
        implicit none
        integer m,n,ia,ja,info
        integer descA(*),ipiv(*)
        complex*16 A(*)
        end subroutine pzgetrf

        subroutine pzgetrs(trans,n,nrhs,A,ia,ja,descA,                        &
     &           ipiv,B,ib,jb,descB,info)
        implicit none
        character trans
        integer n,nrhs,ia,ja,ib,jb,info
        integer descA(*),descB(*),ipiv(*)
        complex*16 A(*),B(*)
        end subroutine pzgetrs

        end interface


!       ---------------------------
!       BLACS communication library
!       ---------------------------
        interface

        subroutine blacs_abort(icontxt,ierrnum)
        implicit none
        integer icontxt,ierrnum
        end subroutine blacs_abort

        subroutine blacs_pinfo(iam,nprocs)
        implicit none
        integer iam,nprocs
        end subroutine blacs_pinfo

        subroutine blacs_gridinit( icontext, order, nprow,npcol)
        implicit none
        character order
        integer icontext,nprow,npcol
        end subroutine blacs_gridinit

        subroutine blacs_gridinfo(icontext,nprow,npcol,myrow,mycol)
        implicit none
        integer icontext,nprow,npcol,myrow,mycol
        end subroutine blacs_gridinfo


        subroutine blacs_gridexit( icontext )
        implicit none
        integer icontext
        end subroutine blacs_gridexit

        subroutine blacs_get( icontext,what,ival)
        implicit none
        integer icontext,what,ival
        end subroutine blacs_get

        end interface

        interface

        subroutine zgebs2d( icontext,scope,top,m,n,A,lda)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        complex*16 A(lda,*)
        end subroutine zgebs2d

        subroutine cgebs2d( icontext,scope,top,m,n,A,lda)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        complex*8 A(lda,*)
        end subroutine cgebs2d

        subroutine dgebs2d( icontext,scope,top,m,n,A,lda)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        real*8 A(lda,*)
        end subroutine dgebs2d

        subroutine sgebs2d( icontext,scope,top,m,n,A,lda)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        real*4 A(lda,*)
        end subroutine sgebs2d

        end interface


        interface

        subroutine sgebr2d( icontext,scope,top,m,n,A,lda,rsrc,csrc)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        real*4 A(lda,*)
        integer rsrc,csrc
        end subroutine sgebr2d

        subroutine dgebr2d( icontext,scope,top,m,n,A,lda,rsrc,csrc)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        real*8 A(lda,*)
        integer rsrc,csrc
        end subroutine dgebr2d

        subroutine cgebr2d( icontext,scope,top,m,n,A,lda,rsrc,csrc)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        complex*8 A(lda,*)
        integer rsrc,csrc
        end subroutine cgebr2d

        subroutine zgebr2d( icontext,scope,top,m,n,A,lda,rsrc,csrc)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        complex*16 A(lda,*)
        integer rsrc,csrc
        end subroutine zgebr2d

        end interface


        contains


        subroutine pzgecopy(m,n, A,ia,ja,descA, B,ib,jb,descB)
        implicit none
        integer, intent(in) :: m,n,ia,ja, ib,jb
        integer, dimension(dlen_),intent(in) :: descA,descB
        complex*16, dimension(*),intent(in) :: A
        complex*16, dimension(*),intent(inout) :: B

        integer :: j,jja,jjb,incx,incy
        logical isok
        complex*16 :: alpha, beta

        logical, parameter :: use_geadd = .true.


        call profstart('pzgecopy')

        call assert( m.ge.1,'pzgecopy:invalid m ',m)
        call assert( n.ge.1,'pzgecopy:invalid n ',n)

        call assert( descA(ctxt_).eq.descB(ctxt_),                            &
     &      'pzgecopy: should have same context ',descB(ctxt_))

        isok = descA(m_).ge.(ia+m-1)
        call assert( isok,'pzgecopy: invalid descA(m_) ',descA(m_))

        isok = descA(n_).ge.(ja+n-1)
        if (.not.isok) then
          write(*,*) 'm,n ', m,n
          write(*,*) 'ia,ja,descA(M_),descA(N_) ',                            &
     &        ia,ja,descA(M_),descA(N_)
          write(*,*) 'ib,jb,descB(M_),descB(N_) ',                            &
     &        ib,jb,descB(M_),descB(N_)
        endif
        call assert( isok,'pzgecopy: invalid descA(n_) ',descA(n_))

        isok = descB(m_).ge.(ib+m-1)
        call assert( isok,'pzgecopy: invalid descB(m_) ',descB(m_))

        isok = descB(n_).ge.(jb+n-1)
        call assert( isok,'pzgecopy: invalid descB(n_) ',descB(n_))

        if (use_geadd) then
!       -----------------------------
!       pzgeadd is in PBLAS Version 2 
!       -----------------------------
        alpha = 1.0
        beta = 0.0
        call pzgeadd( 'N', m,n, alpha, A,ia,ja,descA,                        &
     &       beta, B,ib,jb,descB)

        else

        incx = 1
        incy = 1
        do j=1,n
          jjb = jb + (j-1)
          jja = ja + (j-1)
          call pzcopy( m, A,ia,jja,descA,incx, B,ib,jjb,descB,incy)
        enddo

        endif

        call profend('pzgecopy')
        return
        end subroutine pzgecopy

        subroutine pzgecopy_2to3(k,nrhs, xk,descxk, xvec,descxvec)
        implicit none
        integer, intent(in) :: k,nrhs
        integer, dimension(dlen_),intent(in) :: descxk,descxvec
        complex*16, dimension(:), intent(in) :: xk
        complex*16, dimension(:,:),intent(inout) :: xvec

        integer :: m, ix,jx,  iix,jjx,  irhs,incx,incy
        logical :: isok

        call profstart('pzgecopy_2to3')

!       --------------
!       error checking
!       --------------
        m = descxk(m_)

        isok  = (m.le.descxvec(m_))
        if (.not.isok) then
           write(*,*) 'pzgecopy_2to3:k,nrhs ',k,nrhs 
           write(*,*) 'descxk(M_),descxk(N_) ',                              &
     &            descxk(M_),descxk(N_)                                       
           write(*,*) 'descxvec(M_),descxvec(N_) ',                          &
     &            descxvec(M_),descxvec(N_)
        endif
           
        call assert(isok,                                                    &
     &    'pzgecopy_2to3: invalid descxk(m_)',descxk(m_))

        isok = (descxk(n_).ge.nrhs)
        call assert(isok,                                                    &
     &    'pzgecopy_2to3: invalid descxk(n_)',descxk(n_))

        isok = (size(xvec,2).ge.nrhs)
        call assert(isok,                                                    &
     &     'pzgecopy_2to3: invalid size(xvec,2) ',size(xvec,2))


        

        incx = 1
        incy = 1
        do irhs=1,nrhs
          ix = 1
          jx = irhs
          iix = 1
          jjx = k
          call pzcopy( m, xk,ix,jx,descxk,incx,                               &
     &         xvec(:,irhs),iix,jjx,descxvec,incy )
        enddo

        call profend('pzgecopy_2to3')
        return

        end subroutine pzgecopy_2to3

        subroutine pzgecopy_3to2(k,nrhs, yvec,descyvec,yk,descyk)
        integer, intent(in) :: k,nrhs
        integer, dimension(dlen_),intent(in) :: descyvec,descyk
        complex*16,dimension(:,:),intent(in) :: yvec
        complex*16,dimension(:),intent(inout) :: yk

        integer :: m,  iy,jy,  iiy,jjy,  irhs,  incx,incy
        logical :: isok


        call profstart('pzgecopy_3to2')

        m = descyk(m_)


!       --------------
!       error checking
!       --------------
        isok = (descyk(m_) .le. descyvec(m_)) 
        call assert(isok,                                                    &
     &      'pzgecopy_3to2:invalid descyk(m_)',descyk(m_))

        isok = (descyk(n_).ge.nrhs)
        call assert(isok,                                                    &
     &      'pzgecopy_3to2: invalid  descyk(n_)',descyk(n_))

        isok = (size(yvec,2).ge.nrhs)
        call assert(isok,                                                    &
     &      'pzgecopy_3to2: invalid size(yvec,2)',size(yvec,2))
      

        incx = 1
        incy = 1

        do irhs=1,nrhs
          iiy = 1
          jjy = k
          iy = 1
          jy = irhs
          call pzcopy( m, yvec(:,irhs),iiy,jjy,descyvec,incx,                 &
     &                    yk,iy,jy,descyk,incy )
        enddo



        call profend('pzgecopy_3to2')
        return
        end subroutine pzgecopy_3to2

!       -----------------------------------------------
!       fill a section of array by constant value alpha
!       -----------------------------------------------
        subroutine pzgefill( m,n, A, ia,ja, descA, alpha )
        implicit none
        integer, intent(in) :: m,n, ia,ja
        integer, dimension(DLEN_) :: descA
        complex*16, dimension(*), intent(inout) :: A
        complex*16, intent(in) :: alpha

        complex*16 :: beta


        call profstart('pzgefill')
        beta = alpha
        call pzlaset( 'F', m,n, alpha, beta, A,ia,ja,descA)
        call profend('pzgefill')

        return
        end subroutine pzgefill

        subroutine pzgefill_org( m,n, A, ia,ja, descA, alpha )
        implicit none
        integer, intent(in) :: m,n, ia,ja
        integer, dimension(DLEN_) :: descA
        complex*16, dimension(*), intent(inout) :: A
        complex*16, intent(in) :: alpha

        complex*16, parameter :: zero = 0.0d0

        integer :: rsrc,csrc,iia,jja
        integer :: lrindx,lcindx,ipos,iproc
        integer :: nprow,npcol,myrow,mycol

        call profstart('pzgefill')
!       ------------
!       special case
!       ------------
        if (alpha.eq.zero) then
           call pzgeadd( 'N', m,n, zero, A,ia,ja,descA,                      &
     &       zero, A,ia,ja,descA)
        else

           call blacs_gridinfo( descA(CTXT_), nprow,npcol,myrow,mycol)

           do jja=ja,ja+n-1
              iproc = indxg2p( jja, descA(NB_),mycol,descA(CSRC_),npcol)
              if (iproc.ne.mycol) cycle

              do iia=ia,ia+m-1
                iproc=indxg2p(iia,descA(MB_),myrow,descA(RSRC_),nprow)
                if (iproc.ne.myrow) cycle


                call infog2l(iia,jja,descA,nprow,npcol,                      &
     &           myrow,mycol,lrindx,lcindx,   rsrc,csrc )
                call assert(rsrc.eq.myrow,                                   &
     &               'pzgefill: invalid rsrc ',rsrc)
                call assert(csrc.eq.mycol,                                   &
     &               'pzgefill: invalid csrc ',csrc)

                ipos = lrindx + (lcindx-1)*descA(LLD_)
                A(ipos) = alpha
              enddo
            enddo
         endif


        call profend('pzgefill')
        return
        end subroutine pzgefill_org

        subroutine pzgeqrs( m,n, nrhs, A,ia,ja,descA, tau,                   &
     &             B,ib,jb,descB, info)
        implicit none
        integer, intent(in) :: m,n,nrhs
        integer, intent(in) :: ia,ja,descA(*)
        integer, intent(in) :: ib,jb,descB(*)
        complex*16, intent(in) :: tau(*)
        integer, intent(inout) :: info
        complex*16, intent(in) :: A(*)
        complex*16, intent(inout) :: B(*)

        integer :: lzwork
        complex*16, allocatable, dimension(:) :: zwork
        complex*16 :: zwork1(1)
        complex*16 :: alpha
        integer :: mm,nn,kk

        call profstart('pzgeqrs')
!       ----------
!       apply Q'*B
!       ----------

        call profstart('pzgeqrs:pzunmqr')
        lzwork = -1
        mm = n
        nn = nrhs
        kk = m

!       --------------------------
!       query amount of work space
!       --------------------------
        call pzunmqr('Left','Conjugate', mm,nn,kk,                          &
     &          A,ia,ja,descA,tau,                                          &
     &          B,ib,jb,descB,   zwork1, lzwork, info )

        lzwork = abs(zwork1(1)) + 1
        allocate( zwork( lzwork ) )

!       ------------
!       perform Q'*B
!       ------------
        call pzunmqr('Left','Conjugate', mm,nn,kk,                           &
     &          A,ia,ja,descA,tau,                                           &
     &          B,ib,jb,descB,   zwork, lzwork, info )

        deallocate( zwork )
        call profend('pzgeqrs:pzunmqr')


!       ----------------
!       triangular solve
!       R\(Q'B)
!       ----------------
        call profstart('pzgeqrs:pztrsm')
        alpha = 1.0d0
        mm = n
        nn = nrhs
        call pztrsm( 'Left', 'Upper', 'No-transpose', 'Non-unit',             &
     &                 mm, nn, alpha, A,ia,ja,descA,                          &
     &                 B,ib,jb,descB )

        call profend('pzgeqrs:pztrsm')




        call profend('pzgeqrs')
        return

        end subroutine pzgeqrs

        subroutine pgeadd_z(trans,m,n,alpha, A,ia,ja,descA,              &
     &    beta, C,ic,jc,descC)
        character trans
        integer m,n,  ia,ja,descA(*),   ic,jc,descC(*)
        complex*16 alpha,beta
        complex*16 A(*), C(*)

        logical :: isok, no_trans

        isok = (1.le.ic).and.((ic+m-1).le.descC(M_)) .and.                  &
     &         (1.le.jc).and.((jc+n-1).le.descC(N_)) 
        no_trans = (trans.eq.'N').or.(trans.eq.'n')
        if (no_trans) then
           isok = isok.and.                                               &
     &            (1.le.ia).and.((ia+m-1).le.descA(M_)) .and.             &
     &            (1.le.ja).and.((ja+n-1).le.descA(N_))
        else
           isok = isok.and.                                               &
     &            (1.le.ia).and.((ia+n-1).le.descA(M_)) .and.             &
     &            (1.le.ja).and.((ja+m-1).le.descA(N_))
        endif
        if (.not.isok) then
          write(*,*) 'pzgeadd_z: trans,m,n ',trans,m,n
          write(*,*) 'alpha, ia,ja ',alpha,ia,ja
          write(*,*) 'beta, ic,jc ', beta,ic,jc
          write(*,*) 'descA(M_) ',descA(M_),' descC(M_) ',descC(M_)
          write(*,*) 'descA(N_) ',descA(N_),' descC(N_) ',descC(N_)
          write(*,*) 'descA(MB_) ',descA(MB_),' descC(MB_) ',descC(MB_)
          write(*,*) 'descA(NB_) ',descA(NB_),' descC(NB_) ',descC(NB_)
          write(*,*)'descA(LLD_)',descA(LLD_),' descC(LLD_)',descC(LLD_)
!          stop '** error ** '
          call blacs_gridexit( descA(CTXT_) )
        endif

        call pzgeadd(trans,m,n,alpha,A,ia,ja,descA,beta,C,ic,jc,descC)
        return

        end subroutine pgeadd_z

        subroutine pgeadd_c(trans,m,n,alpha, A,ia,ja,descA,              &
     &    beta, C,ic,jc,descC)
        character trans
        integer m,n,  ia,ja,descA(*),   ic,jc,descC(*)
        complex*8 alpha,beta
        complex*8 A(*), C(*)

        call pcgeadd(trans,m,n,alpha,A,ia,ja,descA,beta,C,ic,jc,descC)
        return

        end subroutine pgeadd_c

        subroutine pgeadd_d(trans,m,n,alpha, A,ia,ja,descA,              &
     &    beta, C,ic,jc,descC)
        character trans
        integer m,n,  ia,ja,descA(*),   ic,jc,descC(*)
        real*8 alpha,beta
        real*8 A(*), C(*)

        call pdgeadd(trans,m,n,alpha,A,ia,ja,descA,beta,C,ic,jc,descC)
        return

        end subroutine pgeadd_d

        subroutine pgeadd_s(trans,m,n,alpha, A,ia,ja,descA,              &
     &    beta, C,ic,jc,descC)
        character trans
        integer m,n,  ia,ja,descA(*),   ic,jc,descC(*)
        real*4 alpha,beta
        real*4 A(*), C(*)

        call psgeadd(trans,m,n,alpha,A,ia,ja,descA,beta,C,ic,jc,descC)
        return

        end subroutine pgeadd_s


        subroutine pgeadd_z22(trans,m,n,alpha, A,ia,ja,descA,              &
     &    beta, C,ic,jc,descC)
        character trans
        integer m,n,  ia,ja,descA(*),   ic,jc,descC(*)
        complex*16 alpha,beta
        complex*16 A(:,:), C(:,:)
        logical :: isok

        isok = (1.le.jc).and.(jc+n-1.le.descC(N_))
        if (.not.isok) then
          write(*,*) 'pgeadd_z22:m,n,ic,jc, descC(N_) ',                   &
     &         m,n,ic,jc,descC(N_)
        endif
        call pzgeadd(trans,m,n,alpha,A,ia,ja,descA,beta,C,ic,jc,descC)
        return

        end subroutine pgeadd_z22

        subroutine pgeadd_c22(trans,m,n,alpha, A,ia,ja,descA,              &
     &    beta, C,ic,jc,descC)
        character trans
        integer m,n,  ia,ja,descA(*),   ic,jc,descC(*)
        complex*8 alpha,beta
        complex*8 A(:,:), C(:,:)

        call pcgeadd(trans,m,n,alpha,A,ia,ja,descA,beta,C,ic,jc,descC)
        return

        end subroutine pgeadd_c22

        subroutine pgeadd_d22(trans,m,n,alpha, A,ia,ja,descA,              &
     &    beta, C,ic,jc,descC)
        character trans
        integer m,n,  ia,ja,descA(*),   ic,jc,descC(*)
        real*8 alpha,beta
        real*8 A(:,:), C(:,:)

        call pdgeadd(trans,m,n,alpha,A,ia,ja,descA,beta,C,ic,jc,descC)
        return

        end subroutine pgeadd_d22

        subroutine pgeadd_s22(trans,m,n,alpha, A,ia,ja,descA,              &
     &    beta, C,ic,jc,descC)
        character trans
        integer m,n,  ia,ja,descA(*),   ic,jc,descC(*)
        real*4 alpha,beta
        real*4 A(:,:), C(:,:)

        call psgeadd(trans,m,n,alpha,A,ia,ja,descA,beta,C,ic,jc,descC)
        return

        end subroutine pgeadd_s22

        subroutine pgeadd_z12(trans,m,n,alpha, A,ia,ja,descA,              &
     &    beta, C,ic,jc,descC)
        character trans
        integer m,n,  ia,ja,descA(*),   ic,jc,descC(*)
        complex*16 alpha,beta
        complex*16 A(*), C(:,:)
        logical :: isok

        isok = (1.le.jc).and.(jc+n-1.le.descC(N_))
        if (.not.isok) then
          write(*,*) 'pgeadd_z12:m,n,ic,jc, descC(N_) ',                   &
     &         m,n,ic,jc,descC(N_)
        endif
        call pzgeadd(trans,m,n,alpha,A,ia,ja,descA,beta,C,ic,jc,descC)
        return

        end subroutine pgeadd_z12

        subroutine pgeadd_c12(trans,m,n,alpha, A,ia,ja,descA,              &
     &    beta, C,ic,jc,descC)
        character trans
        integer m,n,  ia,ja,descA(*),   ic,jc,descC(*)
        complex*8 alpha,beta
        complex*8 A(*), C(:,:)

        call pcgeadd(trans,m,n,alpha,A,ia,ja,descA,beta,C,ic,jc,descC)
        return

        end subroutine pgeadd_c12

        subroutine pgeadd_d12(trans,m,n,alpha, A,ia,ja,descA,              &
     &    beta, C,ic,jc,descC)
        character trans
        integer m,n,  ia,ja,descA(*),   ic,jc,descC(*)
        real*8 alpha,beta
        real*8 A(*), C(:,:)

        call pdgeadd(trans,m,n,alpha,A,ia,ja,descA,beta,C,ic,jc,descC)
        return

        end subroutine pgeadd_d12

        subroutine pgeadd_s12(trans,m,n,alpha, A,ia,ja,descA,              &
     &    beta, C,ic,jc,descC)
        character trans
        integer m,n,  ia,ja,descA(*),   ic,jc,descC(*)
        real*4 alpha,beta
        real*4 A(*), C(:,:)

        call psgeadd(trans,m,n,alpha,A,ia,ja,descA,beta,C,ic,jc,descC)
        return

        end subroutine pgeadd_s12

        subroutine pgeadd_z21(trans,m,n,alpha, A,ia,ja,descA,              &
     &    beta, C,ic,jc,descC)
        character trans
        integer m,n,  ia,ja,descA(*),   ic,jc,descC(*)
        complex*16 alpha,beta
        complex*16 A(:,:), C(*)
        logical :: isok

        isok = (1.le.jc).and.(jc+n-1.le.descC(N_))
        if (.not.isok) then
          write(*,*) 'pgeadd_z21:m,n,ic,jc, descC(N_) ',                   &
     &         m,n,ic,jc,descC(N_)
        endif
        call pzgeadd(trans,m,n,alpha,A,ia,ja,descA,beta,C,ic,jc,descC)
        return

        end subroutine pgeadd_z21

        subroutine pgeadd_c21(trans,m,n,alpha, A,ia,ja,descA,              &
     &    beta, C,ic,jc,descC)
        character trans
        integer m,n,  ia,ja,descA(*),   ic,jc,descC(*)
        complex*8 alpha,beta
        complex*8 A(:,:), C(*)

        call pcgeadd(trans,m,n,alpha,A,ia,ja,descA,beta,C,ic,jc,descC)
        return

        end subroutine pgeadd_c21

        subroutine pgeadd_d21(trans,m,n,alpha, A,ia,ja,descA,              &
     &    beta, C,ic,jc,descC)
        character trans
        integer m,n,  ia,ja,descA(*),   ic,jc,descC(*)
        real*8 alpha,beta
        real*8 A(:,:), C(*)

        call pdgeadd(trans,m,n,alpha,A,ia,ja,descA,beta,C,ic,jc,descC)
        return

        end subroutine pgeadd_d21

        subroutine pgeadd_s21(trans,m,n,alpha, A,ia,ja,descA,              &
     &    beta, C,ic,jc,descC)
        character trans
        integer m,n,  ia,ja,descA(*),   ic,jc,descC(*)
        real*4 alpha,beta
        real*4 A(:,:), C(*)

        call psgeadd(trans,m,n,alpha,A,ia,ja,descA,beta,C,ic,jc,descC)
        return

        end subroutine pgeadd_s21

        subroutine pgemr2d_s(m,n,a,ia,ja,adesc,b,ib,jb,bdesc,ctxt)
        integer m,n,ia,ja,ib,jb,ctxt,adesc(*),bdesc(*)
        real*4 a(:,:),b(:,:)
        call psgemr2d(m,n,a,ia,ja,adesc,b,ib,jb,bdesc,ctxt)
        end subroutine pgemr2d_s

        subroutine pgemr2d_d(m,n,a,ia,ja,adesc,b,ib,jb,bdesc,ctxt)
        integer m,n,ia,ja,ib,jb,ctxt,adesc(*),bdesc(*)
        real*8 a(:,:),b(:,:)
        call pdgemr2d(m,n,a,ia,ja,adesc,b,ib,jb,bdesc,ctxt)
        end subroutine pgemr2d_d

        subroutine pgemr2d_c(m,n,a,ia,ja,adesc,b,ib,jb,bdesc,ctxt)
        integer m,n,ia,ja,ib,jb,ctxt,adesc(*),bdesc(*)
        complex*8 a(:,:),b(:,:)
        call pcgemr2d(m,n,a,ia,ja,adesc,b,ib,jb,bdesc,ctxt)
        end subroutine pgemr2d_c

        subroutine pgemr2d_z(m,n,a,ia,ja,adesc,b,ib,jb,bdesc,ctxt)
        integer m,n,ia,ja,ib,jb,ctxt,adesc(*),bdesc(*)
        complex*16 a(:,:),b(:,:)
        call pzgemr2d(m,n,a,ia,ja,adesc,b,ib,jb,bdesc,ctxt)
        end subroutine pgemr2d_z

        subroutine pgemr2d_s1(m,n,a,ia,ja,adesc,b,ib,jb,bdesc,ctxt)
        integer m,n,ia,ja,ib,jb,ctxt,adesc(*),bdesc(*)
        real*4 a(:),b(:)
        call psgemr2d(m,n,a,ia,ja,adesc,b,ib,jb,bdesc,ctxt)
        end subroutine pgemr2d_s1

        subroutine pgemr2d_d1(m,n,a,ia,ja,adesc,b,ib,jb,bdesc,ctxt)
        integer m,n,ia,ja,ib,jb,ctxt,adesc(*),bdesc(*)
        real*8 a(:),b(:)
        call pdgemr2d(m,n,a,ia,ja,adesc,b,ib,jb,bdesc,ctxt)
        end subroutine pgemr2d_d1

        subroutine pgemr2d_c1(m,n,a,ia,ja,adesc,b,ib,jb,bdesc,ctxt)
        integer m,n,ia,ja,ib,jb,ctxt,adesc(*),bdesc(*)
        complex*8 a(:),b(:)
        call pcgemr2d(m,n,a,ia,ja,adesc,b,ib,jb,bdesc,ctxt)
        end subroutine pgemr2d_c1

        subroutine pgemr2d_z1(m,n,a,ia,ja,adesc,b,ib,jb,bdesc,ctxt)
        integer m,n,ia,ja,ib,jb,ctxt,adesc(*),bdesc(*)
        complex*16 a(:),b(:)
        call pzgemr2d(m,n,a,ia,ja,adesc,b,ib,jb,bdesc,ctxt)
        end subroutine pgemr2d_z1


        subroutine gsum2d_i(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        character scope,top
        integer icontxt,m,n,rdest,cdest,lda
        integer a(:,:)
        call igsum2d(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        end subroutine gsum2d_i

        subroutine gsum2d_s(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        character scope,top
        integer icontxt,m,n,rdest,cdest,lda
        real*4 a(:,:)
        call sgsum2d(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        end subroutine gsum2d_s

        subroutine gsum2d_d(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        character scope,top
        integer icontxt,m,n,rdest,cdest,lda
        real*8 a(:,:)
        call dgsum2d(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        end subroutine gsum2d_d

        subroutine gsum2d_c(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        character scope,top
        integer icontxt,m,n,rdest,cdest,lda
        complex*8 a(:,:)
        call cgsum2d(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        end subroutine gsum2d_c

        subroutine gsum2d_z(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        character scope,top
        integer icontxt,m,n,rdest,cdest,lda
        complex*16 a(:,:)
        call zgsum2d(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        end subroutine gsum2d_z

        subroutine gsum2d_i1(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        character scope,top
        integer icontxt,m,n,rdest,cdest,lda
        integer a(:)
        call igsum2d(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        end subroutine gsum2d_i1

        subroutine gsum2d_s1(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        character scope,top
        integer icontxt,m,n,rdest,cdest,lda
        real*4 a(:)
        call sgsum2d(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        end subroutine gsum2d_s1

        subroutine gsum2d_d1(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        character scope,top
        integer icontxt,m,n,rdest,cdest,lda
        real*8 a(:)
        call dgsum2d(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        end subroutine gsum2d_d1

        subroutine gsum2d_c1(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        character scope,top
        integer icontxt,m,n,rdest,cdest,lda
        complex*8 a(:)
        call cgsum2d(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        end subroutine gsum2d_c1

        subroutine gsum2d_z1(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        character scope,top
        integer icontxt,m,n,rdest,cdest,lda
        complex*16 a(:)
        call zgsum2d(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        end subroutine gsum2d_z1

        subroutine gsum2d_i0(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        character scope,top
        integer icontxt,m,n,rdest,cdest,lda
        integer a
        call igsum2d(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        end subroutine gsum2d_i0

        subroutine gsum2d_s0(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        character scope,top
        integer icontxt,m,n,rdest,cdest,lda
        real*4 a
        call sgsum2d(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        end subroutine gsum2d_s0

        subroutine gsum2d_d0(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        character scope,top
        integer icontxt,m,n,rdest,cdest,lda
        real*8 a
        call dgsum2d(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        end subroutine gsum2d_d0

        subroutine gsum2d_c0(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        character scope,top
        integer icontxt,m,n,rdest,cdest,lda
        complex*8 a
        call cgsum2d(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        end subroutine gsum2d_c0

        subroutine gsum2d_z0(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        character scope,top
        integer icontxt,m,n,rdest,cdest,lda
        complex*16 a
        call zgsum2d(icontxt,scope,top,m,n,a,lda,rdest,cdest)
        end subroutine gsum2d_z0


        subroutine pelget_z2( scope, top, alpha, A, ia, ja, descA )
        character scope
        character top
        integer ia,ja,descA(*)
        complex*16 alpha
        complex*16 A(:,:)

        call pzelget(scope,top,alpha,A,ia,ja,descA)
        return

        end subroutine pelget_z2




        subroutine pelget_c2( scope, top, alpha, A, ia, ja, descA )
        character scope
        character top
        integer ia,ja,descA(*)
        complex*8 alpha
        complex*8 A(:,:)

        call pcelget(scope,top,alpha,A,ia,ja,descA)
        return

        end subroutine pelget_c2

        subroutine pelget_d2( scope, top, alpha, A, ia, ja, descA )
        character scope
        character top
        integer ia,ja,descA(*)
        real*8 alpha
        real*8 A(:,:)

        call pdelget(scope,top,alpha,A,ia,ja,descA)
        return
        end subroutine pelget_d2

        subroutine pelget_s2( scope, top, alpha, A, ia, ja, descA )
        character scope
        character top
        integer ia,ja,descA(*)
        real*4 alpha
        real*4 A(:,:)

        call pselget(scope,top,alpha,A,ia,ja,descA)
        return
        end subroutine pelget_s2

        subroutine pelget_z( scope, top, alpha, A, ia, ja, descA )
        character scope
        character top
        integer ia,ja,descA(*)
        complex*16 alpha
        complex*16 A(*)

        call pzelget(scope,top,alpha,A,ia,ja,descA)
        return

        end subroutine pelget_z




        subroutine pelget_c( scope, top, alpha, A, ia, ja, descA )
        character scope
        character top
        integer ia,ja,descA(*)
        complex*8 alpha
        complex*8 A(*)

        call pcelget(scope,top,alpha,A,ia,ja,descA)
        return

        end subroutine pelget_c

        subroutine pelget_d( scope, top, alpha, A, ia, ja, descA )
        character scope
        character top
        integer ia,ja,descA(*)
        real*8 alpha
        real*8 A(*)

        call pdelget(scope,top,alpha,A,ia,ja,descA)
        return
        end subroutine pelget_d

        subroutine pelget_s( scope, top, alpha, A, ia, ja, descA )
        character scope
        character top
        integer ia,ja,descA(*)
        real*4 alpha
        real*4 A(*)

        call pselget(scope,top,alpha,A,ia,ja,descA)
        return
        end subroutine pelget_s


        subroutine gebr2d_z( icontext,scope,top,m,n,A,lda,rsrc,csrc)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        complex*16 A(lda,*)
        integer rsrc,csrc

        call zgebr2d(icontext,scope,top,m,n,A,lda,rsrc,csrc)
        return
        end subroutine gebr2d_z



        subroutine gebr2d_c( icontext,scope,top,m,n,A,lda,rsrc,csrc)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        complex*8 A(lda,*)
        integer rsrc,csrc

        call cgebr2d(icontext,scope,top,m,n,A,lda,rsrc,csrc)
        return
        end subroutine gebr2d_c

        subroutine gebr2d_d( icontext,scope,top,m,n,A,lda,rsrc,csrc)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        real*8 A(lda,*)
        integer rsrc,csrc

        call dgebr2d(icontext,scope,top,m,n,A,lda,rsrc,csrc)
        return
        end subroutine gebr2d_d

        subroutine gebr2d_s( icontext,scope,top,m,n,A,lda,rsrc,csrc)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        real*4 A(lda,*)
        integer rsrc,csrc

        call sgebr2d(icontext,scope,top,m,n,A,lda,rsrc,csrc)
        return
        end subroutine gebr2d_s


        subroutine gebr2d_z1( icontext,scope,top,m,n,A,lda,rsrc,csrc)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        complex*16 A(*)
        integer rsrc,csrc

        call zgebr2d(icontext,scope,top,m,n,A,lda,rsrc,csrc)
        return
        end subroutine gebr2d_z1



        subroutine gebr2d_c1( icontext,scope,top,m,n,A,lda,rsrc,csrc)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        complex*8 A(*)
        integer rsrc,csrc

        call cgebr2d(icontext,scope,top,m,n,A,lda,rsrc,csrc)
        return
        end subroutine gebr2d_c1

        subroutine gebr2d_d1( icontext,scope,top,m,n,A,lda,rsrc,csrc)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        real*8 A(*)
        integer rsrc,csrc

        call dgebr2d(icontext,scope,top,m,n,A,lda,rsrc,csrc)
        return
        end subroutine gebr2d_d1

        subroutine gebr2d_s1( icontext,scope,top,m,n,A,lda,rsrc,csrc)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        real*4 A(*)
        integer rsrc,csrc

        call sgebr2d(icontext,scope,top,m,n,A,lda,rsrc,csrc)
        return
        end subroutine gebr2d_s1

        subroutine gebr2d_s0( icontext,scope,top,m,n,A,lda,rsrc,csrc)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        real*4  A
        real*4 A1(1,1)
        integer rsrc,csrc

        call sgebr2d(icontext,scope,top,m,n,A1,lda,rsrc,csrc)

        A = A1(1,1)
        return
        end subroutine gebr2d_s0


        subroutine gebr2d_d0( icontext,scope,top,m,n,A,lda,rsrc,csrc)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        real*8  A
        real*8 A1(1,1)
        integer rsrc,csrc

        call dgebr2d(icontext,scope,top,m,n,A1,lda,rsrc,csrc)

        A = A1(1,1)
        return
        end subroutine gebr2d_d0


        subroutine gebr2d_c0( icontext,scope,top,m,n,A,lda,rsrc,csrc)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        complex*8  A
        complex*8 A1(1,1)
        integer rsrc,csrc

        call cgebr2d(icontext,scope,top,m,n,A1,lda,rsrc,csrc)

        A = A1(1,1)
        return
        end subroutine gebr2d_c0



        subroutine gebr2d_z0( icontext,scope,top,m,n,A,lda,rsrc,csrc)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        complex*16  A
        complex*16 A1(1,1)
        integer rsrc,csrc

        call zgebr2d(icontext,scope,top,m,n,A1,lda,rsrc,csrc)

        A = A1(1,1)
        return
        end subroutine gebr2d_z0

        subroutine gebs2d_z( icontext,scope,top,m,n,A,lda)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        complex*16 A(lda,*)

        call zgebs2d(icontext,scope,top,m,n,A,lda)
        return
        end subroutine gebs2d_z

        subroutine gebs2d_c( icontext,scope,top,m,n,A,lda)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        complex*8 A(lda,*)

        call cgebs2d(icontext,scope,top,m,n,A,lda)
        return
        end subroutine gebs2d_c

        subroutine gebs2d_d( icontext,scope,top,m,n,A,lda)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        real*8 A(lda,*)

        call dgebs2d(icontext,scope,top,m,n,A,lda)
        return
        end subroutine gebs2d_d

        subroutine gebs2d_s( icontext,scope,top,m,n,A,lda)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        real*4 A(lda,*)

        call sgebs2d(icontext,scope,top,m,n,A,lda)
        return
        end subroutine gebs2d_s

        subroutine gebs2d_z1( icontext,scope,top,m,n,A,lda)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        complex*16 A(*)

        call zgebs2d(icontext,scope,top,m,n,A,lda)
        return
        end subroutine gebs2d_z1

        subroutine gebs2d_c1( icontext,scope,top,m,n,A,lda)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        complex*8 A(*)

        call cgebs2d(icontext,scope,top,m,n,A,lda)
        return
        end subroutine gebs2d_c1

        subroutine gebs2d_d1( icontext,scope,top,m,n,A,lda)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        real*8 A(*)

        call dgebs2d(icontext,scope,top,m,n,A,lda)
        return
        end subroutine gebs2d_d1

        subroutine gebs2d_s1( icontext,scope,top,m,n,A,lda)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        real*4 A(*)

        call sgebs2d(icontext,scope,top,m,n,A,lda)
        return
        end subroutine gebs2d_s1

        subroutine gesd2d_z( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        complex*16 A(lda,*)

        call zgesd2d(icontext,m,n,A,lda)
        return
        end subroutine gesd2d_z

        subroutine gesd2d_c( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        complex*8 A(lda,*)

        call cgesd2d(icontext,m,n,A,lda)
        return
        end subroutine gesd2d_c

        subroutine gesd2d_d( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        real*8 A(lda,*)

        call dgesd2d(icontext,m,n,A,lda)
        return
        end subroutine gesd2d_d

        subroutine gesd2d_s( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        real*4 A(lda,*)

        call sgesd2d(icontext,m,n,A,lda)
        return
        end subroutine gesd2d_s

        subroutine gesd2d_z1( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        complex*16 A(*)

        call zgesd2d(icontext,m,n,A,lda)
        return
        end subroutine gesd2d_z1

        subroutine gesd2d_c1( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        complex*8 A(*)

        call cgesd2d(icontext,m,n,A,lda)
        return
        end subroutine gesd2d_c1

        subroutine gesd2d_d1( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        real*8 A(*)

        call dgesd2d(icontext,m,n,A,lda)
        return
        end subroutine gesd2d_d1

        subroutine gesd2d_s1( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        real*4 A(*)

        call sgesd2d(icontext,m,n,A,lda)
        return
        end subroutine gesd2d_s1

        subroutine gesd2d_z0( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        complex*16 A

        call zgesd2d(icontext,m,n,A,lda)
        return
        end subroutine gesd2d_z0

        subroutine gesd2d_c0( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        complex*8 A

        call cgesd2d(icontext,m,n,A,lda)
        return
        end subroutine gesd2d_c0

        subroutine gesd2d_d0( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        real*8 A

        call dgesd2d(icontext,m,n,A,lda)
        return
        end subroutine gesd2d_d0

        subroutine gesd2d_s0( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        real*4 A

        call sgesd2d(icontext,m,n,A,lda)
        return
        end subroutine gesd2d_s0


        subroutine gerv2d_z( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        complex*16 A(lda,*)

        call zgerv2d(icontext,m,n,A,lda)
        return
        end subroutine gerv2d_z

        subroutine gerv2d_c( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        complex*8 A(lda,*)

        call cgerv2d(icontext,m,n,A,lda)
        return
        end subroutine gerv2d_c

        subroutine gerv2d_d( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        real*8 A(lda,*)

        call dgerv2d(icontext,m,n,A,lda)
        return
        end subroutine gerv2d_d

        subroutine gerv2d_s( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        real*4 A(lda,*)

        call sgerv2d(icontext,m,n,A,lda)
        return
        end subroutine gerv2d_s

        subroutine gerv2d_z1( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        complex*16 A(*)

        call zgerv2d(icontext,m,n,A,lda)
        return
        end subroutine gerv2d_z1

        subroutine gerv2d_c1( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        complex*8 A(*)

        call cgerv2d(icontext,m,n,A,lda)
        return
        end subroutine gerv2d_c1

        subroutine gerv2d_d1( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        real*8 A(*)

        call dgerv2d(icontext,m,n,A,lda)
        return
        end subroutine gerv2d_d1

        subroutine gerv2d_s1( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        real*4 A(*)

        call sgerv2d(icontext,m,n,A,lda)
        return
        end subroutine gerv2d_s1

        subroutine gerv2d_z0( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        complex*16 A

        call zgerv2d(icontext,m,n,A,lda)
        return
        end subroutine gerv2d_z0

        subroutine gerv2d_c0( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        complex*8 A

        call cgerv2d(icontext,m,n,A,lda)
        return
        end subroutine gerv2d_c0

        subroutine gerv2d_d0( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        real*8 A

        call dgerv2d(icontext,m,n,A,lda)
        return
        end subroutine gerv2d_d0

        subroutine gerv2d_s0( icontext,m,n,A,lda,rdest,cdest)
        implicit none
        integer icontext
        integer m,n,lda,rdest,cdest
        real*4 A

        call sgerv2d(icontext,m,n,A,lda)
        return
        end subroutine gerv2d_s0

        subroutine gebs2d_s0( icontext,scope,top,m,n,A,lda)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        real*4 A
        real*4 A1(1,1)

        A1(1,1) = A
        call sgebs2d(icontext,scope,top,m,n,A1,lda)
        return
        end subroutine gebs2d_s0





        subroutine gebs2d_d0( icontext,scope,top,m,n,A,lda)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        real*8 A
        real*8 A1(1,1)

        A1(1,1) = A
        call dgebs2d(icontext,scope,top,m,n,A1,lda)
        return
        end subroutine gebs2d_d0



        subroutine gebs2d_c0( icontext,scope,top,m,n,A,lda)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        complex*8 A
        complex*8 A1(1,1)

        A1(1,1) = A
        call cgebs2d(icontext,scope,top,m,n,A1,lda)
        return
        end subroutine gebs2d_c0



        subroutine gebs2d_z0( icontext,scope,top,m,n,A,lda)
        implicit none
        integer icontext
        character scope
        character top
        integer m,n,lda
        complex*16 A
        complex*16 A1(1,1)

        A1(1,1) = A
        call zgebs2d(icontext,scope,top,m,n,A1,lda)
        return
        end subroutine gebs2d_z0

        subroutine pgemm_d( transA, transB,                                  &
     &       m,n,k,                                                          &
     &       alpha,   A,ia,ja,descA,                                         &
     &                B,ib,jb,descB,                                         &
     &       beta,    C,ic,jc,descC )

        implicit none
        character transA,transB
        integer m,n,k
        integer ia,ja,  ib,jb,  ic,jc
        integer descA(*),descB(*),descC(*)
        real*8 alpha,beta
        real*8  A(*),B(*),C(*)

        call pdgemm( transA, transB,                                         &
     &       m,n,k,                                                          &
     &       alpha,   A,ia,ja,descA,                                         &
     &                B,ib,jb,descB,                                         &
     &       beta,    C,ic,jc,descC )
        return
        end subroutine pgemm_d



        subroutine pgemm_z( transA, transB,                                  &
     &       m,n,k,                                                          &
     &       alpha,   A,ia,ja,descA,                                         &
     &                B,ib,jb,descB,                                         &
     &       beta,    C,ic,jc,descC )

        implicit none
        character transA,transB
        integer m,n,k
        integer ia,ja,  ib,jb,  ic,jc
        integer descA(*),descB(*),descC(*)
        complex*16 alpha,beta
        complex*16  A(*),B(*),C(*)

        call pzgemm( transA, transB,                                          &
     &       m,n,k,                                                           &
     &       alpha,   A,ia,ja,descA,                                          &
     &                B,ib,jb,descB,                                          &
     &       beta,    C,ic,jc,descC )
        return
        end subroutine pgemm_z


        subroutine pgemm_c( transA, transB,                                  &
     &       m,n,k,                                                          &
     &       alpha,   A,ia,ja,descA,                                         &
     &                B,ib,jb,descB,                                         &
     &       beta,    C,ic,jc,descC )

        implicit none
        character transA,transB
        integer m,n,k
        integer ia,ja,  ib,jb,  ic,jc
        integer descA(*),descB(*),descC(*)
        complex*8 alpha,beta
        complex*8  A(*),B(*),C(*)

        call pcgemm( transA, transB,                                         &
     &       m,n,k,                                                          &
     &       alpha,   A,ia,ja,descA,                                         &
     &                B,ib,jb,descB,                                         &
     &       beta,    C,ic,jc,descC )
        return
        end subroutine pgemm_c

        
        subroutine pgemm_s( transA, transB,                                  &
     &       m,n,k,                                                          &
     &       alpha,   A,ia,ja,descA,                                         &
     &                B,ib,jb,descB,                                         &
     &       beta,    C,ic,jc,descC )

        implicit none
        character transA,transB
        integer m,n,k
        integer ia,ja,  ib,jb,  ic,jc
        integer descA(*),descB(*),descC(*)
        real*4 alpha,beta
        real*4  A(*),B(*),C(*)

        call psgemm( transA, transB,                                         &
     &       m,n,k,                                                          &
     &       alpha,   A,ia,ja,descA,                                         &
     &                B,ib,jb,descB,                                         &
     &       beta,    C,ic,jc,descC )
        return
        end subroutine pgemm_s

        subroutine pgemm_d2( transA, transB,                                 &
     &       m,n,k,                                                          &
     &       alpha,   A,ia,ja,descA,                                         &
     &                B,ib,jb,descB,                                         &
     &       beta,    C,ic,jc,descC ) 

        implicit none
        character transA,transB
        integer m,n,k
        integer ia,ja,  ib,jb,  ic,jc
        integer descA(*),descB(*),descC(*)
        real*8 alpha,beta
        real*8  A(:,:),B(:,:),C(:,:)

        call pdgemm( transA, transB,                                         &
     &       m,n,k,                                                          &
     &       alpha,   A,ia,ja,descA,                                         &
     &                B,ib,jb,descB,                                         &
     &       beta,    C,ic,jc,descC )
        return
        end subroutine pgemm_d2



        subroutine pgemm_z2( transA, transB,                                 &
     &       m,n,k,                                                          &
     &       alpha,   A,ia,ja,descA,                                         &
     &                B,ib,jb,descB,                                         &
     &       beta,    C,ic,jc,descC )

        implicit none
        character transA,transB
        integer m,n,k
        integer ia,ja,  ib,jb,  ic,jc
        integer descA(*),descB(*),descC(*)
        complex*16 alpha,beta
        complex*16  A(:,:),B(:,:),C(:,:)

        call pzgemm( transA, transB,                                         &
     &       m,n,k,                                                          &
     &       alpha,   A,ia,ja,descA,                                         &
     &                B,ib,jb,descB,                                         &
     &       beta,    C,ic,jc,descC )
        return
        end subroutine pgemm_z2


        subroutine pgemm_c2( transA, transB,                                &
     &       m,n,k,                                                         &
     &       alpha,   A,ia,ja,descA,                                        &
     &                B,ib,jb,descB,                                        &
     &       beta,    C,ic,jc,descC )

        implicit none
        character transA,transB
        integer m,n,k
        integer ia,ja,  ib,jb,  ic,jc
        integer descA(*),descB(*),descC(*)
        complex*8 alpha,beta
        complex*8  A(:,:),B(:,:),C(:,:)

        call pcgemm( transA, transB,                                         &
     &       m,n,k,                                                          &
     &       alpha,   A,ia,ja,descA,                                         &
     &                B,ib,jb,descB,                                         &
     &       beta,    C,ic,jc,descC )
        return
        end subroutine pgemm_c2

        
        subroutine pgemm_s2( transA, transB,                                 &
     &       m,n,k,                                                          &
     &       alpha,   A,ia,ja,descA,                                         &
     &                B,ib,jb,descB,                                         &
     &       beta,    C,ic,jc,descC )

        implicit none
        character transA,transB
        integer m,n,k
        integer ia,ja,  ib,jb,  ic,jc
        integer descA(*),descB(*),descC(*)
        real*4 alpha,beta
        real*4  A(:,:),B(:,:),C(:,:)

        call psgemm( transA, transB,                                           &
     &       m,n,k,                                                            &
     &       alpha,   A,ia,ja,descA,                                           &
     &                B,ib,jb,descB,                                           &
     &       beta,    C,ic,jc,descC )
        return
        end subroutine pgemm_s2


        subroutine pger_d( m,n,                                                &
     &       alpha,   X,ix,jx,descX,incX,                                      &
     &                Y,iy,jy,descY,incY,                                      &
     &                A,ia,ja,descA )

        implicit none
        integer m,n
        integer ix,jx,incX,  iy,jy,incY,  ia,ja
        integer descX(*),descY(*),descA(*)
        real*8 alpha
        real*8  X(*),Y(*),A(*)

        call pdger( m, n,                                                      &
     &       alpha,   X,ix,jx,descX,incX,                                      &
     &                Y,iy,jy,descY,incY,                                      &
     &                A,ia,ja,descA )
        return
        end subroutine pger_d

        subroutine pger_s( m,n,                                                &
     &       alpha,   X,ix,jx,descX,incX,                                      &
     &                Y,iy,jy,descY,incY,                                      &
     &                A,ia,ja,descA )

        implicit none
        integer m,n
        integer ix,jx,incX,  iy,jy,incY,  ia,ja
        integer descX(*),descY(*),descA(*)
        real*4 alpha
        real*4  X(*),Y(*),A(*)

        call psger( m, n,                                                      &
     &       alpha,   X,ix,jx,descX,incX,                                      &
     &                Y,iy,jy,descY,incY,                                      &
     &                A,ia,ja,descA )
        return
        end subroutine pger_s

        subroutine pger_z( m,n,                                                &
     &       alpha,   X,ix,jx,descX,incX,                                      &
     &                Y,iy,jy,descY,incY,                                      &
     &                A,ia,ja,descA )

        implicit none
        integer m,n
        integer ix,jx,incX,  iy,jy,incY,  ia,ja
        integer descX(*),descY(*),descA(*)
        complex*16 alpha
        complex*16  X(*),Y(*),A(*)

        call pzgeru( m, n,                                                     &
     &       alpha,   X,ix,jx,descX,incX,                                      &
     &                Y,iy,jy,descY,incY,                                      &
     &                A,ia,ja,descA )
        return
        end subroutine pger_z

        subroutine pger_c( m,n,                                                &
     &       alpha,   X,ix,jx,descX,incX,                                      &
     &                Y,iy,jy,descY,incY,                                      &
     &                A,ia,ja,descA )

        implicit none
        integer m,n
        integer ix,jx,incX,  iy,jy,incY,  ia,ja
        integer descX(*),descY(*),descA(*)
        complex*8 alpha
        complex*8  X(*),Y(*),A(*)

        call pcgeru( m, n,                                                     &
     &       alpha,   X,ix,jx,descX,incX,                                      &
     &                Y,iy,jy,descY,incY,                                      &
     &                A,ia,ja,descA )
        return
        end subroutine pger_c

        subroutine pger_d2( m,n,                                               &
     &       alpha,   X,ix,jx,descX,incX,                                      &
     &                Y,iy,jy,descY,incY,                                      &
     &                A,ia,ja,descA )

        implicit none
        integer m,n
        integer ix,jx,incX,  iy,jy,incY,  ia,ja
        integer descX(*),descY(*),descA(*)
        real*8 alpha
        real*8  X(:,:),Y(:,:),A(:,:)

        call pdger( m, n,                                                      &
     &       alpha,   X,ix,jx,descX,incX,                                      &
     &                Y,iy,jy,descY,incY,                                      &
     &                A,ia,ja,descA )
        return
        end subroutine pger_d2

        subroutine pger_s2( m,n,                                               &
     &       alpha,   X,ix,jx,descX,incX,                                      &
     &                Y,iy,jy,descY,incY,                                      &
     &                A,ia,ja,descA )

        implicit none
        integer m,n
        integer ix,jx,incX,  iy,jy,incY,  ia,ja
        integer descX(*),descY(*),descA(*)
        real*4 alpha
        real*4  X(:,:),Y(:,:),A(:,:)

        call psger( m, n,                                                      &
     &       alpha,   X,ix,jx,descX,incX,                                      &
     &                Y,iy,jy,descY,incY,                                      &
     &                A,ia,ja,descA )
        return
        end subroutine pger_s2

        subroutine pger_z2( m,n,                                               &
     &       alpha,   X,ix,jx,descX,incX,                                      &
     &                Y,iy,jy,descY,incY,                                      &
     &                A,ia,ja,descA )

        implicit none
        integer m,n
        integer ix,jx,incX,  iy,jy,incY,  ia,ja
        integer descX(*),descY(*),descA(*)
        complex*16 alpha
        complex*16  X(:,:),Y(:,:),A(:,:)

        call pzgeru( m, n,                                                     &
     &       alpha,   X,ix,jx,descX,incX,                                      &
     &                Y,iy,jy,descY,incY,                                      &
     &                A,ia,ja,descA )
        return
        end subroutine pger_z2

        subroutine pger_c2( m,n,                                               &
     &       alpha,   X,ix,jx,descX,incX,                                      &
     &                Y,iy,jy,descY,incY,                                      &
     &                A,ia,ja,descA )

        implicit none
        integer m,n
        integer ix,jx,incX,  iy,jy,incY,  ia,ja
        integer descX(*),descY(*),descA(*)
        complex*8 alpha
        complex*8  X(:,:),Y(:,:),A(:,:)

        call pcgeru( m, n,                                                     &
     &       alpha,   X,ix,jx,descX,incX,                                      &
     &                Y,iy,jy,descY,incY,                                      &
     &                A,ia,ja,descA )
        return
        end subroutine pger_c2

        real*4 function plange_c( norm,m,n,A,ia,ja,descA,work)
        implicit none
        character norm
        integer m,n,ia,ja,descA(*)
        complex*8 A(*)
        real*4 work(*)

        plange_c = pclange(norm,m,n,A,ia,ja,descA,work)
        return
        end function plange_c



        real*8 function plange_z( norm,m,n,A,ia,ja,descA,work)
        implicit none
        character norm
        integer m,n,ia,ja,descA(*)
        complex*16 A(*)
        real*8 work(*)

        plange_z = pzlange(norm,m,n,A,ia,ja,descA,work)
        return
        end function plange_z


        real*4 function plange_s( norm,m,n,A,ia,ja,descA,work)
        implicit none
        character norm
        integer m,n,ia,ja,descA(*)
        real*4 A(*)
        real*4 work(*)

        plange_s = pslange(norm,m,n,A,ia,ja,descA,work)
        return
        end function plange_s


        real*8 function plange_d( norm,m,n,A,ia,ja,descA,work)
        implicit none
        character norm
        integer m,n,ia,ja,descA(*)
        real*8 A(*)
        real*8 work(*)

        plange_d = pdlange(norm,m,n,A,ia,ja,descA,work)
        return
        end function plange_d

        real*4 function plange_c2( norm,m,n,A,ia,ja,descA,work)
        implicit none
        character norm
        integer m,n,ia,ja,descA(*)
        complex*8 A(:,:)
        real*4 work(*)

        plange_c2 = pclange(norm,m,n,A,ia,ja,descA,work)
        return
        end function plange_c2



        real*8 function plange_z2( norm,m,n,A,ia,ja,descA,work)
        implicit none
        character norm
        integer m,n,ia,ja,descA(*)
        complex*16 A(:,:)
        real*8 work(*)

        plange_z2 = pzlange(norm,m,n,A,ia,ja,descA,work)
        return
        end function plange_z2


        real*4 function plange_s2( norm,m,n,A,ia,ja,descA,work)
        implicit none
        character norm
        integer m,n,ia,ja,descA(*)
        real*4 A(:,:)
        real*4 work(*)

        plange_s2 = pslange(norm,m,n,A,ia,ja,descA,work)
        return
        end function plange_s2


        real*8 function plange_d2( norm,m,n,A,ia,ja,descA,work)
        implicit none
        character norm
        integer m,n,ia,ja,descA(*)
        real*8 A(:,:)
        real*8 work(*)

        plange_d2 = pdlange(norm,m,n,A,ia,ja,descA,work)
        return
        end function plange_d2

        subroutine pdot_d22(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        implicit none
        integer n,ix,jx,descx(*),incx,iy,jy,descy(*),incy
        real*8 dot,x(:,:),y(:,:)

        call pddot(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        return
        end subroutine pdot_d22

        subroutine pdot_s22(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        implicit none
        integer n,ix,jx,descx(*),incx,iy,jy,descy(*),incy
        real*4 dot,x(:,:),y(:,:)

        call psdot(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        return
        end subroutine pdot_s22

        subroutine pdot_d(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        implicit none
        integer n,ix,jx,descx(*),incx,iy,jy,descy(*),incy
        real*8 dot,x(*),y(*)

        call pddot(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        return
        end subroutine pdot_d

        subroutine pdot_s(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        implicit none
        integer n,ix,jx,descx(*),incx,iy,jy,descy(*),incy
        real*4 dot,x(*),y(*)

        call psdot(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        return
        end subroutine pdot_s


        subroutine pdotc_c(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        implicit none
        integer n,ix,jx,descx(*),incx,iy,jy,descy(*),incy
        complex*8 dot,x(*),y(*)

        call pcdotc(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        return
        end subroutine pdotc_c

        subroutine pdotu_c(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        implicit none
        integer n,ix,jx,descx(*),incx,iy,jy,descy(*),incy
        complex*8 dot,x(*),y(*)

        call pcdotu(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        return
        end subroutine pdotu_c




        subroutine pdotc_z(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        implicit none
        integer n,ix,jx,descx(*),incx,iy,jy,descy(*),incy
        complex*16 dot,x(*),y(*)

        call pzdotc(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        return
        end subroutine pdotc_z

        subroutine pdotu_z(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        implicit none
        integer n,ix,jx,descx(*),incx,iy,jy,descy(*),incy
        complex*16 dot,x(*),y(*)

        call pzdotu(n,dot,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        return
        end subroutine pdotu_z




      subroutine ptran_s(m,n,alpha,a,ia,ja,desca,beta,c,ic,jc,descc)
      implicit none
      integer m,n,ia,ja,desca(*),ic,jc,descc(*)
      real*4 a(:,:),c(:,:),alpha,beta

      call pstran(m,n,alpha,a,ia,ja,desca,beta,c,ic,jc,descc)
      return
      end subroutine ptran_s
      
      subroutine ptran_d(m,n,alpha,a,ia,ja,desca,beta,c,ic,jc,descc)
      implicit none
      integer m,n,ia,ja,desca(*),ic,jc,descc(*)
      real*8 a(:,:),c(:,:),alpha,beta

      call pdtran(m,n,alpha,a,ia,ja,desca,beta,c,ic,jc,descc)
      return
      end subroutine ptran_d
      

        subroutine placp3_c( M, I, A, DESCA, B, LDB, II, JJ, REV)
        implicit none
        INTEGER            I, II, JJ, LDB, M, REV
        INTEGER            DESCA( * )
        COMPLEX   A( * ), B( LDB, * )

        call pclacp3( M, I, A, DESCA, B, LDB, II, JJ, REV)
        return
        end subroutine placp3_c

        subroutine placp3_z( M, I, A, DESCA, B, LDB, II, JJ, REV)
        implicit none
        INTEGER            I, II, JJ, LDB, M, REV
        INTEGER            DESCA( * )
        COMPLEX*16   A( * ), B( LDB, * )

        call pzlacp3( M, I, A, DESCA, B, LDB, II, JJ, REV)
        return
        end subroutine placp3_z

        subroutine placp3_s( M, I, A, DESCA, B, LDB, II, JJ, REV)
        implicit none
        INTEGER            I, II, JJ, LDB, M, REV
        INTEGER            DESCA( * )
        REAL   A( * ), B( LDB, * )

        call pslacp3( M, I, A, DESCA, B, LDB, II, JJ, REV)
        return
        end subroutine placp3_s

        subroutine placp3_d( M, I, A, DESCA, B, LDB, II, JJ, REV)
        implicit none
        INTEGER            I, II, JJ, LDB, M, REV
        INTEGER            DESCA( * )
        DOUBLE PRECISION   A( * ), B( LDB, * )

        call pdlacp3( M, I, A, DESCA, B, LDB, II, JJ, REV)
        return
        end subroutine placp3_d



        subroutine plaprnt_z(m,n,A,ia,ja,descA,irprnt,icprnt,                  &
     &       cmatnm,nout,work)
        implicit none
        integer m,n,ia,ja,descA(*),irprnt,icprnt,nout
        complex*16 A(*),work(*)
        character*(*) cmatnm

        call pzlaprnt(m,n,A,ia,ja,descA,irprnt,icprnt,cmatnm,nout,work)
        return
        end subroutine plaprnt_z

        subroutine plaprnt_c(m,n,A,ia,ja,descA,irprnt,icprnt,                  &
     &       cmatnm,nout,work)
        implicit none
        integer m,n,ia,ja,descA(*),irprnt,icprnt,nout
        complex*8 A(*),work(*)
        character*(*) cmatnm

        call pclaprnt(m,n,A,ia,ja,descA,irprnt,icprnt,cmatnm,nout,work)
        return
        end subroutine plaprnt_c

        subroutine plaprnt_d(m,n,A,ia,ja,descA,irprnt,icprnt,                  &
     &       cmatnm,nout,work)
        implicit none
        integer m,n,ia,ja,descA(*),irprnt,icprnt,nout
        real*8 A(*),work(*)
        character*(*) cmatnm

        call pdlaprnt(m,n,A,ia,ja,descA,irprnt,icprnt,cmatnm,nout,work)
        return
        end subroutine plaprnt_d

        subroutine plaprnt_s(m,n,A,ia,ja,descA,irprnt,icprnt,                  &
     &       cmatnm,nout,work)
        implicit none
        integer m,n,ia,ja,descA(*),irprnt,icprnt,nout
        real*4 A(*),work(*)
        character*(*) cmatnm

        call pslaprnt(m,n,A,ia,ja,descA,irprnt,icprnt,cmatnm,nout,work)
        return
        end subroutine plaprnt_s


        subroutine plaprnt_z2(m,n,A,ia,ja,descA,irprnt,icprnt,                 &
     &       cmatnm,nout,work)
        implicit none
        integer m,n,ia,ja,descA(*),irprnt,icprnt,nout
        complex*16 A(:,:),work(*)
        character*(*) cmatnm

        call pzlaprnt(m,n,A,ia,ja,descA,irprnt,icprnt,cmatnm,nout,work)
        return
        end subroutine plaprnt_z2

        subroutine plaprnt_c2(m,n,A,ia,ja,descA,irprnt,icprnt,                 &
     &       cmatnm,nout,work)
        implicit none
        integer m,n,ia,ja,descA(*),irprnt,icprnt,nout
        complex*8 A(:,:),work(*)
        character*(*) cmatnm

        call pclaprnt(m,n,A,ia,ja,descA,irprnt,icprnt,cmatnm,nout,work)
        return
        end subroutine plaprnt_c2

        subroutine plaprnt_d2(m,n,A,ia,ja,descA,irprnt,icprnt,                 &
     &       cmatnm,nout,work)
        implicit none
        integer m,n,ia,ja,descA(*),irprnt,icprnt,nout
        real*8 A(:,:),work(*)
        character*(*) cmatnm

        call pdlaprnt(m,n,A,ia,ja,descA,irprnt,icprnt,cmatnm,nout,work)
        return
        end subroutine plaprnt_d2

        subroutine plaprnt_s2(m,n,A,ia,ja,descA,irprnt,icprnt,                 &
     &       cmatnm,nout,work)
        implicit none
        integer m,n,ia,ja,descA(*),irprnt,icprnt,nout
        real*4 A(:,:),work(*)
        character*(*) cmatnm

        call pslaprnt(m,n,A,ia,ja,descA,irprnt,icprnt,cmatnm,nout,work)
        return
        end subroutine plaprnt_s2




        subroutine pelset_z(A,ia,ja,descA,alpha)
        integer ia,ja,descA(*)
        complex*16 alpha
        complex*16 A(*)

        call pzelset(A,ia,ja,descA,alpha)
        return
        end subroutine pelset_z

        subroutine pelset_c(A,ia,ja,descA,alpha)
        integer ia,ja,descA(*)
        complex*8 alpha
        complex*8 A(*)

        call pcelset(A,ia,ja,descA,alpha)
        return
        end subroutine pelset_c

        subroutine pelset_d(A,ia,ja,descA,alpha)
        integer ia,ja,descA(*)
        real*8 alpha
        real*8 A(*)

        call pdelset(A,ia,ja,descA,alpha)
        return
        end subroutine pelset_d

        subroutine pelset_s(A,ia,ja,descA,alpha)
        integer ia,ja,descA(*)
        real*4 alpha
        real*4 A(*)

        call pselset(A,ia,ja,descA,alpha)
        return
        end subroutine pelset_s

        subroutine pelset_z2(A,ia,ja,descA,alpha)
        integer ia,ja,descA(*)
        complex*16 alpha
        complex*16 A(:,:)

        call pzelset(A,ia,ja,descA,alpha)
        return
        end subroutine pelset_z2

        subroutine pelset_c2(A,ia,ja,descA,alpha)
        integer ia,ja,descA(*)
        complex*8 alpha
        complex*8 A(:,:)

        call pcelset(A,ia,ja,descA,alpha)
        return
        end subroutine pelset_c2

        subroutine pelset_d2(A,ia,ja,descA,alpha)
        integer ia,ja,descA(*)
        real*8 alpha
        real*8 A(:,:)

        call pdelset(A,ia,ja,descA,alpha)
        return
        end subroutine pelset_d2

        subroutine pelset_s2(A,ia,ja,descA,alpha)
        integer ia,ja,descA(*)
        real*4 alpha
        real*4 A(:,:)

        call pselset(A,ia,ja,descA,alpha)
        return
        end subroutine pelset_s2


        subroutine pcopy_z(n,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        implicit none
        integer n,ix,jx,descx(*),incx,iy,jy,descy(*),incy
        complex*16 x(:,:), y(:,:)

        call pzcopy(n,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        return

        end subroutine pcopy_z

        subroutine pcopy_c(n,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        implicit none
        integer n,ix,jx,descx(*),incx,iy,jy,descy(*),incy
        complex*8 x(:,:), y(:,:)

        call pccopy(n,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        return

        end subroutine pcopy_c


        subroutine pcopy_s(n,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        implicit none
        integer n,ix,jx,descx(*),incx,iy,jy,descy(*),incy
        real*4 x(:,:), y(:,:)

        call pscopy(n,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        return

        end subroutine pcopy_s

        subroutine pcopy_d(n,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        implicit none
        integer n,ix,jx,descx(*),incx,iy,jy,descy(*),incy
        real*8 x(:,:), y(:,:)

        call pdcopy(n,x,ix,jx,descx,incx,y,iy,jy,descy,incy)
        return

        end subroutine pcopy_d

        subroutine pcopy2d_z(m,n,A,ia,ja,descA,B,ib,jb,descB)
        implicit none
        integer m,n,ia,ja,descA(*),ib,jb,descB(*)
        complex*16 A(*), B(*)
        complex*16 alpha,beta

        alpha = 1.0d0
        beta = 0.0d0
        call pgeadd('N',m,n,alpha,A,ia,ja,descA,beta,B,ib,jb,descB)
        return

        end subroutine pcopy2d_z


        subroutine pcopy2d_c(m,n,A,ia,ja,descA,B,ib,jb,descB)
        implicit none
        integer m,n,ia,ja,descA(*),ib,jb,descB(*)
        complex*8 A(*), B(*)
        complex*8 alpha,beta

        alpha = 1.0d0
        beta = 0.0d0
        call pgeadd('N',m,n,alpha,A,ia,ja,descA,beta,B,ib,jb,descB)
        return

        end subroutine pcopy2d_c


        subroutine pcopy2d_d(m,n,A,ia,ja,descA,B,ib,jb,descB)
        implicit none
        integer m,n,ia,ja,descA(*),ib,jb,descB(*)
        real*8 A(*), B(*)
        real*8 alpha,beta

        alpha = 1.0d0
        beta = 0.0d0
        call pgeadd('N',m,n,alpha,A,ia,ja,descA,beta,B,ib,jb,descB)
        return

        end subroutine pcopy2d_d

        subroutine pcopy2d_s(m,n,A,ia,ja,descA,B,ib,jb,descB)
        implicit none
        integer m,n,ia,ja,descA(*),ib,jb,descB(*)
        real*4 A(*), B(*)
        real*4 alpha,beta

        alpha = 1.0d0
        beta = 0.0d0
        call pgeadd('N',m,n,alpha,A,ia,ja,descA,beta,B,ib,jb,descB)
        return

        end subroutine pcopy2d_s

        subroutine pcopy2d_z2(m,n,A,ia,ja,descA,B,ib,jb,descB)
        implicit none
        integer m,n,ia,ja,descA(*),ib,jb,descB(*)
        complex*16 A(:,:), B(:,:)
        complex*16 alpha,beta

        alpha = 1.0d0
        beta = 0.0d0
        call pgeadd('N',m,n,alpha,A,ia,ja,descA,beta,B,ib,jb,descB)
        return

        end subroutine pcopy2d_z2


        subroutine pcopy2d_c2(m,n,A,ia,ja,descA,B,ib,jb,descB)
        implicit none
        integer m,n,ia,ja,descA(*),ib,jb,descB(*)
        complex*8 A(:,:), B(:,:)
        complex*8 alpha,beta

        alpha = 1.0d0
        beta = 0.0d0
        call pgeadd('N',m,n,alpha,A,ia,ja,descA,beta,B,ib,jb,descB)
        return

        end subroutine pcopy2d_c2


        subroutine pcopy2d_d2(m,n,A,ia,ja,descA,B,ib,jb,descB)
        implicit none
        integer m,n,ia,ja,descA(*),ib,jb,descB(*)
        real*8 A(:,:), B(:,:)
        real*8 alpha,beta

        alpha = 1.0d0
        beta = 0.0d0
        call pgeadd('N',m,n,alpha,A,ia,ja,descA,beta,B,ib,jb,descB)
        return

        end subroutine pcopy2d_d2

        subroutine pcopy2d_s2(m,n,A,ia,ja,descA,B,ib,jb,descB)
        implicit none
        integer m,n,ia,ja,descA(*),ib,jb,descB(*)
        real*4 A(:,:), B(:,:)
        real*4 alpha,beta

        alpha = 1.0d0
        beta = 0.0d0
        call pgeadd('N',m,n,alpha,A,ia,ja,descA,beta,B,ib,jb,descB)
        return

        end subroutine pcopy2d_s2


        subroutine pgemv_z111(trans, m,n, alpha, A,ia,ja,descA,                &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)
        character trans
        integer m,n,  ia,ja,descA(*) 
        integer ix,jx,descX(*),incx
        integer iy,jy,descY(*),incy
        complex*16 alpha,beta
        complex*16 A(*),X(*),Y(*)

        call  pzgemv(trans, m,n, alpha, A,ia,ja,descA,                         &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)

        return
        end subroutine pgemv_z111

        subroutine pgemv_z221(trans, m,n, alpha, A,ia,ja,descA,                &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)
        character trans
        integer m,n,  ia,ja,descA(*) 
        integer ix,jx,descX(*),incx
        integer iy,jy,descY(*),incy
        complex*16 alpha,beta
        complex*16 A(:,:),X(:,:),Y(*)

        call  pzgemv(trans, m,n, alpha, A,ia,ja,descA,                         &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)
        return
        end subroutine pgemv_z221

        subroutine pgemv_z211(trans, m,n, alpha, A,ia,ja,descA,                &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)
        character trans
        integer m,n,  ia,ja,descA(*) 
        integer ix,jx,descX(*),incx
        integer iy,jy,descY(*),incy
        complex*16 alpha,beta
        complex*16 A(:,:),X(*),Y(*)

        call  pzgemv(trans, m,n, alpha, A,ia,ja,descA,                         &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)
        return
        end subroutine pgemv_z211


        subroutine pgemv_c111(trans, m,n, alpha, A,ia,ja,descA,                &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)
        character trans
        integer m,n,  ia,ja,descA(*) 
        integer ix,jx,descX(*),incx
        integer iy,jy,descY(*),incy
        complex*8 alpha,beta
        complex*8 A(*),X(*),Y(*)

        call  pcgemv(trans, m,n, alpha, A,ia,ja,descA,                         &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)

        return
        end subroutine pgemv_c111

        subroutine pgemv_c221(trans, m,n, alpha, A,ia,ja,descA,                &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)
        character trans
        integer m,n,  ia,ja,descA(*) 
        integer ix,jx,descX(*),incx
        integer iy,jy,descY(*),incy
        complex*8 alpha,beta
        complex*8 A(:,:),X(:,:),Y(*)

        call  pcgemv(trans, m,n, alpha, A,ia,ja,descA,                         &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)
        return
        end subroutine pgemv_c221

        subroutine pgemv_c211(trans, m,n, alpha, A,ia,ja,descA,                &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)
        character trans
        integer m,n,  ia,ja,descA(*) 
        integer ix,jx,descX(*),incx
        integer iy,jy,descY(*),incy
        complex*8 alpha,beta
        complex*8 A(:,:),X(*),Y(*)

        call  pcgemv(trans, m,n, alpha, A,ia,ja,descA,                         &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)
        return
        end subroutine pgemv_c211



        subroutine pgemv_s111(trans, m,n, alpha, A,ia,ja,descA,                &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)
        character trans
        integer m,n,  ia,ja,descA(*) 
        integer ix,jx,descX(*),incx
        integer iy,jy,descY(*),incy
        real*4 alpha,beta
        real*4 A(*),X(*),Y(*)

        call  psgemv(trans, m,n, alpha, A,ia,ja,descA,                         &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)

        return
        end subroutine pgemv_s111

        subroutine pgemv_s221(trans, m,n, alpha, A,ia,ja,descA,                &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)
        character trans
        integer m,n,  ia,ja,descA(*) 
        integer ix,jx,descX(*),incx
        integer iy,jy,descY(*),incy
        real*4 alpha,beta
        real*4 A(:,:),X(:,:),Y(*)

        call  psgemv(trans, m,n, alpha, A,ia,ja,descA,                         &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)
        return
        end subroutine pgemv_s221

        subroutine pgemv_s211(trans, m,n, alpha, A,ia,ja,descA,                &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)
        character trans
        integer m,n,  ia,ja,descA(*) 
        integer ix,jx,descX(*),incx
        integer iy,jy,descY(*),incy
        real*4 alpha,beta
        real*4 A(:,:),X(*),Y(*)

        call  psgemv(trans, m,n, alpha, A,ia,ja,descA,                         &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)
        return
        end subroutine pgemv_s211

        subroutine pgemv_d111(trans, m,n, alpha, A,ia,ja,descA,                &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)
        character trans
        integer m,n,  ia,ja,descA(*) 
        integer ix,jx,descX(*),incx
        integer iy,jy,descY(*),incy
        real*8 alpha,beta
        real*8 A(*),X(*),Y(*)

        call  pdgemv(trans, m,n, alpha, A,ia,ja,descA,                         &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)

        return
        end subroutine pgemv_d111

        subroutine pgemv_d211(trans, m,n, alpha, A,ia,ja,descA,                &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)
        character trans
        integer m,n,  ia,ja,descA(*) 
        integer ix,jx,descX(*),incx
        integer iy,jy,descY(*),incy
        real*8 alpha,beta
        real*8 A(:,:),X(*),Y(*)

        call  pdgemv(trans, m,n, alpha, A,ia,ja,descA,                         &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)
        return
        end subroutine pgemv_d211

        subroutine pgemv_d221(trans, m,n, alpha, A,ia,ja,descA,                &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)
        character trans
        integer m,n,  ia,ja,descA(*) 
        integer ix,jx,descX(*),incx
        integer iy,jy,descY(*),incy
        real*8 alpha,beta
        real*8 A(:,:),X(:,:),Y(*)

        call  pdgemv(trans, m,n, alpha, A,ia,ja,descA,                         &
     &                                       X,ix,jx,descX,incx,               &
     &                                beta,  Y,iy,jy,descY,incy)
        return
        end subroutine pgemv_d221


        subroutine pcopy2d_gl_c(m,n,A,ia,ja,descA,B,ldB)
        implicit none
!
!       copy from global matrix A(ia:(ia+m-1),ja:(ja+n-1))
!       into local matrix B(1:m,1:n)
!
        integer m,n,ia,ja,descA(DLEN_),ldB
        complex*8 A(*), B(ldB,*)

        complex*8 alpha,beta
        integer descB(DLEN_),info

        call descinit( descB, m,n,m,n,0,0,descA(CTXT_),ldB,info)
        descB(RSRC_) = -1
        descB(CSRC_) = -1
        if (info.ne.0) then
           call BLACS_ABORT(0,info)
        end if

        alpha = 1.0d0
        beta = 0.0d0
        call pcgeadd('N',m,n,alpha,A,ia,ja,descA,beta,B,1,1,descB)
        return
        end subroutine pcopy2d_gl_c

        subroutine pcopy2d_gl_c2(m,n,A,ia,ja,descA,B,ldB)
        implicit none
!
!       copy from global matrix A(ia:(ia+m-1),ja:(ja+n-1))
!       into local matrix B(1:m,1:n)
!
        integer m,n,ia,ja,descA(DLEN_),ldB
        complex*8 A(:,:), B(ldB,*)

        complex*8 alpha,beta
        integer descB(DLEN_),info

        call descinit( descB, m,n,m,n,0,0,descA(CTXT_),ldB,info)
        if (info.ne.0) call BLACS_ABORT(0,info)
        descB(RSRC_) = -1
        descB(CSRC_) = -1

        alpha = 1.0d0
        beta = 0.0d0
        call pcgeadd('N',m,n,alpha,A,ia,ja,descA,beta,B,1,1,descB)
        return
        end subroutine pcopy2d_gl_c2


        subroutine pcopy2d_lg_c(m,n,B,ldB,A,ia,ja,descA)
        implicit none
!
!       copy from local matrix B(1:m,1:n) into 
!       global matrix A(ia:(ia+m-1),ja:(ja+n-1))
!
        integer m,n,ia,ja,descA(DLEN_),ldB
        complex*8 A(:,:), B(ldB,*)

        complex*8 alpha,beta
        integer descB(DLEN_),info

        call descinit( descB, m,n,m,n,0,0,descA(CTXT_),ldB,info)
        if (info.ne.0) call BLACS_ABORT(0,info)
        descB(RSRC_) = -1
        descB(CSRC_) = -1

        alpha = 1.0d0
        beta = 0.0d0
        call pcgeadd('N',m,n,alpha,B,1,1,descB,beta,A,ia,ja,descA)
        return
        end subroutine pcopy2d_lg_c

        
        subroutine pcopy2d_gl_z(m,n,A,ia,ja,descA,B,ldB)
        implicit none
!
!       copy from global matrix A(ia:(ia+m-1),ja:(ja+n-1))
!       into local matrix B(1:m,1:n)
!
        integer m,n,ia,ja,descA(DLEN_),ldB
        complex*16 A(*), B(ldB,*)

        complex*16 alpha,beta
        integer descB(DLEN_),info

        call descinit( descB, m,n,m,n,0,0,descA(CTXT_),ldB,info)
        if (info.ne.0) call BLACS_ABORT(0,info)
        descB(RSRC_) = -1
        descB(CSRC_) = -1

        alpha = 1.0d0
        beta = 0.0d0
        call pzgeadd('N',m,n,alpha,A,ia,ja,descA,beta,B,1,1,descB)
        return
        end subroutine pcopy2d_gl_z

        subroutine pcopy2d_gl_z2(m,n,A,ia,ja,descA,B,ldB)
        implicit none
!
!       copy from global matrix A(ia:(ia+m-1),ja:(ja+n-1))
!       into local matrix B(1:m,1:n)
!
        integer m,n,ia,ja,descA(DLEN_),ldB
        complex*16 A(:,:), B(ldB,*)

        complex*16 alpha,beta
        integer descB(DLEN_),info

        call descinit( descB, m,n,m,n,0,0,descA(CTXT_),ldB,info)
        if (info.ne.0) call BLACS_ABORT(0,info)
        descB(RSRC_) = -1
        descB(CSRC_) = -1

        alpha = 1.0d0
        beta = 0.0d0
        call pzgeadd('N',m,n,alpha,A,ia,ja,descA,beta,B,1,1,descB)
        return
        end subroutine pcopy2d_gl_z2


        subroutine pcopy2d_lg_z(m,n,B,ldB,A,ia,ja,descA)
        implicit none
!
!       copy from local matrix B(1:m,1:n) into 
!       global matrix A(ia:(ia+m-1),ja:(ja+n-1))
!
        integer m,n,ia,ja,descA(DLEN_),ldB
        complex*16 A(:,:), B(ldB,*)

        complex*16 alpha,beta
        integer descB(DLEN_),info

        call descinit( descB, m,n,m,n,0,0,descA(CTXT_),ldB,info)
        if (info.ne.0) call BLACS_ABORT(0,info)
        descB(RSRC_) = -1
        descB(CSRC_) = -1

        alpha = 1.0d0
        beta = 0.0d0
        call pzgeadd('N',m,n,alpha,B,1,1,descB,beta,A,ia,ja,descA)
        return
        end subroutine pcopy2d_lg_z

        
        subroutine pcopy2d_gl_s(m,n,A,ia,ja,descA,B,ldB)
        implicit none
!
!       copy from global matrix A(ia:(ia+m-1),ja:(ja+n-1))
!       into local matrix B(1:m,1:n)
!
        integer m,n,ia,ja,descA(DLEN_),ldB
        real*4 A(*), B(ldB,*)

        real*4 alpha,beta
        integer descB(DLEN_),info

        call descinit( descB, m,n,m,n,0,0,descA(CTXT_),ldB,info)
        if (info.ne.0) call BLACS_ABORT(0,info)
        descB(RSRC_) = -1
        descB(CSRC_) = -1

        alpha = 1.0d0
        beta = 0.0d0
        call psgeadd('N',m,n,alpha,A,ia,ja,descA,beta,B,1,1,descB)
        return
        end subroutine pcopy2d_gl_s

        subroutine pcopy2d_gl_s2(m,n,A,ia,ja,descA,B,ldB)
        implicit none
!
!       copy from global matrix A(ia:(ia+m-1),ja:(ja+n-1))
!       into local matrix B(1:m,1:n)
!
        integer m,n,ia,ja,descA(DLEN_),ldB
        real*4 A(:,:), B(ldB,*)

        real*4 alpha,beta
        integer descB(DLEN_),info

        call descinit( descB, m,n,m,n,0,0,descA(CTXT_),ldB,info)
        if (info.ne.0) call BLACS_ABORT(0,info)
        descB(RSRC_) = -1
        descB(CSRC_) = -1

        alpha = 1.0d0
        beta = 0.0d0
        call psgeadd('N',m,n,alpha,A,ia,ja,descA,beta,B,1,1,descB)
        return
        end subroutine pcopy2d_gl_s2


        subroutine pcopy2d_lg_s(m,n,B,ldB,A,ia,ja,descA)
        implicit none
!
!       copy from local matrix B(1:m,1:n) into 
!       global matrix A(ia:(ia+m-1),ja:(ja+n-1))
!
        integer m,n,ia,ja,descA(DLEN_),ldB
        real*4 A(:,:), B(ldB,*)

        real*4 alpha,beta
        integer descB(DLEN_),info

        call descinit( descB, m,n,m,n,0,0,descA(CTXT_),ldB,info)
        if (info.ne.0) call BLACS_ABORT(0,info)
        descB(RSRC_) = -1
        descB(CSRC_) = -1

        alpha = 1.0d0
        beta = 0.0d0
        call psgeadd('N',m,n,alpha,B,1,1,descB,beta,A,ia,ja,descA)
        return
        end subroutine pcopy2d_lg_s

        
        subroutine pcopy2d_gl_d(m,n,A,ia,ja,descA,B,ldB)
        implicit none
!
!       copy from global matrix A(ia:(ia+m-1),ja:(ja+n-1))
!       into local matrix B(1:m,1:n)
!
        integer m,n,ia,ja,descA(DLEN_),ldB
        real*8 A(*), B(ldB,*)

        real*8 alpha,beta
        integer descB(DLEN_),info

        call descinit( descB, m,n,m,n,0,0,descA(CTXT_),ldB,info)
        if (info.ne.0) call BLACS_ABORT(0,info)
        descB(RSRC_) = -1
        descB(CSRC_) = -1

        alpha = 1.0d0
        beta = 0.0d0
        call pdgeadd('N',m,n,alpha,A,ia,ja,descA,beta,B,1,1,descB)
        return
        end subroutine pcopy2d_gl_d

        subroutine pcopy2d_gl_d2(m,n,A,ia,ja,descA,B,ldB)
        implicit none
!
!       copy from global matrix A(ia:(ia+m-1),ja:(ja+n-1))
!       into local matrix B(1:m,1:n)
!
        integer m,n,ia,ja,descA(DLEN_),ldB
        real*8 A(:,:), B(ldB,*)

        real*8 alpha,beta
        integer descB(DLEN_),info

        call descinit( descB, m,n,m,n,0,0,descA(CTXT_),ldB,info)
        if (info.ne.0) call BLACS_ABORT(0,info)
        descB(RSRC_) = -1
        descB(CSRC_) = -1

        alpha = 1.0d0
        beta = 0.0d0
        call pdgeadd('N',m,n,alpha,A,ia,ja,descA,beta,B,1,1,descB)
        return
        end subroutine pcopy2d_gl_d2


        subroutine pcopy2d_lg_d(m,n,B,ldB,A,ia,ja,descA)
        implicit none
!
!       copy from local matrix B(1:m,1:n) into 
!       global matrix A(ia:(ia+m-1),ja:(ja+n-1))
!
        integer m,n,ia,ja,descA(DLEN_),ldB
        real*8 A(:,:), B(ldB,*)

        real*8 alpha,beta
        integer descB(DLEN_),info

        call descinit( descB, m,n,m,n,0,0,descA(CTXT_),ldB,info)
        if (info.ne.0) call BLACS_ABORT(0,info)
        descB(RSRC_) = -1
        descB(CSRC_) = -1

        alpha = 1.0d0
        beta = 0.0d0
        call pdgeadd('N',m,n,alpha,B,1,1,descB,beta,A,ia,ja,descA)
        return
        end subroutine pcopy2d_lg_d

        

        end module scalapack_mod
