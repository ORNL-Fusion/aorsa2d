program test_sigma

use constants
use sigma_mod
use aorsa2din_mod, &
only: nPhi, damping, delta0, lMax, nzfun, xnuomg, upshift, &
    read_nameList, xkPerp_cutOff
use netcdf
use check_mod

implicit none

integer :: nPts, i, j, nK, s
real :: rLeft, rRigh, dr, mSpec(2), qSpec(2), b0
real, allocatable :: r(:), b(:), omgc(:)
real :: kx, ky, kVec(3)
real :: k_cutOff, freq, omgrf, k0, density, omgp2
real :: ktSpec
complex :: sigma_tmp(3,3)
real :: U_cyl(3,3), U_xyz(3,3)
integer :: zSpec(2)
character(len=100) :: fileName, outFName
complex, allocatable, dimension(:,:) :: sig11, sig12, sig13, &
    sig21, sig22, sig23, &
    sig31, sig32, sig33
integer :: nc_id, n_id, nK_id, &
    r_id, k_id, &
    sig11_re_id, sig11_im_id,&
    sig12_re_id, sig12_im_id,&
    sig13_re_id, sig13_im_id,&
    sig21_re_id, sig21_im_id,&
    sig22_re_id, sig22_im_id,&
    sig23_re_id, sig23_im_id,&
    sig31_re_id, sig31_im_id,&
    sig32_re_id, sig32_im_id,&
    sig33_re_id, sig33_im_id

real, allocatable :: kArr(:)
real :: kPhi

fileName = 'test_sigma.in'

call read_nameList ( fileName = fileName )

nPts    = 2
rLeft   = 0.4
rRigh   = 1.6
dr      = (rRigh-rLeft)/(nPts-1)

allocate ( r(nPts),b(nPts),omgc(nPts) )

r       = [ (i*dr,i=1,nPts) ] + rLeft

b0  = 0.55
b   = b0 / r
mSpec   = (/ xme, xmh /)
ktSpec  = 0.003e3 * q
zSpec   = (/-1, 1/)
qSpec   = zSpec * q
density = 4e19
kx  = 20.0
ky  = 2.0
k_cutOff =  150.0
freq   = 30e6
omgrf = 2.0 * pi * freq
k0 = omgrf / clight
U_xyz   = real( reshape ( (/1,0,0,0,1,0,0,0,1/) , (/3,3/) ) )
U_cyl   = real( reshape ( (/1,0,0,0,1,0,0,0,1/) , (/3,3/) ) )
xNuOmg = 0.00

nK = 1

allocate (  sig11(nPts,nK),sig21(nPts,nK),sig31(nPts,nK), &
            sig12(nPts,nK),sig22(nPts,nK),sig32(nPts,nK), & 
            sig13(nPts,nK),sig23(nPts,nK),sig33(nPts,nK), &
            kArr(nK) )

sig11 = 0 
sig12 = 0 
sig13 = 0

sig21 = 0
sig22 = 0
sig23 = 0

sig31 = 0
sig32 = 0
sig33 = 0

kArr    = 0

do s=1,2

omgc    = qSpec(s) * b / mSpec(s)
omgp2   = density * qSpec(s)**2 / ( eps0 * mSpec(s) )
!write(*,*) omgc
!write(*,*) omgp2

do i=1,nPts

    do j=1,nK

    nPhi = j
    kPhi = nPhi/r(i)

    kArr(j) = j
    kVec = (/ kx, kPhi, ky /)

    sigma_tmp = sigmaHot_maxwellian&
        ( mSpec(s), &
        ktSpec, omgc(i), omgp2, &
        kVec, r(i), &
        omgrf, k0, &
        k_cutoff, s, &
        0.0, 0.0, b(i), 0.0, &
        xnuOmg )

    !sigma_tmp = sigmaCold_stix &
    !    ( omgc(i), omgp2, omgrf, xnuOmg )

        sig11(i,j) = sig11(i,j) + sigma_tmp(1,1)
        sig12(i,j) = sig12(i,j) + sigma_tmp(1,2)
        sig13(i,j) = sig13(i,j) + sigma_tmp(1,3)

        sig21(i,j) = sig21(i,j) + sigma_tmp(2,1)
        sig22(i,j) = sig22(i,j) + sigma_tmp(2,2)
        sig23(i,j) = sig23(i,j) + sigma_tmp(2,3)

        sig31(i,j) = sig31(i,j) + sigma_tmp(3,1)
        sig32(i,j) = sig32(i,j) + sigma_tmp(3,2)
        sig33(i,j) = sig33(i,j) + sigma_tmp(3,3)

    enddo
enddo                              

write(*,*) sig11(1,1),sig12(1,1),sig13(1,1)
write(*,*) sig21(1,1),sig22(1,1),sig23(1,1)
write(*,*) sig31(1,1),sig32(1,1),sig33(1,1)


enddo

! write value to file for plotting in IDL
! ---------------------------------------


    outFName = 'test_sigma.nc'
    call check ( nf90_create ( outFName, nf90_clobber, nc_id ) )
    call check ( nf90_def_dim ( nc_id, "nPts", nPts, n_id ) )
    call check ( nf90_def_dim ( nc_id, "nK", nK, nK_id ) )

    call check ( nf90_def_var ( nc_id, "sig31_re", nf90_real, &
        (/n_id,nK_id/), sig31_re_id ) ) 
    call check ( nf90_def_var ( nc_id, "sig31_im", nf90_real, &
        (/n_id,nK_id/), sig31_im_id ) ) 
 
    call check ( nf90_def_var ( nc_id, "sig32_re", nf90_real, &
        (/n_id,nK_id/), sig32_re_id ) ) 
    call check ( nf90_def_var ( nc_id, "sig32_im", nf90_real, &
        (/n_id,nK_id/), sig32_im_id ) ) 
 
    call check ( nf90_def_var ( nc_id, "sig33_re", nf90_real, &
        (/n_id,nK_id/), sig33_re_id ) ) 
    call check ( nf90_def_var ( nc_id, "sig33_im", nf90_real, &
        (/n_id,nK_id/), sig33_im_id ) ) 

    call check ( nf90_def_var ( nc_id, "sig21_re", nf90_real, &
        (/n_id,nK_id/), sig21_re_id ) ) 
    call check ( nf90_def_var ( nc_id, "sig21_im", nf90_real, &
        (/n_id,nK_id/), sig21_im_id ) ) 
 
    call check ( nf90_def_var ( nc_id, "sig22_re", nf90_real, &
        (/n_id,nK_id/), sig22_re_id ) ) 
    call check ( nf90_def_var ( nc_id, "sig22_im", nf90_real, &
        (/n_id,nK_id/), sig22_im_id ) ) 
 
    call check ( nf90_def_var ( nc_id, "sig23_re", nf90_real, &
        (/n_id,nK_id/), sig23_re_id ) ) 
    call check ( nf90_def_var ( nc_id, "sig23_im", nf90_real, &
        (/n_id,nK_id/), sig23_im_id ) ) 

    call check ( nf90_def_var ( nc_id, "sig11_re", nf90_real, &
        (/n_id,nK_id/), sig11_re_id ) ) 
    call check ( nf90_def_var ( nc_id, "sig11_im", nf90_real, &
        (/n_id,nK_id/), sig11_im_id ) ) 
 
    call check ( nf90_def_var ( nc_id, "sig12_re", nf90_real, &
        (/n_id,nK_id/), sig12_re_id ) ) 
    call check ( nf90_def_var ( nc_id, "sig12_im", nf90_real, &
        (/n_id,nK_id/), sig12_im_id ) ) 
 
    call check ( nf90_def_var ( nc_id, "sig13_re", nf90_real, &
        (/n_id,nK_id/), sig13_re_id ) ) 
    call check ( nf90_def_var ( nc_id, "sig13_im", nf90_real, &
        (/n_id,nK_id/), sig13_im_id ) ) 


    call check ( nf90_def_var ( nc_id, "r", nf90_real, &
        (/n_id/), r_id ) ) 
    call check ( nf90_def_var ( nc_id, "k", nf90_real, &
        (/nk_id/), k_id ) ) 
 
    call check ( nf90_enddef ( nc_id ) )

    call check ( nf90_put_var ( nc_id, sig11_re_id, real  ( sig11 ) ) )
    call check ( nf90_put_var ( nc_id, sig11_im_id, aimag ( sig11 ) ) )
    call check ( nf90_put_var ( nc_id, sig12_re_id, real  ( sig12 ) ) )
    call check ( nf90_put_var ( nc_id, sig12_im_id, aimag ( sig12 ) ) )
    call check ( nf90_put_var ( nc_id, sig13_re_id, real  ( sig13 ) ) )
    call check ( nf90_put_var ( nc_id, sig13_im_id, aimag ( sig13 ) ) )

    call check ( nf90_put_var ( nc_id, sig21_re_id, real  ( sig21 ) ) )
    call check ( nf90_put_var ( nc_id, sig21_im_id, aimag ( sig21 ) ) )
    call check ( nf90_put_var ( nc_id, sig22_re_id, real  ( sig22 ) ) )
    call check ( nf90_put_var ( nc_id, sig22_im_id, aimag ( sig22 ) ) )
    call check ( nf90_put_var ( nc_id, sig23_re_id, real  ( sig23 ) ) )
    call check ( nf90_put_var ( nc_id, sig23_im_id, aimag ( sig23 ) ) )

    call check ( nf90_put_var ( nc_id, sig31_re_id, real  ( sig31 ) ) )
    call check ( nf90_put_var ( nc_id, sig31_im_id, aimag ( sig31 ) ) )
    call check ( nf90_put_var ( nc_id, sig32_re_id, real  ( sig32 ) ) )
    call check ( nf90_put_var ( nc_id, sig32_im_id, aimag ( sig32 ) ) )
    call check ( nf90_put_var ( nc_id, sig33_re_id, real  ( sig33 ) ) )
    call check ( nf90_put_var ( nc_id, sig33_im_id, aimag ( sig33 ) ) )

    call check ( nf90_put_var ( nc_id, r_id, r ) )
    call check ( nf90_put_var ( nc_id, k_id, kArr ) )

    call check ( nf90_close ( nc_id ) )

end program test_sigma
