module hammett

interface Zfun
    module procedure Zfun_zd
    module procedure Zfun_rd
end interface

contains

complex function Zfun_zd(zeta)

  use constants
  use z_erf

  implicit none

  complex(kind=dbl), intent(in) :: zeta
  logical :: flag
  real(kind=dbl), parameter :: sqrt_pi = 1.772453850905516
  complex(kind=dbl) :: y
  
  call wofz_f90(zeta, y, flag)

  if(flag)then 
    write(*,*) 'CRAP: wofz_f90 failed'
    stop
  endif
 
  Zfun_zd = sqrt_pi*zi*y

end function Zfun_zd

complex function Zfun_rd(zeta)

  use constants
  use z_erf

  implicit none

  real(kind=dbl), intent(in) :: zeta
  logical :: flag
  real(kind=dbl), parameter :: sqrt_pi = 1.772453850905516
  complex(kind=dbl) :: y
 
  call wofz_f90(cmplx(zeta,0,dbl),y, flag)

  if(flag)then 
    write(*,*) 'CRAP: wofz_f90 failed'
    stop
  endif

  Zfun_rd = sqrt_pi*zi*y
  
end function Zfun_rd

end module hammett
