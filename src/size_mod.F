      module size_mod
      implicit none

      integer nmodesmax, mmodesmax
	
!     ---------------------------------------------------
!     Maximum number of modes allowed:
!     IMPORTANT NOTE: when changing these parameters, you 
!     must touch all files and recomile everything!!!
!
!	DLG: Altering a parameter (which is set at compile time)
!	requires a "make clean" and then "make". 
!     ---------------------------------------------------
#if DLGHP 
      parameter (nmodesmax = 256)
      parameter (mmodesmax = 256)
#else
      parameter (nmodesmax = 256)
      parameter (mmodesmax = 256)
#endif

      end module size_mod

