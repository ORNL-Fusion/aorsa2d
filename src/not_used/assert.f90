module assert_mod

contains

      subroutine assert( lcond, msg, ivalue )
      logical lcond
      character*(*) msg
      integer ivalue

      if (.not.lcond) then
          write(*,*) msg,ivalue
          stop '** error in assertion ** '
      endif

      return
      end subroutine assert

end module assert_mod
