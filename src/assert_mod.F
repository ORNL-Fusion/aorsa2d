        module assert_mod

        implicit none
        private
        public :: assert

        contains

        subroutine assert( lcond, mesg, ivalue )
        implicit none
        logical, intent(in) :: lcond
        character*(*),intent(in) :: mesg
        integer, intent(in) :: ivalue
        integer, parameter :: stderr=0

        if (.not.lcond) then
          write(stderr,*) mesg,ivalue
          stop '** assertion error ** '
        endif
        return
        end subroutine assert

        end module assert_mod
