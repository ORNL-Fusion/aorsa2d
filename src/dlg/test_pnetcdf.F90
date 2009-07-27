module checkf
    
contains
    subroutine check ( status )
#       include <pnetcdf.inc>
        integer, intent ( in) :: status
     
        if(status /= nf_noerr) then 
            print *, trim(nfmpi_strerror(status))
            stop "Stopped"
        end if

    end subroutine check 
end module checkf

!   Compilation :
!   /home/dg6/code/openmpi/bin/mpif90 test_pnetcdf.F90 -I
!   /home/dg6/code/pNetCdf/include/ /home/dg6/code/pNetCdf/lib/libpnetcdf.a
!   -fdefault-real-8
!
!   The -fdefault-real-8 is what AORSA uses so I have to get it working with 
!   that. This is why we have the double nf type and call.
 
program test_pnetcdf
    use checkf
    use mpi
    
#   include <pnetcdf.inc>

    integer :: pnc_id, pnR_id, pnz_id, &
        pnuper_id, pnupar_id, pql_b_id, mpi_iErr
    character(len=100) :: pncFileName
    integer(KIND=MPI_OFFSET_KIND) :: nnodex, nnodey, nuper, nupar
    integer(KIND=MPI_OFFSET_KIND) :: start(4), cnt(4)
    real, allocatable :: bql_store (:,:,:)

    nnodex  = 4 
    nnodey  = 1 
    nuper   = 2 
    nupar   = 2 

    pncFileName = 'test_pncdf.nc'

    call mpi_init ( mpi_iErr )
    call mpi_comm_size ( MPI_COMM_WORLD, mpi_nProcs, mpi_iErr )
    call mpi_comm_rank ( MPI_COMM_WORLD, mpi_pId, mpi_iErr )

    allocate ( bql_store (nnodex,nuper,nupar) )

    bql_store   = 0.0
    bql_store(:,1,1) = 10.0

    call check ( nfmpi_create ( MPI_COMM_WORLD, pncFileName, &
        NF_CLOBBER, MPI_INFO_NULL, pnc_id ) )
    
    call check ( nfmpi_def_dim ( pnc_id, "nR", nnodex, pnR_id ) )
    call check ( nfmpi_def_dim ( pnc_id, "nz", nnodey, pnz_id ) )
    call check ( nfmpi_def_dim ( pnc_id, "nuper", nuper, pnuper_id ) )
    call check ( nfmpi_def_dim ( pnc_id, "nupar", nupar, pnupar_id ) )

    call check ( nfmpi_def_var ( pnc_id, "ql_b", NF_DOBULE, 4, &
        (/ pnR_id, pnz_id, pnuper_id, pnupar_id /), pql_b_id ) )

    call check (nfmpi_enddef ( pnc_id ) )

    start(1)    = nnodex / mpi_nProcs * mpi_pId + 1
    start(2)    = 1
    start(3)    = 1
    start(4)    = 1

    cnt(1)  = nnodex / mpi_nProcs
    cnt(2)  = nnodey
    cnt(3)  = nuper
    cnt(4)  = nupar

    write(*,'(i1.1,4x,2(i2.2,1x,i2.2,1x,i2.2,1x,i2.2,3x))') mpi_pId, start, cnt

    call check ( &
        nfmpi_put_vara_double_all ( pnc_id, pql_b_id, &
        start, cnt, &
        !bql_store( start(1):start(1)+cnt(1)-1,start(2):start(2)+cnt(2)-1, &
        !     start(3):start(3)+cnt(3)-1,start(4):start(4)+cnt(4)-1 ) ) )
        bql_store( start(1):start(1)+cnt(1)-1,:,:) ) )


    call check ( nfmpi_close ( pnc_id ) )
    call mpi_finalize ( mpi_iErr )

end program test_pnetcdf

