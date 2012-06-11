module Performance

use constants

implicit none

type, public :: RunPerfData 

        integer(kind=long) :: nProcs
        integer(kind=long) :: nSpatialPoints
        integer(kind=long) :: nRowLocal
        integer(kind=long) :: nColLocal
        integer(kind=long) :: nRowGlobal
        integer(kind=long) :: nColGlobal
        real :: MatSizeLocal_GB
        real :: MatSizeGlobal_GB
        real :: TimeWorkList
        real :: TimeFill
        real :: TimeSolve
        real :: TimeCurrent
        real :: TimeTotal
        real :: GflopsFillLocal
        real :: GflopsFillGlobal
        real :: GflopsSolveLocal
        real :: GflopsSolveGlobal
        real :: GflopsCurrentLocal
        real :: GflopsCurrentGlobal

end type RunPerfData

type, public :: MemoryUsage

        integer :: nAllocations
        real :: MBytes

end type MemoryUsage

type(MemoryUsage) :: Mem
type(RunPerfData) :: Perf

end module Performance
