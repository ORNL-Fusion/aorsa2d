module Performance

implicit none

type, public :: RunPerfData 

        integer :: nProcs
        integer :: nSpatialPoints
        integer :: nRowLocal
        integer :: nColLocal
        integer :: nRowGlobal
        integer :: nColGlobal
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
