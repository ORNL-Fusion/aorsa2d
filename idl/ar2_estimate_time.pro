pro ar2_estimate_time

nR = 257 
nZ = 257

npRow = 48 
npCol = 48 

; Hopper parameters
MemPerCPU_GB = 32.0/24.0*0.8
Scalapack_GFlops = 6.0

; Memory required

N = 3.0 * nR * nZ

DoubleComplexBytes = 16.0

Mem_B = N^2 * DoubleComplexBytes
Mem_KB = Mem_B/1024d0
Mem_MB = Mem_B/(1024d0^2)
Mem_GB = Mem_B/(1024d0^3)
Mem_TB = Mem_B/(1024d0^4)

print, 'Total Mem Required: [Bytes]', Mem_B
print, 'Total Mem Required: [MBytes]', Mem_MB
print, 'Total Mem Required: [GBytes]', Mem_GB
print, 'Total Mem Required: [TBytes]', Mem_TB

MinNProcs = Mem_GB / MemPerCPU_GB

print, 'Minimum number of procs: ', fix(MinNProcs)


TotalNProcs = npRow*npCol*1.0

print, 'Actual number of procs: ', TotalNProcs

if(MinNProcs gt TotalNProcs ) then begin
		print, 'PROBLEM: Not enough memory, use more procs.'
		stop
endif

print, 'Mem Per CPU [MBytes]: ', Mem_MB / TotalNProcs

NFOperations = 2.67 * N^3

print, 'Number of Scalapack operations: ', NFOperations

Time_s = NFOperations / ( Scalapack_GFlops * 1024d0^3 * TotalNProcs )
Time_m = Time_s / 60d0
Time_h = Time_m / 60d0

print, '-----------------------------------------------'
print, 'For a ',string(nR,format='(i3.3)'),'x',string(nZ,format='(i3.3)'),$
		' grid with ',string(npRow,format='(i3.3)'),'x',string(npCol,format='(i3.3)'),' procs: '
print, 'Estimates solve time: ', fix(Time_s), ' seconds'
print, 'Estimates solve time: ', fix(Time_m), ' minutes'
print, 'Estimates solve time: ', fix(Time_h), ' hours'

stop

end
