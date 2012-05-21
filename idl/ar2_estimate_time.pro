pro ar2_estimate_time

nR = 513d0
nZ = 257d0 

npRow = 96d0
npCol = 48d0

; Hopper parameters
MemPerCPUAvail_GB = 32.0/24.0*0.8
Scalapack_GFlops = 5.0

; Memory required

N = 3d0 * nR * nZ

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

MinNProcs = Mem_GB / MemPerCPUAvail_GB

print, 'Minimum number of procs: ', string(MinNProcs,format='(i)')


TotalNProcs = npRow*npCol*1.0

print, 'Actual number of procs: ', TotalNProcs

if(MinNProcs gt TotalNProcs ) then begin
		print, 'PROBLEM: Not enough memory, use more procs.'
		stop
endif

MemPerCPU_MB = Mem_MB / TotalNProcs

print, 'Mem Per CPU [MBytes]: ', MemPerCPU_MB
print, 'Mem Per CPU Avail [MBytes]: ', MemPerCPUAvail_GB*1024.0

NFOperations = 2.67d0 * N^3

print, 'Number of Scalapack operations: ', NFOperations

Time_s = NFOperations / ( Scalapack_GFlops * 1024d0^3 * TotalNProcs )
Time_m = Time_s / 60d0
Time_h = Time_m / 60d0

print, '-----------------------------------------------'
print, 'For a ',string(nR,format='(i3.3)'),'x',string(nZ,format='(i3.3)'),$
		' grid with ',string(npRow,format='(i3.3)'),'x',string(npCol,format='(i3.3)'),' procs: '
print, 'Estimates solve time: ', string(Time_s,format='(i)'), ' seconds'
print, 'Estimates solve time: ', string(Time_m,format='(i)'), ' minutes'
print, 'Estimates solve time: ', Time_h, ' hours'

FillTime_s =  MemPerCPU_MB*0.1*1.3
CurrentTime_s = MemPerCPU_MB*0.1*1.15
WorkListTime_s = MemPerCPU_MB*0.1*1.4 
TotalTime_s = WorkListTime_s+$
		FillTime_s+$
		CurrentTime_s+$
		WorkListTime_s+$
		Time_s
TotalTime_m = TotalTime_s / 60.0
TotalTime_h = TotalTime_m / 60.0

print, 'Estimated fill time: ', FillTime_s, ' seconds'
print, 'Estimated current time: ', CurrentTime_s, ' seconds'
print, 'Estimated workList time: ', WorkListTime_s, ' seconds'
print, 'Estimated total time: ', TotalTime_s, ' seconds' 
print, 'Estimated total time: ', TotalTime_m, ' minutes' 
print, 'Estimated total time: ', TotalTime_h, ' hours' 


stop

end
