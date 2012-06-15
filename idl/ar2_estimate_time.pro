pro ar2_estimate_time

nR = 257d0
nZ = 257d0 

npRow = 32d0
npCol = 64d0

nSpec = 2
lMax = 2 

;; Hopper parameters
;; -----------------
;print, 'HOPPER'
;CoresPerNode = 24
;MemPerNode_GB = 32.0
;CPU_GFlops = 5.0

; Jaguarpf parameters
; -------------------
print, 'JAGUARPF'
CoresPerNode = 16 
MemPerNode_GB = 32.0
CPU_GFlops = 5.5
CPU_PGESVR_GFlops = 9.8
GPU_GFlops = 170.0
GPU_PGESVR_GFlops = 

; Memory required

MemPerCPUAvail_GB = MemPerNode_GB/CoresPerNode*0.8
N = 3d0 * nR * nZ
nNodes = npRow*npCol/CoresPerNode

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

Time_s_CPU = NFOperations / ( CPU_GFlops * 1024d0^3 * TotalNProcs )
Time_m_CPU = Time_s / 60d0
Time_h_CPU = Time_m / 60d0

Time_s_CPU_PGESVR = NFOperations / ( CPU_PGESVR_GFlops * 1024d0^3 * TotalNProcs )
Time_m_CPU_PGESVR = Time_s / 60d0
Time_h_CPU_PGESVR = Time_m / 60d0

Time_s_GPU = NFOperations / ( GPU_GFlops * 1024d0^3 * nNodes )
Time_m_GPU = Time_s_GPU / 60d0
Time_h_GPU = Time_m_GPU / 60d0

Time_s_GPU_PGESVR = NFOperations / ( GPU_PGESVR_GFlops * 1024d0^3 * nNodes )
Time_m_GPU_PGESVR = Time_s_GPU / 60d0
Time_h_GPU_PGESVR = Time_m_GPU / 60d0


print, '-----------------------------------------------'
print, 'For a ',string(nR,format='(i3.3)'),'x',string(nZ,format='(i3.3)'),$
		' grid with ',string(npRow,format='(i3.3)'),'x',string(npCol,format='(i3.3)'),' procs: '

print, 'Estimates solve time (CPU): ', string(Time_s_CPU,format='(i)'), ' seconds'
print, 'Estimates solve time (CPU): ', string(Time_m_CPU,format='(i)'), ' minutes'
print, 'Estimates solve time (CPU): ', Time_h_CPU, ' hours'

print, 'Estimates solve time (GPU): ', string(Time_s_GPU,format='(i)'), ' seconds'
print, 'Estimates solve time (GPU): ', string(Time_m_GPU,format='(i)'), ' minutes'
print, 'Estimates solve time (GPU): ', Time_h_GPU, ' hours'

print, 'Estimates solve time (CPU PGESVR): ', string(Time_s_CPU_PGESVR,format='(i)'), ' seconds'
print, 'Estimates solve time (CPU PGESVR): ', string(Time_m_CPU_PGESVR,format='(i)'), ' minutes'
print, 'Estimates solve time (CPU PGESVR): ', Time_h_CPU_PGESVR, ' hours'

print, 'Estimates solve time (GPU PGESVR): ', string(Time_s_GPU_PGESVR,format='(i)'), ' seconds'
print, 'Estimates solve time (GPU PGESVR): ', string(Time_m_GPU_PGESVR,format='(i)'), ' minutes'
print, 'Estimates solve time (GPU PGESVR): ', Time_h_GPU_PGESVR, ' hours'


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
print, 'Estimated CPU hours: ', TotalTime_h*TotalNProcs
print, 'nNodes: ', nNodes

stop

end
