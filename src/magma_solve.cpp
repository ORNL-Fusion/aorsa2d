#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cuda.h"
#include "cuda_runtime_api.h"
#include "cublas.h"
#include "magma.h"

extern "C"
int magma_solve ( int *dA_dim, int *lWork, double2 *A, int *ipiv, int *N ){

	// Check inputs
	//
	fprintf (stderr, "Using MAGMA solve\n" );
	fprintf (stderr, "	dA_dim: %i\n", *dA_dim );
	fprintf (stderr, "	N: %i\n", *N );
	fprintf (stderr, "	lWork: %i\n", *lWork );

	cuInit(0);
	cublasInit();
	printout_devices();

	cublasStatus status;

	double2 *d_A, *work;
	status = cublasAlloc ( *dA_dim, sizeof(double2), (void**)&d_A );

	if ( status != CUBLAS_STATUS_SUCCESS ){
			fprintf (stderr, "ERROR: device memory allocation error (d_A)\n" );
			fprintf (stderr, "ERROR: dA_dim: %i\n", dA_dim );
	}

	cudaError_t err;
	err = cudaMallocHost ( (void**)&work, *lWork * sizeof(double2) );

	if(err != cudaSuccess){
		fprintf (stderr, "ERROR: cudaMallocHost error (work)\n" );
	}

	int info[1];
	TimeStruct start, end;

	start = get_current_time ();
	magma_zgetrf ( N, N, A, N, ipiv, work, d_A, info );
	end = get_current_time ();

	double gpu_perf;
	gpu_perf = 4.*2.*(*N)*(*N)*(*N)/(3.*1000000*GetTimerValue(start,end));

	if ( info[0] != 0 ){
			fprintf (stderr, "ERROR: magma_zgetrf failed\n" );
	}

	printf("	GPU performance: %6.2f GFlop/s\n", gpu_perf);

	int stat = 0;
	return stat;

}
