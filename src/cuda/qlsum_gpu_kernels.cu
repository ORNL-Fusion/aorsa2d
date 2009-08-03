#include <stdio.h>
#include "cuda.h"

#include "ccmplx.h"
#include "zcmplx.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

#define INT int

//These are the types used for computation on the GPU
#define REAL float
#define CMPLX ccmplx

//These are the types used for referencing host memory
#define REAL_H double
#define CMPLX_H zcmplx

//This includes must follow the type declaration. I apologize, but this is the easiest way
#include "qlsum_gpu_kernels.cuh"

#define NUMBLOCKS 60

//Dimension constants
INT nkx1, nkx2,
  nky1, nky2,
  nzeta,
  nuper, nupar,
  nkdim1, nkdim2,
  mkdim1, mkdim2,
  lmaxdim;

//Device pointers
CMPLX
//Accumulator temporaries
*sum_p,
  *sum2_p, *sumkx2_p, *sumky2_p,
  *sumwdot_p, *sumwdotkx_p, *sumwdotky_p,
  //Various other temporaries
  *eps_t_p, *sumb_11_nm_t_p, *sumb_31_nm_t_p,
  *zbeta_p, *zbeta_iharm_p,
  //Output
  *sum_wdot_p, *sum_fx0_p, *sum_fy0_p,
  *b_sum_p, *c_sum_p, *e_sum_p, *f_sum_p,
  //Read only
  *xx_p, *yy_p,
  *ealphak_p, *ebetak_p, *ebk_p;

//All these pointers are read only
REAL *uper_p, *upar_p,
  *xkxsav_p, *xkysav_p,
  *dfduper_p, *dfdupar_p,
  *npara_sav_p,
  *xkperpn_tmp_p, *zetai_p, *Jni_p,
  *zetamin_tmp_p, *dzetai_tmp_p;

INT *nres_p, *mres_p;

//Single values
REAL sqmut0, dui, wcw;

//Various useful indexing constants
#define F 0
#define E 2
#define C 4
#define B 6

#define BUFF_SIZE 2048

float float_tmp[BUFF_SIZE];

void cudaMemcpyD2SH2D(void *targ_, void *source_, unsigned int n) {
  int i = 0, j = 0, blocks = n / BUFF_SIZE, remainder = n % BUFF_SIZE;
  float *targ;
  double *source;

  targ = (float *)targ_;
  source = (double *)source_;

  for(i = 0; i < blocks; i++) {
    for(j = 0; j < BUFF_SIZE; j++) {
      float_tmp[j] = (float)source[i * BUFF_SIZE + j];
    }

    cudaMemcpy(&targ[i * BUFF_SIZE], &float_tmp[0], sizeof(float) * BUFF_SIZE, cudaMemcpyHostToDevice);
  }

  if(remainder > 0) {
    for(j = 0; j < remainder; j++) {
      float_tmp[j] = (float)source[blocks * BUFF_SIZE + j];
    }

    cudaMemcpy(&targ[blocks * BUFF_SIZE], &float_tmp[0], sizeof(float) * remainder, cudaMemcpyHostToDevice);
  }
}

void cudaMemcpyS2DD2H(void *targ_, void *source_, unsigned int n) {
  int i = 0, j = 0, blocks = n / BUFF_SIZE, remainder = n % BUFF_SIZE;
  double *targ;
  float *source;

  targ = (double *)targ_;
  source = (float *)source_;

  for(i = 0; i < blocks; i++) {
    cudaMemcpy(&float_tmp[0], &source[i * BUFF_SIZE], sizeof(float) * BUFF_SIZE, cudaMemcpyDeviceToHost);

    for(j = 0; j < BUFF_SIZE; j++) {
      targ[i * BUFF_SIZE + j] = (double)float_tmp[j];
    }
  }

  if(remainder > 0) {
    cudaMemcpy(&float_tmp[0], &source[i * BUFF_SIZE], sizeof(float) * remainder, cudaMemcpyDeviceToHost);

    for(j = 0; j < remainder; j++) {
      targ[blocks * BUFF_SIZE + j] = (double)float_tmp[j];
    }
  }
}

//Description:
//Allocates all necessary GPU-side arrays and sets global dimension constants
extern "C" void qlsum_gpu_initialize_(INT *nuper_, INT *nupar_, INT *nzeta_,
				      INT *nkdim1_, INT *nkdim2_, INT *mkdim1_, INT *mkdim2_,
				      INT *lmaxdim_, INT *nkx1_, INT *nkx2_, INT *nky1_, INT *nky2_,
				      REAL_H *xkxsav_, REAL_H *xkysav_,
				      CMPLX_H *ealphak_, CMPLX_H *ebetak_, CMPLX_H *ebk_) {
  nkx1 = *nkx1_; nkx2 = *nkx2_;
  nky1 = *nky1_; nky2 = *nky2_;
  nzeta = *nzeta_;
  nuper = *nuper_; nupar = *nupar_;
  nkdim1 = *nkdim1_; nkdim2 = *nkdim2_;
  mkdim1 = *mkdim1_; mkdim2 = *mkdim2_;
  lmaxdim = *lmaxdim_;

  //Exercise the cudaMalloc function
  cudaMalloc((void **)&sum_p, 8 * sizeof(CMPLX) * nuper * nupar * NUMBLOCKS);

  cudaMalloc((void **)&xx_p, sizeof(CMPLX) * (nkdim2 - nkdim1 + 1));
  cudaMalloc((void **)&yy_p, sizeof(CMPLX) * (mkdim2 - mkdim1 + 1));
  cudaMalloc((void **)&ealphak_p, sizeof(CMPLX) * (nkdim2 - nkdim1 + 1) * (mkdim2 - mkdim1 + 1));
  cudaMalloc((void **)&ebetak_p, sizeof(CMPLX) * (nkdim2 - nkdim1 + 1) * (mkdim2 - mkdim1 + 1));
  cudaMalloc((void **)&ebk_p, sizeof(CMPLX) * (nkdim2 - nkdim1 + 1) * (mkdim2 - mkdim1 + 1));

  cudaMalloc((void **)&sum_wdot_p, sizeof(CMPLX) * nuper * NUMBLOCKS);
  cudaMalloc((void **)&sum_fx0_p, sizeof(CMPLX) * nuper * NUMBLOCKS);
  cudaMalloc((void **)&sum_fy0_p, sizeof(CMPLX) * nuper * NUMBLOCKS);

  cudaMalloc((void **)&b_sum_p, sizeof(CMPLX) * nuper * nupar * NUMBLOCKS);
  cudaMalloc((void **)&c_sum_p, sizeof(CMPLX) * nuper * nupar * NUMBLOCKS);
  cudaMalloc((void **)&e_sum_p, sizeof(CMPLX) * nuper * nupar * NUMBLOCKS);
  cudaMalloc((void **)&f_sum_p, sizeof(CMPLX) * nuper * nupar * NUMBLOCKS);

  cudaMalloc((void **)&nres_p, sizeof(INT) * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1));
  cudaMalloc((void **)&mres_p, sizeof(INT) * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1));
  cudaMalloc((void **)&zbeta_p, sizeof(CMPLX) * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1));
  cudaMalloc((void **)&zbeta_iharm_p, sizeof(CMPLX) * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1));

  cudaMalloc((void **)&uper_p, sizeof(REAL) * nuper);
  cudaMalloc((void **)&upar_p, sizeof(REAL) * nupar);
  cudaMalloc((void **)&xkxsav_p, sizeof(REAL) * (nkdim2 - nkdim1 + 1));
  cudaMalloc((void **)&xkysav_p, sizeof(REAL) * (mkdim2 - mkdim1 + 1));
  cudaMalloc((void **)&dfduper_p, sizeof(REAL) * (nuper * nupar));
  cudaMalloc((void **)&dfdupar_p, sizeof(REAL) * (nuper * nupar));
  cudaMalloc((void **)&npara_sav_p, sizeof(REAL) * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1));

  cudaMalloc((void **)&sum2_p, 3 * sizeof(CMPLX) * nuper * NUMBLOCKS);
  cudaMalloc((void **)&sumkx2_p, 3 * sizeof(CMPLX) * nuper * NUMBLOCKS);
  cudaMalloc((void **)&sumky2_p, 3 * sizeof(CMPLX) * nuper * NUMBLOCKS);
  cudaMalloc((void **)&sumwdot_p, 2 * sizeof(CMPLX) * nuper * NUMBLOCKS);
  cudaMalloc((void **)&sumwdotkx_p, 2 * sizeof(CMPLX) * nuper * NUMBLOCKS);
  cudaMalloc((void **)&sumwdotky_p, 2 * sizeof(CMPLX) * nuper * NUMBLOCKS);

  cudaMalloc((void **)&xkperpn_tmp_p, sizeof(REAL) * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1));
  cudaMalloc((void **)&zetai_p, sizeof(REAL) * (nzeta + 1));
  cudaMalloc((void **)&Jni_p, sizeof(REAL) * nuper * (2 * lmaxdim + 1) * (nzeta + 1));
  cudaMalloc((void **)&zetamin_tmp_p, sizeof(REAL) * nuper);
  cudaMalloc((void **)&dzetai_tmp_p, sizeof(REAL) * nuper);

  cudaMalloc((void **)&eps_t_p, 3 * sizeof(CMPLX) * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1));

  cudaMalloc((void **)&sumb_11_nm_t_p, sizeof(CMPLX) * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1) * nuper);
  cudaMalloc((void **)&sumb_31_nm_t_p, sizeof(CMPLX) * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1) * nuper);

  //Copy over useful arrays
  cudaMemcpyD2SH2D(ealphak_p, ealphak_, 2 * (nkdim2 - nkdim1 + 1) * (mkdim2 - mkdim1 + 1));
  cudaMemcpyD2SH2D(ebetak_p, ebetak_, 2 * (nkdim2 - nkdim1 + 1) * (mkdim2 - mkdim1 + 1));
  cudaMemcpyD2SH2D(ebk_p, ebk_, 2 * (nkdim2 - nkdim1 + 1) * (mkdim2 - mkdim1 + 1));

  cudaMemcpyD2SH2D(xkxsav_p, xkxsav_, (nkdim2 - nkdim1 + 1));
  cudaMemcpyD2SH2D(xkysav_p, xkysav_, (mkdim2 - mkdim1 + 1));
}

//Description:
//Frees all dynamic arrays
extern "C" void qlsum_gpu_cleanup_() {
  cudaFree(sum_p);

  cudaFree(xx_p);
  cudaFree(yy_p);
  cudaFree(ealphak_p);
  cudaFree(ebetak_p);
  cudaFree(ebk_p);

  cudaFree(sum_wdot_p);
  cudaFree(sum_fx0_p);
  cudaFree(sum_fy0_p);

  cudaFree(b_sum_p);
  cudaFree(c_sum_p);
  cudaFree(e_sum_p);
  cudaFree(f_sum_p);

  cudaFree(nres_p);
  cudaFree(mres_p);
  cudaFree(zbeta_p);
  cudaFree(zbeta_iharm_p);

  cudaFree(uper_p);
  cudaFree(upar_p);
  cudaFree(xkxsav_p);
  cudaFree(xkysav_p);
  cudaFree(dfduper_p);
  cudaFree(dfdupar_p);
  cudaFree(npara_sav_p);

  cudaFree(sum2_p);
  cudaFree(sumkx2_p);
  cudaFree(sumky2_p);
  cudaFree(sumwdot_p);
  cudaFree(sumwdotkx_p);
  cudaFree(sumwdotky_p);

  cudaFree(xkperpn_tmp_p);
  cudaFree(zetai_p);
  cudaFree(Jni_p);
  cudaFree(zetamin_tmp_p);
  cudaFree(dzetai_tmp_p);

  cudaFree(eps_t_p);

  cudaFree(sumb_11_nm_t_p);
  cudaFree(sumb_31_nm_t_p);
}	

/*
  Description:
  Computes the shared epsx/epsy/epsz values

  There are iresmax epsx/epsy/epsz values generated
  To save pointer space, they are saved in one large array (eps_t)

  Because iresmax is at most (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1), we allocate
  an array of 3 * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1) elements, and say that
  the first (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1) are for epsx values,
  the second (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1) are for epsy values, etc.

  Output:
  eps_t

  Parallelism:
  There is one loop which operates across the resonant points

  Each block of threads takes a chunk of the 0 -> iresmax - 1 loop
  For each of these chunks, each thread in a block takes a consecutive
  element to operate on. There is no need for the chunksize to be a multiple
  of the threadsize.
*/
__global__ void qlsum_gpu_iharm_shared(INT nkx1, INT nkx2,
				       INT nky1, INT nky2,
				       INT nkdim1, INT nkdim2,
				       INT mkdim1,
				       INT iresmax, INT *nres, INT *mres,
				       CMPLX *zbeta, CMPLX *zbeta_iharm,
				       CMPLX *xx, CMPLX *yy,
				       CMPLX *ealphak, CMPLX *ebetak, CMPLX *ebk,
				       CMPLX *eps_t) {
  INT k_uper = threadIdx.x, block = blockIdx.x,
    ires, iresstart, iresfinish;

  INT xy_ind, dim_ind,
    n, m;

  CMPLX epsx, epsy, epsz,
    cexp1, cexp2, cexp0;

  iresstart = min((iresmax + NUMBLOCKS - 1) / NUMBLOCKS * block, iresmax);
  iresfinish = min(((iresmax + NUMBLOCKS - 1) / NUMBLOCKS) * (block + 1), iresmax);

  for(ires = iresstart + k_uper; ires < iresfinish; ires += blockDim.x) {
    n = nres[ires];
    m = mres[ires];

    xy_ind = (n - nkx1) + (m - nky1) * (nkx2 - nkx1 + 1);
    dim_ind = (n - nkdim1) + (m - mkdim1) * (nkdim2 - nkdim1 + 1);

    cexp2 = xx[(n - nkx1)] * yy[(m - nky1)] * zbeta_iharm[xy_ind];
    cexp0 = cexp2 * zbeta[xy_ind];
    cexp1 = cexp0 * zbeta[xy_ind];
	
    epsx = (ealphak[dim_ind] - zcmplx(0.0, 1.0) * ebetak[dim_ind]) * cexp1 / sqrtf(2.0f);
    epsy = (ealphak[dim_ind] + zcmplx(0.0, 1.0) * ebetak[dim_ind]) * cexp2 / sqrtf(2.0f);
    epsz = ebk[dim_ind] * cexp0;
      
    eps_t[ires + 0 * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1)] = epsx;
    eps_t[ires + 1 * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1)] = epsy;
    eps_t[ires + 2 * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1)] = epsz;
  }
}

/*
  Description:
  Computes various sums used in final reduction/generation of output
  Stores temporary copy of values derived from Bessel function interpolation
  (used by qlsum_gpu_iharm_second)

  Multiple arrays are stored within one array pointer to save argument space
  sum2 expands to sum2_1, sum2_2, and sum2_3 in the Fortran code
  sumkx2 and sumky2 are similar
  sum2_2 is accessed nuper * NUMBLOCKS elements from &sum[0]

  Output:
  sumb_11_nm_t, sumb_31_nm_t
  sum2, sumkx2, sumky2

  Parallelism:
  There is one explicit loop which operates across the resonant points

  Each block takes a chunk of this loop, which they work through in a +1 pattern
  Each thread within a block is assigned a k_uper value from 0 to nuper - 1
  Effectively, a nuper wide SIMD unit operates across the resonant points

  Because values must be accumulated across the resonant points, there must
  be a reduction across each of the sum2, sumkx2, and sumky2 data
*/
__global__ void qlsum_gpu_iharm_first(INT nkx1, INT nkx2,
				      INT nky1, INT nky2,
				      INT nkdim1, INT mkdim1,
				      INT nuper, INT nupar,
				      INT iharm, INT iresmax,
				      INT lmaxdim,
				      REAL nwcw,
				      REAL sqmut0, REAL dui, CMPLX *eps_t,
				      CMPLX *sumb_11_nm_t, CMPLX *sumb_31_nm_t,
				      INT *nres, INT *mres,
				      CMPLX *sum2, CMPLX *sumkx2, CMPLX *sumky2,
				      REAL *uper, REAL *upar,
				      REAL *xkxsav, REAL *xkysav,
				      REAL *dfduper, REAL *npara_sav,
				      REAL *xkperpn_tmp, REAL *zetai, REAL *Jni,
				      REAL *zetamin_tmp, REAL *dzetai_tmp) {
  INT k_uper = threadIdx.x, block = blockIdx.x,
    ires = 0, i2,
    idx2, idx3,
    iresstart, iresfinish;
  __shared__ INT xy_ind, n, m;

  REAL uper_kuper, zeta0, p2, zetamin, dzetai,
    Jnp_t, Jnm_t, Jnn_t;
  __shared__ REAL upar0, factor,  xkxsav_n, xkysav_m;

  CMPLX sum2_1 = 0.0, sum2_2 = 0.0, sum2_3 = 0.0,
    sumkx2_1 = 0.0, sumkx2_2, sumkx2_3 = 0.0,
    sumky2_1 = 0.0, sumky2_2, sumky2_3 = 0.0,
    sumb_11_nm = 0.0, sumb_31_nm = 0.0,
    sumwdot_11 = 0.0, sumwdot_31 = 0.0,
    sumwdotkx_11 = 0.0, sumwdotkx_31 = 0.0,
    sumwdotky_11 = 0.0, sumwdotky_31 = 0.0;
  __shared__ CMPLX epsx, epsy, epsz;

  __syncthreads(); //Every thread calls the __shared__ CMPLX constructor

  iresstart = min((iresmax + NUMBLOCKS - 1) / NUMBLOCKS * block, iresmax);
  iresfinish = min(((iresmax + NUMBLOCKS - 1) / NUMBLOCKS) * (block + 1), iresmax);

  uper_kuper = uper[k_uper];
  zetamin = zetamin_tmp[k_uper];
  dzetai = dzetai_tmp[k_uper];

  for(ires = iresstart; ires < iresfinish; ires++) {
    if(k_uper == 0) {
      n = nres[ires];
      m = mres[ires];

      xy_ind = (n - nkx1) + (m - nky1) * (nkx2 - nkx1 + 1);

      epsx = eps_t[ires + 0 * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1)];
      epsy = eps_t[ires + 1 * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1)];
      epsz = eps_t[ires + 2 * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1)];

      upar0 = sqmut0 / npara_sav[xy_ind] * (1.0 - nwcw);
      factor = M_PI * sqmut0 / abs(npara_sav[xy_ind]);

      xkxsav_n = xkxsav[n - nkdim1];
      xkysav_m = xkysav[m - mkdim1];
    }
    __syncthreads();

    zeta0 = xkperpn_tmp[(n - nkx1) + (m - nky1) * (nkx2 - nkx1 + 1)] * uper_kuper;
  
    i2 = int((zeta0 - zetamin) * dzetai);
    p2 = (zeta0 - (zetamin + (REAL)i2 / dzetai)) * dzetai;

    idx2 = k_uper + (iharm + lmaxdim) * nuper + i2 * (2 * lmaxdim + 1) * nuper;
    idx3 = k_uper + (iharm + lmaxdim) * nuper + (i2 + 1) * (2 * lmaxdim + 1) * nuper;
      
    Jnm_t = Jni[idx2 - nuper] + p2 * (Jni[idx3 - nuper] - Jni[idx2 - nuper]);
    Jnn_t = Jni[idx2] + p2 * (Jni[idx3] - Jni[idx2]);
    Jnp_t = Jni[idx2 + nuper] + p2 * (Jni[idx3 + nuper] - Jni[idx2 + nuper]);

    sum2_1 += ~epsx * Jnp_t; 
    sum2_2 += ~epsy * Jnm_t;
    sum2_3 += ~epsz * Jnn_t;

    sumkx2_1 += ~epsx * Jnp_t * xkxsav_n;
    sumkx2_2 += ~epsy * Jnm_t * xkxsav_n;
    sumkx2_3 += ~epsz * Jnn_t * xkxsav_n;

    sumky2_1 += ~epsx * Jnp_t * xkysav_m;
    sumky2_2 += ~epsy * Jnm_t * xkysav_m;
    sumky2_3 += ~epsz * Jnn_t * xkysav_m;

    sumb_11_nm = (epsx * uper_kuper * uper_kuper * Jnp_t + epsy * uper_kuper * uper_kuper  * Jnm_t + epsz * sqrtf(2.0f) * uper_kuper * upar0 * Jnn_t) * factor;
    sumb_31_nm = (epsx * sqrtf(2.0f) * uper_kuper * upar0 * Jnp_t + epsy * sqrtf(2.0f) * uper_kuper * upar0 * Jnm_t + epsz * 2.0 * upar0 * upar0 * Jnn_t) * factor;

    sumb_11_nm_t[k_uper + ires * nuper] = sumb_11_nm;
    sumb_31_nm_t[k_uper + ires * nuper] = sumb_31_nm;

    __syncthreads();
  }

  sum2[k_uper + block * nuper + 0 * nuper * NUMBLOCKS] += sum2_1; sum2[k_uper + block * nuper + 1 * nuper * NUMBLOCKS] += sum2_2; sum2[k_uper + block * nuper + 2 * nuper * NUMBLOCKS] += sum2_3;
  sumkx2[k_uper + block * nuper + 0 * nuper * NUMBLOCKS] += sumkx2_1; sumkx2[k_uper + block * nuper + 1 * nuper * NUMBLOCKS] += sumkx2_2; sumkx2[k_uper + block * nuper + 2 * nuper * NUMBLOCKS] += sumkx2_3;
  sumky2[k_uper + block * nuper + 0 * nuper * NUMBLOCKS] += sumky2_1; sumky2[k_uper + block * nuper + 1 * nuper * NUMBLOCKS] += sumky2_2; sumky2[k_uper + block * nuper + 2 * nuper * NUMBLOCKS] += sumky2_3;
}

/*
  Description:
  Computes various sums used in final reduction/generation of output

  Multiple arrays are stored within one array pointer to save argument space
  sum expands to sumb_11, sumb_31, sumc_11, sumc_31, sume_11, sume_31,
  sumf_11, and sumf_31 in the Fortran code
  The '31' elements are offset from the '11' elements by nuper * nupar * NUMBLOCKS
  The letters are offset from the base pointer by [B, C, E, or F] * nuper * nupar * NUMBLOCKS

  sumwdot expands to sumwdot_11, and sumwdot_31 in the Fortran code
  sumwdotkx and sumwdotky are similar
  The '31' elements are offset from the '11' elements by nuper * NUMBLOCKS

  Output:
  sum, sumwdot, sumwdotkx, sumwdotky

  Parallelism:
  There is one explicit loop which operates across the resonant points

  Each block takes a chunk of this loop, which they work through in a +1 pattern
  Each thread within a block is assigned a k_uper value from 0 to nuper - 1
  Effectively, a nuper wide SIMD unit operates across the resonant points

  Because values must be accumulated across the resonant points, there must
  be a reduction across each of effective arrays in the the sum, sumwdot,
  sumwdotkx, and sumwdotky data
*/
__global__ void qlsum_gpu_iharm_second(INT nkx1, INT nkx2,
				       INT nky1, INT nky2,
				       INT nkdim1, INT mkdim1,
				       INT nuper, INT nupar,
				       INT iharm, INT iresmax,
				       REAL nwcw,
				       REAL sqmut0, REAL dui, CMPLX *eps_t,
				       CMPLX *sumb_11_nm_t, CMPLX *sumb_31_nm_t,
				       INT *nres, INT *mres,
				       CMPLX *sum,
				       CMPLX *sumwdot, CMPLX *sumwdotkx, CMPLX *sumwdotky,
				       REAL *uper, REAL *upar,
				       REAL *xkxsav, REAL *xkysav,
				       REAL *dfduper, REAL *dfdupar,
				       REAL *npara_sav,
				       REAL *xkperpn_tmp, REAL *zetai, REAL *Jni,
				       REAL *zetamin_tmp, REAL *dzetai_tmp) {

  INT k_uper = threadIdx.x, block = blockIdx.x, ires = 0, iresstart, iresfinish;
  __shared__ INT i, xy_ind, n, m;

  REAL u, sinth, facte, u0, uper_kuper,
    dfduper0, dfdupar0,
    dfactpar, dfactper;
  __shared__ REAL npara1, upar0, p, xkxsav_n, xkysav_m;

  CMPLX sumf_11_t = 0.0, sumf_31_t = 0.0,
    sume_11_t = 0.0, sume_31_t = 0.0,
    sumc_11_t = 0.0, sumc_31_t = 0.0,
    sumb_11_t = 0.0, sumb_31_t = 0.0,
    sumb_11_nm = 0.0, sumb_31_nm = 0.0,
    sumwdot_11 = 0.0, sumwdot_31 = 0.0,
    sumwdotkx_11 = 0.0, sumwdotkx_31 = 0.0,
    sumwdotky_11 = 0.0, sumwdotky_31 = 0.0;
  __shared__ CMPLX epsx, epsy, epsz;

  __syncthreads(); //Every thread calls the __shared__ CMPLX constructor

  iresstart = min((iresmax + NUMBLOCKS - 1) / NUMBLOCKS * block, iresmax);
  iresfinish = min(((iresmax + NUMBLOCKS - 1) / NUMBLOCKS) * (block + 1), iresmax);

  uper_kuper = uper[k_uper];

  if(k_uper == 0) {
    n = nres[0];
    m = mres[0];

    xy_ind = (n - nkx1) + (m - nky1) * (nkx2 - nkx1 + 1);

    upar0 = sqmut0 / npara_sav[xy_ind] * (1.0 - nwcw);
    i = (int)floor((upar0 - upar[0]) * dui);
  }

  __syncthreads();

  for(ires = iresstart; ires < iresfinish; ires++) {
    if(k_uper == 0) {
      n = nres[ires];
      m = mres[ires];

      xy_ind = (n - nkx1) + (m - nky1) * (nkx2 - nkx1 + 1);
      npara1 = npara_sav[xy_ind];

      epsx = eps_t[ires + 0 * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1)];
      epsy = eps_t[ires + 1 * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1)];
      epsz = eps_t[ires + 2 * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1)];

      upar0 = sqmut0 / npara1 * (1.0 - nwcw);
      p = (upar0 - upar[i]) * dui;

      xkxsav_n = xkxsav[n - nkdim1];
      xkysav_m = xkysav[m - mkdim1];
    }
    __syncthreads();

    u = sqrt(upar0 * upar0 + uper_kuper * uper_kuper) + 0.00000001;
    sinth = uper_kuper * __frcp_rn((float)u) + 0.00000001;      
    facte = (nwcw - sinth * sinth) * __frcp_rn((float)upar0);
    
    //// This is the non-Maxwellian version (which covers the Maxwellian case)
    // Replacing the u0 = ... statement with
    // u0 = dfduper0 gives the Maxwellian qlsum
    //
    if(i <= (nupar - 1)) {
      dfduper0 = dfduper[k_uper + i * nuper] + (dfduper[k_uper + (i + 1) * nuper] - dfduper[k_uper + i * nuper]) * p;
      dfdupar0 = dfdupar[k_uper + i * nuper] + (dfdupar[k_uper + (i + 1) * nuper] - dfdupar[k_uper + i * nuper]) * p;
    } else {
      dfduper0 = dfduper[k_uper + i * nuper];
      dfdupar0 = dfdupar[k_uper + i * nuper];
    }

    dfactpar = npara1 * upar0 * __frcp_rn((float)sqmut0);
    dfactper = npara1 * uper_kuper * __frcp_rn((float)sqmut0);
    
    u0 = (1.0 - dfactpar) * dfduper0 + dfactper * dfdupar0;
    ////

    REAL sinth_inv = __frcp_rn((float)sinth);

    sumb_11_nm = sumb_11_nm_t[k_uper + ires * nuper];
    sumb_11_t += sumb_11_nm;
    sume_11_t += sumb_11_nm * facte;
    sumc_11_t += sumb_11_nm * facte * sinth_inv;
    sumf_11_t += sumb_11_nm * facte * facte * sinth_inv;
    sumwdot_11 += sumb_11_nm * u0;
    sumwdotkx_11 += sumb_11_nm * xkxsav_n * u0;
    sumwdotky_11 += sumb_11_nm * xkysav_m * u0;

    sumb_31_nm = sumb_31_nm_t[k_uper + ires * nuper];
    sumb_31_t += sumb_31_nm;
    sume_31_t += sumb_31_nm * facte;
    sumc_31_t += sumb_31_nm * facte * sinth_inv;
    sumf_31_t += sumb_31_nm * facte * facte * sinth_inv;
    sumwdot_31 += sumb_31_nm * u0;
    sumwdotkx_31 += sumb_31_nm * xkxsav_n * u0;
    sumwdotky_31 += sumb_31_nm * xkysav_m * u0;

    __syncthreads();
  }

  sum[k_uper + i * nuper + block * nuper * nupar + (F + 0) * nuper * nupar * NUMBLOCKS] += sumf_11_t;
  sum[k_uper + i * nuper + block * nuper * nupar + (F + 1) * nuper * nupar * NUMBLOCKS] += sumf_31_t;

  sum[k_uper + i * nuper + block * nuper * nupar + (E + 0) * nuper * nupar * NUMBLOCKS] += sume_11_t;
  sum[k_uper + i * nuper + block * nuper * nupar + (E + 1) * nuper * nupar * NUMBLOCKS] += sume_31_t;

  sum[k_uper + i * nuper + block * nuper * nupar + (C + 0) * nuper * nupar * NUMBLOCKS] += sumc_11_t;
  sum[k_uper + i * nuper + block * nuper * nupar + (C + 1) * nuper * nupar * NUMBLOCKS] += sumc_31_t;

  sum[k_uper + i * nuper + block * nuper * nupar + (B + 0) * nuper * nupar * NUMBLOCKS] += sumb_11_t;
  sum[k_uper + i * nuper + block * nuper * nupar + (B + 1) * nuper * nupar * NUMBLOCKS] += sumb_31_t;

  sumwdot[k_uper + block * nuper + 0 * nuper * NUMBLOCKS] += sumwdot_11;
  sumwdot[k_uper + block * nuper + 1 * nuper * NUMBLOCKS] += sumwdot_31;

  sumwdotkx[k_uper + block * nuper + 0 * nuper * NUMBLOCKS] += sumwdotkx_11;
  sumwdotkx[k_uper + block * nuper + 1 * nuper * NUMBLOCKS] += sumwdotkx_31;

  sumwdotky[k_uper + block * nuper + 0 * nuper * NUMBLOCKS] += sumwdotky_11;
  sumwdotky[k_uper + block * nuper + 1 * nuper * NUMBLOCKS] += sumwdotky_31;
}

__device__ CMPLX reduceCMPLX(CMPLX *array, int stride, int n) {
  INT i = 0;
  CMPLX acc = 0.0;

  for(i = 0; i < n; i++) {
    acc += array[i * stride];
  }

  return acc;
}

/*
  Description:
  Performs the reduction and calculates the qlsum output values

  Output:
  sum_wdot, sum_fx0, sum_fy0, b_sum, c_sum, e_sum, f_sum

  Parallelism:
  The reduction operates as a single nuper wide SIMD unit which runs along,
  reduces variables as necessary, and produces output

  This code only uses one multiprocessor
*/
__global__ void qlsum_gpu_iharm_reduction(INT nuper, INT nupar,
					  CMPLX *sum_wdot, CMPLX *sum_fx0, CMPLX *sum_fy0,
					  CMPLX *b_sum, CMPLX *c_sum, CMPLX *e_sum, CMPLX *f_sum,
					  CMPLX *sum,
					  CMPLX *sum2, CMPLX *sumkx2, CMPLX *sumky2,
					  CMPLX *sumwdot, CMPLX *sumwdotkx, CMPLX *sumwdotky) {
  INT k_uper = threadIdx.x, i = 0, j = 0, i_uprl = 0;

  CMPLX sumwdot_11 = 0.0, sumwdot_31 = 0.0,
    sumwdotkx_11 = 0.0, sumwdotkx_31 = 0.0,
    sumwdotky_11 = 0.0, sumwdotky_31 = 0.0,
    sum2_1 = 0.0, sum2_2 = 0.0, sum2_3 = 0.0,
    sumkx2_1 = 0.0, sumkx2_2 = 0.0, sumkx2_3 = 0.0,
    sumky2_1 = 0.0, sumky2_2 = 0.0, sumky2_3 = 0.0;
  
  //Perform reduction
  sumwdot_11 = reduceCMPLX(&sumwdot[k_uper + 0 * nuper * NUMBLOCKS], nuper, NUMBLOCKS);
  sumwdot_31 = reduceCMPLX(&sumwdot[k_uper + 1 * nuper * NUMBLOCKS], nuper, NUMBLOCKS);

  sumwdotkx_11 = reduceCMPLX(&sumwdotkx[k_uper + 0 * nuper * NUMBLOCKS], nuper, NUMBLOCKS);
  sumwdotky_31 = reduceCMPLX(&sumwdotkx[k_uper + 1 * nuper * NUMBLOCKS], nuper, NUMBLOCKS);

  sumwdotky_11 = reduceCMPLX(&sumwdotky[k_uper + 0 * nuper * NUMBLOCKS], nuper, NUMBLOCKS);
  sumwdotky_31 = reduceCMPLX(&sumwdotky[k_uper + 1 * nuper * NUMBLOCKS], nuper, NUMBLOCKS);

  sum2_1 = reduceCMPLX(&sum2[k_uper + 0 * nuper * NUMBLOCKS], nuper, NUMBLOCKS);
  sum2_2 = reduceCMPLX(&sum2[k_uper + 1 * nuper * NUMBLOCKS], nuper, NUMBLOCKS);
  sum2_3 = reduceCMPLX(&sum2[k_uper + 2 * nuper * NUMBLOCKS], nuper, NUMBLOCKS);

  sumkx2_1 = reduceCMPLX(&sumkx2[k_uper + 0 * nuper * NUMBLOCKS], nuper, NUMBLOCKS);
  sumkx2_2 = reduceCMPLX(&sumkx2[k_uper + 1 * nuper * NUMBLOCKS], nuper, NUMBLOCKS);
  sumkx2_3 = reduceCMPLX(&sumkx2[k_uper + 2 * nuper * NUMBLOCKS], nuper, NUMBLOCKS);

  sumky2_1 = reduceCMPLX(&sumky2[k_uper + 0 * nuper * NUMBLOCKS], nuper, NUMBLOCKS);
  sumky2_2 = reduceCMPLX(&sumky2[k_uper + 1 * nuper * NUMBLOCKS], nuper, NUMBLOCKS);
  sumky2_3 = reduceCMPLX(&sumky2[k_uper + 2 * nuper * NUMBLOCKS], nuper, NUMBLOCKS);
  
  //Calculate output
  sum_wdot[k_uper] +=
    sum2_1 * sumwdot_11 +
    sum2_2 * sumwdot_11 +
    sum2_3 * sumwdot_31;

  sum_fx0[k_uper] +=
    sumkx2_1 * sumwdot_11 +
    sumkx2_2 * sumwdot_11 +
    sumkx2_3 * sumwdot_31 +
    sum2_1 * sumwdotkx_11 +
    sum2_2 * sumwdotkx_11 +
    sum2_3 * sumwdotkx_31;

  sum_fy0[k_uper] +=
    sumky2_1 * sumwdot_11 +
    sumky2_2 * sumwdot_11 +
    sumky2_3 * sumwdot_31 +
    sum2_1 * sumwdotky_11 +
    sum2_2 * sumwdotky_11 +
    sum2_3 * sumwdotky_31;

  //Perform reduction
  for(j = 0; j < nupar; j++) {
    sum[k_uper + j * nuper + (B + 0) * nuper * nupar * NUMBLOCKS] = reduceCMPLX(&sum[k_uper + j * nuper + (B + 0) * nuper * nupar * NUMBLOCKS], nupar * nuper, NUMBLOCKS);
    sum[k_uper + j * nuper + (B + 1) * nuper * nupar * NUMBLOCKS] = reduceCMPLX(&sum[k_uper + j * nuper + (B + 1) * nuper * nupar * NUMBLOCKS], nupar * nuper, NUMBLOCKS);
      
    sum[k_uper + j * nuper + (C + 0) * nuper * nupar * NUMBLOCKS] = reduceCMPLX(&sum[k_uper + j * nuper + (C + 0) * nuper * nupar * NUMBLOCKS], nupar * nuper, NUMBLOCKS);
    sum[k_uper + j * nuper + (C + 1) * nuper * nupar * NUMBLOCKS] = reduceCMPLX(&sum[k_uper + j * nuper + (C + 1) * nuper * nupar * NUMBLOCKS], nupar * nuper, NUMBLOCKS);
      
    sum[k_uper + j * nuper + (E + 0) * nuper * nupar * NUMBLOCKS] = reduceCMPLX(&sum[k_uper + j * nuper + (E + 0) * nuper * nupar * NUMBLOCKS], nupar * nuper, NUMBLOCKS);
    sum[k_uper + j * nuper + (E + 1) * nuper * nupar * NUMBLOCKS] = reduceCMPLX(&sum[k_uper + j * nuper + (E + 1) * nuper * nupar * NUMBLOCKS], nupar * nuper, NUMBLOCKS);
      
    sum[k_uper + j * nuper + (F + 0) * nuper * nupar * NUMBLOCKS] = reduceCMPLX(&sum[k_uper + j * nuper + (F + 0) * nuper * nupar * NUMBLOCKS], nupar * nuper, NUMBLOCKS);
    sum[k_uper + j * nuper + (F + 1) * nuper * nupar * NUMBLOCKS] = reduceCMPLX(&sum[k_uper + j * nuper + (F + 1) * nuper * nupar * NUMBLOCKS], nupar * nuper, NUMBLOCKS);
  }

  //Calculate output
  for(i_uprl = 0; i_uprl < nupar; i_uprl++) {
    b_sum[k_uper + i_uprl * nuper] +=
      sum2_1 * sum[k_uper + i_uprl * nuper + (B + 0) * nuper * nupar * NUMBLOCKS] +
      sum2_2 * sum[k_uper + i_uprl * nuper + (B + 0) * nuper * nupar * NUMBLOCKS] +
      sum2_3 * sum[k_uper + i_uprl * nuper + (B + 1) * nuper * nupar * NUMBLOCKS];

    c_sum[k_uper + i_uprl * nuper] +=
      sum2_1 * sum[k_uper + i_uprl * nuper + (C + 0) * nuper * nupar * NUMBLOCKS] +
      sum2_2 * sum[k_uper + i_uprl * nuper + (C + 0) * nuper * nupar * NUMBLOCKS] +
      sum2_3 * sum[k_uper + i_uprl * nuper + (C + 1) * nuper * nupar * NUMBLOCKS];

    e_sum[k_uper + i_uprl * nuper] +=
      sum2_1 * sum[k_uper + i_uprl * nuper + (E + 0) * nuper * nupar * NUMBLOCKS] +
      sum2_2 * sum[k_uper + i_uprl * nuper + (E + 0) * nuper * nupar * NUMBLOCKS] +
      sum2_3 * sum[k_uper + i_uprl * nuper + (E + 1) * nuper * nupar * NUMBLOCKS];

    f_sum[k_uper + i_uprl * nuper] +=
      sum2_1 * sum[k_uper + i_uprl * nuper + (F + 0) * nuper * nupar * NUMBLOCKS] +
      sum2_2 * sum[k_uper + i_uprl * nuper + (F + 0) * nuper * nupar * NUMBLOCKS] +
      sum2_3 * sum[k_uper + i_uprl * nuper + (F + 1) * nuper * nupar * NUMBLOCKS];

  }
}

/*
  Description:
  Updates the zbeta values by multiplying in respective elements of zbeta_iharm

  Output:
  zbeta

  Parallelism:
  There is one loop which operates across the grid points (zbeta and zbeta_iharm
  really point to 2D arrays)

  Each block of threads takes a chunk of the 0 -> n - 1 loop
  For each of these chunks, each thread in a block takes a consecutive
  element to operate on. There is no need for the chunksize to be a multiple
  of the threadsize.
*/
__global__ void qlsum_gpu_iharm_zbeta_update(INT n, CMPLX *zbeta, CMPLX *zbeta_iharm) {
  INT thread = threadIdx.x, block = blockIdx.x, numthreads = blockDim.x, start, finish, i;

  start = min((n + NUMBLOCKS - 1) / NUMBLOCKS * block, n);
  finish = min(((n + NUMBLOCKS - 1) / NUMBLOCKS) * (block + 1), n);

  for(i = start + thread; i < finish; i += numthreads) {
    zbeta_iharm[i] *= zbeta[i];
  }
}

//Description:
//Initializes all the reduction temporaries
void qlsum_cpu_iharm_setup() {
  cudaMemset(sum_p, 0, 8 * sizeof(CMPLX) * nuper * nupar * NUMBLOCKS);
  cudaMemset(sum2_p, 0, 3 * sizeof(CMPLX) * nuper * NUMBLOCKS);
  cudaMemset(sumkx2_p, 0, 3 * sizeof(CMPLX) * nuper * NUMBLOCKS);
  cudaMemset(sumky2_p, 0, 3 * sizeof(CMPLX) * nuper * NUMBLOCKS);
  cudaMemset(sumwdot_p, 0, 2 * sizeof(CMPLX) * nuper * NUMBLOCKS);
  cudaMemset(sumwdotkx_p, 0, 2 * sizeof(CMPLX) * nuper * NUMBLOCKS);
  cudaMemset(sumwdotky_p, 0, 2 * sizeof(CMPLX) * nuper * NUMBLOCKS);
}

//Description:
//Provides a wrapper around the GPU grunt work kernel calls which process
// a worklist for a given iharm number
void qlsum_cpu_iharm_accumulate(INT iharm, INT iresmax,
				INT *nres, INT *mres) {
  dim3 threads(nuper);
  dim3 blocks(NUMBLOCKS);

  if(iresmax > 0) {
    cudaMemcpy(nres_p, nres, sizeof(INT) * iresmax, cudaMemcpyHostToDevice);
    cudaMemcpy(mres_p, mres, sizeof(INT) * iresmax, cudaMemcpyHostToDevice);

    //cudaError_t error;
    //error = cudaGetLastError();
    //if(error != cudaSuccess) { printf("Before memcopies, cuError: %s\n", cudaGetErrorString(error)); }

    qlsum_gpu_iharm_shared<<<blocks, 64>>>(nkx1, nkx2,
					   nky1, nky2,
					   nkdim1, nkdim2,
					   mkdim1,
					   iresmax, nres_p, mres_p,
					   zbeta_p, zbeta_iharm_p,
					   xx_p, yy_p,
					   ealphak_p, ebetak_p, ebk_p,
					   eps_t_p);

    qlsum_gpu_iharm_first<<<blocks, threads>>>(nkx1, nkx2,
					       nky1, nky2,
					       nkdim1, mkdim1,
					       nuper, nupar,
					       iharm, iresmax,
					       lmaxdim,
					       iharm * wcw,
					       sqmut0, dui, eps_t_p,
					       sumb_11_nm_t_p, sumb_31_nm_t_p,
					       nres_p, mres_p,
					       sum2_p, sumkx2_p, sumky2_p,
					       uper_p, upar_p,
					       xkxsav_p, xkysav_p,
					       dfduper_p, npara_sav_p,
					       xkperpn_tmp_p, zetai_p, Jni_p,
					       zetamin_tmp_p, dzetai_tmp_p);

    qlsum_gpu_iharm_second<<<blocks, threads>>>(nkx1, nkx2,
						nky1, nky2,
						nkdim1, mkdim1,
						nuper, nupar,
						iharm, iresmax,
						iharm * wcw,
						sqmut0, dui, eps_t_p,
						sumb_11_nm_t_p, sumb_31_nm_t_p,
						nres_p, mres_p,
						sum_p,
						sumwdot_p, sumwdotkx_p, sumwdotky_p,
						uper_p, upar_p,
						xkxsav_p, xkysav_p,
						dfduper_p, dfdupar_p,
						npara_sav_p,
						xkperpn_tmp_p, zetai_p, Jni_p,
						zetamin_tmp_p, dzetai_tmp_p);
  }

  cudaThreadSynchronize();
}

//Description:
//Provides a wrapper around the GPU reduction
void qlsum_cpu_iharm_finish() {
  dim3 threads(nuper);

  qlsum_gpu_iharm_reduction<<<1, threads>>>(nuper, nupar,
					    sum_wdot_p, sum_fx0_p, sum_fy0_p,
					    b_sum_p, c_sum_p, e_sum_p, f_sum_p,
					    sum_p,
					    sum2_p, sumkx2_p, sumky2_p,
					    sumwdot_p, sumwdotkx_p, sumwdotky_p);

  cudaThreadSynchronize();
}

//Description:
//Provides a wrapper around the GPU zbeta update
void qlsum_cpu_iharm_zbeta_update() {
  dim3 blocks(NUMBLOCKS);

  qlsum_gpu_iharm_zbeta_update<<<blocks, 64>>>((nkx2 - nkx1 + 1) * (nky2 - nky1 + 1), zbeta_p, zbeta_iharm_p);

  cudaThreadSynchronize();
}

//Description:
//Zeroes all the output values and copies over all the CPU generated arrays
void qlsum_cpu_nharm_setup(REAL_H *uper_, REAL_H *upar_,
			   REAL_H *dfduper_, REAL_H *dfdupar_,
			   REAL_H *dui_, REAL_H *wcw_,
			   REAL_H *sqmut0_,
			   REAL_H *npara_sav_, REAL_H *Jni_,
			   REAL_H *xkperpn_tmp_, REAL_H *zetamin_tmp_,
			   REAL_H *zetai_, REAL_H *dzetai_tmp_,
			   CMPLX_H *xx_, CMPLX_H *yy_,
			   CMPLX_H *zbeta, CMPLX_H *zbeta_iharm) {
  //Zero accumulators
  cudaMemset(sum_wdot_p, 0, sizeof(CMPLX) * nuper * NUMBLOCKS);
  cudaMemset(sum_fx0_p, 0, sizeof(CMPLX) * nuper * NUMBLOCKS);
  cudaMemset(sum_fy0_p, 0, sizeof(CMPLX) * nuper * NUMBLOCKS);

  cudaMemset(b_sum_p, 0, sizeof(CMPLX) * nuper * nupar * NUMBLOCKS);
  cudaMemset(c_sum_p, 0, sizeof(CMPLX) * nuper * nupar * NUMBLOCKS);
  cudaMemset(e_sum_p, 0, sizeof(CMPLX) * nuper * nupar * NUMBLOCKS);
  cudaMemset(f_sum_p, 0, sizeof(CMPLX) * nuper * nupar * NUMBLOCKS);

  //Copy over useful arrays
  cudaMemcpyD2SH2D(xx_p, xx_, 2 * (nkx2 - nkx1 + 1));
  cudaMemcpyD2SH2D(yy_p, yy_, 2 * (nky2 - nky1 + 1));
  cudaMemcpyD2SH2D(uper_p, uper_, nuper);
  cudaMemcpyD2SH2D(upar_p, upar_, nupar);
  cudaMemcpyD2SH2D(zbeta_p, zbeta, 2 * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1));
  cudaMemcpyD2SH2D(zbeta_iharm_p, zbeta_iharm, 2 * (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1));

  cudaMemcpyD2SH2D(xkperpn_tmp_p, xkperpn_tmp_, (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1));
  cudaMemcpyD2SH2D(zetai_p, zetai_, (nzeta + 1));
  cudaMemcpyD2SH2D(Jni_p, Jni_, nuper * (2 * lmaxdim + 1) * (nzeta + 1));
  cudaMemcpyD2SH2D(zetamin_tmp_p, zetamin_tmp_, nuper);
  cudaMemcpyD2SH2D(dzetai_tmp_p, dzetai_tmp_, nuper);

  cudaMemcpyD2SH2D(dfduper_p, dfduper_, nuper * nupar);
  cudaMemcpyD2SH2D(dfdupar_p, dfdupar_, nuper * nupar);
  cudaMemcpyD2SH2D(npara_sav_p, npara_sav_, (nkx2 - nkx1 + 1) * (nky2 - nky1 + 1));

  //Set single values
  sqmut0 = *sqmut0_;
  dui = *dui_;
  wcw = *wcw_;
}

//Description:
//Copies the output values back to the CPU
void qlsum_cpu_nharm_finish(CMPLX_H *sum_wdot_, CMPLX_H *sum_fx0_, CMPLX_H *sum_fy0_,
			    CMPLX_H *b_sum_, CMPLX_H *c_sum_, CMPLX_H *e_sum_, CMPLX_H *f_sum_) {
  cudaMemcpyS2DD2H(sum_wdot_, sum_wdot_p, 2 * nuper);
  cudaMemcpyS2DD2H(sum_fx0_, sum_fx0_p, 2 * nuper);
  cudaMemcpyS2DD2H(sum_fy0_, sum_fy0_p, 2 * nuper);
    
  cudaMemcpyS2DD2H(b_sum_, b_sum_p, 2 * nuper * nupar);
  cudaMemcpyS2DD2H(c_sum_, c_sum_p, 2 * nuper * nupar);
  cudaMemcpyS2DD2H(e_sum_, e_sum_p, 2 * nuper * nupar);
  cudaMemcpyS2DD2H(f_sum_, f_sum_p, 2 * nuper * nupar);
}
