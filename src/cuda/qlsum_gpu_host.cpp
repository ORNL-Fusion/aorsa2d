#include <stdio.h>
#include <complex>
#include <vector>
#include <omp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

#define INT int
#define REAL double
#define CMPLX std::complex<double>

#define REAL_H double
#define CMPLX_H std::complex<double>

//This includes must follow the type declaration. I apologize, but this is the easiest way
#include "qlsum_gpu_kernels.cuh"

extern "C" void qlsum_gpu_(INT *nkx1_, INT *nkx2_,
			   INT *nky1_, INT *nky2_,
			   INT *nharm,
			   INT *nuper_, INT *nupar_,
			   INT *nres, INT *mres,
			   REAL *wcw,
			   REAL *uparmax, REAL *uparmin,
			   REAL *uper, REAL *upar,
			   REAL *dfduper, REAL *dfdupar,
			   REAL *dui, REAL *sqmut0,
			   REAL *xkxsav, REAL *xkysav,
			   REAL *npara_sav, REAL *Jni,
			   REAL *xkperpn_tmp, REAL *zetamin_tmp,
			   REAL *zetai, REAL *dzetai_tmp,
			   CMPLX *zbeta, CMPLX *zbeta_iharm,
			   CMPLX *xx, CMPLX *yy,
			   CMPLX *ealphak, CMPLX *ebetak, CMPLX *ebk,
			   CMPLX *sum_wdot, CMPLX *sum_fx0, CMPLX *sum_fy0,
			   CMPLX *b_sum, CMPLX *c_sum, CMPLX *e_sum, CMPLX *f_sum) {
  int nkx1 = *nkx1_, nkx2 = *nkx2_, nky1 = *nky1_, nky2 = *nky2_,
    nuper = *nuper_, nupar = *nupar_;

  int i = 0, n = 0, m = 0, iharm = 0, ires = 0, iresmax = 0, nm_ind, irestotal = 0;

  double time = 0.0, copyin_time = 0.0, copyout_time = 0.0, worklist_time = 0.0, gpu_setup_time = 0.0, gpu_acc_time = 0.0, gpu_finish_time = 0.0, zbeta_time = 0.0;

  std::vector<std::vector<INT> > nres_vecs;
  std::vector<std::vector<INT> > mres_vecs;

  nres_vecs.resize(nupar);
  mres_vecs.resize(nupar);

  time = omp_get_wtime();
  qlsum_cpu_nharm_setup(uper, upar,
			dfduper, dfdupar,
			dui, wcw,
			sqmut0,
			npara_sav, Jni,
			xkperpn_tmp, zetamin_tmp,
			zetai, dzetai_tmp,
			xx, yy,
			zbeta, zbeta_iharm);
  copyin_time = omp_get_wtime() - time;

  for(iharm = -(*nharm); iharm <= *nharm; iharm++) {
    REAL nwcw = iharm * *wcw;

    time = omp_get_wtime();
    for(i = 0; i < nupar; i++) {
      nres_vecs[i].clear();
      mres_vecs[i].clear();
    }

    ires = 0;

    for(n = nkx1; n <= nkx2; n++) {
      for(m = nky1; m <= nky2; m++) {
	nm_ind = (n - nkx1) + (m - nky1) * (nkx2 - nkx1 + 1);

	REAL rrp = 1.0 - nwcw - npara_sav[nm_ind] * *uparmax / *sqmut0;
	REAL rrm = 1.0 - nwcw - npara_sav[nm_ind] * *uparmin / *sqmut0;

	if(rrp * rrm < 0.0) {
	  REAL upar0 = *sqmut0 / npara_sav[nm_ind] * (1.0 - nwcw);
	  INT it = (int)floor((upar0 - upar[0]) * *dui) + 1;

	  nres_vecs[it].push_back(n);
	  mres_vecs[it].push_back(m);

	  ires++;
	}
      }
    }

    iresmax = ires;

    irestotal += iresmax;

    worklist_time += (omp_get_wtime() - time);

    INT numblocks = 0, usedblocks = 0;

    if(iresmax > 0) {
      time = omp_get_wtime();
      qlsum_cpu_iharm_setup();
      gpu_setup_time += (omp_get_wtime() - time);

      for(i = 0; i < nupar; i++) {
	INT worklength = nres_vecs[i].size();
	numblocks = std::min(worklength / 16 + 1, MAXBLOCKS);
	usedblocks = std::max(numblocks, usedblocks);

	time = omp_get_wtime();
	qlsum_cpu_iharm_accumulate(numblocks, iharm, nres_vecs[i].size(),
				   &nres_vecs[i].front(), &mres_vecs[i].front());
	gpu_acc_time += (omp_get_wtime() - time);
      }

      time = omp_get_wtime();
      qlsum_cpu_iharm_finish(usedblocks);
      gpu_finish_time += (omp_get_wtime() - time);
    }

    time = omp_get_wtime();
    qlsum_cpu_iharm_zbeta_update();
    zbeta_time += (omp_get_wtime() - time);
  }

  //printf("Worklist time %f, zbeta time %f\n", worklist_time, zbeta_time);
  //printf("gpu setup time %f, gpu acc time %f, gpu finish time %f\n", gpu_setup_time, gpu_acc_time, gpu_finish_time);

  time = omp_get_wtime();
  qlsum_cpu_nharm_finish(sum_wdot, sum_fx0, sum_fy0,
			 b_sum, c_sum, e_sum, f_sum);
  copyout_time = omp_get_wtime() - time;

  printf("%d %f %f %f %f %f %f %f\n", irestotal, copyin_time, copyout_time, worklist_time, gpu_setup_time, gpu_acc_time, gpu_finish_time, zbeta_time);
}

