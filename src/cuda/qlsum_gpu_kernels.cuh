#ifndef QLSUM_GPU_KERNELS_CUH_
#define QLSUM_GPU_KERNELS_CUH_

extern "C" void qlsum_cpu_iharm_setup();

extern "C" void qlsum_cpu_iharm_accumulate(INT iharm, INT iresmax,
					   INT *nres, INT *mres);

extern "C" void qlsum_cpu_iharm_finish();

extern "C" void qlsum_cpu_iharm_zbeta_update();

extern "C" void qlsum_cpu_nharm_setup(REAL_H *nuper_, REAL_H *nupar_,
				      REAL_H *dfduper_, REAL_H *dfdupar_,
				      REAL_H *dui_, REAL_H *wcw_,
				      REAL_H *sqmut0_,
				      REAL_H *npara_sav_, REAL_H *Jni_,
				      REAL_H *xkperpn_tmp_, REAL_H *zetamin_tmp_,
				      REAL_H *zetai_, REAL_H *dzetai_tmp_,
				      CMPLX_H *xx_, CMPLX_H *yy_,
				      CMPLX_H *zbeta, CMPLX_H *zbeta_iharm);

extern "C" void qlsum_cpu_nharm_finish(CMPLX_H *sum_wdot_, CMPLX_H *sum_fx0_, CMPLX_H *sum_fy0_,
				       CMPLX_H *b_sum_, CMPLX_H *c_sum_, CMPLX_H *e_sum_, CMPLX_H *f_sum_);

#endif
