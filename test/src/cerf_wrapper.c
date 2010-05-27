#include <complex.h>

int cerf_wrap_ (double z_re, double z_im, double out_re, double out_im){

	double complex z, result;

	z = ( z_re, z_im );

	result = cerf ( z );

	int stat;	
	return stat;

}
