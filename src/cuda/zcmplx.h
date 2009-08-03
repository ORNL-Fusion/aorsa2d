#include <cuComplex.h>

#ifndef ZCMPLX_H_
#define ZCMPLX_H_

// Depending on whether we're running inside the CUDA compiler, define the __host_
// and __device__ intrinsics, otherwise just make the functions static to prevent
// linkage issues (duplicate symbols and such)
#ifdef __CUDACC__
#define HOST __host__
#define DEVICE __device__
#define HOSTDEVICE __host__ __device__
#else
#define HOST
#define DEVICE
#define HOSTDEVICE
#endif

// Struct alignment is handled differently between the CUDA compiler and other
// compilers (e.g. GCC, MS Visual C++ .NET)
#ifdef __CUDACC__
#define ALIGN(x)  __align__(x)
#else
#if defined(__GNUC__)
#define ALIGN(x)  __attribute__ ((aligned (x)))
#endif
#endif

class zcmplx : public cuDoubleComplex {
 public:
  //Regular constructors/destructor
  HOSTDEVICE zcmplx() { x = 0.0; y = 0.0; }
  HOSTDEVICE zcmplx(const double &_real, const double &_imag) { x = _real; y = _imag; }
  HOSTDEVICE ~zcmplx() { }

  //Various copy constructors
  HOSTDEVICE zcmplx(const cuFloatComplex &nval) { x = (double)nval.x; y = (double)nval.y; }
  HOSTDEVICE zcmplx(const cuDoubleComplex &nval) { x = nval.x; y = nval.y; }
  HOSTDEVICE zcmplx(const zcmplx &nval) { x = nval.real(); y = nval.imag(); }
  HOSTDEVICE zcmplx(const int &nval) { x = (double)nval; y = 0; }
  HOSTDEVICE zcmplx(const float &nval) { x = nval; y = 0; }
  HOSTDEVICE zcmplx(const double &nval) { x = (double)nval; y = 0; }

  //Assignment operators
  HOSTDEVICE zcmplx& operator=(const cuFloatComplex &nval) { x = (double)nval.x; y = (double)nval.y; return *this; }
  HOSTDEVICE zcmplx& operator=(const cuDoubleComplex &nval) { x = nval.x; y = nval.y; return *this; }
  HOSTDEVICE zcmplx& operator=(const zcmplx &nval) { x = nval.x; y = nval.y; return *this; }
  HOSTDEVICE zcmplx& operator=(const int &nval) { x = (double)nval; y = 0.0; return *this; }
  HOSTDEVICE zcmplx& operator=(const float &nval) { x = nval; y = 0.0; return *this; }
  HOSTDEVICE zcmplx& operator=(const double &nval) { x = (double)nval; y = 0.0; return *this; }

  //Addition
  HOSTDEVICE zcmplx operator+(const zcmplx &nval) const;
  HOSTDEVICE zcmplx operator+(const double &nval) const;
  HOSTDEVICE zcmplx& operator+=(const zcmplx &nval);
  //Subtraction
  HOSTDEVICE zcmplx operator-(const zcmplx &nval) const;
  HOSTDEVICE zcmplx operator-(const double &nval) const;
  HOSTDEVICE zcmplx& operator-=(const zcmplx &nval);
  //Multiplication
  HOSTDEVICE zcmplx operator*(const zcmplx &nval) const;
  HOSTDEVICE zcmplx operator*(const double &nval) const;
  HOSTDEVICE zcmplx& operator*=(const zcmplx &nval);
  //Division
  HOSTDEVICE zcmplx operator/(const zcmplx &nval) const;
  HOSTDEVICE zcmplx operator/(const double &nval) const;
  HOSTDEVICE zcmplx& operator/=(const zcmplx &nval);

  //Conjugate
  HOSTDEVICE zcmplx operator~();

  HOSTDEVICE double real() const { return cuCreal(*this); }
  HOSTDEVICE double imag() const { return cuCimag(*this); }
};

//Addition
HOSTDEVICE zcmplx zcmplx::operator+(const zcmplx &nval) const {
  return cuCadd(*this, nval);
}

/*HOSTDEVICE zcmplx operator+(const double &sval, const zcmplx &nval) {
  return nval + sval;
  }*/

HOSTDEVICE zcmplx zcmplx::operator+(const double &nval) const {
  return zcmplx(x + nval, y);
}

HOSTDEVICE zcmplx& zcmplx::operator+=(const zcmplx &nval) {
  return *this = *this + nval;
}

//Subtraction
HOSTDEVICE zcmplx zcmplx::operator-(const zcmplx &nval) const {
  return cuCsub(*this, nval);
}

/*HOSTDEVICE zcmplx operator-(const double &sval, const zcmplx &nval) {
  return zcmplx(sval - nval.x, nval.y);
  }*/

HOSTDEVICE zcmplx zcmplx::operator-(const double &nval) const {
  return zcmplx(x - nval, y);
}

HOSTDEVICE zcmplx& zcmplx::operator-=(const zcmplx &nval) {
  return *this = *this - nval;
}

//Multiplication
HOSTDEVICE zcmplx zcmplx::operator*(const zcmplx &nval) const {
  return cuCmul(*this, nval);
}

/*HOSTDEVICE zcmplx operator*(const double &sval, const zcmplx &nval) {
  return nval * sval;
  }*/

HOSTDEVICE zcmplx zcmplx::operator*(const double &nval) const {
  return zcmplx(x * nval, y * nval);
}

HOSTDEVICE zcmplx& zcmplx::operator*=(const zcmplx &nval) {
  return *this = *this * nval;
}

//Division
HOSTDEVICE zcmplx zcmplx::operator/(const zcmplx &nval) const {
  return cuCdiv(*this, nval);
}

/*HOSTDEVICE zcmplx operator/(const double &sval, const zcmplx &nval) {
  return ((zcmplx)sval) / nval;
  }*/

HOSTDEVICE zcmplx zcmplx::operator/(const double &nval) const {
  return zcmplx(x / nval, y / nval);
}

HOSTDEVICE zcmplx& zcmplx::operator/=(const zcmplx &nval) {
  return *this = *this / nval;
}

//Conjugate
HOSTDEVICE zcmplx zcmplx::operator~() {
  return cuConj(*this);
}

#endif // #ifndef ZCMPLX_H_
