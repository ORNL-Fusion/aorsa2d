#include <cuComplex.h>

#ifndef CCMPLX_H_
#define CCMPLX_H_

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

class ccmplx : public cuFloatComplex {
 public:
  //Regular constructors/destructor
  HOSTDEVICE ccmplx() { x = 0.0; y = 0.0; }
  HOSTDEVICE ccmplx(const float &_real, const float &_imag) { x = _real; y = _imag; }
  HOSTDEVICE ~ccmplx() { }

  //Various copy constructors
  HOSTDEVICE ccmplx(const cuFloatComplex &nval) { x = nval.x; y = nval.y; }
  HOSTDEVICE ccmplx(const cuDoubleComplex &nval) { x = (float)nval.x; y = (float)nval.y; }
  HOSTDEVICE ccmplx(const ccmplx &nval) { x = nval.real(); y = nval.imag(); }
  HOSTDEVICE ccmplx(const int &nval) { x = (float)nval; y = 0.0; }
  HOSTDEVICE ccmplx(const float &nval) { x = nval; y = 0.0; }
  HOSTDEVICE ccmplx(const double &nval) { x = (float)nval; y = 0.0; }

  //Assignment operators
  HOSTDEVICE ccmplx& operator=(const cuFloatComplex &nval) { x = nval.x; y = nval.y; return *this; }
  HOSTDEVICE ccmplx& operator=(const cuDoubleComplex &nval) { x = (float)nval.x; y = (float)nval.y; return *this; }
  HOSTDEVICE ccmplx& operator=(const ccmplx &nval) { x = nval.x; y = nval.y; return *this; }
  HOSTDEVICE ccmplx& operator=(const int &nval) { x = (float)nval; y = 0.0; return *this; }
  HOSTDEVICE ccmplx& operator=(const float &nval) { x = nval; y = 0.0; return *this; }
  HOSTDEVICE ccmplx& operator=(const double &nval) { x = (float)nval; y = 0.0; return *this; }

  //Addition
  HOSTDEVICE ccmplx operator+(const ccmplx &nval) const;
  HOSTDEVICE ccmplx operator+(const float &nval) const;
  HOSTDEVICE ccmplx& operator+=(const ccmplx &nval);
  //Subtraction
  HOSTDEVICE ccmplx operator-(const ccmplx &nval) const;
  HOSTDEVICE ccmplx operator-(const float &nval) const;
  HOSTDEVICE ccmplx& operator-=(const ccmplx &nval);
  //Multiplication
  HOSTDEVICE ccmplx operator*(const ccmplx &nval) const;
  HOSTDEVICE ccmplx operator*(const float &nval) const;
  HOSTDEVICE ccmplx& operator*=(const ccmplx &nval);
  //Division
  HOSTDEVICE ccmplx operator/(const ccmplx &nval) const;
  HOSTDEVICE ccmplx operator/(const float &nval) const;
  HOSTDEVICE ccmplx& operator/=(const ccmplx &nval);

  //Conjugate
  HOSTDEVICE ccmplx operator~();

  HOSTDEVICE float real() const { return cuCrealf(*this); }
  HOSTDEVICE float imag() const { return cuCimagf(*this); }
};

//Addition
HOSTDEVICE ccmplx ccmplx::operator+(const ccmplx &nval) const {
  return cuCaddf(*this, nval);
}

HOSTDEVICE ccmplx operator+(const float &sval, const ccmplx &nval) {
  return nval + sval;
}

HOSTDEVICE ccmplx ccmplx::operator+(const float &nval) const {
  return ccmplx(x + nval, y);
}

HOSTDEVICE ccmplx& ccmplx::operator+=(const ccmplx &nval) {
  return *this = *this + nval;
}

//Subtraction
HOSTDEVICE ccmplx ccmplx::operator-(const ccmplx &nval) const {
  return cuCsubf(*this, nval);
}

HOSTDEVICE ccmplx operator-(const float &sval, const ccmplx &nval) {
  return ccmplx(sval - nval.x, nval.y);
}

HOSTDEVICE ccmplx ccmplx::operator-(const float &nval) const {
  return ccmplx(x - nval, y);
}

HOSTDEVICE ccmplx& ccmplx::operator-=(const ccmplx &nval) {
  return *this = *this - nval;
}

//Multiplication
HOSTDEVICE ccmplx ccmplx::operator*(const ccmplx &nval) const {
  return cuCmulf(*this, nval);
}

HOSTDEVICE ccmplx operator*(const float &sval, const ccmplx &nval) {
  return nval * sval;
}

HOSTDEVICE ccmplx ccmplx::operator*(const float &nval) const {
  return ccmplx(x * nval, y * nval);
}

HOSTDEVICE ccmplx& ccmplx::operator*=(const ccmplx &nval) {
  return *this = *this * nval;
}

//Division
HOSTDEVICE ccmplx ccmplx::operator/(const ccmplx &nval) const {
  return cuCdivf(*this, nval);
}

HOSTDEVICE ccmplx operator/(const float &sval, const ccmplx &nval) {
  return ((ccmplx)sval) / nval;
}

HOSTDEVICE ccmplx ccmplx::operator/(const float &nval) const {
  return ccmplx(x / nval, y / nval);
}

HOSTDEVICE ccmplx& ccmplx::operator/=(const ccmplx &nval) {
  return *this = *this / nval;
}

//Conjugate
HOSTDEVICE ccmplx ccmplx::operator~() {
  return cuConjf(*this);
}

#endif // #ifndef CCMPLX_H_
