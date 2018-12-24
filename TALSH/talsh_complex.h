/** ExaTensor::TAL-SH: Complex arithmetic header.
REVISION: 2018/12/24

Copyright (C) 2014-2018 Dmitry I. Lyakh (Liakh)
Copyright (C) 2014-2018 Oak Ridge National Laboratory (UT-Battelle)

This file is part of ExaTensor.

ExaTensor is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ExaTensor is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with ExaTensor. If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------
**/

#ifndef TALSH_COMPLEX_H_
#define TALSH_COMPLEX_H_

#include <math.h>

#ifdef __cplusplus
#include <complex>
#endif

#ifndef NO_GPU
#include <cuComplex.h>
#endif

//DECLARATIONS:
// Complex number:
#ifndef NO_GPU
typedef cuFloatComplex talshComplex4;
typedef cuDoubleComplex talshComplex8;
#else
#ifdef __cplusplus
typedef std::complex<float> talshComplex4;
typedef std::complex<double> talshComplex8;
#else
typedef struct{float real; float imag;} talshComplex4;
typedef struct{double real; double imag;} talshComplex8;
#endif
#endif /*NO_GPU*/

/* TAL-SH complex arithmetic headers:
inline talshComplex4 talshComplex4Set(float real, float imag);
inline talshComplex8 talshComplex8Set(double real, double imag);
inline float talshComplex4Real(talshComplex4 cmplx);
inline double talshComplex8Real(talshComplex8 cmplx);
inline float talshComplex4Imag(talshComplex4 cmplx);
inline double talshComplex8Imag(talshComplex8 cmplx);
inline talshComplex4 talshComplex4Conjg(talshComplex4 cmplx);
inline talshComplex8 talshComplex8Conjg(talshComplex8 cmplx);
inline float talshComplex4Abs(talshComplex4 cmplx);
inline double talshComplex8Abs(talshComplex8 cmplx);
inline talshComplex4 talshComplex4Add(talshComplex4 x, talshComplex4 y);
inline talshComplex8 talshComplex8Add(talshComplex8 x, talshComplex8 y);
inline talshComplex4 talshComplex4Sub(talshComplex4 x, talshComplex4 y);
inline talshComplex8 talshComplex8Sub(talshComplex8 x, talshComplex8 y);
inline talshComplex4 talshComplex4Mul(talshComplex4 x, talshComplex4 y);
inline talshComplex8 talshComplex8Mul(talshComplex8 x, talshComplex8 y);
inline talshComplex4 talshComplex4Div(talshComplex4 x, talshComplex4 y);
inline talshComplex8 talshComplex8Div(talshComplex8 x, talshComplex8 y);
*/

//DEFINITIONS:
// Complex arithmetic:
#ifndef NO_GPU
__host__ __device__ inline talshComplex4 talshComplex4Set(float real, float imag)
#else
inline talshComplex4 talshComplex4Set(float real, float imag)
#endif
{
#ifndef NO_GPU
 talshComplex4 result = make_cuFloatComplex(real,imag);
#else
#ifdef __cplusplus
 talshComplex4 result(real,imag);
#else
 talshComplex4 result = {real,imag};
#endif
#endif
 return result;
}

#ifndef NO_GPU
__host__ __device__ inline talshComplex8 talshComplex8Set(double real, double imag)
#else
inline talshComplex8 talshComplex8Set(double real, double imag)
#endif
{
#ifndef NO_GPU
 talshComplex8 result = make_cuDoubleComplex(real,imag);
#else
#ifdef __cplusplus
 talshComplex8 result(real,imag);
#else
 talshComplex8 result = {real,imag};
#endif
#endif
 return result;
}

#ifndef NO_GPU
__host__ __device__ inline float talshComplex4Real(talshComplex4 cmplx)
#else
inline float talshComplex4Real(talshComplex4 cmplx)
#endif
{
#ifndef NO_GPU
 return cuCrealf(cmplx);
#else
#ifdef __cplusplus
 return cmplx.real();
#else
 return cmplx.real;
#endif
#endif
}

#ifndef NO_GPU
__host__ __device__ inline double talshComplex8Real(talshComplex8 cmplx)
#else
inline double talshComplex8Real(talshComplex8 cmplx)
#endif
{
#ifndef NO_GPU
 return cuCreal(cmplx);
#else
#ifdef __cplusplus
 return cmplx.real();
#else
 return cmplx.real;
#endif
#endif
}

#ifndef NO_GPU
__host__ __device__ inline float talshComplex4Imag(talshComplex4 cmplx)
#else
inline float talshComplex4Imag(talshComplex4 cmplx)
#endif
{
#ifndef NO_GPU
 return cuCimagf(cmplx);
#else
#ifdef __cplusplus
 return cmplx.imag();
#else
 return cmplx.imag;
#endif
#endif
}

#ifndef NO_GPU
__host__ __device__ inline double talshComplex8Imag(talshComplex8 cmplx)
#else
inline double talshComplex8Imag(talshComplex8 cmplx)
#endif
{
#ifndef NO_GPU
 return cuCimag(cmplx);
#else
#ifdef __cplusplus
 return cmplx.imag();
#else
 return cmplx.imag;
#endif
#endif
}

#ifndef NO_GPU
__host__ __device__ inline talshComplex4 talshComplex4Conjg(talshComplex4 cmplx)
#else
inline talshComplex4 talshComplex4Conjg(talshComplex4 cmplx)
#endif
{
#ifndef NO_GPU
 return cuConjf(cmplx);
#else
#ifdef __cplusplus
 return std::conj(cmplx);
#else
 talshComplex4 result = {cmplx.real,-cmplx.imag};
 return result;
#endif
#endif
}

#ifndef NO_GPU
__host__ __device__ inline talshComplex8 talshComplex8Conjg(talshComplex8 cmplx)
#else
inline talshComplex8 talshComplex8Conjg(talshComplex8 cmplx)
#endif
{
#ifndef NO_GPU
 return cuConj(cmplx);
#else
#ifdef __cplusplus
 return std::conj(cmplx);
#else
 talshComplex8 result = {cmplx.real,-cmplx.imag};
 return result;
#endif
#endif
}

#ifndef NO_GPU
__host__ __device__ inline float talshComplex4Abs(talshComplex4 cmplx)
#else
inline float talshComplex4Abs(talshComplex4 cmplx)
#endif
{
#ifndef NO_GPU
 return cuCabsf(cmplx);
#else
#ifdef __cplusplus
 return std::abs(cmplx);
#else
 return (float)sqrt((double)((cmplx.real)*(cmplx.real)) + (double)((cmplx.imag)*(cmplx.imag)));
#endif
#endif
}

#ifndef NO_GPU
__host__ __device__ inline double talshComplex8Abs(talshComplex8 cmplx)
#else
inline double talshComplex8Abs(talshComplex8 cmplx)
#endif
{
#ifndef NO_GPU
 return cuCabs(cmplx);
#else
#ifdef __cplusplus
 return std::abs(cmplx);
#else
 return sqrt(((cmplx.real)*(cmplx.real)) + ((cmplx.imag)*(cmplx.imag)));
#endif
#endif
}

#ifndef NO_GPU
__host__ __device__ inline talshComplex4 talshComplex4Add(talshComplex4 x, talshComplex4 y)
#else
inline talshComplex4 talshComplex4Add(talshComplex4 x, talshComplex4 y)
#endif
{
#ifndef NO_GPU
 return cuCaddf(x,y);
#else
#ifdef __cplusplus
 return x+y;
#else
 return talshComplex4Set(x.real+y.real,x.imag+y.imag);
#endif
#endif
}

#ifndef NO_GPU
__host__ __device__ inline talshComplex8 talshComplex8Add(talshComplex8 x, talshComplex8 y)
#else
inline talshComplex8 talshComplex8Add(talshComplex8 x, talshComplex8 y)
#endif
{
#ifndef NO_GPU
 return cuCadd(x,y);
#else
#ifdef __cplusplus
 return x+y;
#else
 return talshComplex8Set(x.real+y.real,x.imag+y.imag);
#endif
#endif
}

#ifndef NO_GPU
__host__ __device__ inline talshComplex4 talshComplex4Sub(talshComplex4 x, talshComplex4 y)
#else
inline talshComplex4 talshComplex4Sub(talshComplex4 x, talshComplex4 y)
#endif
{
#ifndef NO_GPU
 return cuCsubf(x,y);
#else
#ifdef __cplusplus
 return x-y;
#else
 return talshComplex4Set(x.real-y.real,x.imag-y.imag);
#endif
#endif
}

#ifndef NO_GPU
__host__ __device__ inline talshComplex8 talshComplex8Sub(talshComplex8 x, talshComplex8 y)
#else
inline talshComplex8 talshComplex8Sub(talshComplex8 x, talshComplex8 y)
#endif
{
#ifndef NO_GPU
 return cuCsub(x,y);
#else
#ifdef __cplusplus
 return x-y;
#else
 return talshComplex8Set(x.real-y.real,x.imag-y.imag);
#endif
#endif
}

#ifndef NO_GPU
__host__ __device__ inline talshComplex4 talshComplex4Mul(talshComplex4 x, talshComplex4 y)
#else
inline talshComplex4 talshComplex4Mul(talshComplex4 x, talshComplex4 y)
#endif
{
#ifndef NO_GPU
 return cuCmulf(x,y);
#else
#ifdef __cplusplus
 return x*y;
#else
 float rlx = x.real; float imx = x.imag;
 float rly = y.real; float imy = y.imag;
 return talshComplex4Set(rlx*rly-imx*imy,rlx*imy+imx*rly);
#endif
#endif
}

#ifndef NO_GPU
__host__ __device__ inline talshComplex8 talshComplex8Mul(talshComplex8 x, talshComplex8 y)
#else
inline talshComplex8 talshComplex8Mul(talshComplex8 x, talshComplex8 y)
#endif
{
#ifndef NO_GPU
 return cuCmul(x,y);
#else
#ifdef __cplusplus
 return x*y;
#else
 double rlx = x.real; double imx = x.imag;
 double rly = y.real; double imy = y.imag;
 return talshComplex8Set(rlx*rly-imx*imy,rlx*imy+imx*rly);
#endif
#endif
}

#ifndef NO_GPU
__host__ __device__ inline talshComplex4 talshComplex4Div(talshComplex4 x, talshComplex4 y)
#else
inline talshComplex4 talshComplex4Div(talshComplex4 x, talshComplex4 y)
#endif
{
#ifndef NO_GPU
 return cuCdivf(x,y);
#else
#ifdef __cplusplus
 return x/y;
#else
 float rlx = x.real; float imx = x.imag;
 float rly = y.real; float imy = y.imag;
 float dny = 1.0f/(rly*rly + imy*imy);
 return talshComplex4Set((rlx*rly+imx*imy)*dny,(imx*rly-rlx*imy)*dny);
#endif
#endif
}

#ifndef NO_GPU
__host__ __device__ inline talshComplex8 talshComplex8Div(talshComplex8 x, talshComplex8 y)
#else
inline talshComplex8 talshComplex8Div(talshComplex8 x, talshComplex8 y)
#endif
{
#ifndef NO_GPU
 return cuCdiv(x,y);
#else
#ifdef __cplusplus
 return x/y;
#else
 double rlx = x.real; double imx = x.imag;
 double rly = y.real; double imy = y.imag;
 double dny = 1.0/(rly*rly + imy*imy);
 return talshComplex8Set((rlx*rly+imx*imy)*dny,(imx*rly-rlx*imy)*dny);
#endif
#endif
}

#endif /*TALSH_COMPLEX_H_*/
