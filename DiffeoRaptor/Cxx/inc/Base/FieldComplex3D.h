#ifndef __FIELDCOMPLEX3D_H__
#define __FIELDCOMPLEX3D_H__

#include <complex>
#include <vector>
#include <iostream>
#include <math.h>
#include <fftw3.h>
#include "cstring"
#include "stdlib.h"

using namespace std;

// deltaX, deltaY, deltaZ are assumed to be 1.0
typedef vector< complex<float> > Vec3DComplex;

struct FieldComplex3D
{  
  int xDim, yDim, zDim;
  complex<float>* data;

  FieldComplex3D()
  {
    xDim = 0;
    yDim = 0;
    zDim = 0;
    data = NULL;
  }

  FieldComplex3D(int _xDim, int _yDim, int _zDim)
  {
    xDim = _xDim;
    yDim = _yDim;
    zDim = _zDim;

    data = new complex<float>[xDim * yDim * zDim * 3];
    for(int i = 0; i < xDim * yDim * zDim * 3; i++)
      data[i] = complex<float>(0.0, 0.0);
  }

  ~FieldComplex3D()
  {
    delete [] data;
  }
  Vec3DComplex getVal(int x, int y, int z) const;
  Vec3DComplex getValSafe(int x, int y, int z) const;
  void setVal(FieldComplex3D& h_o, Vec3DComplex val, int x, int y, int z) const;
  float Norm() const;
  void randGauss(float mean, float sigma) const;
  void initVal(complex<float> val) const;
};

//Returns a = a1 + a2
Vec3DComplex Add_Vec(Vec3DComplex a1, Vec3DComplex a2);

//Returns a = a1 * a2
Vec3DComplex Mul_Vec(Vec3DComplex a1, Vec3DComplex a2);

// copy Field
void Copy_FieldComplex(FieldComplex3D& out, const FieldComplex3D& in);

// dot product
complex<float> Dotprod(FieldComplex3D& v1, FieldComplex3D& v2);

//Returns v = v + c*v1
void AddI_FieldComplex(FieldComplex3D& v,
		       const FieldComplex3D& v1,
		       float c);

//Returns v = v1 + c*v2
void Add_FieldComplex(FieldComplex3D& v,
		      const FieldComplex3D& v1,
		      const FieldComplex3D& v2,
		      float c);

//Returns v = v + c*(v1+v2)
void AddIMul_FieldComplex(FieldComplex3D& v,
			  const FieldComplex3D& v1,
			  const FieldComplex3D& v2,
			  float c);

//Returns v = v1*v2
void Mul_FieldComplex(FieldComplex3D& v,
		      const FieldComplex3D& v1,
		      const FieldComplex3D& v2);

//Returns v = v*v1
void MulI_FieldComplex(FieldComplex3D& v,
		       const FieldComplex3D& v1);

//Returns v1 = v2 * a
void MulC_FieldComplex(FieldComplex3D& v1,
                       const FieldComplex3D& v2,
		       float c);

//Returns v1 *= a
void MulCI_FieldComplex(FieldComplex3D& v1,
			float c);

// calculate transpose jacobian matrix on 3D grid
void Jacobian(FieldComplex3D& JacX,
	      FieldComplex3D& JacY,
	      FieldComplex3D& JacZ,
	      const FieldComplex3D& CDcoeff,
	      const FieldComplex3D& Mat);

void JacobianT(FieldComplex3D& JacX,
	       FieldComplex3D& JacY,
	       FieldComplex3D& JacZ,
	       const FieldComplex3D& CDcoeff,
	       const FieldComplex3D& Mat);	 
  
#endif // __FIELDCOMPLEX3D_H__
