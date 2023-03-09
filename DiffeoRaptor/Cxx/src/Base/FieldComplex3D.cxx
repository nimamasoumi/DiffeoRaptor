#include "FieldComplex3D.h"

// Returns vector value for each voxel
Vec3DComplex FieldComplex3D::getVal(int x, int y, int z)  const
{
  Vec3DComplex ret(3, complex<float>(0.0, 0.0));
  int index = 3*(xDim*(z*yDim + y) + x);
  ret[0] = data[index];
  ret[1] = data[index+1];
  ret[2] = data[index+2];
  return ret; // reshape according to col
}

// Check boundary when returns value for each voxel
Vec3DComplex FieldComplex3D::getValSafe(int x, int y, int z) const
{
  Vec3DComplex ret(3, complex<float>(0.0, 0.0));
  if (x<0 || x>=xDim || y<0 || y>=yDim || z<0 || z>=zDim)
    return ret;
  else 
    return getVal(x, y, z);
}

//Set value at (x, y, z) location
void FieldComplex3D::setVal(FieldComplex3D& h_o, Vec3DComplex val, int x, int y, int z) const
{
  int index = 3*(xDim*(z*yDim + y) + x);
  h_o.data[index] = val[0];
  h_o.data[index+1] = val[1];
  h_o.data[index+2] = val[2];
}

// get Norm of FieldComplex3D
float FieldComplex3D::Norm() const
{
  float ret=0.0;
  for(int i = 0; i < xDim * yDim * zDim * 3; i++)
       ret += norm(data[i]);
  return (sqrt(ret));
}

// generate FieldComplex3D from random Gauss distribution
void FieldComplex3D::randGauss(float mean, float sigma) const
{  
  float r2, x, y;
  srand (time(NULL));
  for(int i = 0; i < xDim * yDim * zDim * 3; i++)
    {   
      r2 = 0;
      while(r2 > 1.0 || r2 == 0)
	{
	  x = 2.0 * (float)rand() / (RAND_MAX + 1.0) - 1.0;
	  y = 2.0 * (float)rand() / (RAND_MAX + 1.0) - 1.0;
	  r2 = x * x + y * y;
	}
      x *= sigma * sqrt(-2.0 * log(r2) / r2);
      y *= sigma * sqrt(-2.0 * log(r2) / r2);
      data[i] = complex<float>(x+mean, y+mean);
    } 
}

// set the same value for whole field
void FieldComplex3D::initVal(complex<float> val) const
{
  for(int i = 0; i < xDim * yDim * zDim * 3; i++)
    data[i] = val;
}

//Returns a = a1 + a2
Vec3DComplex Add_Vec(Vec3DComplex a1, Vec3DComplex a2)
{
  Vec3DComplex ret(3, complex<float>(0.0, 0.0));
  for (int i = 0; i < 3; i++)
    ret[i] = a1[i] + a2[i];
  return ret;
}

//Returns a = a1 * a2
Vec3DComplex Mul_Vec(Vec3DComplex a1, Vec3DComplex a2)
{
  Vec3DComplex ret(3, complex<float>(0.0, 0.0));
  for (int i = 0; i < 3; i++)
       ret[i] = a1[i] * a2[i];
  return ret;
}

// copy Field value
void Copy_FieldComplex(FieldComplex3D& out, const FieldComplex3D& in)
{
  memcpy(&out.data[0], &in.data[0], out.xDim * out.yDim * out.zDim * 3 * sizeof(complex<float>));
}

// dot product
complex<float> Dotprod(FieldComplex3D& v1, FieldComplex3D& v2)
{
     complex<float> ret(0.0, 0.0);
     for(int i = 0; i < v1.xDim * v1.yDim * v1.zDim * 3; i++)
	  ret += conj(v1.data[i]) * v2.data[i]; 
     return (ret);
}

//Returns v = v + c*v1
void AddI_FieldComplex(FieldComplex3D& v,
		       const FieldComplex3D& v1,
		       float c)
{
     int size = v.xDim * v.yDim * v.zDim * 3;

     for(int i = 0; i < size; i++)
	  v.data[i] += c * v1.data[i];
}

//Returns v = v + c*(v1+v2)
void AddIMul_FieldComplex(FieldComplex3D& v,
			  const FieldComplex3D& v1,
			  const FieldComplex3D& v2,
			  float c)
{
  int size = v.xDim * v.yDim * v.zDim * 3;
  
  for(int i = 0; i < size; i++)
    v.data[i] += c * (v1.data[i] + v2.data[i]);
}

//Returns v = v1 + c*v2
void Add_FieldComplex(FieldComplex3D& v,
		      const FieldComplex3D& v1,
		      const FieldComplex3D& v2,
		      float c)
{
  int size = v.xDim * v.yDim * v.zDim * 3;

  for(int i = 0; i < size; i++)
    v.data[i] = v1.data[i] + c * v2.data[i];
}

//Returns v = v1*v2
void Mul_FieldComplex(FieldComplex3D& v,
		      const FieldComplex3D& v1,
		      const FieldComplex3D& v2)
{
  int size = v.xDim * v.yDim * v.zDim * 3;

  for(int i = 0; i < size; i++)
       v.data[i] = v1.data[i] * v2.data[i];
}

//Returns v = v*v1
void MulI_FieldComplex(FieldComplex3D& v,
		       const FieldComplex3D& v1)
{
  int size = v.xDim * v.yDim * v.zDim * 3;

  for(int i = 0; i < size; i++)
       v.data[i] = v.data[i] * v1.data[i];
}

//Returns v1 = v2 * a
void MulC_FieldComplex(FieldComplex3D& v1,
                       const FieldComplex3D& v2,
		       float c)
{
  int size = v1.xDim * v1.yDim * v1.zDim * 3;

  for(int i = 0; i < size; i++)
    v1.data[i] = c * v2.data[i];
}

//Returns v1 *= a
void MulCI_FieldComplex(FieldComplex3D& v1,
			float c)
{
  int size = v1.xDim * v1.yDim * v1.zDim * 3;

  for(int i = 0; i < size; i++)
    v1.data[i] *= c;
}

// calculate Dv on 3D grid
void Jacobian(FieldComplex3D& JacX,
	      FieldComplex3D& JacY,
	      FieldComplex3D& JacZ,
	      const FieldComplex3D& CDcoeff,
	      const FieldComplex3D& Mat)
{
  int size = CDcoeff.xDim * CDcoeff.yDim * CDcoeff.zDim;
  int index;

  for(int i = 0; i < size; i++)
    {
      index = 3*i;
      JacX.data[index] = CDcoeff.data[index]*Mat.data[index];
      JacX.data[index+1] = CDcoeff.data[index]*Mat.data[index+1];
      JacX.data[index+2] = CDcoeff.data[index]*Mat.data[index+2];

      JacY.data[index] = CDcoeff.data[index+1]*Mat.data[index];
      JacY.data[index+1] = CDcoeff.data[index+1]*Mat.data[index+1];
      JacY.data[index+2] = CDcoeff.data[index+1]*Mat.data[index+2];

      JacZ.data[index] = CDcoeff.data[index+2]*Mat.data[index];
      JacZ.data[index+1] = CDcoeff.data[index+2]*Mat.data[index+1];
      JacZ.data[index+2] = CDcoeff.data[index+2]*Mat.data[index+2];
    }
}

// calculate Dv on 3D grid
void JacobianT(FieldComplex3D& JacX,
	       FieldComplex3D& JacY,
	       FieldComplex3D& JacZ,
	       const FieldComplex3D& CDcoeff,
	       const FieldComplex3D& Mat)
{
  int size = CDcoeff.xDim * CDcoeff.yDim * CDcoeff.zDim;
  int index;
  
  for(int i = 0; i < size; i++)
    {
      index = 3*i;
      JacX.data[index] = CDcoeff.data[index]*Mat.data[index];
      JacX.data[index+1] = CDcoeff.data[index+1]*Mat.data[index];
      JacX.data[index+2] = CDcoeff.data[index+2]*Mat.data[index];
      
      JacY.data[index] = CDcoeff.data[index]*Mat.data[index+1];
      JacY.data[index+1] = CDcoeff.data[index+1]*Mat.data[index+1];
      JacY.data[index+2] = CDcoeff.data[index+2]*Mat.data[index+1];
      
      JacZ.data[index] = CDcoeff.data[index]*Mat.data[index+2];
      JacZ.data[index+1] = CDcoeff.data[index+1]*Mat.data[index+2];
      JacZ.data[index+2] = CDcoeff.data[index+2]*Mat.data[index+2];
    }
}
