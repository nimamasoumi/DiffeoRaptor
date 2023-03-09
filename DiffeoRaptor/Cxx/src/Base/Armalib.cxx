#include "Armalib.h"

// convert from Field3D to Matrix type (nx1)
void Field2Mat(Mat<float>& mat, const Field3D& v)
{
  int nVox = v.nVox();
  memcpy(&mat[0], v.getX(), nVox * sizeof(float));
  memcpy(&mat[nVox], v.getY(), nVox * sizeof(float));
  memcpy(&mat[2*nVox], v.getZ(), nVox * sizeof(float));
}

void Field2Mat(Mat<float>& mat, Field3D**& v, int nv)
{
  for (int i = 0; i < nv; i++)
    {
      Mat<float> temp(mat.colptr(i), mat.n_rows, 1, false, true);
      Field2Mat(temp, *v[i]);
    }
}

// convert from Matrix type to Field3D
void Mat2Field(Field3D& v, const Mat<float>& mat)
{
  int nVox = v.nVox();
  memcpy(v.getX(), &mat[0], nVox * sizeof(float));
  memcpy(v.getY(), &mat[nVox], nVox * sizeof(float));
  memcpy(v.getZ(), &mat[2*nVox], nVox * sizeof(float));
}

void Mat2Field(Field3D**& v, const Mat<float>& mat, int nv)
{
  for (int i = 0; i < nv; i++)
    Mat2Field(*v[i], mat.col(i));
}
