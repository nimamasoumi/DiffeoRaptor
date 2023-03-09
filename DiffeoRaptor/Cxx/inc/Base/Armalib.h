#ifndef __Armalib_h
#define __Armalib_h

#include <armadillo>
#include "Field3D.h"

using namespace PyCA;
using namespace arma;

// convert from Field3D to Matrix type (nx1)
void Field2Mat(Mat<float>& mat, Field3D**& v, int nv);

void Field2Mat(Mat<float>& mat, const Field3D& v);

// convert from Matrix type to Field3D
void Mat2Field(Field3D**& v, const Mat<float>& mat, int nv);

void Mat2Field(Field3D& v, const Mat<float>& mat);

#endif
