#ifndef __MPIlib_h
#define __MPIlib_h

#include <mpi.h>
#include "Image3D.h"
#include "Field3D.h"

using namespace PyCA;

// convert Field3D to 1D array type
// for MPI
void field3D2arr(float* arr, const Field3D& v, int len);

// convert 1D array to Field3D type
// for MPI
void arr2Field3D(Field3D& v, float* arr, int len);

#endif

