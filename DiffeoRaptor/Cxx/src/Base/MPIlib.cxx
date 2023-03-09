#include "MPIlib.h"

// convert Field3D to 1D array type
// for MPI
void field3D2arr(float* arr, const Field3D& v, int len) 
{
  for (int i = 0; i < len; i++)
    {
      arr[3*i] = v.x[i];
      arr[3*i+1] = v.y[i];
      arr[3*i+2] = v.z[i];
    }
}

// convert 1D array to Field3D type
// for MPI
void arr2Field3D(Field3D& v, float* arr, int len) 
{
  for (int i = 0; i < len; i++)
    {
      v.x[i] = arr[3*i];
      v.y[i] = arr[3*i+1];
      v.z[i] = arr[3*i+2];
    }
}
