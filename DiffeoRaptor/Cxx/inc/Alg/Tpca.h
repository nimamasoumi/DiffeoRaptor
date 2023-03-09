/*==============================================================================
  File: Tpca.h

  ============================================================================= */
#ifndef __Tpca_h
#define __Tpca_h

#include "IOpers.h"
#include "FOpers.h"
#include "IFOpers.h"
#include "HFOpers.h"
#include "Armalib.h"
#include "FluidKernelFFT.h"

class Tpca
{
 protected:

  Field3D** v0;
  Field3D* sumV0;
  int nV0, nVox, lpower;
  float alpha, gamma;
  FluidKernelFFT<EXEC_CPU> diffOper;

 public:

  Tpca(Field3D** _v0,
       const int _nV0,
       const int _lpower,
       const float _alpha,
       const float _gamma);
  ~Tpca();

  void tangentPCA(Field3D** princomp);
};

#endif
