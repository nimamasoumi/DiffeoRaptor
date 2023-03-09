// Tangent principal component analysis on initial velocity
// This is non-MPI version

#include "Tpca.h"

Tpca::Tpca(Field3D** _v0,
	   const int _nV0,
	   const int _lpower,
	   const float _alpha,
	   const float _gamma)
{
  v0 = _v0;
  nV0 = _nV0;
  lpower = _lpower;
  alpha = _alpha;
  gamma = _gamma;
  nVox = (*v0)->nVox();
  sumV0 = new Field3D((*v0)->grid(), (*v0)->memType());

  // create DiffOper
  diffOper.setAlpha(alpha);
  diffOper.setBeta(0.0);
  diffOper.setGamma(gamma);
  diffOper.setLPow(lpower);
  diffOper.setGrid((*v0)->grid());
}

Tpca::~Tpca()
{
  delete sumV0;
}

// main algorithm for tanget PCA
void Tpca::tangentPCA(Field3D** princomp)
{
  Opers::SetMem(*sumV0, 0.0);
  Mat<float> *storeV0, *LstoreV0, U, V, newU;
  Col<float> s;
  storeV0 = new Mat<float>(nVox*3, nV0);
  LstoreV0 = new Mat<float>(nVox*3, nV0);
  int i;
  // center to mean 0 
  for (i = 0; i < nV0; i++)
    Opers::Add_I(*sumV0, *v0[i]);
  Opers::MulC_I(*sumV0, 1.0/static_cast<float>(nV0));

  for (i = 0; i < nV0; i++)
    {
      Opers::Add_MulC_I(*v0[i], *sumV0, -1.0);
      Field2Mat(*storeV0, v0, nV0);
      diffOper.applyOperator(*v0[i]);
      Field2Mat(*LstoreV0, v0, nV0);
    }
 
  svd(U, s, V, (*storeV0).t()*(*LstoreV0)); // transpose trick
  newU = (*storeV0)*U*sqrt(diagmat(1.0/s));
  Mat2Field(princomp, newU, nV0);

  delete storeV0;
  delete LstoreV0;
}
