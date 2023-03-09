/*==============================================================================
  File: GeodesicShooting.h

  ============================================================================= */
#ifndef __GeodesicShooting_h
#define __GeodesicShooting_h

#include "FftOper.h"
#include "FieldComplex3D.h"
#include "ITKFileIO.h"
#include "IOpers.h"
#include "FOpers.h"
#include "IFOpers.h"
#include "HFOpers.h"
#include "Reduction.h"
#include "FluidKernelFFT.h"

using namespace PyCA;
using namespace std;

class GeodesicShooting
{

public:

     GeodesicShooting(FftOper* _fftOper,
		      const MemoryType _mType,
          const int _numTimeSteps,
          bool _doRK4 = false);

     ~GeodesicShooting();

     // ad operator 
     void ad(FieldComplex3D& advw,
	     const FieldComplex3D& v, 
	     const FieldComplex3D& w);

     // adTranspose
     // K*(ConvolveComplex(CD*v, L*w)+CD*ConvolveComplex(L*w, v))
     void adTranspose(FieldComplex3D& adTransvw,
		      const FieldComplex3D& v, 
		      const FieldComplex3D& w);

     // forward integration of v
     void fwdIntegrateV(const FieldComplex3D& v0);
     
     // forward integration
     void fwdIntegration(FieldComplex3D& v0,
			 FieldComplex3D& fwd_gradvfft);

     // update adjoint variable
     // for backward integration
     void bwdUpdate(FieldComplex3D& dvadj,
		    const FieldComplex3D& vadj,
		    const FieldComplex3D& ad, 
		    const FieldComplex3D& adTrans);
     
     // backward integration, reduced adjoint jacobi fields
     void bwdIntegration(FieldComplex3D& dvadj,
			 FieldComplex3D& vadj); // forward gradient

     // parallel translation helper function
     void parallelTranslateF(FieldComplex3D& fvw,
        const FieldComplex3D& v,
        const FieldComplex3D& w);
     
     // parallel translation
     void parallelTranslate(FieldComplex3D& parTransvw,
        const FieldComplex3D& v,
        const FieldComplex3D& w);

     // get gradient term
     void Gradient(FieldComplex3D& v0,
		   FieldComplex3D& gradv);

     // ImageMatching includes extra work for atlas building
     // fixed stepsize for gradient descent
     void ImageMatching(FieldComplex3D& v0,
			FieldComplex3D& gradv,
			int maxIter,
			float stepSizeGD);

     // HMC sampling for initial velocities
     void HmcV(FieldComplex3D** vSamples,
	       FieldComplex3D* current_q, // initial velocity
	       int num_samples,
	       int burn_samples,
	       int LeapfrogSteps,
	       float HmcStepsInitial);

     Image3D *I0, *I1, *deformIm, *splatI, *splatOnes;
     float VEnergy, IEnergy, TotalEnergy, sigma;
     Field3D *phiinv;

protected:
     FftOper* fftOper;
     FieldComplex3D** storeV;
     FieldComplex3D *imMatchGradient;
     FieldComplex3D *fwd_gradvfft;
     Field3D *identity;

     Image3D *residualIm;
     int numTimeSteps;

     int xDim, yDim, zDim; // truncated dimension
     float dt;
     float *idxf, *idyf, *idzf;
     FieldComplex3D *scratch1, *scratch2, *scratch3; // scratch variables
     FieldComplex3D *adScratch1, *adScratch2; // scratch variables
     FieldComplex3D *JacX, *JacY, *JacZ; // scratch variables
     Field3D *scratchV1, *scratchV2;
     GridInfo grid;
     MemoryType mType;
     bool doRK4;
};

#endif
