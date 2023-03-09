/*==============================================================================
  File: DiffeoRaptor.h

  ============================================================================= */
#ifndef __DiffeoRaptor_h
#define __DiffeoRaptor_h

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

class DiffeoRaptor
{

public:

     DiffeoRaptor(FftOper* _fftOper,
                  const MemoryType _mType,
                  const int _numTimeSteps,
                  bool _doRK4 = false);

     ~DiffeoRaptor();

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

     // get gradient term
     void Gradient(FieldComplex3D& v0,
		   FieldComplex3D& gradv);

     // ImageMatching includes extra work for atlas building
     // fixed stepsize for gradient descent
     void ImageMatching(FieldComplex3D& v0,
			FieldComplex3D& gradv,
			int maxIter,
			float stepSizeGD);

     // normalizing images
     void US_normalize(Image3D& In,
                       const Image3D* I);

     // Functions manager for calculating the gradient of CR and CR
     void H_CRc(double& CR,
                 Image3D &GR,
                 const Image3D *Iref,
                 const Image3D *Imov,
                 const int bins,
                 const int n_pts,
                 double& deNUM);

     // Total CR and gradient calculator
     void CRtGrad(double& CR,
                    Image3D &dd,
                    const Image3D *Iref,
                    const Image3D *Imov);
     // Masking function
     void MR_grad_US_mask_v3(std::vector<std::vector<double>>& idx2yxz,
                             const Image3D* If,
                             const Image3D* Im);
     // Ratio Image
     void calc_ratio(Image3D& Irat,
                     const Image3D* Irefg,
                     const Image3D* Imovg,
                     const double rr);

     void outlierThGen(double &outlier_TH,
                       const std::vector<double> pointTarget,
                       const Image3D* Iref,
                       const Image3D* Imov,
                       const double rr);

     Image3D *I0, *I1, *deformIm, *splatI, *splatOnes;
     float VEnergy, IEnergy, TotalEnergy, sigma, corrFactor, outlier_TH;
     Field3D *phiinv;
     std::vector<double> collectedCR, collectedE, collectedMSE;
     std::vector<std::vector<double>> Ifg_idx2yxz;
     int TotalPatch, iterCount, rrPatchSize;

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
