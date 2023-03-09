#ifndef __FFTOPER_H__
#define __FFTOPER_H__

#include "FieldComplex3D.h"
#include <iostream>
#include <math.h>
#include <fftw3.h>
#include <vector>
#include "Field3D.h" // including functions from pyca
#include "Image3D.h"

using namespace PyCA;

class FftOper
{
public:

     FftOper(const float _alpha,
	     const float _gamma,
	     const int _lpow,
	     const GridInfo _grid,
	     const int _truncX,
	     const int _truncY,
	     const int _truncZ);

     ~FftOper();
    
     void FourierCoefficient();

     void spatial2fourier(FieldComplex3D& vtrunc,
			  const Field3D& vf);

     void fourier2spatial(Field3D& vf,
			  const FieldComplex3D& vtrunc);

     void spatial2fourier_F(float* xf,
			    float* yf,
			    float* zf,
			    const Field3D& vf);
       
     void fourier2spatial_addH(Field3D& vf,
			       const FieldComplex3D& vtrunc,
			       const float* idx,
			       const float* idy,
			       const float* idz);
       
     void ConvolveComplexFFT(FieldComplex3D& h_o_ori,
			     const int flag,
			     const FieldComplex3D& h_ix_ori, 
			     const FieldComplex3D& h_iy_ori,
			     const FieldComplex3D& h_iz_ori, 
			     const FieldComplex3D& h_kernel_ori);
     
     void CorrComplexFFT(FieldComplex3D& h_o_ori,
			 const FieldComplex3D& h_i_ori, 
			 const FieldComplex3D& h_kernel_ori,
			 const FieldComplex3D& CDcoeff);
     
     FieldComplex3D *Lcoeff, *Kcoeff, *CDcoeff;
     GridInfo grid;
     float alpha, gamma;
     int lpow;
     int truncX, truncY, truncZ;
     int fsx, fsy, fsz, fsxFFT;
     
protected:
     void fftshiftMat();
     void PointwiseMultiply(float *im, float *in, int flag);
     int size, fsize, padsize;
     int *fftLoc, *corrLoc, *conjLoc;
     float *scratch, *data;
     double spx, spy, spz;
     int beginX, beginY, beginZ, endX, endY, endZ;
     fftwf_plan fftwForwardPlan, fftwBackwardPlan;

     int padX, padY, padZ;
     float *scratchI, *scratchK, *dataI, *dataK;
     float *scratchX, *scratchY, *scratchZ;
     fftwf_plan fftwFwdPlanI, fftwFwdPlanK, fftwBwdPlan;
     FieldComplex3D *h_ix, *h_iy, *h_iz, *h_kernel, *h_i;
};

#endif // __FFTOPER_H__
