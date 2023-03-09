#include "GeodesicShooting.h"
#include <time.h>

GeodesicShooting::GeodesicShooting(FftOper* _fftOper,
				   const MemoryType _mType,
           const int _numTimeSteps,
           bool _doRK4)
{
     fftOper = _fftOper;
     xDim = fftOper->truncX;
     yDim = fftOper->truncY;
     zDim = fftOper->truncZ;
     mType = _mType;
     grid = fftOper->grid;
     numTimeSteps = _numTimeSteps;
     dt = 1.0 / numTimeSteps;
     doRK4 = _doRK4;

     storeV = new FieldComplex3D* [numTimeSteps+1];
     for (int i = 0; i <= numTimeSteps; i++)
	  storeV[i] = new FieldComplex3D(xDim, yDim, zDim);

     adScratch1 = new FieldComplex3D(xDim, yDim, zDim);
     adScratch2 = new FieldComplex3D(xDim, yDim, zDim);
     scratch1 = new FieldComplex3D(xDim, yDim, zDim);
     scratch2 = new FieldComplex3D(xDim, yDim, zDim);
     scratch3 = new FieldComplex3D(xDim, yDim, zDim);
     JacX = new FieldComplex3D(xDim, yDim, zDim);
     JacY = new FieldComplex3D(xDim, yDim, zDim);
     JacZ = new FieldComplex3D(xDim, yDim, zDim);

     imMatchGradient = new FieldComplex3D(xDim, yDim, zDim);
     fwd_gradvfft = new FieldComplex3D(xDim, yDim, zDim);

     scratchV1 = new Field3D(grid, mType);
     scratchV2 = new Field3D(grid,  mType);
     phiinv = new Field3D(grid, mType);

     I0 = new Image3D(grid, mType);
     I1 = new Image3D(grid, mType);
     deformIm = new Image3D(grid, mType);
     splatI = new Image3D(grid, mType);
     splatOnes = new Image3D(grid, mType);
     residualIm = new Image3D(grid, mType);

     // id matrix in Fourier domain can be preset
     identity = new Field3D(grid, mType);
     Opers::SetToIdentity(*identity);
     idxf = new float[2 * fftOper->fsxFFT * fftOper->fsy * fftOper->fsz];
     idyf = new float[2 * fftOper->fsxFFT * fftOper->fsy * fftOper->fsz];
     idzf = new float[2 * fftOper->fsxFFT * fftOper->fsy * fftOper->fsz];
     fftOper->spatial2fourier_F(idxf, idyf, idzf, *identity); 
}

GeodesicShooting::~GeodesicShooting()
{
     delete I0; I0 = NULL;
     delete I1; I1 = NULL;
     delete residualIm; residualIm = NULL;
     delete deformIm; deformIm = NULL;
     delete splatI; splatI = NULL;
     delete splatOnes; splatOnes = NULL;
     delete scratch1; scratch1 = NULL;
     delete scratch2; scratch2 = NULL;
     delete scratch3; scratch3 = NULL;
     delete adScratch1; adScratch1 = NULL;
     delete adScratch2; adScratch2 = NULL;
     delete scratchV1; scratchV1 = NULL;
     delete scratchV2; scratchV2 = NULL;
     delete phiinv; phiinv = NULL;
     delete imMatchGradient; imMatchGradient = NULL;
     delete fwd_gradvfft; fwd_gradvfft = NULL;
     delete JacX; JacX = NULL;
     delete JacY; JacY = NULL;
     delete JacZ; JacZ = NULL;
     delete identity; identity = NULL;
     delete idxf; idxf = NULL;
     delete idyf; idyf = NULL;
     delete idzf; idzf = NULL;
     
     for (int i = 0; i <= numTimeSteps; i++)
     {delete storeV[i]; storeV[i] = NULL;}
	       delete [] storeV; storeV = NULL;
}

// ad operator
// CD: central difference
// ConvolveComplexFFT(CD*v, w)-ConvolveComplexFFT(CD*w, v))
void GeodesicShooting::ad(FieldComplex3D& advw,
			  const FieldComplex3D& v,
			  const FieldComplex3D& w)
{
     Jacobian(*JacX, *JacY, *JacZ, *(fftOper->CDcoeff), v); 
     fftOper->ConvolveComplexFFT(advw, 0, *JacX, *JacY, *JacZ, w);

     Jacobian(*JacX, *JacY, *JacZ, *(fftOper->CDcoeff), w); 
     fftOper->ConvolveComplexFFT(*adScratch1, 0, *JacX, *JacY, *JacZ, v);
    
     AddI_FieldComplex(advw, *adScratch1, -1.0);
}

// adTranspose
// spatial domain: K(Dv^T Lw + div(L*w x v))
// K*(CorrComplexFFT(CD*v^T, L*w) + TensorCorr(L*w, v) * D)
void GeodesicShooting::adTranspose(FieldComplex3D& adTransvw,
				   const FieldComplex3D& v,
				   const FieldComplex3D& w)
{
     Mul_FieldComplex(*adScratch1, *(fftOper->Lcoeff), w);

     JacobianT(*JacX, *JacY, *JacZ, *(fftOper->CDcoeff), v); 
     fftOper->ConvolveComplexFFT(adTransvw, 1, *JacX, *JacY, *JacZ, *adScratch1);

     fftOper->CorrComplexFFT(*adScratch2, v, *adScratch1, *(fftOper->CDcoeff));
     AddI_FieldComplex(adTransvw, *adScratch2, 1.0);

     MulI_FieldComplex(adTransvw, *(fftOper->Kcoeff));
} 


void GeodesicShooting::fwdIntegrateV(const FieldComplex3D& v0)
{
   Copy_FieldComplex(*storeV[0], v0);

   if (doRK4)
   {
     for (int i = 1; i <= numTimeSteps; i++)
     {
       // v1 = v0 - (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
       // k1
       adTranspose(*scratch1, *storeV[i-1], *storeV[i-1]);
       // partially update v1 = v0 - (dt/6)*k1
       Copy_FieldComplex(*storeV[i], *storeV[i-1]);
       AddI_FieldComplex(*storeV[i], *scratch1, -dt / 6.0);

       // k2
       Copy_FieldComplex(*scratch3, *storeV[i-1]);
       AddI_FieldComplex(*scratch3, *scratch1, -0.5*dt);
       adTranspose(*scratch2, *scratch3, *scratch3);
       // partially update v1 = v1 - (dt/3)*k2
       AddI_FieldComplex(*storeV[i], *scratch2, -dt / 3.0);

       // k3 (stored in scratch1)
       Copy_FieldComplex(*scratch3, *storeV[i-1]);
       AddI_FieldComplex(*scratch3, *scratch2, -0.5*dt);
       adTranspose(*scratch1, *scratch3, *scratch3);
       // partially update v1 = v1 - (dt/3)*k3
       AddI_FieldComplex(*storeV[i], *scratch1, -dt / 3.0);

       // k4 (stored in scratch2)
       Copy_FieldComplex(*scratch3, *storeV[i-1]);
       AddI_FieldComplex(*scratch3, *scratch1, -dt);
       adTranspose(*scratch2, *scratch3, *scratch3);
       // finish updating v1 = v1 - (dt/6)*k4
       AddI_FieldComplex(*storeV[i], *scratch2, -dt / 6.0);
     }
   }
   else
   {
     for (int i = 1; i <= numTimeSteps; i++)
     {
	    // v0 = v0 - dt * adTranspose(v0, v0)
      Copy_FieldComplex(*storeV[i], *storeV[i-1]);
	    adTranspose(*scratch1, *storeV[i], *storeV[i]);
	    AddI_FieldComplex(*storeV[i], *scratch1, -dt);
     }
   }

}
// forward integration
// Euler integration
// calculate phiinv by integrating backward
void GeodesicShooting::fwdIntegration(FieldComplex3D& v0,
				      FieldComplex3D& fwd_gradvfft)
{
     fwdIntegrateV(v0);
     Copy_FieldComplex(v0, *storeV[numTimeSteps]); // done because prev version was updating v0

     // Opers::SetToIdentity(*phiinv);

     // for (int i = 0; i < numTimeSteps; i++) // generate phi^{-1} under left invariant metric
     // {
     // 	  fftOper->fourier2spatial(*scratchV1, *storeV[i]);
     // 	  Opers::JacobianXY(*scratchV2, *phiinv, *scratchV1); // scratchV2: D(phiinv) v_t
     // 	  Opers::Add_MulC_I(*phiinv, *scratchV2, -dt);
     // }

     scratch1->initVal(complex<float>(0.0, 0.0)); // displacement field
     for (int i = 0; i < numTimeSteps; i++) // generate phi^{-1} under left invariant metric
     {
       Jacobian(*JacX, *JacY, *JacZ, *(fftOper->CDcoeff), *scratch1); 
       fftOper->ConvolveComplexFFT(*scratch2, 0, *JacX, *JacY, *JacZ, *storeV[i]);
       AddIMul_FieldComplex(*scratch1, *scratch2, *storeV[i], -dt);
     }

     fftOper->fourier2spatial_addH(*phiinv, *scratch1, idxf, idyf, idzf);

     // calculate gradient
     Opers::ApplyH(*deformIm, *I0, *phiinv, BACKGROUND_STRATEGY_WRAP);
     Opers::Gradient(*scratchV1, *deformIm, DIFF_CENTRAL, BC_WRAP);
     Opers::Sub(*residualIm, *deformIm, *I1);
     Opers::MulMulC_I(*scratchV1, *residualIm, -1.0);
     fftOper->spatial2fourier(fwd_gradvfft, *scratchV1);

     MulI_FieldComplex(fwd_gradvfft, *(fftOper->Kcoeff));
}

// update adjoint variable
// for backward integration
void GeodesicShooting::bwdUpdate(FieldComplex3D& dvadj,
				 const FieldComplex3D& vadj,
				 const FieldComplex3D& ad, 
				 const FieldComplex3D& adTrans)
{
     for (int i = 0; i < xDim * yDim * zDim * 3; i++)
	  dvadj.data[i] += dt * (vadj.data[i] - ad.data[i] + adTrans.data[i]);
}

// backward integration, reduced adjoint jacobi fields
// Euler integration
void GeodesicShooting::bwdIntegration(FieldComplex3D& dvadj,
				      FieldComplex3D& vadj) // forward gradient
{
     dvadj.initVal(complex<float>(0.0, 0.0)); // backward to t=0
     for (int i = numTimeSteps; i > 0; i--) // reduced adjoint jacobi fields
     {
       // dvadj = dvadj - dt * (-vadj + ad(v, dvadj) - adTranspose(dvadj, v));
       // dvadj = dvadj + dt * (vadj - ad(v, dvadj) + adTranspose(dvadj, v));
       // vadj = vadj + dt * adTranspose(v, vadj);

       ad(*scratch1, *storeV[i], dvadj);
       adTranspose(*scratch2, dvadj, *storeV[i]);
       bwdUpdate(dvadj, vadj, *scratch1, *scratch2);
       
       adTranspose(*scratch1, *storeV[i], vadj);
       AddI_FieldComplex(vadj, *scratch1, dt);
     }
}

// parallel translation function (compute f for either RK4 or Euler)
void GeodesicShooting::parallelTranslateF(FieldComplex3D& fvw,
                                          const FieldComplex3D& v,
                                          const FieldComplex3D& w)
{
  // compute 1/2{adTranspose(v,w)+adTranspose(w,v)-ad(v,w)}
  adTranspose(fvw, v, w);
  MulCI_FieldComplex(fvw, 0.5f);
  adTranspose(*scratch3, w, v);
  AddI_FieldComplex(fvw, *scratch3, 0.5f);
  ad(*scratch3, v, w);
  AddI_FieldComplex(fvw, *scratch3, -0.5f);
}

// parallel translation
void GeodesicShooting::parallelTranslate(FieldComplex3D& parTransvw,
                                         const FieldComplex3D& v,
                                         const FieldComplex3D& w)
{
  fwdIntegrateV(v);
  Copy_FieldComplex(parTransvw, w);

  if (doRK4)
  {
    FieldComplex3D* curstep = new FieldComplex3D(xDim, yDim, zDim);
    Copy_FieldComplex(*curstep, parTransvw);
    for (int i = 0; i < numTimeSteps; i++)
    {
      // w1 = w0 - (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
      // dw/dt=-1/2{adTranspose(v,w)+adTranspose(w,v)-ad(v,w)}
      // k1
      parallelTranslateF(*scratch1, *storeV[i], parTransvw);
      // partially update w1 = w0-(dt/6)*k1
      AddI_FieldComplex(*curstep, *scratch1, -dt / 6.0);

      // k2
      Copy_FieldComplex(*scratch2, parTransvw);
      AddI_FieldComplex(*scratch2, *scratch1, -0.5*dt);
      parallelTranslateF(*scratch1, *storeV[i], *scratch2);
      // partially update w1 = w1-(dt/3)*k2
      AddI_FieldComplex(*curstep, *scratch1, -dt / 3.0);

      // k3
      Copy_FieldComplex(*scratch2, parTransvw);
      AddI_FieldComplex(*scratch2, *scratch1, -0.5*dt);
      parallelTranslateF(*scratch1, *storeV[i], *scratch2);
      // partially update w1 = w1-(dt/3)*k3
      AddI_FieldComplex(*curstep, *scratch1, -dt / 3.0);

      // k4
      Copy_FieldComplex(*scratch2, parTransvw);
      AddI_FieldComplex(*scratch2, *scratch1, -dt);
      parallelTranslateF(*scratch1, *storeV[i], *scratch2);
      // partially update w1 = w1-(dt/6)*k4
      AddI_FieldComplex(*curstep, *scratch1, -dt / 6.0);

      Copy_FieldComplex(parTransvw, *curstep);

      // // // display norms for debug purposes
      // // // TODO should be able to compute the norm in the fourier domain
      // // float normv, normw, normvw;
      // // Mul_FieldComplex(*scratch1, *(fftOper->Lcoeff), *storeV[i]);
      // // fftOper->fourier2spatial(*scratchV1, *scratch1);
      // // fftOper->fourier2spatial(*scratchV2, *storeV[i]);
      // // Opers::Dot(normv, *scratchV1, *scratchV2);

      // // fftOper->fourier2spatial(*scratchV2, parTransvw);
      // // Opers::Dot(normvw, *scratchV1, *scratchV2);

      // // Mul_FieldComplex(*scratch1, *(fftOper->Lcoeff), parTransvw);
      // // fftOper->fourier2spatial(*scratchV1, *scratch1);
      // // Opers::Dot(normw, *scratchV1, *scratchV2);

      // // typedef std::numeric_limits< float > flt;
      // // std::cout.precision(flt::max_digits10);
      // // std::cout << "Step " << i
      // //           << ", Lvv: " << fixed << normv
      // //           << ", Lvw: " << fixed << normvw
      // //           << ", Lww: " << fixed << normw << std::endl;
    }
    delete curstep;
  }
  else
  {
    for (int i = 0; i < numTimeSteps; i++)
    {
      parallelTranslateF(*scratch1, *storeV[i], parTransvw);
      AddI_FieldComplex(parTransvw, *scratch1, -dt);

      // // // display norms for debug purposes
      // // // TODO should be able to compute the norm in the fourier domain
      // // float normv, normw, normvw;
      // // Mul_FieldComplex(*scratch1, *(fftOper->Lcoeff), *storeV[i]);
      // // fftOper->fourier2spatial(*scratchV1, *scratch1);
      // // fftOper->fourier2spatial(*scratchV2, *storeV[i]);
      // // Opers::Dot(normv, *scratchV1, *scratchV2);

      // // fftOper->fourier2spatial(*scratchV2, parTransvw);
      // // Opers::Dot(normvw, *scratchV1, *scratchV2);

      // // Mul_FieldComplex(*scratch1, *(fftOper->Lcoeff), parTransvw);
      // // fftOper->fourier2spatial(*scratchV1, *scratch1);
      // // Opers::Dot(normw, *scratchV1, *scratchV2);

      // // std::cout << "Step " << i
      // //           << ", Lvv: " << normv
      // //           << ", Lvw: " << normvw
      // //           << ", Lww: " << normw << std::endl;
    }
  }
}

// gradient term for image matching
void GeodesicShooting::Gradient(FieldComplex3D& v0,
				FieldComplex3D& gradv)
{
     Copy_FieldComplex(gradv, v0);

     // calculate Venergy
     Mul_FieldComplex(*scratch1, *(fftOper->Lcoeff), v0);
     fftOper->fourier2spatial(*scratchV1, *scratch1);
     fftOper->fourier2spatial(*scratchV2, v0);
     Opers::Dot(VEnergy, *scratchV1, *scratchV2);
     VEnergy *= 0.5;

     // forward integration
     fwdIntegration(v0, *fwd_gradvfft);
 
     // backward integration
     // imMatchGradient should be initialized as zero before passing
     bwdIntegration(*imMatchGradient, *fwd_gradvfft);

     // calculate Ienergy
     Opers::Sum2(IEnergy, *residualIm);
     TotalEnergy = VEnergy + 0.5/(sigma*sigma)*IEnergy;

     // gradient term
     // gradv = v0 + \tilde{v};
     AddI_FieldComplex(gradv, *imMatchGradient, 1.0/(sigma*sigma));
}

// gradient descent strategy for image matching
// adaptive stepsize
void GeodesicShooting::ImageMatching(FieldComplex3D& v0,
				     FieldComplex3D& gradv,
				     int maxIter,
				     float stepSizeGD)
{
     float epsilon = 1.0e-10;
     float prevEnergy = 1e20;

     FieldComplex3D *currV0 = new FieldComplex3D(xDim, yDim, zDim);
     FieldComplex3D *prevV0 = new FieldComplex3D(xDim, yDim, zDim);
     FieldComplex3D *newV0 = new FieldComplex3D(xDim, yDim, zDim);
     Copy_FieldComplex(*prevV0, v0);

     for (int i = 0; i < maxIter; i++)
     {
      printf("\n %i th iteration!\n",i);
	  Copy_FieldComplex(*currV0, v0);

	  Gradient(v0, gradv);

	  // newV0 = v0 - stepSizeGD*gradv;

	  Add_FieldComplex(*newV0, *currV0, gradv, -stepSizeGD);

	  if (TotalEnergy > prevEnergy)
	  {
	       stepSizeGD *= 0.8;
	       Copy_FieldComplex(v0, *prevV0);
	  }
	  else
	  {
	       Copy_FieldComplex(*prevV0, *currV0);
	       Copy_FieldComplex(v0, *newV0);
	       prevEnergy = TotalEnergy;
	  }
	      
	  if (gradv.Norm() < epsilon)
	       break;
     }

     Opers::SetMem(*residualIm, 1.0);
     Opers::Splat(*splatI, *phiinv, *I1, BACKGROUND_STRATEGY_WRAP);
     Opers::Splat(*splatOnes, *phiinv, *residualIm, BACKGROUND_STRATEGY_WRAP);

     delete currV0; currV0 = NULL;
     delete prevV0; prevV0 = NULL;
     delete newV0; newV0 = NULL;
}

// HMC sampling for initial velocities
// p: auxiliary random variable
// q: sample of initial velocity
// current_U: potential energy
// current_K: kinetic energy
void GeodesicShooting::HmcV(FieldComplex3D** vSamples,
			    FieldComplex3D* current_q, // initial velocity
			    int num_samples,
			    int burn_samples,
			    int LeapfrogSteps,
			    float HmcStepSize)
{
     float current_U, current_K, proposed_U, proposed_K;
     FieldComplex3D* q = new FieldComplex3D(xDim, yDim, zDim);
     FieldComplex3D* p = new FieldComplex3D(xDim, yDim, zDim);
     FieldComplex3D* gradvTemp = new FieldComplex3D(xDim, yDim, zDim);
     FieldComplex3D* qTemp = new FieldComplex3D(xDim, yDim, zDim);

     srand(time(NULL));
     int acc_sample=0, acc=0;

     while (acc_sample < num_samples)
     {
	  Copy_FieldComplex(*q, *current_q);

	  p->randGauss(0.0, 1.0); // generate random Gaussian
	  MulI_FieldComplex(*p, *(fftOper->Kcoeff)); // smooth p in frenquency domain
	  current_K = 0.5 * pow(p->Norm(), 2);

	  // first half step for p
	  // p = p - 0.5*HmcStepSize*gradv
	  Copy_FieldComplex(*qTemp, *q);
	  Gradient(*qTemp, *gradvTemp);
	  AddI_FieldComplex(*p, *gradvTemp, -0.5*HmcStepSize);
	  current_U = TotalEnergy;

	  // full step for p and q
	  // q = q + HmcStepSize*gradp
	  for (int i = 0; i < LeapfrogSteps-1; i++)
	  {
	       AddI_FieldComplex(*q, *p, HmcStepSize);

	       Copy_FieldComplex(*qTemp, *q);
	       Gradient(*qTemp, *gradvTemp);
	       AddI_FieldComplex(*p, *gradvTemp, -1.0*HmcStepSize);
	  }
	  AddI_FieldComplex(*q, *p, HmcStepSize);

	  // final half step for p
	  Copy_FieldComplex(*qTemp, *q);
	  Gradient(*qTemp, *gradvTemp);
	  AddI_FieldComplex(*p, *gradvTemp, -0.5*HmcStepSize);
	  proposed_U = TotalEnergy;
	  proposed_K = 0.5 * pow(p->Norm(), 2);

	  if ((float)rand() / RAND_MAX < exp(current_U-proposed_U+current_K-proposed_K))
	  {
	       Copy_FieldComplex(*current_q, *q);
	       // std::cout << "current_U, " << current_U << ",";
	       // std::cout << "current_K, " << current_K << ",";
	       // std::cout << "proposed_U, " << proposed_U << ",";
	       // std::cout << "proposed_K, " << proposed_K << ",";
	       // std::cout <<  exp((current_U-proposed_U+current_K-proposed_K)) << std::endl;
	  }

	  acc += 1;
	  if (acc > burn_samples) // save samples after burn in
	  {
	       if ((float)rand() / RAND_MAX < exp(current_U-proposed_U+current_K-proposed_K))
		    Copy_FieldComplex(*current_q, *q);

	       Copy_FieldComplex(*vSamples[acc_sample], *current_q);
	       acc_sample +=1;
	  }
     }

     delete q;
     delete p;
     delete gradvTemp;
     delete qTemp;
}
