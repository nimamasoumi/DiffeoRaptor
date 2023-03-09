# include "FftOper.h"

FftOper::FftOper(const float _alpha,
		 const float _gamma,
		 const int _lpow,
		 const GridInfo _grid,
		 const int _truncX,
		 const int _truncY, 
		 const int _truncZ)
{
     alpha = _alpha;
     gamma = _gamma;
     lpow = _lpow;
     grid = _grid;
     fsx = grid.size().x;
     fsy = grid.size().y;
     fsz = grid.size().z;

     truncX = _truncX;
     truncY = _truncY;
     truncZ = _truncZ;

     size = truncX * truncY * truncZ;
     fsize = fsx * fsy * fsz;

     Lcoeff = new FieldComplex3D(truncX, truncY, truncZ);
     Kcoeff = new FieldComplex3D(truncX, truncY, truncZ);
     CDcoeff = new FieldComplex3D(truncX, truncY, truncZ);
     fftLoc = new int[3 * size];
     corrLoc = new int[3 * size];
     conjLoc = new int[3 * size];

     // for converting data between fourier and spatial domain
     fsxFFT = fsx / 2 + 1;
     scratch = new float[2 * fsxFFT * fsy * fsz];
     data = new float[fsize];

     fftwForwardPlan = fftwf_plan_dft_r2c_3d(fsz, fsy, fsx, data,
					     (fftwf_complex*)scratch, FFTW_MEASURE);

     fftwBackwardPlan = fftwf_plan_dft_c2r_3d(fsz, fsy, fsx, (fftwf_complex*)scratch,
					      data, FFTW_MEASURE);
     beginX = -(truncX-1) / 2.0;
     beginY = -(truncY-1) / 2.0;
     beginZ = -(truncZ-1) / 2.0;
     endX = beginX + truncX;
     endY = beginY + truncY;
     endZ = beginZ + truncZ;

     // for convolution in Fourier domain
     padX = 2 * truncX - 1; 
     padY = 2 * truncY - 1;
     padZ = 2 * truncZ - 1;
     padsize = padX * padY * padZ;

     scratchI = new float[2 * padsize];  // complex: real + imaginary
     scratchK = new float[2 * padsize];
     dataI = new float[2 * padsize];
     dataK = new float[2 * padsize];

     scratchX = new float[2 * padsize];
     scratchY = new float[2 * padsize];
     scratchZ = new float[2 * padsize];

     h_ix = new FieldComplex3D(padX, padY, padZ);
     h_iy = new FieldComplex3D(padX, padY, padZ);
     h_iz = new FieldComplex3D(padX, padY, padZ);
     h_kernel = new FieldComplex3D(padX, padY, padZ);
     h_i = new FieldComplex3D(padX, padY, padZ);

     fftwFwdPlanI = fftwf_plan_dft_3d(padZ, padY, padX, (fftwf_complex*)dataI,
				      (fftwf_complex*)scratchI, -1, FFTW_MEASURE);
     fftwFwdPlanK = fftwf_plan_dft_3d(padZ, padY, padX, (fftwf_complex*)dataK,
				      (fftwf_complex*)scratchK, -1, FFTW_MEASURE);
     fftwBwdPlan = fftwf_plan_dft_3d(padZ, padY, padX, (fftwf_complex*)scratchI,
				     (fftwf_complex*)dataI, +1, FFTW_MEASURE);
}

FftOper::~FftOper()
{
     delete [] scratch; scratch = NULL;
     delete [] data; data = NULL;
     delete [] dataI; dataI = NULL;
     delete [] dataK; dataK = NULL;
     delete [] scratchX; scratchX = NULL;
     delete [] scratchY; scratchY = NULL;
     delete [] scratchZ; scratchZ = NULL;
     delete h_ix; h_ix = NULL;
     delete h_iy; h_iy = NULL;
     delete h_iz; h_iz = NULL;
     delete h_i; h_i = NULL;
     delete h_kernel; h_kernel = NULL;
     delete Lcoeff; Lcoeff = NULL;
     delete Kcoeff; Kcoeff = NULL;
     delete CDcoeff; CDcoeff = NULL;
     delete [] fftLoc; fftLoc = NULL;
     delete [] corrLoc; corrLoc = NULL;
     delete [] conjLoc; conjLoc = NULL;

     fftwf_destroy_plan (fftwForwardPlan);
     fftwf_destroy_plan (fftwBackwardPlan);

     fftwf_destroy_plan (fftwFwdPlanI);
     fftwf_destroy_plan (fftwFwdPlanK);
     fftwf_destroy_plan (fftwBwdPlan);
}

// extract low frequency value by fftshift
// precalculate the location matrix
// truncDimX, truncDimY, truncDimZ: truncated Dim of XYZ
void FftOper::fftshiftMat()
{
     int id = 0;
     int i, j, k;
     int fftX, fftY, fftZ, corrX, corrY, corrZ, conjX, conjY, conjZ;

     for(k = beginZ; k < endZ; k++) {
	  for(j = beginY; j < endY; j++) { 
	       for(i = beginX; i < endX; i++, id++) {
		    fftZ = k; corrZ = k; conjZ = k;
		    fftY = j; corrY = j; conjY = j; 
		    fftX = i; corrX = i; conjX = i;

		    if(k < 0) {fftZ = k + fsz; corrZ = k + padZ; conjZ = fftZ;}
		    if(j < 0) {fftY = j + fsy; corrY = j + padY; conjY = fftY;}
		    if(i < 0) {fftX = i + fsx; corrX = i + padX;
			       conjX = -i ; conjY = -j; conjZ = -k;
                               if (j > 0) conjY = -j + fsy;
			       if (k > 0) conjZ = -k + fsz;}
    
		    fftLoc[3*id] = fftX; corrLoc[3*id] = corrX; conjLoc[3*id] = conjX;
		    fftLoc[3*id+1] = fftY; corrLoc[3*id+1] = corrY; conjLoc[3*id+1] = conjY;
		    fftLoc[3*id+2] = fftZ; corrLoc[3*id+2] = corrZ; conjLoc[3*id+2] = conjZ;
	       }
	  }
     }
}


// Returns Fourier coefficient
// L, K, central difference, and divergence operator
// L, K, and CD need to be initialized as zero
void FftOper::FourierCoefficient()
{
     fftshiftMat(); // generate fftLoc

     double sX = 2.0 * M_PI / fsx;
     double sY = 2.0 * M_PI / fsy;
     double sZ = 2.0 * M_PI / fsz;

     double xcoeff, ycoeff, zcoeff;
     int index, id;
     double val;

     double spx = (double)grid.spacing().x;
     double spy = (double)grid.spacing().y;
     double spz = (double)grid.spacing().z;

     for (id = 0; id < size; id++)   
     {
	  index = 3*id;
	  xcoeff = (-2.0*cos(sX*fftLoc[index]) + 2.0)/(spx*spx);
	  ycoeff = (-2.0*cos(sY*fftLoc[index+1]) + 2.0)/(spy*spy);
	  zcoeff = (-2.0*cos(sZ*fftLoc[index+2]) + 2.0)/(spz*spz);

	  val = pow(alpha*(xcoeff + ycoeff + zcoeff)+gamma, lpow);

	  Lcoeff->data[index].real( val ); 
	  Lcoeff->data[index+1].real( val );
	  Lcoeff->data[index+2].real( val );

	  Kcoeff->data[index].real( 1.0/val );
	  Kcoeff->data[index+1].real( 1.0/val );
	  Kcoeff->data[index+2].real( 1.0/val );

	  CDcoeff->data[index].imag( sin(sX*fftLoc[index])/spx );
	  CDcoeff->data[index+1].imag( sin(sY*fftLoc[index+1])/spy );
	  CDcoeff->data[index+2].imag( sin(sZ*fftLoc[index+2])/spz );
     }
}

// convert to truncated fourier domain
void FftOper::spatial2fourier(FieldComplex3D& vtrunc,
			      const Field3D& vf)
{
     int i, j, k, index, findex, id; // findex: index in full spatial domain
     int findexArr[size];

     // for x component
     for(i = 0; i < fsize; i++)
	  data[i] = vf.x[i] / static_cast<float>(fsize);

     fftwf_execute(fftwForwardPlan);

     // extract low frequency
     id = 0;
     for(k = beginZ; k < endZ; k++) {
	  for(j = beginY; j < endY; j++) { 
	       for(i = beginX; i < endX; i++, id++) {
		    index = 3*id;
		    findex = 2*(fsxFFT*(conjLoc[index+2]*fsy + conjLoc[index+1]) + conjLoc[index]);
		    findexArr[id] = findex;
		    if (i < 0)
			 vtrunc.data[index] = complex<float>(scratch[findex], 
							     -1.0*scratch[findex+1]);
		    else
			 vtrunc.data[index] = complex<float>(scratch[findex],
							     scratch[findex+1]);
	       }
	  }
     }
 
     //============================================================================
     // for y component
     for(i = 0; i < fsize; i++)
	  data[i] = vf.y[i] / static_cast<float>(fsize);

     fftwf_execute(fftwForwardPlan);

     id = 0;
     for(k = beginZ; k < endZ; k++) {
	  for(j = beginY; j < endY; j++) { 
	       for(i = beginX; i < endX; i++, id++) {
		    if (i < 0)
			 vtrunc.data[3*id+1] = complex<float>(scratch[findexArr[id]],
							      -1.0*scratch[findexArr[id]+1]);
		    else
			 vtrunc.data[3*id+1] = complex<float>(scratch[findexArr[id]],
							      scratch[findexArr[id]+1]);
	       } 
	  }
     }
      
     //============================================================================
     // for z component
     for(i = 0; i < fsize; i++)
	  data[i] = vf.z[i] / static_cast<float>(fsize);

     fftwf_execute(fftwForwardPlan);

     // extract low frequency
     id = 0;
     for(k = beginZ; k < endZ; k++) {
	  for(j = beginY; j < endY; j++) { 
	       for(i = beginX; i < endX; i++, id++) {
		    if (i < 0)
			 vtrunc.data[3*id+2] = complex<float>(scratch[findexArr[id]],
							      -1.0*scratch[findexArr[id]+1]);
		    else
			 vtrunc.data[3*id+2] = complex<float>(scratch[findexArr[id]],
							      scratch[findexArr[id]+1]);
	       } 
	  }
     }
}

// convert to full image domain
void FftOper::spatial2fourier_F(float* xf,
				float* yf,
				float* zf,
				const Field3D& vf)
{    
  // for x component
  for(int i = 0; i < fsize; i++)
    data[i] = vf.x[i] / static_cast<float>(fsize);
  
  fftwf_execute(fftwForwardPlan);
  memcpy(xf, scratch, 2 * fsxFFT * fsy * fsz * sizeof(float)); 
  
  //============================================================================
  // for y component
  for(int i = 0; i < fsize; i++)
    data[i] = vf.y[i] / static_cast<float>(fsize);
  
  fftwf_execute(fftwForwardPlan);
  memcpy(yf, scratch, 2 * fsxFFT * fsy * fsz * sizeof(float));     
  
  //============================================================================
  // for z component
  for(int i = 0; i < fsize; i++)
    data[i] = vf.z[i] / static_cast<float>(fsize);

  fftwf_execute(fftwForwardPlan);
  memcpy(zf, scratch, 2 * fsxFFT * fsy * fsz * sizeof(float));  
}

// convert to full spatial domain
void FftOper::fourier2spatial(Field3D& vf,
			      const FieldComplex3D& vtrunc)
{
     int i, j, k, id, index, findex;
     int findexArr[size];

     // for x component
     memset(scratch, 0.0, 2 * fsxFFT * fsy * fsz * sizeof(float));

     id = 0;
     for(k = beginZ; k < endZ; k++) {
	  for(j = beginY; j < endY; j++) { 
	       for(i = beginX; i < endX; i++, id++) {
		    index = 3 * id;
		    findex = 2*(fsxFFT*(conjLoc[index+2]*fsy + conjLoc[index+1]) + conjLoc[index]);
		    findexArr[id] = findex;
		    if (i >= 0) {
			 scratch[findex] = vtrunc.data[index].real();
			 scratch[findex+1] = vtrunc.data[index].imag();
		    } 
	       }
	  }
     }
  
     fftwf_execute(fftwBackwardPlan);
     memcpy(vf.getX(), data, fsize * sizeof(float)); 

     //============================================================================
     // for y component
     memset(scratch, 0.0, 2 * fsxFFT * fsy * fsz * sizeof(float));

     id = 0;
     for(k = beginZ; k < endZ; k++) {
	  for(j = beginY; j < endY; j++) { 
	       for(i = beginX; i < endX; i++, id++) {
		    index = 3 * id;
		    if (i >= 0) {
			 scratch[findexArr[id]] = vtrunc.data[index+1].real();
			 scratch[findexArr[id]+1] = vtrunc.data[index+1].imag();
		    } 
	       }
	  }
     }

     fftwf_execute(fftwBackwardPlan);
 
     memcpy(vf.getY(), data, fsize * sizeof(float));

     //============================================================================
     // for z component
     memset(scratch, 0.0, 2 * fsxFFT * fsy * fsz * sizeof(float));

     id = 0;
     for(k = beginZ; k < endZ; k++) {
	  for(j = beginY; j < endY; j++) { 
	       for(i = beginX; i < endX; i++, id++) {
		    index = 3 * id;
		    if (i >= 0) {
			 scratch[findexArr[id]] = vtrunc.data[index+2].real();
			 scratch[findexArr[id]+1] = vtrunc.data[index+2].imag();
		    } 
	       }
	  }
     }

     fftwf_execute(fftwBackwardPlan);
     memcpy(vf.getZ(), data, fsize * sizeof(float));
}

// add high frequencies and then ifft
void FftOper::fourier2spatial_addH(Field3D& vf,
				 const FieldComplex3D& vtrunc,
				 const float* idx,
				 const float* idy,
				 const float* idz)
{
     int i, j, k, id, index, findex;
     int findexArr[size];

     // for x component
     memcpy(scratch, idx, 2 * fsxFFT * fsy * fsz * sizeof(float));

     id = 0;
     for(k = beginZ; k < endZ; k++) {
	  for(j = beginY; j < endY; j++) { 
	       for(i = beginX; i < endX; i++, id++) {
		    index = 3 * id;
		    findex = 2*(fsxFFT*(conjLoc[index+2]*fsy + conjLoc[index+1]) + conjLoc[index]);
		    findexArr[id] = findex;
		    if (i >= 0) {
			 scratch[findex] += vtrunc.data[index].real();
			 scratch[findex+1] += vtrunc.data[index].imag();
		    } 
	       }
	  }
     }
  
     fftwf_execute(fftwBackwardPlan);
     memcpy(vf.getX(), data, fsize * sizeof(float)); 

     //============================================================================
     // for y component
     memcpy(scratch, idy, 2 * fsxFFT * fsy * fsz * sizeof(float));

     id = 0;
     for(k = beginZ; k < endZ; k++) {
	  for(j = beginY; j < endY; j++) { 
	       for(i = beginX; i < endX; i++, id++) {
		    index = 3 * id;
		    if (i >= 0) {
		      scratch[findexArr[id]] += vtrunc.data[index+1].real();
		      scratch[findexArr[id]+1] += vtrunc.data[index+1].imag();
		    } 
	       }
	  }
     }

     fftwf_execute(fftwBackwardPlan);
 
     memcpy(vf.getY(), data, fsize * sizeof(float));

     //============================================================================
     // for z component
     memcpy(scratch, idz, 2 * fsxFFT * fsy * fsz * sizeof(float));

     id = 0;
     for(k = beginZ; k < endZ; k++) {
	  for(j = beginY; j < endY; j++) { 
	       for(i = beginX; i < endX; i++, id++) {
		    index = 3 * id;
		    if (i >= 0) {
		      scratch[findexArr[id]] += vtrunc.data[index+2].real();
		      scratch[findexArr[id]+1] += vtrunc.data[index+2].imag();
		    } 
	       }
	  }
     }

     fftwf_execute(fftwBackwardPlan);
     memcpy(vf.getZ(), data, fsize * sizeof(float));
}

// convolution : 0 / correlation : 1, multiplication in fourier domain
void FftOper::PointwiseMultiply(float *im, float *in, int flag)
{
     complex<float> tempI;
     for(int i = 0; i < padsize; i++)
     { 
	  if (flag == 0)
	  {
	       tempI = complex<float>(im[2*i], im[2*i+1])*complex<float>(in[2*i], in[2*i+1]);
	  }
	  else if (flag == 1)
	  {
	       tempI = complex<float>(im[2*i], -1.0*im[2*i+1])*complex<float>(in[2*i], in[2*i+1]);
	  }
	  scratchI[2*i] = tempI.real() / static_cast<float>(padsize);
	  scratchI[2*i+1] = tempI.imag() / static_cast<float>(padsize);
     }
     fftwf_execute(fftwBwdPlan);
}

//Return convolution by fft
//Matrix multiplication conv(Dv, w)
void FftOper::ConvolveComplexFFT(FieldComplex3D& h_o_ori,
				 const int flag,
				 const FieldComplex3D& h_ix_ori, 
				 const FieldComplex3D& h_iy_ori,
				 const FieldComplex3D& h_iz_ori, 
				 const FieldComplex3D& h_kernel_ori)
{
     int i, x, y, z, index, index_ori, id;
     h_ix->initVal(complex<float>(0.0, 0.0));
     h_iy->initVal(complex<float>(0.0, 0.0));
     h_iz->initVal(complex<float>(0.0, 0.0));
     h_kernel->initVal(complex<float>(0.0, 0.0));

     // pad zeros
     for (z = 0; z < truncZ; z++)
	  for (y = 0; y < truncY; y++)
	  {
	       index_ori = 3*truncX*(z*truncY + y);
	       index = 3*padX*(z*padY + y);
	       
	       memcpy(&h_ix->data[index], &h_ix_ori.data[index_ori], truncX*3*sizeof(complex<float>));
	       memcpy(&h_iy->data[index], &h_iy_ori.data[index_ori], truncX*3*sizeof(complex<float>));
	       memcpy(&h_iz->data[index], &h_iz_ori.data[index_ori], truncX*3*sizeof(complex<float>));
	       memcpy(&h_kernel->data[index], &h_kernel_ori.data[index_ori], truncX*3*sizeof(complex<float>));
	  }
     
     // fourier transform
     int indexArr[size];
     
     //==============================================
     // for x component
     // assign data for kernel and jacobian matrix x
     for(i = 0; i < padsize; i++)
     { 
	  index = 3*i;
	  dataI[2*i] = h_ix->data[index].real();
	  dataI[2*i+1] = h_ix->data[index].imag();

	  dataK[2*i] = h_kernel->data[index].real();
	  dataK[2*i+1] = h_kernel->data[index].imag();
     }
     fftwf_execute(fftwFwdPlanI);
     fftwf_execute(fftwFwdPlanK);

     PointwiseMultiply(scratchI, scratchK, flag);

     // put data back to output
     int cx =  truncX / 2;
     int cy =  truncY / 2;
     int cz =  truncZ / 2;
     id = 0;
     for (z = 0; z < truncZ; z++)
	  for (y = 0; y < truncY; y++)
	       for (x = 0; x < truncX; x++, id++)
	       {
		    if (flag == 0)
			 index = 2*(padX*((z+cz)*padY + (y+cy)) + (x+cx));
		    else if (flag == 1)
		     	 index = 2*(padX*(corrLoc[3*id+2]*padY + corrLoc[3*id+1]) + corrLoc[3*id]); 
		    h_o_ori.data[3*id] = complex<float>(dataI[index], dataI[index+1]);
		    indexArr[id] = index;
	       }

     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     //assign data for jacobian matrix y
     for(i = 0; i < padsize; i++)
     { 
	  index = 3*i;
	  dataI[2*i] = h_ix->data[index+1].real();
	  dataI[2*i+1] = h_ix->data[index+1].imag();
     }
 
     fftwf_execute(fftwFwdPlanI);

     PointwiseMultiply(scratchI, scratchK, flag);

     // put data back to output
     for(id = 0; id < size; id++)
	  h_o_ori.data[3*id+1] = complex<float>(dataI[indexArr[id]], dataI[indexArr[id]+1]);

     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     //assign data for jacobian matrix z
     for(i = 0; i < padsize; i++)
     { 
	  index = 3*i;
	  dataI[2*i] = h_ix->data[index+2].real();
	  dataI[2*i+1] = h_ix->data[index+2].imag();
     }
 
     fftwf_execute(fftwFwdPlanI);

     PointwiseMultiply(scratchI, scratchK, flag);

     // put data back to output
     for(id = 0; id < size; id++)
	  h_o_ori.data[3*id+2] = complex<float>(dataI[indexArr[id]], dataI[indexArr[id]+1]);
  
     //================================================
     // for y component
     //assign data
     for(i = 0; i < padsize; i++)
     { 
	  index = 3*i;
	  dataI[2*i] = h_iy->data[index].real();
	  dataI[2*i+1] = h_iy->data[index].imag();
      
	  dataK[2*i] = h_kernel->data[index+1].real();
	  dataK[2*i+1] = h_kernel->data[index+1].imag();
     }
     fftwf_execute(fftwFwdPlanI);
     fftwf_execute(fftwFwdPlanK);

     PointwiseMultiply(scratchI, scratchK, flag);

     // put data back to output
     for(id = 0; id < size; id++)
	  h_o_ori.data[3*id] += complex<float>(dataI[indexArr[id]], dataI[indexArr[id]+1]);

     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     //assign data for jacobian matrix y
     for(i = 0; i < padsize; i++)
     { 
	  index = 3*i;
	  dataI[2*i] = h_iy->data[index+1].real();
	  dataI[2*i+1] = h_iy->data[index+1].imag();
     }
 
     fftwf_execute(fftwFwdPlanI);

     PointwiseMultiply(scratchI, scratchK, flag);

     // put data back to output
     for(id = 0; id < size; id++)
	  h_o_ori.data[3*id+1] += complex<float>(dataI[indexArr[id]], dataI[indexArr[id]+1]);

     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     //assign data for jacobian matrix z
     for(i = 0; i < padsize; i++)
     { 
	  index = 3*i;
	  dataI[2*i] = h_iy->data[index+2].real();
	  dataI[2*i+1] = h_iy->data[index+2].imag();
     }
 
     fftwf_execute(fftwFwdPlanI);

     PointwiseMultiply(scratchI, scratchK, flag);

     // put data back to output
     for(id = 0; id < size; id++)
	  h_o_ori.data[3*id+2] += complex<float>(dataI[indexArr[id]], dataI[indexArr[id]+1]);

     //==================================================
     // for z component
     //assign data
     for(i = 0; i < padsize; i++)
     { 
	  index = 3*i;
	  dataI[2*i] = h_iz->data[index].real();
	  dataI[2*i+1] = h_iz->data[index].imag();
      
	  dataK[2*i] = h_kernel->data[index+2].real();
	  dataK[2*i+1] = h_kernel->data[index+2].imag();
     }
     fftwf_execute(fftwFwdPlanI);
     fftwf_execute(fftwFwdPlanK);

     PointwiseMultiply(scratchI, scratchK, flag);

     // put data back to output
     for(id = 0; id < size; id++)
	  h_o_ori.data[3*id] += complex<float>(dataI[indexArr[id]], dataI[indexArr[id]+1]);

     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     //assign data for jacobian matrix y
     for(i = 0; i < padsize; i++)
     { 
	  index = 3*i;
	  dataI[2*i] = h_iz->data[index+1].real();
	  dataI[2*i+1] = h_iz->data[index+1].imag();
     }
 
     fftwf_execute(fftwFwdPlanI);

     PointwiseMultiply(scratchI, scratchK, flag);

     // put data back to output
     for(id = 0; id < size; id++)
	  h_o_ori.data[3*id+1] += complex<float>(dataI[indexArr[id]], dataI[indexArr[id]+1]);

     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     //assign data for jacobian matrix z
     for(i = 0; i < padsize; i++)
     { 
	  index = 3*i;
	  dataI[2*i] = h_iz->data[index+2].real();
	  dataI[2*i+1] = h_iz->data[index+2].imag();
     }
 
     fftwf_execute(fftwFwdPlanI);

     PointwiseMultiply(scratchI, scratchK, flag);
   
     // put data back to output
     for(id = 0; id < size; id++)
	  h_o_ori.data[3*id+2] += complex<float>(dataI[indexArr[id]], dataI[indexArr[id]+1]);
}

//Return convolution by fft
// Div(Lw, v)
void FftOper::CorrComplexFFT(FieldComplex3D& h_o_ori,
			     const FieldComplex3D& h_i_ori, //v
			     const FieldComplex3D& h_kernel_ori, //lw
			     const FieldComplex3D& CDcoeff) 
{
     if (h_i_ori.Norm() < 1.0e-20 || h_kernel_ori.Norm() < 1.0e-20)
     {
     	  h_o_ori.initVal(complex<float>(0.0, 0.0));
     	  return;
     }

     int i, x, y, z, index, index_ori, id;
     h_i->initVal(complex<float>(0.0, 0.0));
     h_kernel->initVal(complex<float>(0.0, 0.0));

     // pad zeros
     for (z = 0; z < truncZ; z++)
	  for (y = 0; y < truncY; y++)
	  {
	       index_ori = 3*truncX*(z*truncY + y);
	       index = 3*padX*(z*padY + y);
	       memcpy(&h_i->data[index], &h_i_ori.data[index_ori], truncX*3*sizeof(complex<float>));
	       memcpy(&h_kernel->data[index], &h_kernel_ori.data[index_ori], truncX*3*sizeof(complex<float>));
	  }

     // fourier transform
     int indexArr[size];
    
     // fwd fft on x 
     for(i = 0; i < padsize; i++)
     { 
	  index = 3*i;
	  dataI[2*i] = h_i->data[index].real();
	  dataI[2*i+1] = h_i->data[index].imag();
     }
     fftwf_execute(fftwFwdPlanI);
     memcpy(scratchX, scratchI, 2 * padsize * sizeof(float));

     // fwd fft on y
     for(i = 0; i < padsize; i++)
     { 
	  index = 3*i;
	  dataI[2*i] = h_i->data[index+1].real();
	  dataI[2*i+1] = h_i->data[index+1].imag();
     }
     fftwf_execute(fftwFwdPlanI);
     memcpy(scratchY, scratchI, 2 * padsize * sizeof(float));

     // fwd fft on z
     for(i = 0; i < padsize; i++)
     { 
	  index = 3*i;
	  dataI[2*i] = h_i->data[index+2].real();
	  dataI[2*i+1] = h_i->data[index+2].imag();
     }
     fftwf_execute(fftwFwdPlanI);
     memcpy(scratchZ, scratchI, 2 * padsize * sizeof(float));

     //==============================================
     // for x component
     // assign data for kernel
     for(i = 0; i < padsize; i++)
     { 
	  index = 3*i;
	  dataK[2*i] = h_kernel->data[index].real();
	  dataK[2*i+1] = h_kernel->data[index].imag();
     }
     fftwf_execute(fftwFwdPlanK);

     PointwiseMultiply(scratchX, scratchK, 1);
    
     // put data back to output
     id = 0;
     for (z = 0; z < truncZ; z++)
	  for (y = 0; y < truncY; y++)
	       for (x = 0; x < truncX; x++, id++)
	       {
		    index = 2*(padX*(corrLoc[3*id+2]*padY + corrLoc[3*id+1]) + corrLoc[3*id]); 
		    h_o_ori.data[3*id] = complex<float>(dataI[index], dataI[index+1])
			 * CDcoeff.data[3*id];
		    indexArr[id] = index;
	       }

     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     PointwiseMultiply(scratchY, scratchK, 1);
     
     // put data back to output
     for(id = 0; id < size; id++)
	  h_o_ori.data[3*id] += complex<float>(dataI[indexArr[id]], dataI[indexArr[id]+1])
	       * CDcoeff.data[3*id+1];

     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     PointwiseMultiply(scratchZ, scratchK, 1);
     
     // put data back to output
     for(id = 0; id < size; id++)
	  h_o_ori.data[3*id] += complex<float>(dataI[indexArr[id]], dataI[indexArr[id]+1])
	       * CDcoeff.data[3*id+2];

     //================================================
     // for y component
     for(i = 0; i < padsize; i++)
     { 
	  index = 3*i;
	  dataK[2*i] = h_kernel->data[index+1].real();
	  dataK[2*i+1] = h_kernel->data[index+1].imag();
     }
     fftwf_execute(fftwFwdPlanK);

     PointwiseMultiply(scratchX, scratchK, 1);
     
     // put data back to output
     for(id = 0; id < size; id++)
	  h_o_ori.data[3*id+1] = complex<float>(dataI[indexArr[id]], dataI[indexArr[id]+1])
	       * CDcoeff.data[3*id];

     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     PointwiseMultiply(scratchY, scratchK, 1);
     
     // put data back to output
     for(id = 0; id < size; id++)
	  h_o_ori.data[3*id+1] += complex<float>(dataI[indexArr[id]], dataI[indexArr[id]+1])
	       * CDcoeff.data[3*id+1];

     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     PointwiseMultiply(scratchZ, scratchK, 1);
     
     // put data back to output
     for(id = 0; id < size; id++)
	  h_o_ori.data[3*id+1] += complex<float>(dataI[indexArr[id]], dataI[indexArr[id]+1])
	       * CDcoeff.data[3*id+2];

     //==================================================
     // for z component
     for(i = 0; i < padsize; i++)
     { 
	  index = 3*i;
	  dataK[2*i] = h_kernel->data[index+2].real();
	  dataK[2*i+1] = h_kernel->data[index+2].imag();
     }
     fftwf_execute(fftwFwdPlanK);

     PointwiseMultiply(scratchX, scratchK, 1);
     
     // put data back to output
     for(id = 0; id < size; id++)
	  h_o_ori.data[3*id+2] = complex<float>(dataI[indexArr[id]], dataI[indexArr[id]+1])
	       * CDcoeff.data[3*id];

     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     PointwiseMultiply(scratchY, scratchK, 1);
     
     // put data back to output
     for(id = 0; id < size; id++)
	  h_o_ori.data[3*id+2] += complex<float>(dataI[indexArr[id]], dataI[indexArr[id]+1])
	       * CDcoeff.data[3*id+1];

     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     PointwiseMultiply(scratchZ, scratchK, 1);
     
     // put data back to output
     for(id = 0; id < size; id++)
	  h_o_ori.data[3*id+2] += complex<float>(dataI[indexArr[id]], dataI[indexArr[id]+1])
	       * CDcoeff.data[3*id+2];
}
