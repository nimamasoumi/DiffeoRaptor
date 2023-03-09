#include "DiffeoRaptor.h"
#include "Vec2D.h"
#include "GaussianFilter.h"
#include <time.h>
#include <random>
#include <algorithm>

DiffeoRaptor::DiffeoRaptor(FftOper* _fftOper,
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
     //dd = new Image3D(grid, mType);
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

     //Ifg_idx2yxz =new std::vector<std::vector<double>>(3);
     Ifg_idx2yxz.resize(3);
     TotalPatch=1000;
     iterCount=0;
     outlier_TH=1;
     corrFactor=1;
     rrPatchSize=3;
}

DiffeoRaptor::~DiffeoRaptor()
{
     delete I0; I0 = NULL;
     delete I1; I1 = NULL;
     //delete dd; dd = NULL;
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
void DiffeoRaptor::ad(FieldComplex3D& advw,
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
void DiffeoRaptor::adTranspose(FieldComplex3D& adTransvw,
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


void DiffeoRaptor::fwdIntegrateV(const FieldComplex3D& v0)
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
void DiffeoRaptor::fwdIntegration(FieldComplex3D& v0,
                      FieldComplex3D& fwd_gradvfft)
{
     fwdIntegrateV(v0);
     Copy_FieldComplex(v0, *storeV[numTimeSteps]); // done because prev version was updating v0

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

     double tmpcr=0;
     Image3D* dd=new Image3D(grid, mType);     
     Opers::SetMem(*dd,0.0);
     CRtGrad(tmpcr,*dd,I1,deformIm);          
     Opers::MulC_I(*dd,corrFactor);
     float fcE=0;     
     Opers::Sum2(fcE,*dd);
     std::cout<<"RaPTOR MSE "<<fcE<<std::endl;

     collectedCR.push_back(tmpcr);
     Opers::Copy(*residualIm,*dd);
     Opers::MulMulC_I(*scratchV1, *residualIm, -1.0);
     fftOper->spatial2fourier(fwd_gradvfft, *scratchV1);

     MulI_FieldComplex(fwd_gradvfft, *(fftOper->Kcoeff));
}

// update adjoint variable
// for backward integration
void DiffeoRaptor::bwdUpdate(FieldComplex3D& dvadj,
                 const FieldComplex3D& vadj,
                 const FieldComplex3D& ad,
                 const FieldComplex3D& adTrans)
{
     for (int i = 0; i < xDim * yDim * zDim * 3; i++)
      dvadj.data[i] += dt * (vadj.data[i] - ad.data[i] + adTrans.data[i]);
}

// backward integration, reduced adjoint jacobi fields
// Euler integration
void DiffeoRaptor::bwdIntegration(FieldComplex3D& dvadj,
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

// gradient term for image matching
void DiffeoRaptor::Gradient(FieldComplex3D& v0,
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
     fwdIntegration(v0, *fwd_gradvfft); //specifically this one
     // backward integration
     // imMatchGradient should be initialized as zero before passing
     bwdIntegration(*imMatchGradient, *fwd_gradvfft);

     // calculate Ienergy
     Image3D* Itmp = new Image3D(grid,mType);float fE=0;
     Opers::Sub(*Itmp, *deformIm, *I1);
     Opers::Sum2(fE, *Itmp); // This line should be replaced with the value of RaPTOR     

     IEnergy=corrFactor*collectedCR.back();     
     TotalEnergy = VEnergy + 0.5/(sigma*sigma)*IEnergy;
     //TotalEnergy = VEnergy + IEnergy;
     std::cout<<"MSE Venergy Ienergy Tenergy"<<std::endl;
     std::cout<<fE<<" "<<VEnergy<<" "<<0.5/(sigma*sigma)*IEnergy<<" "<<TotalEnergy<<std::endl;
     collectedE.push_back(TotalEnergy);
     collectedMSE.push_back(fE);
     
     // gradient term
     // gradv = v0 + \tilde{v};
     AddI_FieldComplex(gradv, *imMatchGradient, 1.0/(sigma*sigma));

     delete Itmp;Itmp=NULL;
}

// gradient descent strategy for image matching
// adaptive stepsize
void DiffeoRaptor::ImageMatching(FieldComplex3D& v0,
                     FieldComplex3D& gradv,
                     int maxIter,
                     float stepSizeGD)
{
     float epsilon = 1.0e-10;
     float prevEnergy = 1e20;

     //parameters
     float beta1=0.9,beta2=0.999,epsL=1e-8;
     FieldComplex3D *mom1= new FieldComplex3D(xDim, yDim, zDim);
     FieldComplex3D *umom1= new FieldComplex3D(xDim, yDim, zDim);
     FieldComplex3D *prevmom1= new FieldComplex3D(xDim, yDim, zDim);
     FieldComplex3D *currmom1= new FieldComplex3D(xDim, yDim, zDim);
     FieldComplex3D *newmom1= new FieldComplex3D(xDim, yDim, zDim);
     FieldComplex3D *temp1= new FieldComplex3D(xDim, yDim, zDim);
     FieldComplex3D *temp3= new FieldComplex3D(xDim, yDim, zDim);
     FieldComplex3D *temp4= new FieldComplex3D(xDim, yDim, zDim);

     FieldComplex3D *currV0 = new FieldComplex3D(xDim, yDim, zDim);
     FieldComplex3D *prevV0 = new FieldComplex3D(xDim, yDim, zDim);
     FieldComplex3D *newV0 = new FieldComplex3D(xDim, yDim, zDim);     
     Copy_FieldComplex(*prevV0, v0);
     MR_grad_US_mask_v3(Ifg_idx2yxz, I1, I0);     
     MulCI_FieldComplex(*mom1,0);
     Copy_FieldComplex(*prevmom1, *mom1);

     for (int i = 0; i < maxIter; i++)
     {
      iterCount=i;
      printf("\n%i th iteration!\n",i);
      Copy_FieldComplex(*currV0, v0);
      Copy_FieldComplex(*currmom1, *mom1);

      Gradient(v0, gradv);

      //Momentum gradient descent
      MulC_FieldComplex(*temp1,*currmom1,beta1);
      Add_FieldComplex(*newmom1,*temp1,gradv,1-beta1);
      MulC_FieldComplex(*umom1,*newmom1,1/(1-pow(beta1,i+1)));
      Copy_FieldComplex(*prevmom1, *currmom1);
      Copy_FieldComplex(*mom1, *newmom1);

      Add_FieldComplex(*newV0, *currV0, *umom1, -stepSizeGD);

      Copy_FieldComplex(*prevV0, *currV0);
      Copy_FieldComplex(v0, *newV0);
      prevEnergy = TotalEnergy;

     }

     Opers::SetMem(*residualIm, 1.0);
     Opers::Splat(*splatI, *phiinv, *I1, BACKGROUND_STRATEGY_WRAP);
     Opers::Splat(*splatOnes, *phiinv, *residualIm, BACKGROUND_STRATEGY_WRAP);

     delete currV0; currV0 = NULL;
     delete prevV0; prevV0 = NULL;
     delete newV0; newV0 = NULL;
     delete mom1;mom1 =NULL;
     delete umom1;umom1 =NULL;
     delete currmom1;currmom1 =NULL;
     delete prevmom1;prevmom1 =NULL;
     delete newmom1;newmom1 =NULL;
     delete temp1;temp1=NULL;
}

void DiffeoRaptor::US_normalize(Image3D& In,
                                const Image3D *I)
{
    GridInfo gridP = I->grid();
    Image3D* newIref=new Image3D(gridP, mType);
    Vec2Df* aaa=new Vec2Df();
    Opers::MaxMin(*aaa,*I);
    float maxref =aaa->x;
    float minref =aaa->y;    
    Opers::SubC(*newIref,*I,minref);
    Opers::MaxMin(*aaa,*newIref);
    maxref =aaa->x;
    minref =aaa->y;
    Opers::DivC(In,*newIref,maxref);

    delete newIref;newIref=NULL;
}

void DiffeoRaptor::H_CRc(double& CR,
                         Image3D& GR,
                         const Image3D *Iref,
                         const Image3D *Imov,
                         const int bins,
                         const int n_pts,
                         double& deNUM)
{
    // zero min and max 1
    GridInfo gridP = Iref->grid();    
    Image3D* nIref = new Image3D(gridP,mType);Image3D* pIref = new Image3D(gridP,mType);
    Image3D* IxB = new Image3D(gridP, mType);Image3D* b_n = new Image3D(gridP, mType);
    Image3D* b_f = new Image3D(gridP, mType);Image3D* b_f_1 = new Image3D(gridP, mType);

    double CR_n[bins+2]={0},CR_m[bins+2]={0},CR_mean[bins+2]={0};
    float Sum_Imov2=0,mu_Imov=0,DeNum=0,pSum=0;

    // bins determination
    US_normalize(*nIref,Iref);

    //CR calculation
    for (int ii = 0; ii < n_pts; ii++)
    {
        IxB->set(ii,nIref->get(ii)*bins);
        b_n->set(ii,std::round(IxB->get(ii))+1);
        b_f->set(ii,IxB->get(ii)-b_n->get(ii)+1.5);
        b_f_1->set(ii,1-b_f->get(ii));

        CR_n[(int) b_n->get(ii)-1]   += b_f_1->get(ii);
        CR_n[(int) b_n->get(ii)] += b_f->get(ii);

        CR_m[(int) b_n->get(ii)-1] += (b_f_1->get(ii)*Imov->get(ii));
        CR_m[(int) b_n->get(ii)] += (b_f->get(ii)*Imov->get(ii));
    }

    Opers::Sum(mu_Imov, *Imov);
    mu_Imov=mu_Imov/(double)n_pts;
    Opers::Sum2(Sum_Imov2, *Imov);
    for (int jj=0; jj<bins+2; jj++)
    {                
        if(CR_n[jj]==0.0)
        {
            CR_n[jj]=1;
        }
        CR_mean[jj]=CR_m[jj]/CR_n[jj];
        pSum+=CR_mean[jj]*CR_m[jj];
    }
    DeNum = (double)n_pts * ((Sum_Imov2/(double)n_pts) - pow(mu_Imov,2));
    CR = (Sum_Imov2 - pSum)/ DeNum;

    for (int ii = 0; ii < n_pts; ii++)
    {
        GR.set(ii, Imov->get(ii) - CR_mean[(int) b_n->get(ii)-1]*b_f_1->get(ii)
                - CR_mean[(int) b_n->get(ii)]*b_f->get(ii) - CR*(Imov->get(ii) - mu_Imov));
    }
    Opers::DivC_I(GR,DeNum);
    Opers::MulC_I(GR,2.0);
    deNUM=DeNum;
}

void DiffeoRaptor::CRtGrad(double& CR,
                             Image3D& dd,
                             const Image3D *Iref,
                             const Image3D *Imov)
{    
    float tmp_GR_sum=0;
    int bins=16;
    int Total_patches=TotalPatch;
    double rr =rrPatchSize;

    Image3D *Grd_bxx, *Grd_bxy, *Grd_bxz, *Irefg, *Imovg, *Irat;
    Irefg = new Image3D(grid,mType);Imovg = new Image3D(grid,mType);Irat = new Image3D(grid,mType);
    Field3D *GR = new Field3D(grid, mType);
    Opers::Gradient(*GR, *Imov, DIFF_CENTRAL, BC_WRAP);
    Opers::GradientMag(*Irefg,*Iref, DIFF_CENTRAL, BC_WRAP);
    Opers::GradientMag(*Imovg,*Imov, DIFF_CENTRAL, BC_WRAP);

    int Ifg_sum = Ifg_idx2yxz[0].size();
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<int> distribution (1,Ifg_sum-1);
    CR =0;

    Vec3Di mSize = grid.size();
    Vec3Di PatchSize;
    PatchSize.x=(int)(2*rr+1);PatchSize.y=(int)(2*rr+1);PatchSize.z=(int)(2*rr+1);
    if(mSize.z==1)
    {
        PatchSize.z=1;
    }
    double n_pts = (double)(PatchSize.x*PatchSize.y*PatchSize.z);

    double rv[(int)PatchSize.x], tmp_ddx0[3][(int)n_pts],var[3],dd_sum[3],dd_sum2[3], samSum[3], GRD_mag,rg,rg2, dd_tmp,rT;
    int ii=0,mm=0,ooo=0;
    double xx[(int)PatchSize.x],yy[(int)PatchSize.y],zz[(int)PatchSize.z];
    double tmpCR=0,deNUM=0;

    Vec3Di gridSize((int)PatchSize.x,(int)PatchSize.y,(int)PatchSize.z);
    GridInfo patchGrid;
    patchGrid.setSize(gridSize);
    float sum_tmp,sum_tmp2,tmp_ddx0_sum=0;
    Vec3Df tmp_grad;
    int ind =0,safeG=10000;
    if(false)//if(mSize.z!=1)
    {
        calc_ratio(*Irat, Irefg, Imovg, rr);
    }


    Image3D* Iref_bx0 = new Image3D(patchGrid,mType);
    Image3D* Imov_bx =  new Image3D(patchGrid,mType);
    Image3D* tmpGr0 =  new Image3D(patchGrid,mType);
    Grd_bxx =  new Image3D(patchGrid,mType);
    Grd_bxy =  new Image3D(patchGrid,mType);
    Grd_bxz =  new Image3D(patchGrid,mType);
    bool indexFlag=true;
    std::vector<double> dicInd;
    dicInd.push_back(-1);    
    while((ii<Total_patches))
    {
        while(indexFlag)
        {
            ind = distribution(generator);//rand()%Ifg_sum;
            auto it =std::find(dicInd.begin(),dicInd.end(),ind);
            if(it==dicInd.end())
            {
                dicInd.push_back((double)ind);
                indexFlag=false;
                break;
            }
        }
        indexFlag=true;   

        double xxi =Ifg_idx2yxz[0][ind];
        double yyi =Ifg_idx2yxz[1][ind];
        double zzi =Ifg_idx2yxz[2][ind];

        for(int jj=0;jj<(int)PatchSize.x; jj++)
        {
            rv[jj]=-rr+jj;
            yy[jj]=rv[jj]+yyi;
            xx[jj]=rv[jj]+xxi;
            if(mSize.z!=1)
            {
                zz[jj]=rv[jj]+zzi;
            }
            else
            {
                zz[0]=0;
            }
        }        

        for(int jj=0;jj<(int)PatchSize.x;jj++)
        {
            for(int kk=0;kk<(int)PatchSize.y;kk++)
            {
                for(int pp=0;pp<(int)PatchSize.z;pp++)
                {
                    Iref_bx0->set(jj,kk,pp,Iref->get((int)xx[jj],(int)yy[kk],(int)zz[pp]));
                    Imov_bx->set(jj,kk,pp,Imov->get((int)xx[jj],(int)yy[kk],(int)zz[pp]));
                    tmp_grad = GR->get((int)xx[jj],(int)yy[kk],(int)zz[pp]);
                    Grd_bxx->set(jj,kk,pp,tmp_grad.x);
                    Grd_bxy->set(jj,kk,pp,tmp_grad.y);
                    Grd_bxz->set(jj,kk,pp,tmp_grad.z);
                }
            }
        }
        if(false)//if(mSize.z!=1)
        {
            rg2=Irat->get((int)xxi,(int)yyi,(int)zzi);
        }

        Opers::Sum2(sum_tmp,*Iref_bx0);
        Opers::Sum2(sum_tmp2,*Imov_bx);

        if((sum_tmp>0.01)&&(sum_tmp2>0.01))
        {
            tmpCR=0;
            H_CRc(tmpCR, *tmpGr0, Iref_bx0, Imov_bx, bins, (int)n_pts, deNUM);

            dd_sum[0]=0;dd_sum[1]=0;dd_sum[2]=0;
            dd_sum2[0]=0;dd_sum2[1]=0;dd_sum2[2]=0;
            for(int i2=0;i2<(int)n_pts;i2++)
            {
                tmp_ddx0[0][i2]=Grd_bxx->get(i2)*tmpGr0->get(i2);
                tmp_ddx0[1][i2]=Grd_bxy->get(i2)*tmpGr0->get(i2);
                tmp_ddx0[2][i2]=Grd_bxz->get(i2)*tmpGr0->get(i2);
                dd_sum[0]=dd_sum[0]+tmp_ddx0[0][i2];
                dd_sum[1]=dd_sum[1]+tmp_ddx0[1][i2];
                dd_sum[2]=dd_sum[2]+tmp_ddx0[2][i2];
                dd_sum2[0]=dd_sum2[0]+pow(tmp_ddx0[0][i2],2);
                dd_sum2[1]=dd_sum2[1]+pow(tmp_ddx0[1][i2],2);
                dd_sum2[2]=dd_sum2[2]+pow(tmp_ddx0[2][i2],2);
            }
            GRD_mag=sqrt(pow(dd_sum[0],2)+pow(dd_sum[1],2)+pow(dd_sum[2],2));
            tmp_ddx0_sum=sqrt((dd_sum2[0]/n_pts)-pow(dd_sum[0]/n_pts,2))+
                         sqrt((dd_sum2[1]/n_pts)-pow(dd_sum[1]/n_pts,2))+
                         sqrt((dd_sum2[2]/n_pts)-pow(dd_sum[2]/n_pts,2));

            if(GRD_mag==0)
            {
                GRD_mag=1e-7;
            }
            if(false)//if(mSize.z!=1)
            {
                rT = rg2*(tmp_ddx0_sum/GRD_mag);
            }
            else
            {
                rT = (tmp_ddx0_sum/GRD_mag);
            }
            if(!std::isfinite(tmpCR))
            {
                ooo++;
            }

            if((rT<outlier_TH)&&std::isfinite(tmpCR))
            {
                Opers::Sum2(tmp_GR_sum,*tmpGr0);
                CR = CR + tmpCR;
                for(int jj=0;jj<(int)PatchSize.x;jj++)
                {
                    for(int kk=0;kk<(int)PatchSize.y;kk++)
                    {
                        for(int pp=0;pp<(int)PatchSize.z;pp++)
                        {
                            if(isfinite(tmpGr0->get(jj,kk,pp)))
                            {
                                dd_tmp = dd.get((int)xx[jj],(int)yy[kk],(int)zz[pp])+tmpGr0->get(jj,kk,pp);
                                dd.set((int)xx[jj],(int)yy[kk],(int)zz[pp],dd_tmp);
                            }                            
                        }
                    }
                }
                ii=ii+1;
             }
        }
        mm++;
        if(mm==Total_patches+safeG)
        {
            std::cout<<"Reached this point! "<<ii<<endl;
            break;
        }
    }
    std::cout<<"Singular patches: "<<ooo<<std::endl;    
    Opers::DivC_I(dd,(float)Total_patches);
    CR = CR/(float)Total_patches;

    delete Iref_bx0;Iref_bx0=NULL;
    delete Imov_bx;Imov_bx=NULL;
    delete Irefg;Irefg= NULL;
    delete Imovg;Imovg= NULL;
    delete GR;GR= NULL;
}

void DiffeoRaptor::MR_grad_US_mask_v3(std::vector<std::vector<double>>& idx2yxz,
                                      const Image3D *If,
                                      const Image3D *Im)
{
    Image3D* Ifg = new Image3D(grid, mType);
    Image3D* Img = new Image3D(grid, mType);
    Image3D* Im_mas = new Image3D(grid, mType);
    Image3D* Img_mas = new Image3D(grid, mType);
    Image3D* If_mas = new Image3D(grid, mType);
    Image3D* Ifg_mas = new Image3D(grid, mType);
    Image3D* MASK = new Image3D(grid, mType);
    float sum_Img,sum_Ifg;
    Vec3Di mSize = grid.size();
    int fsx = mSize.x;
    int fsy = mSize.y;
    int fsz = mSize.z;

    Opers::GradientMag(*Ifg,*If, DIFF_CENTRAL, BC_WRAP);
    Opers::GradientMag(*Img,*Im, DIFF_CENTRAL, BC_WRAP);
    Opers::GTC(*Im_mas,*Im,0.001);
    Opers::Sum(sum_Img,*Img);
    sum_Img = sum_Img/(fsx*fsy*fsz);
    Opers::GTC(*Img_mas,*Img,sum_Img);
    Opers::GTC(*If_mas,*If,0.001);
    Opers::Sum(sum_Ifg,*Ifg);
    sum_Ifg = sum_Ifg/(fsx*fsy*fsz);
    Opers::GTC(*Ifg_mas,*Ifg,sum_Ifg);
    Opers::Mul(*MASK,*Im_mas,*Img_mas);
    Opers::Mul(*MASK,*MASK,*If_mas);
    Opers::Mul(*MASK,*MASK,*Ifg_mas);


    int boundCut=rrPatchSize,boundCutz=rrPatchSize;
    if(fsz==1)
    {
        boundCutz=0;
    }
    for(int ii=boundCut;ii<fsx-boundCut;ii++)
    {
        for(int jj=boundCut;jj<fsy-boundCut;jj++)
        {
            for(int kk=boundCutz;kk<fsz-boundCutz;kk++)
            {
                if(MASK->get(ii,jj,kk)>0)
                {
                    idx2yxz[0].push_back((double)ii);
                    idx2yxz[1].push_back((double)jj);
                    idx2yxz[2].push_back((double)kk);
                }
            }
        }
    }

    delete Ifg;Ifg=NULL;
    delete Img;Img=NULL;
    delete Im_mas;Im_mas=NULL;
    delete Img_mas;Img_mas=NULL;
    delete If_mas;If_mas=NULL;
    delete Ifg_mas;Ifg_mas=NULL;
    delete MASK;MASK=NULL;
}

void DiffeoRaptor::calc_ratio(Image3D& Irat,
                              const Image3D *Irefg,
                              const Image3D *Imovg,
                              const double rr)
{    

    Image3D* Irefg_sum = new Image3D(grid, mType);
    Image3D* Irefg_sum1 = new Image3D(grid, mType);
    Image3D* Irefg_sum2 = new Image3D(grid, mType);
    Image3D* Imovg_sum = new Image3D(grid, mType);
    Image3D* Imovg_sum1 = new Image3D(grid, mType);
    Image3D* Imovg_sum2 = new Image3D(grid, mType);

    int PatchSize=2*rr+1;   
    double v[PatchSize]={1/7,2/7,3/7,4/7,5/7,6/7,1};

    Vec3Di mSize = grid.size();
    int fsx = mSize.x;
    int fsy = mSize.y;
    int fsz = mSize.z;
    for (int i=rr;i<fsx-rr;i++)
    {
        for (int j=0;j<fsy;j++)
        {
            for(int k=0;k<fsz;k++)
            {

                    Irefg_sum1->set(i,j,k,Irefg->get(i-3,j,k)*v[PatchSize-1]+
                                          Irefg->get(i-2,j,k)*v[PatchSize-2]+
                                          Irefg->get(i-1,j,k)*v[PatchSize-3]+
                                          Irefg->get(i,j,k)*v[PatchSize-4]+
                                          Irefg->get(i+1,j,k)*v[PatchSize-5]+
                                          Irefg->get(i+2,j,k)*v[PatchSize-6]+
                                          Irefg->get(i+3,j,k)*v[PatchSize-7]);
                    Imovg_sum1->set(i,j,k,Imovg->get(i-3,j,k)*v[PatchSize-1]+
                                          Imovg->get(i-2,j,k)*v[PatchSize-2]+
                                          Imovg->get(i-1,j,k)*v[PatchSize-3]+
                                          Imovg->get(i,j,k)*v[PatchSize-4]+
                                          Imovg->get(i+1,j,k)*v[PatchSize-5]+
                                          Imovg->get(i+2,j,k)*v[PatchSize-6]+
                                          Imovg->get(i+3,j,k)*v[PatchSize-7]);

            }
        }
    }
    for (int i=0;i<fsx;i++)
    {
        for (int j=rr;j<fsy-rr;j++)
        {
            for(int k=0;k<fsz;k++)
            {
                Irefg_sum2->set(i,j,k,Irefg_sum1->get(i,j-3,k)*v[PatchSize-1]+
                                      Irefg_sum1->get(i,j-2,k)*v[PatchSize-2]+
                                      Irefg_sum1->get(i,j-1,k)*v[PatchSize-3]+
                                      Irefg_sum1->get(i,j,k)*v[PatchSize-4]+
                                      Irefg_sum1->get(i,j+1,k)*v[PatchSize-5]+
                                      Irefg_sum1->get(i,j+2,k)*v[PatchSize-6]+
                                      Irefg_sum1->get(i,j+3,k)*v[PatchSize-7]);
                Imovg_sum2->set(i,j,k,Imovg_sum1->get(i,j-3,k)*v[PatchSize-1]+
                                      Imovg_sum1->get(i,j-2,k)*v[PatchSize-2]+
                                      Imovg_sum1->get(i,j-1,k)*v[PatchSize-3]+
                                      Imovg_sum1->get(i,j,k)*v[PatchSize-4]+
                                      Imovg_sum1->get(i,j+1,k)*v[PatchSize-5]+
                                      Imovg_sum1->get(i,j+2,k)*v[PatchSize-6]+
                                      Imovg_sum1->get(i,j+3,k)*v[PatchSize-7]);
            }
        }
    }
    for (int i=0;i<fsx;i++)
    {
        for (int j=0;j<fsy;j++)
        {
            for(int k=rr;k<fsz-rr;k++)
            {

                Irefg_sum->set(i,j,k,Irefg_sum2->get(i,j,k-3)*v[PatchSize-1]+
                                     Irefg_sum2->get(i,j,k-2)*v[PatchSize-2]+
                                     Irefg_sum2->get(i,j,k-1)*v[PatchSize-3]+
                                     Irefg_sum2->get(i,j,k)*v[PatchSize-4]+
                                     Irefg_sum2->get(i,j,k+1)*v[PatchSize-5]+
                                     Irefg_sum2->get(i,j,k+2)*v[PatchSize-6]+
                                     Irefg_sum2->get(i,j,k+3)*v[PatchSize-7]);
                Imovg_sum->set(i,j,k,Imovg_sum2->get(i,j,k-3)*v[PatchSize-1]+
                                     Imovg_sum2->get(i,j,k-2)*v[PatchSize-2]+
                                     Imovg_sum2->get(i,j,k-1)*v[PatchSize-3]+
                                     Imovg_sum2->get(i,j,k)*v[PatchSize-4]+
                                     Imovg_sum2->get(i,j,k+1)*v[PatchSize-5]+
                                     Imovg_sum2->get(i,j,k+2)*v[PatchSize-6]+
                                     Imovg_sum2->get(i,j,k+3)*v[PatchSize-7]);
                if(Imovg_sum->get(i,j,k)==0)
                {
                    Imovg_sum->set(i,j,k,1e-7);
                    Irat.set(i,j,k,Irefg_sum->get(i,j,k)/Imovg_sum->get(i,j,k));
                }
                else
                {
                    Irat.set(i,j,k,Irefg_sum->get(i,j,k)/Imovg_sum->get(i,j,k));
                }
            }
        }
    }
}

void DiffeoRaptor::outlierThGen(double& outlier_TH,
                                const std::vector<double> pointTarget,
                                const Image3D* Iref,
                                const Image3D* Imov,
                                const double rr)
{
    double Irefg_sum=0,Imovg_sum=0,h;
    double k[3];
    int PatchSize=2*rr+1;

    for(int ii=0;ii<PatchSize;ii++)
    {
        for(int jj=0;jj<PatchSize;jj++)
        {
            for(int kk=0;kk<PatchSize;kk++)
            {
                h=((double)ii+1)*((double)jj+1)*((double)kk+1)/(double)pow(PatchSize,3);
                k[0]=pointTarget[0]-(double)ii;
                k[1]=pointTarget[1]-(double)jj;
                k[2]=pointTarget[2]-(double)kk;
                Irefg_sum+=(Iref->get(k[0],k[1],k[2]))*h;
                Imovg_sum+=(Imov->get(k[0],k[1],k[2]))*h;
            }
        }
    }
    outlier_TH=(Imovg_sum==0)?(Irefg_sum/1e-7):(Irefg_sum/Imovg_sum);
}






