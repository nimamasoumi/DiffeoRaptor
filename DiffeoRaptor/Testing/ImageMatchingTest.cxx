#include "GeodesicShooting.h"
#include "CalcWarpingIndex.h"
#include "DiffeoRaptor.h"
#include <time.h>

int main(int argc, char** argv)
{
  char I0Path[150];
  char I1Path[150];
  char I0PathSeg[150];
  char I1PathSeg[150];
  char I0PathSegStruct[150];
  char I1PathSegStruct[150];
  //char proTocol[100];
  //double Rpar[6];

  int argi = 1; 
  strcpy(I0Path, argv[argi++]);
  strcpy(I1Path, argv[argi++]);
  strcpy(I0PathSeg, argv[argi++]);
  strcpy(I1PathSeg, argv[argi++]);
  strcpy(I0PathSegStruct, argv[argi++]);
  strcpy(I1PathSegStruct, argv[argi++]);
  int numStep = atoi(argv[argi++]);
  int lpower = atoi(argv[argi++]);
  int truncX = atoi(argv[argi++]);
  int truncY = atoi(argv[argi++]);
  int truncZ = atoi(argv[argi++]);
  float sigma = atof(argv[argi++]);
  float alpha = atof(argv[argi++]);
  float gamma = atof(argv[argi++]);
  int maxIter = atoi(argv[argi++]);
  float stepSizeGD = atof(argv[argi++]);
  //int memType = atoi(argv[argi++]);
  //Rpar[0]=atof(argv[argi++])*(M_PI/180);
  //Rpar[1]=atof(argv[argi++])*(M_PI/180);
  //Rpar[2]=atof(argv[argi++])*(M_PI/180);
  //Rpar[3]=atof(argv[argi++]);
  //Rpar[4]=atof(argv[argi++]);
  //Rpar[5]=atof(argv[argi++]);
  //int caseNo = atoi(argv[argi++]);
  //strcpy(proTocol, argv[argi++]);

  MemoryType mType= MEM_HOST;
  // runs on CPU or GPU
  //if (memType == 0)
    //mType = MEM_HOST;
  //else
    //mType = MEM_DEVICE;
	  
  // read data
  Image3D *I0 = new Image3D(mType);
  Image3D *I1 = new Image3D(mType);
  Image3D *I0seg = new Image3D(mType);
  Image3D *I1seg = new Image3D(mType);
  Image3D* I0segstruct = new Image3D(mType);
  Image3D* I1segstruct = new Image3D(mType);

  ITKFileIO::LoadImage(*I0, I0Path);
  ITKFileIO::LoadImage(*I1, I1Path);
  ITKFileIO::LoadImage(*I0seg, I0PathSeg);
  ITKFileIO::LoadImage(*I1seg, I1PathSeg);
  ITKFileIO::LoadImage(*I0segstruct, I0PathSegStruct);
  ITKFileIO::LoadImage(*I1segstruct, I1PathSegStruct);
  
  // access parameter
  GridInfo grid = I0->grid();
  Vec3Di mSize = grid.size();

  int fsx = mSize.x;
  int fsy = mSize.y;
  int fsz = mSize.z;

  // precalculate low frequency location
  if (truncX % 2 == 0) truncX -= 1; // set last dimension as zero if it is even
  if (truncY % 2 == 0) truncY -= 1;
  if (truncZ % 2 == 0) truncZ -= 1;
  FftOper *fftOper = new FftOper(alpha, gamma, lpower, grid, truncX, truncY, truncZ); 
  fftOper->FourierCoefficient();

  // forward shooting
  FieldComplex3D *v0 = new FieldComplex3D(truncX, truncY, truncZ);
  Field3D *v0Spatial = new Field3D(grid, mType);
  FieldComplex3D *gradv = new FieldComplex3D(truncX, truncY, truncZ);

  Image3D *I0n = new Image3D(grid, mType);
  Image3D *I1n = new Image3D(grid, mType);
  Image3D *deformI0seg = new Image3D(grid, mType);
  Image3D* deformI0segstruct = new Image3D(grid, mType);

//  Field3D* vi0spatial;
//  vi0spatial=new Field3D(grid,mType);
//  Image3D* I0w = new Image3D(grid,mType);
//  CalcWarpingIndex *calcwarpingindex = new CalcWarpingIndex();
//  calcwarpingindex->imageWarp(*vi0spatial,*I0w,Rpar,I1Path,I0Path);
  DiffeoRaptor *diffeoraptor = new DiffeoRaptor(fftOper, mType, numStep, 1e4);
//  diffeoraptor->US_normalize(*I0n,I0w);
  diffeoraptor->US_normalize(*I0n,I0);
  diffeoraptor->US_normalize(*I1n,I1);
  Opers::Copy(*(diffeoraptor->I0), *I0n);
  Opers::Copy(*(diffeoraptor->I1), *I1n);
  diffeoraptor->sigma = sigma;
  diffeoraptor->TotalPatch =1000;

  GeodesicShooting *geodesicshooting = new GeodesicShooting(fftOper, mType, numStep);
  Opers::Copy(*(geodesicshooting->I0), *I0n);
  Opers::Copy(*(geodesicshooting->I1), *I1n);
  geodesicshooting->sigma = sigma;
  
  clock_t t;
  t = clock();

  geodesicshooting->ImageMatching(*v0, *gradv, maxIter, stepSizeGD);

  t = clock() - t;
  printf("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
	
  //ITKFileIO::SaveImage(*(geodesicshooting->deformIm), "deformIm.nii.gz");
//ITKFileIO::SaveImage(*(I0), "deformIf.nii.gz");
  fftOper->fourier2spatial(*v0Spatial, *v0);
  //ITKFileIO::SaveField(*v0Spatial, "v0Spatial.mhd");
  //Opers::ApplyH(*deformI0seg, *I0seg, *(geodesicshooting->phiinv), BACKGROUND_STRATEGY_WRAP);
  Opers::ApplyH(*deformI0seg, *I0seg, *(geodesicshooting->phiinv), BACKGROUND_STRATEGY_WRAP,INTERP_NN);
  Opers::ApplyH(*deformI0segstruct, *I0segstruct, *(geodesicshooting->phiinv), BACKGROUND_STRATEGY_WRAP,INTERP_NN);

//  ITKFileIO::SaveField(*vi0spatial, "vi0Spatial.mhd");

//  double ind0=0,ind1=0;
//  Image3D* mask=new Image3D(grid, mType);
//  Image3D* I0orgn = new Image3D(grid,mType);

//  diffeoraptor->US_normalize(*I0orgn,I0);
//  std::vector<std::vector<double>> Ifg_idx2yxz(3);
//  diffeoraptor->MR_grad_US_mask_v3(Ifg_idx2yxz, I0orgn, I0orgn);
//  //std::cout<<"\nNumber of points in the mask: "<<Ifg_idx2yxz[0].size()<<std::endl;
//  for(int ii=0;ii<Ifg_idx2yxz[0].size();ii++)
//  {
//      mask->set(Ifg_idx2yxz[0][ii],Ifg_idx2yxz[1][ii],Ifg_idx2yxz[2][ii],100);
//  }
//  calcwarpingindex->initialFinalIndex(ind0, ind1, vi0spatial, v0Spatial, mask);
//  std::cout<<"Initial index: "<<ind0<<std::endl;
//  std::cout<<"Final index: "<<ind1<<std::endl;

  float fE=0,fC=0;
  Image3D* Itmp = new Image3D(grid,mType);
  Image3D* Itmp2 = new Image3D(grid,mType);

  Opers::Sub(*Itmp, *I0n, *I1n);
  Opers::Sub(*Itmp2, *(geodesicshooting->deformIm), *I1n);
  Opers::Sum2(fE, *Itmp);
  Opers::Sum2(fC, *Itmp2);

  std::cout<<"\nInitial MSE: "<<fE<<endl;
  std::cout<<"Final MSE: "<<fC<<endl;
  std::cout<<"MSE Reduction Ratio: "<<(fE-fC)/fE<<endl;
  std::cout<<"MSE Ratio: "<<fC/fE<<endl;

  float sumJac=0, msumJac=0;int newJac=0, numJac=0;
  Image3D* modIm=new Image3D(grid,mType);
  Opers::Copy(*modIm,*(geodesicshooting->deformIm));
  Image3D* jacDiff = new Image3D(grid, mType);
  Opers::JacDetH(*jacDiff, *(diffeoraptor->phiinv));
  for (int zz = 0; zz < fsx; zz++)
  {
     for(int jj=0;jj<fsy;jj++)
     {
        for(int kk=0;kk<fsz;kk++)
        {
            if(!isfinite(modIm->get(zz,jj,kk)))
            {
                modIm->set(zz,jj,kk,0.0);
            }
            if(deformI0seg->get(zz,jj,kk)>0)
            {
                numJac++;
                if(jacDiff->get(zz,jj,kk)<0)
                {
                    newJac++;
                }
            }
        }
     }
  }
  ITKFileIO::SaveImage(*modIm, "deformIm.nii.gz");
  ITKFileIO::SaveImage(*jacDiff, "jacFlash.nii.gz");
  ITKFileIO::SaveImage(*deformI0seg, "deformI0seg.nii.gz");
  ITKFileIO::SaveImage(*deformI0segstruct, "deformI0segstruct.nii.gz");

  CalcWarpingIndex* calcwarpingindex=new CalcWarpingIndex();
  std::vector<double> diceScoresStruct, diceScoresClassic;

  calcwarpingindex->calcDice_selective(diceScoresStruct,I1segstruct,I0segstruct,deformI0segstruct);
  std::cout<<"\nStructural initial Dice score: "<<diceScoresStruct[0]<<std::endl;
  std::cout<<"Structural final Dice score: "<<diceScoresStruct[diceScoresStruct.size()-1]<<std::endl;
  for (int kk=1;kk<diceScoresStruct.size()-1;kk++)
  {
      std::cout<<"Label "<<kk<<" Dice score: "<<diceScoresStruct[kk]<<std::endl;
  }

  calcwarpingindex->calcDice_selective(diceScoresClassic,I1seg,I0seg,deformI0seg);
  std::cout<<"\nClassic initial Dice score: "<<diceScoresClassic[0]<<std::endl;
  std::cout<<"Classic final Dice score: "<<diceScoresClassic[diceScoresClassic.size()-1]<<std::endl;
  for (int kk=1;kk<diceScoresClassic.size()-1;kk++)
  {
      std::cout<<"Label "<<kk<<" Dice score: "<<diceScoresClassic[kk]<<std::endl;
  }

  //Opers::Sum(sumJac,*jacDiff);
  msumJac=1-((double)newJac/(double)(numJac));
  std::cout<<"DiffeoRaptor Determinent of Jacobian: "<<msumJac<<std::endl; 

  delete fftOper;
  delete diffeoraptor;
  delete geodesicshooting;
  delete I0;
  delete I1;
  delete v0;
  delete v0Spatial;
  delete gradv;
}
