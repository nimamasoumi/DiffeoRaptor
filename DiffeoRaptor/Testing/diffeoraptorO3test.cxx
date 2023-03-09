#include "DiffeoRaptor.h"
#include "CalcWarpingIndex.h"
#include "Vec2D.h"
#include "MultiscaleManager.h"
#include "MultiscaleResampler.h"
#include "GaussianFilter.h"
#include "IdentityFilter.h"
#include "IOpers.h"
#include "FOpers.h"
#include "IFOpers.h"
#include "HFOpers.h"
#include "VFOpers.h"
#include "VFieldOpers.h"
#include "ImageFieldOpers.h"
#include <fftw3.h>
#include <time.h>

int main(int argc, char** argv)
{
  char I0Path[150];
  char I1Path[150];
  char I0PathSegStruct[150];
  char I1PathSegStruct[150];
  char I0PathSegClassic[150];
  char I1PathSegClassic[150];
  int maxIter[3],totalPatch[3],rrPatchSize[3];
  float stepSizeGD[3],corrFactor[3];

  int argi = 1;
  strcpy(I0Path, argv[argi++]);
  strcpy(I1Path, argv[argi++]);
  strcpy(I0PathSegStruct, argv[argi++]);
  strcpy(I1PathSegStruct, argv[argi++]);
  strcpy(I0PathSegClassic, argv[argi++]);
  strcpy(I1PathSegClassic, argv[argi++]);
  int numStep = atoi(argv[argi++]);
  int lpower = atoi(argv[argi++]);
  int truncX = atoi(argv[argi++]);
  int truncY = atoi(argv[argi++]);
  int truncZ = atoi(argv[argi++]);
  float sigma = atof(argv[argi++]);
  float alpha = atof(argv[argi++]);
  float gamma = atof(argv[argi++]);
  maxIter[0] = atoi(argv[argi++]);
  maxIter[1] = atoi(argv[argi++]);
  maxIter[2] = atoi(argv[argi++]);
  stepSizeGD[0] = atof(argv[argi++]);
  stepSizeGD[1] = atof(argv[argi++]);
  stepSizeGD[2] = atof(argv[argi++]);
  totalPatch[0]=atoi(argv[argi++]);
  totalPatch[1]=atoi(argv[argi++]);
  totalPatch[2]=atoi(argv[argi++]);
  corrFactor[0]=atof(argv[argi++]);
  corrFactor[1]=atof(argv[argi++]);
  corrFactor[2]=atof(argv[argi++]);
  double outlier_TH=atof(argv[argi++]);
  rrPatchSize[0]=atof(argv[argi++]);
  rrPatchSize[1]=atof(argv[argi++]);
  rrPatchSize[2]=atof(argv[argi++]);
  int numLevel=atoi(argv[argi++]);

  MemoryType mType= MEM_HOST;

  // read data
  Image3D* I0 = new Image3D(mType);
  Image3D* I1 = new Image3D(mType);
  Image3D* I0segstruct = new Image3D(mType);
  Image3D* I1segstruct = new Image3D(mType);
  Image3D* I0segclassic = new Image3D(mType);
  Image3D* I1segclassic = new Image3D(mType);

  ITKFileIO::LoadImage(*I0, I0Path);
  ITKFileIO::LoadImage(*I1, I1Path);
  ITKFileIO::LoadImage(*I0segstruct, I0PathSegStruct);
  ITKFileIO::LoadImage(*I1segstruct, I1PathSegStruct);
  ITKFileIO::LoadImage(*I0segclassic, I0PathSegClassic);
  ITKFileIO::LoadImage(*I1segclassic, I1PathSegClassic);

  GridInfo grid = I0->grid();
  Field3D *H0tmp = new Field3D(grid, mType);
  Field3D *H0Spatial = new Field3D(grid, mType);
  Image3D* I0n = new Image3D(grid, mType);
  Image3D* I1n = new Image3D(grid, mType);
  Image3D* I0d = new Image3D(grid, mType);
  Image3D* I1d = new Image3D(grid, mType);
  Image3D* Iest = new Image3D(grid, mType);
  Image3D* deformI0segstruct = new Image3D(grid, mType);
  Image3D* deformI0segclassic = new Image3D(grid, mType);

  MultiscaleManager* scaleManager = new MultiscaleManager(grid);
  switch(numLevel)
  {
    case 3:
      scaleManager->addScaleLevel(4);
      scaleManager->addScaleLevel(2);
      scaleManager->addScaleLevel(1);
      break;
    case 2:
      scaleManager->addScaleLevel(2);
      scaleManager->addScaleLevel(1);
      break;
    case 1:
      scaleManager->addScaleLevel(1);
      break;
  }

  typedef IdentityFilter<0> identityfilter;
  MultiscaleResampler<identityfilter>* resampler = new MultiscaleResampler<identityfilter>(grid);

  // access parameter
  for(int ii=0;ii<numLevel;ii++)
  {
      scaleManager->set(ii);
      resampler->setScaleLevel(*scaleManager);
      GridInfo curGrid = scaleManager->getCurGrid();

      I0n->setGrid(curGrid);
      I1n->setGrid(curGrid);
      I0d->setGrid(curGrid);
      I1d->setGrid(curGrid);
      Iest->setGrid(curGrid);
      H0tmp->setGrid(curGrid);

      // precalculate low frequency location
      if (truncX % 2 == 0) truncX -= 1; // set last dimension as zero if it is even
      if (truncY % 2 == 0) truncY -= 1;
      if (truncZ % 2 == 0) truncZ -= 1;
      FftOper *fftOper = new FftOper(alpha, gamma, lpower, curGrid, truncX, truncY, truncZ);
      fftOper->FourierCoefficient();

      // forward shooting
      FieldComplex3D *v0 = new FieldComplex3D(truncX, truncY, truncZ);
      FieldComplex3D *gradv = new FieldComplex3D(truncX, truncY, truncZ);

      // /////////////////////////////////////////////////
      DiffeoRaptor *diffeoraptor = new DiffeoRaptor(fftOper, mType, numStep, corrFactor[ii]);

      if(scaleManager->isLastScale())
      {
          Opers::Copy(*I0d,*I0);
          Opers::Copy(*I1d,*I1);
      }else
      {
          resampler->downsampleImage(*I0d,*I0);
          resampler->downsampleImage(*I1d,*I1);
      }

      if(scaleManager->isFirstScale())
      {
          H0Spatial->setGrid(curGrid);
          diffeoraptor->US_normalize(*I0n,I0d);
      }else
      {
          resampler->updateHField(*H0Spatial);
          Opers::ApplyH(*Iest, *I0d, *H0Spatial, BACKGROUND_STRATEGY_WRAP);
          diffeoraptor->US_normalize(*I0n,Iest);
      }

      diffeoraptor->US_normalize(*I1n,I1d);
      Opers::Copy(*(diffeoraptor->I0), *I0n);
      Opers::Copy(*(diffeoraptor->I1), *I1n);
      diffeoraptor->sigma = sigma;
      diffeoraptor->TotalPatch =totalPatch[ii];
      diffeoraptor->corrFactor=corrFactor[ii];
      diffeoraptor->outlier_TH=outlier_TH;
      diffeoraptor->rrPatchSize=rrPatchSize[ii];

      clock_t t;
      t = clock();

      diffeoraptor->ImageMatching(*v0, *gradv, maxIter[ii], stepSizeGD[ii]);

      t = clock() - t;
      printf("\n It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);

      Opers::Copy(*Iest,*(diffeoraptor->deformIm));
      Opers::Copy(*H0tmp,*(diffeoraptor->phiinv));

      if(scaleManager->isFirstScale())
      {
        Opers::Copy(*H0Spatial,*H0tmp);
      }else
      {
        Field3D *h0tmp = new Field3D(curGrid, mType);
        Opers::Copy(*h0tmp,*H0Spatial);
        Opers::ComposeHH(*H0Spatial,*h0tmp,*H0tmp);
      }

      std::ostringstream ofileDir;
      ofileDir << "collectedInfoLev"<<ii<<".txt";
      std::ofstream collectedInfo(ofileDir.str());
      std::vector<double> cIE, cE, fE;
      cIE = diffeoraptor->collectedCR;
      cE = diffeoraptor->collectedE;
      fE = diffeoraptor->collectedMSE;
      int jj=0;
      for(auto it = std::begin(cE); it != std::end(cE); ++it)
      {
        collectedInfo<<cIE[jj]<<", "<<cE[jj]<<", "<<fE[jj]<<endl;
        jj++;
      }
      collectedInfo.close();

      if(scaleManager->isLastScale())
      {
          Opers::ApplyH(*deformI0segstruct, *I0segstruct, *H0Spatial, BACKGROUND_STRATEGY_WRAP,INTERP_NN);
          Opers::ApplyH(*deformI0segclassic, *I0segclassic, *H0Spatial, BACKGROUND_STRATEGY_WRAP,INTERP_NN);
          // post-processing no multi-scale
          // Calc DICE score
          GridInfo griD = I0->grid();
          Vec3Di mSize = griD.size();
          int fsx = mSize.x;int fsy = mSize.y;int fsz = mSize.z;
          float sumJac=0, msumJac=0;int newJac=0, numJac=0;
          Image3D* modIm=new Image3D(curGrid,mType);
          Opers::Copy(*modIm,*(diffeoraptor->deformIm));
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
                    if(deformI0segstruct->get(zz,jj,kk)>0)
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
          ITKFileIO::SaveImage(*jacDiff, "jacDiff.nii.gz");
          ITKFileIO::SaveImage(*deformI0segstruct, "deformI0segstruct.nii.gz");
          ITKFileIO::SaveImage(*deformI0segclassic, "deformI0segclassic.nii.gz");

          CalcWarpingIndex* calcwarpingindex=new CalcWarpingIndex();
          std::vector<double> diceScoresStruct, diceScoresClassic;

          calcwarpingindex->calcDice_selective(diceScoresStruct,I1segstruct,I0segstruct,deformI0segstruct);
          std::cout<<"\nStructural initial Dice score: "<<diceScoresStruct[0]<<std::endl;
          std::cout<<"Structural final Dice score: "<<diceScoresStruct[diceScoresStruct.size()-1]<<std::endl;
          for (int kk=1;kk<diceScoresStruct.size()-1;kk++)
          {
              std::cout<<"Label "<<kk<<" Dice score: "<<diceScoresStruct[kk]<<std::endl;
          }

          calcwarpingindex->calcDice_selective(diceScoresClassic,I1segclassic,I0segclassic,deformI0segclassic);
          std::cout<<"\nClassic initial Dice score: "<<diceScoresClassic[0]<<std::endl;
          std::cout<<"Classic final Dice score: "<<diceScoresClassic[diceScoresClassic.size()-1]<<std::endl;
          for (int kk=1;kk<diceScoresClassic.size()-1;kk++)
          {
              std::cout<<"Label "<<kk<<" Dice score: "<<diceScoresClassic[kk]<<std::endl;
          }

          //Opers::Sum(sumJac,*jacDiff);
          msumJac=1-((double)newJac/(double)(numJac));
          std::cout<<"DiffeoRaptor Determinent of Jacobian: "<<msumJac<<std::endl;
      }
  }
}
