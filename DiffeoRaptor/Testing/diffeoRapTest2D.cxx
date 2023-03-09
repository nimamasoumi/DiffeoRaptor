#include "DiffeoRaptor.h"
#include "CalcWarpingIndex.h"
#include "Vec2D.h"
#include "ImageFieldOpers.h"
#include <fftw3.h>
#include <time.h>

int main(int argc, char** argv)
{
  char I0Path[100];
  char I1Path[100];
  char I0PathSeg[100];
  char I1PathSeg[100];

  int argi = 1;
  strcpy(I0Path, argv[argi++]);
  strcpy(I1Path, argv[argi++]);
  strcpy(I0PathSeg, argv[argi++]);
  strcpy(I1PathSeg, argv[argi++]);
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
  int memType = atoi(argv[argi++]);
  int totalPatch=atoi(argv[argi++]);
  double corrFactor=atof(argv[argi++]);
  double outlier_TH=atof(argv[argi++]);
  double rrPatchSize=atof(argv[argi++]);

  MemoryType mType;
  // runs on CPU or GPU
  if (memType == 0)
    mType = MEM_HOST;
  else
    mType = MEM_DEVICE;

  // read data
  Image3D *I0, *I1, *I0n, *I1n, *I1seg, *I0seg, *deformI0seg;
  I0 = new Image3D(mType);
  I1 = new Image3D(mType);
  I0seg = new Image3D(mType);
  I1seg = new Image3D(mType);

  ITKFileIO::LoadImage(*I0, I0Path);
  ITKFileIO::LoadImage(*I1, I1Path);
  ITKFileIO::LoadImage(*I0seg, I0PathSeg);
  ITKFileIO::LoadImage(*I1seg, I1PathSeg);

  // access parameter
  GridInfo grid = I0->grid();
  Vec3Di mSize = grid.size();

  I0n = new Image3D(grid, mType);
  I1n = new Image3D(grid, mType);
  deformI0seg = new Image3D(grid, mType);

  int fsx = mSize.x;
  int fsy = mSize.y;
  int fsz = mSize.z;
  std::cout<<fsx<<", "<<fsy<<", "<<fsz<<std::endl;

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

  // /////////////////////////////////////////////////
  DiffeoRaptor *diffeoraptor = new DiffeoRaptor(fftOper, mType, numStep);
  diffeoraptor->US_normalize(*I0n,I0);
  diffeoraptor->US_normalize(*I1n,I1);
  Opers::Copy(*(diffeoraptor->I0), *I0n);
  Opers::Copy(*(diffeoraptor->I1), *I1n);
  diffeoraptor->sigma = sigma;
  diffeoraptor->TotalPatch =totalPatch;
  diffeoraptor->corrFactor=corrFactor;
  diffeoraptor->outlier_TH=outlier_TH;
  diffeoraptor->rrPatchSize=rrPatchSize;

  clock_t t;
  t = clock();

  diffeoraptor->ImageMatching(*v0, *gradv, maxIter, stepSizeGD);

  t = clock() - t;
  printf("\nIt took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);

  ITKFileIO::SaveImage(*(diffeoraptor->deformIm), "deformIm.nii.gz");
  fftOper->fourier2spatial(*v0Spatial, *v0);
  ITKFileIO::SaveField(*v0Spatial, "v0Spatial.mhd");
  Opers::ApplyH(*deformI0seg, *I0seg, *(diffeoraptor->phiinv), BACKGROUND_STRATEGY_WRAP);

  std::ofstream collectedInfo("collectedInfo.txt");
  std::vector<double> cIE, cE, fE;
  cIE = diffeoraptor->collectedCR;
  cE = diffeoraptor->collectedE;
  fE = diffeoraptor->collectedMSE;
  int ii=0;
  for(auto it = std::begin(cE); it != std::end(cE); ++it)
  {
    collectedInfo<<cIE[ii]<<", "<<cE[ii]<<", "<<fE[ii]<<endl;
    ii++;
  }
  collectedInfo.close();

  // Calc DICE score
  int initInt=0, finalInt=0;
  double initDice=0, finalDice=0;
  for (int ii = 0; ii < fsx; ii++)
  {
     for(int jj=0;jj<fsy;jj++)
     {
        for(int kk=0;kk<fsz;kk++)
        {
            if(std::round(I1seg->get(ii,jj,kk))==std::round(I0seg->get(ii,jj,kk)))
            {
                initInt++;
            }
            if(std::round(I1seg->get(ii,jj,kk))==std::round(deformI0seg->get(ii,jj,kk)))
            {
                finalInt++;
            }
        }
     }
  }
  initDice=((double)initInt)/(fsx*fsy*fsz);
  finalDice=((double)finalInt)/(fsx*fsy*fsz);
  std::cout<<"\nInitial Dice score: "<<initDice<<std::endl;
  std::cout<<"Final Dice score: "<<finalDice<<"\n"<<std::endl;    

  Image3D* jacDiff = new Image3D(grid, mType);
  Opers::JacDetH(*jacDiff, *(diffeoraptor->phiinv));
  ITKFileIO::SaveImage(*jacDiff, "jacDiff.nii.gz");

  delete fftOper;
  delete diffeoraptor;
  delete I0;
  delete I1;
  delete v0;
  delete v0Spatial;
  delete gradv;
}
