#include "DiffeoRaptor.h"
#include "CalcWarpingIndex.h"
#include "Vec2D.h"
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
  char I0PathSegStruct[150];
  char I1PathSegStruct[150];
  char I0PathSegStructDeform[150];
  char I0PathSegClassic[150];
  char I1PathSegClassic[150];
  char I0PathSegClassicDeform[150];
  MemoryType mType = MEM_HOST;
  Image3D* I0segstruct = new Image3D(mType);
  Image3D* I1segstruct = new Image3D(mType);
  Image3D* I0segstructdeform = new Image3D(mType);
  Image3D* I0segclassic = new Image3D(mType);
  Image3D* I1segclassic = new Image3D(mType);
  Image3D* I0segclassicdeform = new Image3D(mType);

  CalcWarpingIndex* calcwarpingindex=new CalcWarpingIndex();
  std::vector<double> diceScoresStruct, diceScoresClassic;  

  int argi = 1;
  // read arguments
  strcpy(I0PathSegStruct, argv[argi++]);
  strcpy(I1PathSegStruct, argv[argi++]);
  strcpy(I0PathSegStructDeform, argv[argi++]);
  strcpy(I0PathSegClassic, argv[argi++]);
  strcpy(I1PathSegClassic, argv[argi++]);
  strcpy(I0PathSegClassicDeform, argv[argi++]);  

  // read data  
  ITKFileIO::LoadImage(*I0segstruct, I0PathSegStruct);
  ITKFileIO::LoadImage(*I1segstruct, I1PathSegStruct);
  ITKFileIO::LoadImage(*I0segstructdeform, I0PathSegStructDeform);
  ITKFileIO::LoadImage(*I0segclassic, I0PathSegClassic);
  ITKFileIO::LoadImage(*I1segclassic, I1PathSegClassic);
  ITKFileIO::LoadImage(*I0segclassicdeform, I0PathSegClassicDeform);

  // calculating the Dice scores
  calcwarpingindex->calcDice_selective(diceScoresStruct,I1segstruct,I0segstruct,I0segstructdeform);
  std::cout<<"\nStructural initial Dice score: "<<diceScoresStruct[0]<<std::endl;
  std::cout<<"Structural final Dice score: "<<diceScoresStruct[diceScoresStruct.size()-1]<<std::endl;
  for (int kk=1;kk<diceScoresStruct.size()-1;kk++)
  {
      std::cout<<"Label "<<kk<<" Dice score: "<<diceScoresStruct[kk]<<std::endl;
  }

  calcwarpingindex->calcDice_selective(diceScoresClassic,I1segclassic,I0segclassic,I0segclassicdeform);
  std::cout<<"\nClassical initial Dice score: "<<diceScoresClassic[0]<<std::endl;
  std::cout<<"Classical final Dice score: "<<diceScoresClassic[diceScoresClassic.size()-1]<<std::endl;
  for (int kk=1;kk<diceScoresClassic.size()-1;kk++)
  {
      std::cout<<"Label "<<kk<<" Dice score: "<<diceScoresClassic[kk]<<std::endl;
  }

  //std::cout<<"\n"<<argc<<std::endl;
  if(argc>7)
  {
      char jacDiffPath[150];
      char methodC[10];
      char outputFileDirClassic[150];
      char outputFileDirStruct[150];
      char caseNo[10];
      Image3D* jacDiff = new Image3D(mType);
      float sumJac=0, msumJac=0, mlogJac=0;int newJac=0, numJac=0;
      strcpy(jacDiffPath, argv[argi++]);
      strcpy(methodC, argv[argi++]);
      strcpy(outputFileDirStruct, argv[argi++]);
      strcpy(outputFileDirClassic, argv[argi++]);
      strcpy(caseNo, argv[argi++]);
      ITKFileIO::LoadImage(*jacDiff,jacDiffPath);

      // calculating the determinent of Jacobian statistics
      GridInfo griD = I0segstruct->grid();
      Vec3Di mSize = griD.size();
      int fsx = mSize.x;int fsy = mSize.y;int fsz = mSize.z;
      for (int zz = 0; zz < fsx; zz++)
      {
         for(int jj=0;jj<fsy;jj++)
         {
            for(int kk=0;kk<fsz;kk++)
            {
                if(I0segstructdeform->get(zz,jj,kk)>0)
                {
                    numJac++;
                    if(jacDiff->get(zz,jj,kk)<0)
                    {
                        newJac++;
                    }
                    if(jacDiff->get(zz,jj,kk)!=0)
                    {
                        sumJac=sumJac+std::log10(jacDiff->get(zz,jj,kk));
                    }
                }
            }
         }
      }
      //Opers::Sum(sumJac,*jacDiff);
      msumJac=1-((double)newJac/(double)(numJac));
      mlogJac=((double)sumJac/(double)(numJac));
      std::cout<<"Determinent of Jacobian: "<<msumJac<<" "<<mlogJac<<" "<<newJac<<" "<<sumJac<<std::endl;
      if(!std::strcmp(methodC,"diff"))
      {
          std::cout<<"\n\nWriting files for DiffRaptor case: "<<caseNo<<std::endl;

      }else if(!std::strcmp(methodC,"ants"))
      {
        std::cout<<"\n\nWriting files for ANTs case: "<<caseNo<<std::endl;

      }else if(!std::strcmp(methodC,"flash"))
      {
        std::cout<<"\n\nWriting files for FLASH case: "<<caseNo<<std::endl;
      }
      std::ostringstream ofileDir1;
      ofileDir1<<outputFileDirStruct;
      std::ofstream diffeoRaptorData1(ofileDir1.str(),std::ios::app);
      diffeoRaptorData1<<caseNo<<" ";
      for(auto it=diceScoresStruct.begin();it!=diceScoresStruct.end();it++)
      {
        diffeoRaptorData1<<*it<<" ";
      }
      diffeoRaptorData1<<std::endl;
      diffeoRaptorData1.close();

      std::ostringstream ofileDir2;
      ofileDir2<<outputFileDirClassic;
      std::ofstream diffeoRaptorData2(ofileDir2.str(),std::ios::app);
      diffeoRaptorData2<<caseNo<<" ";
      for(auto it=diceScoresClassic.begin();it!=diceScoresClassic.end();it++)
      {
        diffeoRaptorData2<<*it<<" ";
      }
      diffeoRaptorData2<<std::endl;
      diffeoRaptorData2.close();
  }
}
