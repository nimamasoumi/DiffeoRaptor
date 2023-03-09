#ifndef CALCWARPINGINDEX_H
#define CALCWARPINGINDEX_H

#include "FftOper.h"
#include "FieldComplex3D.h"
#include "ITKFileIO.h"
#include "Vec2D.h"
#include "IOpers.h"
#include "FOpers.h"
#include "IFOpers.h"
#include "HFOpers.h"
#include "Reduction.h"
#include "FluidKernelFFT.h"
#include <vector>
#include <math.h>
#include <fstream>

#include "itkImage.h"
#include "itkVersorRigid3DTransform.h"
#include "itkEuler3DTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkTransformToDisplacementFieldFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkChangeInformationImageFilter.h"

using namespace PyCA;
using namespace std;

class CalcWarpingIndex
{
public:

    CalcWarpingIndex();
    ~CalcWarpingIndex();

    void initialFinalIndex(double& ind0,
                           double& ind1,
                           const Field3D* vi0,
                           const Field3D* v0,
                           const Image3D* mask);

    void imageWarp(Field3D& vi0spatial,
                   Image3D& I0warped,
                   const double Rpar[],
                   const char I1Path[],
                   const char I0Path[]);

    void initialFinalmTRE(double& ind0,
                          double& ind1,
                          const Field3D* v0,
                          const std::vector<vector<double>> tag_mr,
                          const std::vector<vector<double>> tag_us,
                          const bool modeA);

    void initialFinalmTREV2(double& ind0,
                          double& ind1,
                          const Field3D* phiinv,
                          const std::vector<vector<double>> tag_mr,
                          const std::vector<vector<double>> tag_us);

    void initialFinalmTREV3(double& ind0,
                          double& ind1,
                          const Field3D* H0Spatial,
                          const std::vector<vector<double>> tag_mr,
                          const std::vector<vector<double>> tag_us);

    void trilinearInterpolate(std::vector<vector<double>>& Ttag_us,
                              const Field3D* v0,
                              const std::vector<vector<double>> tag_us,
                              const bool modeA);

    void trilinearInterpolateV2(std::vector<vector<double>>& Ttag_us,
                              const Field3D* phiinv,
                              const std::vector<vector<double>> tag_us);

    void provideInitialGuess(Field3D& v0,
                             const std::vector<vector<double>> tag_mr,
                             const std::vector<vector<double>> tag_us,
                             const int modeG);

    void calcDeterminant(double& deTer,
                         const std::vector<std::vector<double>> maTrix,
                         const int dim);

    void calcInverse(std::vector<std::vector<double>>& invmaTrix,
                     const std::vector<std::vector<double>> maTrix);

    void initialEventualmTRE(double& ind0,
                          double& ind1,
                          const Field3D* v0,
                          const Field3D* vest,
                          const std::vector<vector<double>> tag_mr,
                          const std::vector<vector<double>> tag_us);

    void initialEventualmTREV2(double& ind0,
                          double& ind1,
                          const Field3D* v0,
                          const Field3D* phiinv,
                          const std::vector<vector<double>> tag_mr,
                          const std::vector<vector<double>> tag_us);

    void initialEventualmTREV3(double& ind0,
                          double& ind1,
                          const Field3D* v0,
                          const Field3D* vest,
                          const std::vector<vector<double>> tag_mr,
                          const std::vector<vector<double>> tag_us);

    void calcDice_classic(std::vector<double>& diceScores,
                            const Image3D* I1seg,
                            const Image3D* I0seg,
                            const Image3D* deformI0seg);

    void calcDice_struct(std::vector<double>& diceScores,
                            const Image3D* I1segstruct,
                            const Image3D* I0segstruct,
                            const Image3D* deformI0segstruct,
                            const int Nstruct);

    void calcDice_selective(std::vector<double>& diceScores,
                            const Image3D* I1seg,
                            const Image3D* I0seg,
                            const Image3D* deformI0seg);
};


#endif // CALCWARPINGINDEX_H
