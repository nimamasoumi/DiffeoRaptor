#include "CalcWarpingIndex.h"

CalcWarpingIndex::CalcWarpingIndex(){}

CalcWarpingIndex::~CalcWarpingIndex(){}

void CalcWarpingIndex::initialFinalIndex(double& ind0,
                                         double& ind1,
                                         const Field3D *vi0,
                                         const Field3D *v0,
                                         const Image3D *mask)
{
    //The threshold sets the precision of warping index
    float indsi =0, indsf=0, nz=0, tmp=0, TH=100;
    Vec3Di aa;
    Vec3Df cc, dd;
    aa = vi0->size();
    //Field3D* clonevi0;
    //ITKFileIO::LoadField(*clonevi0,"vi0Spatial.mhd");
    Field3D* composedV = new Field3D(vi0->grid(),vi0->memType());
    Opers::ApplyVInv(*composedV,*v0,*vi0,BACKGROUND_STRATEGY_PARTIAL_ZERO);

    for(unsigned int i=0;i<aa[0];i++)
    {
        for(unsigned int j=0;j<aa[1];j++)
        {
            for(unsigned int k=0;k<aa[2];k++)
            {
                cc = vi0->get(i,j,k);
                dd = composedV->get(i,j,k);

                //tmp = abs(dd[0])+abs(dd[1])+abs(dd[2]);
                tmp=mask->get(i,j,k);

                if(tmp==TH)
                {
                    nz=nz+1;
                    indsi = indsi + sqrt(pow(cc[0],2)+pow(cc[1],2)+pow(cc[2],2));
                    indsf = indsf + sqrt(pow(cc[0]+dd[0],2)+pow(cc[1]+dd[1],2)+pow(cc[2]+dd[2],2));
                }
            }
        }
    }

    if(nz!=0)
    {
        ind0=indsi/nz;
        ind1 = indsf/nz;
    }
    else
    {
        std::cout<<"Resolving index ..."<<std::endl;
        indsi=0;
        Opers::Sum(indsi,*vi0);
        ind0=indsi/(aa[0]*aa[1]*aa[2]);
        ind1=ind0;
    }
}

void CalcWarpingIndex::imageWarp(Field3D& vi0spatial,
                                 Image3D& I0warped,
                                 const double Rpar[],
                                 const char I1Path[],
                                 const char I0Path[])
{
    //   Read the Fixed and Moving images.
    typedef itk::Image< float, 3 > ImageType;
    typedef itk::ImageFileReader<ImageType> FixedImageReaderType;
    typedef itk::ImageFileReader<ImageType> MovingImageReaderType;
    FixedImageReaderType::Pointer  fixedImageReader = FixedImageReaderType::New();
    MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

    fixedImageReader->SetFileName(I1Path);
    movingImageReader->SetFileName(I0Path);

    try
    {
      fixedImageReader->Update();
      movingImageReader->Update();
    }
    catch (const char* err)
    {
      std::cerr << "it is me !\n" << std::endl;
      std::cerr << err << std::endl;
      //return EXIT_FAILURE;
    }

    ImageType::Pointer fixedImage = ImageType::New();
    fixedImage = fixedImageReader->GetOutput();
    ImageType::Pointer movingImage = ImageType::New();
    movingImage = movingImageReader->GetOutput();

    //typedef itk::Euler3DTransform<double> RigidTransformType;
    typedef itk::VersorRigid3DTransform<double> RigidTransformType;
    RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
    RigidTransformType::ParametersType par(6);
    typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleImageFilter;
    ResampleImageFilter::Pointer resampleI = ResampleImageFilter::New();
    resampleI->SetSize(movingImageReader->GetOutput()->GetLargestPossibleRegion().GetSize());
    resampleI->SetOutputOrigin(movingImage->GetOrigin());
    resampleI->SetOutputSpacing(movingImage->GetSpacing());
    resampleI->SetOutputDirection(movingImage->GetDirection());
    resampleI->SetDefaultPixelValue(0);
    resampleI->SetInput(movingImage);
    typedef itk::ImageFileWriter<ImageType> ImageFileWriterType;
    ImageFileWriterType::Pointer writerI = ImageFileWriterType::New();
    char ItempPath[100];
    strcpy(ItempPath, "/home/nima/Projects/Deterministic/Code/DiffeoRaptor_RESECT/Itemp.nii.gz");

    // Setting the parameters
    par[0]=Rpar[0];
    par[1]=Rpar[1];
    par[2]=Rpar[2];
    par[3]=Rpar[3];
    par[4]=Rpar[4];
    par[5]=Rpar[5];
    rigidTransform->SetParameters(par);
    //RigidTransformType::ParametersType fpar;
    //fpar[0]=0;fpar[1]=0;fpar[2]=0;
    //rigidTransform->SetFixedParameters(fpar);
    //RigidTransformType::MatrixType aaa= rigidTransform->GetMatrix();
//    std::cout<<Rpar[3]<<Rpar[4]<<Rpar[5]<<endl;
//    for(int ii=0;ii<aaa.RowDimensions+1;ii++)
//    {
// //        for(int jj=0;jj<aaa.ColumnDimensions+1;jj++)
// //        {
// //            std::cout<<aaa[ii][jj]<<endl;
// //        }
//        std::cout<<aaa[ii][0]<<" "<<aaa[ii][1]<<" "<<aaa[ii][2]<<" "<<aaa[ii][3]<<endl;
//    }

    resampleI->SetTransform(rigidTransform);

    writerI->SetInput(resampleI->GetOutput());
    writerI->SetFileName(ItempPath);
    writerI->Update();

    typedef itk::Vector<float, 3> VectorType;
    typedef itk::Image<VectorType, 3> DisplacementFieldImageType;
    typedef itk::TransformToDisplacementFieldFilter<DisplacementFieldImageType, double> DispFilterType;
    DispFilterType::Pointer dispfilter = DispFilterType::New();
//    DisplacementFieldImageType::DirectionType directionV;
//    directionV.SetIdentity();
//    directionV[0][0]=-1;
//    directionV[1][1]=-1;

    dispfilter->UseReferenceImageOn();
    dispfilter->SetReferenceImage(fixedImage);
    dispfilter->SetTransform(rigidTransform);
    //dispfilter->SetOutputDirection(directionV);
    try
    {
        dispfilter->Update();
    }
    catch (const char* err)
    {
        std::cerr << "Exception detected while generating deformation field";
        std::cerr << " : " << err << std::endl;
        //return EXIT_FAILURE;
    }

    DisplacementFieldImageType::Pointer vi0itk = DisplacementFieldImageType::New();
    vi0itk = dispfilter->GetOutput();

//    typedef itk::ChangeInformationImageFilter<DisplacementFieldImageType> changefilterType;
//    changefilterType::Pointer changefilter = changefilterType::New();
//    changefilter->SetInput(vi0itk);
//    changefilter->SetOutputDirection(directionV);
//    changefilter->ChangeDirectionOn();
//    changefilter->UpdateOutputInformation();

    //directionV=changefilter->GetOutputDirection();
    //std::cout<<directionV[1][0]<<std::endl;
    //DisplacementFieldImageType::ConstPointer outputV=changefilter->GetOutput();
    //std::cout<<"\nChanged Image: "<<outputV<<std::endl;

    typedef itk::ImageFileWriter<DisplacementFieldImageType> writerType;
    writerType::Pointer writer=writerType::New();
    writer->SetFileName("/home/nima/Projects/Deterministic/Code/DiffeoRaptor_RESECT/vitkTemp.mhd");
    writer->SetInput(vi0itk);
    writer->Update();
    ITKFileIO::LoadField(vi0spatial,"/home/nima/Projects/Deterministic/Code/DiffeoRaptor_RESECT/vitkTemp.mhd");
    ITKFileIO::LoadImage(I0warped, ItempPath);

    //ITKFileIO::SaveField(vi0spatial,"/home/nima/Projects/Deterministic/Code/DiffeoRaptor_RESECT/vpycaTemp.mhd");

//    typedef itk::ImageFileReader<DisplacementFieldImageType> readerType;
//    readerType::Pointer reader =readerType::New();
//    reader->SetFileName("/home/nima/Projects/Deterministic/Code/DiffeoRaptor_RESECT/v0Spatial.mhd");
//    reader->Update();
//    DisplacementFieldImageType::ConstPointer outputVv=reader->GetOutput();
//    std::cout<<"\nChanged Image: "<<outputVv<<std::endl;
//    //ImageType::IndexType aaa={10, 20, 30};
    //VectorType tmptmp =vi0itk->GetPixel(aaa);
    //std::cout<<tmptmp<<endl;
    //Field3D *vi0Spatial = new Field3D(grid, mType);

    // itk disp values should be assigned into this field
    // please change the following code for safety
//    typedef itk::ImageRegionConstIterator<DisplacementFieldImageType> ImageRegionItrType;
//    ImageRegionItrType it3(vi0itk, vi0itk->GetLargestPossibleRegion());
//    it3.GoToBegin();
//    VectorType tempVec3;
//    Vec3Df cc;
//    unsigned int i=0;

//    while( !it3.IsAtEnd() )
//    {
//        tempVec3 = it3.Get();
//        cc[0]=tempVec3[0];
//        cc[1]=tempVec3[1];
//        cc[2]=tempVec3[2];
//        vi0spatial.set(i, cc);
//        ++it3;
//        i=i+1;
//    }
}

void CalcWarpingIndex::initialFinalmTRE(double& ind0,
                      double& ind1,
                      const Field3D* v0,
                      const std::vector<std::vector<double> > tag_mr,
                      const std::vector<std::vector<double> > tag_us,
                      const bool modeA)
{
    ind0=0;ind1=0;
    int N=tag_mr.size();
    std::vector<double> sum1(N),sum2(N);
    std::vector<std::vector<double> > Ttag_us(N); 
    //std::ofstream fileDir1("tag_us.txt");
    std::ofstream fileDir2("Ttag_us.txt");
    fileDir2<<"x Tx y Ty z Tz di df"<<std::endl;
    trilinearInterpolate(Ttag_us,v0,tag_us,modeA);
    for(int i=0;i<N;i++)
    {
       for(int j=0;j<tag_mr[0].size();j++)
       {            
            sum1[i]+=pow(tag_mr[i][j]-tag_us[i][j],2);
            sum2[i]+=pow(tag_mr[i][j]-Ttag_us[i][j],2);
            //fileDir1<<pow(Ttag_us[i][j]-tag_mr[i][j],2)<<" ";
            fileDir2<<tag_us[i][j]<<" "<<Ttag_us[i][j]<<" ";
       }       
       ind0+=sqrt(sum1[i])/N;
       ind1+=sqrt(sum2[i])/N;
//       fileDir1<<sqrt(sum1[i])<<std::endl;
       fileDir2<<sqrt(sum1[i])<<" "<<sqrt(sum2[i])<<std::endl;
    }
//    fileDir1.close();
    fileDir2.close();
}

void CalcWarpingIndex::initialFinalmTREV2(double& ind0,
                      double& ind1,
                      const Field3D* phiinv,
                      const std::vector<vector<double>> tag_mr,
                      const std::vector<vector<double>> tag_us)
{
    ind0=0;ind1=0;
    int N=tag_mr.size();
    std::vector<double> sum1(N),sum2(N);
    std::vector<std::vector<double> > Ttag_us(N);
    std::ofstream fileDir1("tag_us.txt");
    std::ofstream fileDir2("Ttag_us.txt");
    trilinearInterpolateV2(Ttag_us,phiinv,tag_us);
    for(int i=0;i<N;i++)
    {
       for(int j=0;j<tag_mr[0].size();j++)
       {
            sum1[i]+=pow(tag_mr[i][j]-tag_us[i][j],2);
            sum2[i]+=pow(tag_mr[i][j]-Ttag_us[i][j],2);
       }
       ind0+=sqrt(sum1[i])/N;
       ind1+=sqrt(sum2[i])/N;
       fileDir1<<sqrt(sum1[i])<<std::endl;fileDir2<<sqrt(sum2[i])<<std::endl;
    }
    fileDir1.close();fileDir2.close();
}

void CalcWarpingIndex::initialFinalmTREV3(double& ind0,
                      double& ind1,
                      const Field3D* H0Spatial,
                      const std::vector<vector<double>> tag_mr,
                      const std::vector<vector<double>> tag_us)
{

}

void CalcWarpingIndex::trilinearInterpolate(std::vector<std::vector<double>>& Ttag_us,
                          const Field3D* v0,
                          const std::vector<std::vector<double>> tag_us,
                          const bool modeA)
{
    int dims=3;
    int N=tag_us.size();
    for(int i=0;i<N;i++)
    {
        std::vector<double> p (dims),p0 (dims),p1 (dims),pd (dims);
        std::vector<double> intp (dims);
        Ttag_us[i].resize(3,0);
        p[0]=tag_us[i][0];
        p[1]=tag_us[i][1];
        p[2]=tag_us[i][2];

        p0[0]=std::floor(p[0]);
        p0[1]=std::floor(p[1]);
        p0[2]=std::floor(p[2]);

        p1[0]=std::ceil(p[0]);
        p1[1]=std::ceil(p[1]);
        p1[2]=std::ceil(p[2]);

        if(p0[0]==p1[0])
        {
            p1[0]=p1[0]+1;
        }
        if(p0[1]==p1[1])
        {
            p1[1]=p1[1]+1;
        }
        if(p0[2]==p1[2])
        {
            p1[2]=p1[2]+1;
        }
//        std::cout<<p[0]<<","<<p[1]<<","<<p[2]<<std::endl;
//        std::cout<<p1[0]<<","<<p1[1]<<","<<p1[2]<<std::endl;
//        std::cout<<p0[0]<<","<<p0[1]<<","<<p0[2]<<std::endl;

        pd[0]=p[0]-p0[0];
        pd[1]=p[1]-p0[1];
        pd[2]=p[2]-p0[2];

        for(int j=0;j<intp.size();j++)
        {
            Vec3Df aa,bb;
            double c00,c01,c10,c11,c0,c1;

            aa=v0->get(int(p0[0]),int(p0[1]),int(p0[2]));
            bb=v0->get(int(p1[0]),int(p0[1]),int(p0[2]));
            c00=(aa[j]*(1-pd[0]))+(bb[j]*pd[0]);
            aa=v0->get(int(p0[0]),int(p0[1]),int(p1[2]));
            bb=v0->get(int(p1[0]),int(p0[1]),int(p1[2]));
            c01=(aa[j]*(1-pd[0]))+(bb[j]*pd[0]);
            aa=v0->get(int(p0[0]),int(p1[1]),int(p0[2]));
            bb=v0->get(int(p1[0]),int(p1[1]),int(p0[2]));
            c10=(aa[j]*(1-pd[0]))+(bb[j]*pd[0]);
            aa=v0->get(int(p0[0]),int(p1[1]),int(p1[2]));
            bb=v0->get(int(p1[0]),int(p1[1]),int(p1[2]));
            c11=(aa[j]*(1-pd[0]))+(bb[j]*pd[0]);

            c0=(c00*(1-pd[1]))+(c10*pd[1]);
            c1=(c01*(1-pd[1]))+(c11*pd[1]);

            intp[j]=(c0*(1-pd[2]))+(c1*pd[2]);
            //aa=v0->get(int(p0[0]),int(p0[1]),int(p0[2]));
            //intp[j]=aa[j];
            //std::cout<<intp[j]<<std::endl;
            if(modeA)
            {
                Ttag_us[i][j]=p[j]-intp[j];
            }
            else
            {
                Ttag_us[i][j]=p[j]+intp[j];
            }
        }
        //std::co`ut<<intp[0]<<","<<intp[1]<<","<<intp[2]<<std::endl;
//        Ttag_us[i][0]=p[0]+intp[0];
//        Ttag_us[i][1]=p[1]+intp[1];
//        Ttag_us[i][2]=p[2]+intp[2];
    }
}

void CalcWarpingIndex:: trilinearInterpolateV2(std::vector<vector<double>>& Ttag_us,
                          const Field3D* phiinv,
                          const std::vector<vector<double>> tag_us)
{
    int dims=3;
    int N=tag_us.size();
    for(int i=0;i<N;i++)
    {
        std::vector<double> p (dims),p0 (dims),p1 (dims),pd (dims);
        std::vector<double> intp (dims);
        Ttag_us[i].resize(3,0);
        p[0]=tag_us[i][0];
        p[1]=tag_us[i][1];
        p[2]=tag_us[i][2];

        p0[0]=std::floor(p[0]);
        p0[1]=std::floor(p[1]);
        p0[2]=std::floor(p[2]);

        p1[0]=std::ceil(p[0]);
        p1[1]=std::ceil(p[1]);
        p1[2]=std::ceil(p[2]);

        if(p0[0]==p1[0])
        {
            p1[0]=p1[0]+1;
        }
        if(p0[1]==p1[1])
        {
            p1[1]=p1[1]+1;
        }
        if(p0[2]==p1[2])
        {
            p1[2]=p1[2]+1;
        }

        pd[0]=p[0]-p0[0];
        pd[1]=p[1]-p0[1];
        pd[2]=p[2]-p0[2];

        for(int j=0;j<intp.size();j++)
        {
            Vec3Df aa,bb;
            double c00,c01,c10,c11,c0,c1;

            aa=phiinv->get(int(p0[0]),int(p0[1]),int(p0[2]));
            bb=phiinv->get(int(p1[0]),int(p0[1]),int(p0[2]));
            c00=(aa[j]*(1-pd[0]))+(bb[j]*pd[0]);
            aa=phiinv->get(int(p0[0]),int(p0[1]),int(p1[2]));
            bb=phiinv->get(int(p1[0]),int(p0[1]),int(p1[2]));
            c01=(aa[j]*(1-pd[0]))+(bb[j]*pd[0]);
            aa=phiinv->get(int(p0[0]),int(p1[1]),int(p0[2]));
            bb=phiinv->get(int(p1[0]),int(p1[1]),int(p0[2]));
            c10=(aa[j]*(1-pd[0]))+(bb[j]*pd[0]);
            aa=phiinv->get(int(p0[0]),int(p1[1]),int(p1[2]));
            bb=phiinv->get(int(p1[0]),int(p1[1]),int(p1[2]));
            c11=(aa[j]*(1-pd[0]))+(bb[j]*pd[0]);

            c0=(c00*(1-pd[1]))+(c10*pd[1]);
            c1=(c01*(1-pd[1]))+(c11*pd[1]);

            intp[j]=(c0*(1-pd[2]))+(c1*pd[2]);
            Ttag_us[i][j]=intp[j];

        }
    }
}

void CalcWarpingIndex::provideInitialGuess(Field3D& v0,
                                           const std::vector<vector<double>> tag_mr,
                                           const std::vector<vector<double>> tag_us,
                                           const int modeG)
{
    int dims=3;
    int N=tag_us.size();
    std::vector<std::vector<double>> MtM (dims+1), MtF(dims+1), invMtM(dims+1), affT(dims+1), guessT(dims+1);
    double auUS=0,auUST=0,auMR=0;
    double thetax=0,thetay=0,thetaz=0,tx=0,ty=0,tz=0,gtx=0,gty=0,gtz=0,
            gx=0,gy=0,gz=0,ex=0,ey=0,ez=0,resx=0,resy=0,resz=0,theSum=0,stepSize=0.05,stepSizeT=0.0005;
    int numIter=50;

    for(int k=0;k<dims+1;k++)
    {
        guessT[k].resize(dims+1,0);
    }
    switch(modeG)
    {
    case 0://affine
        std::cout<<"\nAffine Initial Guess (default)"<<std::endl;
        for(int ii=0; ii<dims+1;ii++)
        {
            MtM[ii].resize(dims+1,0);
            MtF[ii].resize(dims+1,0);
            invMtM[ii].resize(dims+1,0);
            affT[ii].resize(dims+1,0);
            //guessT[ii].resize(dims+1,0);
            for(int jj=0;jj<dims+1;jj++)
            {
                for(int kk=0;kk<N;kk++)
                {
                    if((ii==3)&&(jj==3))
                    {
                        auUST=1;auUS=1;auMR=1;
                    }else if((ii==3)&&(jj!=3))
                    {
                        auUST = 1;auUS=tag_us[kk][jj];auMR=tag_mr[kk][jj];
                    }else if((ii!=3)&&(jj==3))
                    {
                        auUST = tag_us[kk][ii];auUS=1;auMR=1;
                    }else
                    {
                        auUST = tag_us[kk][ii];auUS=tag_us[kk][jj];auMR=tag_mr[kk][jj];
                    }
                    MtM[ii][jj]+=auUST*auUS;
                    MtF[ii][jj]+=auUST*auMR;
                }
            }
        }
        calcInverse(invMtM,MtM);
        for(int i=0;i<dims+1;i++)
        {
            for(int j=0;j<dims+1;j++)
            {
                for(int k=0;k<dims+1;k++)
                {
                    affT[i][j]+=invMtM[i][k]*MtF[k][j];
                }
                guessT[i][j]=affT[i][j];
            }
        }
        break;

    case 1:
        std::cout<<"\nRigid Initial Guess"<<std::endl;
        for(int ii=0;ii<numIter;ii++)
        {
            for(int i=0;i<N;i++)
            {
                ex=(cos(thetay)*cos(thetaz))*tag_us[i][0]+
                        (sin(thetax)*sin(thetay)*cos(thetaz)-cos(thetax)*sin(thetaz))*tag_us[i][1]+
                        (cos(thetax)*sin(thetay)*cos(thetaz)+sin(thetax)*sin(thetaz))*tag_us[i][2]+
                        tx-tag_mr[i][0];
                ey=(cos(thetay)*sin(thetaz))*tag_us[i][0]+
                        (sin(thetax)*sin(thetay)*sin(thetaz)+cos(thetax)*cos(thetaz))*tag_us[i][1]+
                        (cos(thetax)*sin(thetay)*sin(thetaz)-sin(thetax)*cos(thetaz))*tag_us[i][2]+
                        ty-tag_mr[i][1];
                ez=(-sin(thetay))*tag_us[i][0]+
                        (sin(thetax)*cos(thetay))*tag_us[i][1]+
                        (cos(thetax)*cos(thetay))*tag_us[i][2]+
                        tz-tag_mr[i][2];
                theSum=pow(ex,2)+pow(ey,2)+pow(ez,2);
                resx=((cos(thetax)*sin(thetay)*cos(thetaz)+sin(thetax)*sin(thetaz))*tag_us[i][1]+
                        (-sin(thetax)*sin(thetay)*cos(thetaz)+cos(thetax)*sin(thetaz))*tag_us[i][2])*ex+
                        ((cos(thetax)*sin(thetay)*sin(thetaz)-sin(thetax)*cos(thetaz))*tag_us[i][1]+
                        (-sin(thetax)*sin(thetay)*sin(thetaz)-cos(thetax)*cos(thetaz))*tag_us[i][2])*ey+
                        ((cos(thetax)*cos(thetay))*tag_us[i][1]+
                        (-sin(thetax)*cos(thetay))*tag_us[i][2])*ez;
                resy=((-sin(thetay)*cos(thetaz))*tag_us[i][0]+
                        (sin(thetax)*cos(thetay)*cos(thetaz))*tag_us[i][1]+
                        (cos(thetax)*cos(thetay)*cos(thetaz))*tag_us[i][2])*ex+
                        ((-sin(thetay)*sin(thetaz))*tag_us[i][0]+
                        (sin(thetax)*cos(thetay)*sin(thetaz))*tag_us[i][1]+
                        (cos(thetax)*cos(thetay)*sin(thetaz))*tag_us[i][2])*ey+
                        ((-cos(thetay))*tag_us[i][0]+
                        (-sin(thetax)*sin(thetay))*tag_us[i][1]+
                        (-cos(thetax)*sin(thetay))*tag_us[i][2])*ez;
                resz=((-cos(thetay)*sin(thetaz))*tag_us[i][0]+
                        (-sin(thetax)*sin(thetay)*sin(thetaz)-cos(thetax)*cos(thetaz))*tag_us[i][1]+
                        (-cos(thetax)*sin(thetay)*sin(thetaz)+sin(thetax)*cos(thetaz))*tag_us[i][2])*ex+
                        ((cos(thetay)*cos(thetaz))*tag_us[i][0]+
                        (sin(thetax)*sin(thetay)*cos(thetaz)-cos(thetax)*sin(thetaz))*tag_us[i][1]+
                        (cos(thetax)*sin(thetay)*cos(thetaz)+sin(thetax)*sin(thetaz))*tag_us[i][2])*ey;

                gtx+=resx/(sqrt(theSum));
                gty+=resy/(sqrt(theSum));
                gtz+=resz/(sqrt(theSum));
                gx+=ex/(sqrt(theSum));
                gy+=ey/(sqrt(theSum));
                gz+=ez/(sqrt(theSum));
            }
            thetax=thetax-(gtx/N)*stepSizeT;
            thetay=thetay-(gty/N)*stepSizeT;
            thetaz=thetaz-(gtz/N)*stepSizeT;
            tx=tx-(gx/N)*stepSize;
            ty=ty-(gy/N)*stepSize;
            tz=tz-(gz/N)*stepSize;
        }
        guessT[0][0]=cos(thetay)*cos(thetaz);
        guessT[0][1]=cos(thetay)*sin(thetaz);
        guessT[0][2]=-sin(thetay);
        guessT[0][3]=0;
        guessT[1][0]=sin(thetax)*sin(thetay)*cos(thetaz)-cos(thetax)*sin(thetaz);
        guessT[1][1]=sin(thetax)*sin(thetay)*sin(thetaz)+cos(thetax)*cos(thetaz);
        guessT[1][2]=sin(thetax)*cos(thetay);
        guessT[1][3]=0;
        guessT[2][0]=cos(thetax)*sin(thetay)*cos(thetaz)+sin(thetax)*sin(thetaz);
        guessT[2][1]=cos(thetax)*sin(thetay)*sin(thetaz)-sin(thetax)*cos(thetaz);
        guessT[2][2]=cos(thetax)*cos(thetay);
        guessT[2][3]=0;
        guessT[3][0]=tx;
        guessT[3][1]=ty;
        guessT[3][2]=tz;
        guessT[3][3]=1;
        break;

    case 2:
        std::cout<<"\nTranslation Initial Guess"<<std::endl;
        for(int ii=0;ii<numIter;ii++)
        {
            for(int i=0;i<N;i++)
            {
                theSum=pow(tag_mr[i][0]-tag_us[i][0]-tx,2)+pow(tag_mr[i][1]-tag_us[i][1]-ty,2);+pow(tag_mr[i][2]-tag_us[i][2]-tz,2);
                gx+=(tag_us[i][0]+tx-tag_mr[i][0])/(sqrt(theSum));
                gy+=(tag_us[i][1]+ty-tag_mr[i][1])/(sqrt(theSum));
                gz+=(tag_us[i][2]+tz-tag_mr[i][2])/(sqrt(theSum));
            }
            tx=tx-(gx/N)*stepSize;
            ty=ty-(gy/N)*stepSize;
            tz=tz-(gz/N)*stepSize;
        }
        guessT[0][0]=1;guessT[1][1]=1;guessT[2][2]=1;guessT[3][0]=tx;
        guessT[0][1]=0;guessT[1][0]=0;guessT[2][0]=0;guessT[3][1]=ty;
        guessT[0][2]=0;guessT[1][2]=0;guessT[2][1]=0;guessT[3][2]=tz;
        guessT[0][3]=0;guessT[1][3]=0;guessT[2][3]=0;guessT[3][3]=1;
        break;

    case 3:
        std::cout<<"\nExternal Initial Guess"<<std::endl;
        break;

    case 4:
        std::cout<<"\nNo Initial Guess"<<std::endl;
        guessT[0][0]=1;guessT[1][1]=1;guessT[2][2]=1;guessT[3][0]=0;
        guessT[0][1]=0;guessT[1][0]=0;guessT[2][0]=0;guessT[3][1]=0;
        guessT[0][2]=0;guessT[1][2]=0;guessT[2][1]=0;guessT[3][2]=0;
        guessT[0][3]=0;guessT[1][3]=0;guessT[2][3]=0;guessT[3][3]=1;
        break;
    }

    Vec3Di aa;
    Vec3Df bb;
    aa = v0.size();
    for(unsigned int i=0;i<aa[0];i++)
    {
        for(unsigned int j=0;j<aa[1];j++)
        {
            for(unsigned int k=0;k<aa[2];k++)
            {
                bb.x = (i*guessT[0][0])+(j*guessT[1][0])+(k*guessT[2][0])+guessT[3][0]-i;
                bb.y = (i*guessT[0][1])+(j*guessT[1][1])+(k*guessT[2][1])+guessT[3][1]-j;
                bb.z = (i*guessT[0][2])+(j*guessT[1][2])+(k*guessT[2][2])+guessT[3][2]-k;
                v0.set(i,j,k,bb);
            }
        }
    }
//    std::cout<<"\n"<<deT<<std::endl;
//    std::cout<<MtM[0][0]<<" "<<MtM[0][1]<<" "<<MtM[0][2]<<" "<<MtM[0][3]<<std::endl;
//    std::cout<<MtM[1][0]<<" "<<MtM[1][1]<<" "<<MtM[1][2]<<" "<<MtM[1][3]<<std::endl;
//    std::cout<<MtM[2][0]<<" "<<MtM[2][1]<<" "<<MtM[2][2]<<" "<<MtM[2][3]<<std::endl;
//    std::cout<<MtM[3][0]<<" "<<MtM[3][1]<<" "<<MtM[3][2]<<" "<<MtM[3][3]<<std::endl;

//    std::cout<<invMtM[0][0]<<" "<<invMtM[0][1]<<" "<<invMtM[0][2]<<" "<<invMtM[0][3]<<std::endl;
//    std::cout<<invMtM[1][0]<<" "<<invMtM[1][1]<<" "<<invMtM[1][2]<<" "<<invMtM[1][3]<<std::endl;
//    std::cout<<invMtM[2][0]<<" "<<invMtM[2][1]<<" "<<invMtM[2][2]<<" "<<invMtM[2][3]<<std::endl;
//    std::cout<<invMtM[3][0]<<" "<<invMtM[3][1]<<" "<<invMtM[3][2]<<" "<<invMtM[3][3]<<std::endl;
//    std::cout<<affT[0][0]<<" "<<affT[0][1]<<" "<<affT[0][2]<<" "<<affT[0][3]<<std::endl;
//    std::cout<<affT[1][0]<<" "<<affT[1][1]<<" "<<affT[1][2]<<" "<<affT[1][3]<<std::endl;
//    std::cout<<affT[2][0]<<" "<<affT[2][1]<<" "<<affT[2][2]<<" "<<affT[2][3]<<std::endl;
//    std::cout<<affT[3][0]<<" "<<affT[3][1]<<" "<<affT[3][2]<<" "<<affT[3][3]<<std::endl;
}

void CalcWarpingIndex::calcDeterminant(double& deTer,
                     const std::vector<std::vector<double>> maTrix,
                     const int dim)
{
    deTer =0;
    if(dim==2)
    {
        deTer=maTrix[0][0]*maTrix[1][1]-maTrix[0][1]*maTrix[1][0];
    }
    else
    {
        std::vector<std::vector<double>> submaTrix (dim-1);
        for(int i=0;i<dim-1;i++)
        {
            submaTrix[i].resize(dim-1,0);
        }
        for(int x=0;x<dim;x++)
        {
            int subi=0;
            for(int i=1;i<dim;i++)
            {
                int subj=0;
                for(int j=0;j<dim;j++)
                {
                    if(j==x)continue;
                    submaTrix[subi][subj]=maTrix[i][j];
                    subj++;
                }
                subi++;
            }
            double tempD=0;
            calcDeterminant(tempD,submaTrix,dim-1);
            deTer+=std::pow(-1,x)*maTrix[0][x]*tempD;
        }
    }
}

void CalcWarpingIndex::calcInverse(std::vector<std::vector<double>>& invmaTrix,
                 const std::vector<std::vector<double>> maTrix)
{
    int dim = maTrix.size();
    double deTer=0;
    std::vector<std::vector<double>> TmaTrix(dim),submaTrix(dim-1);
    calcDeterminant(deTer,maTrix,dim);
    for(int i=0;i<dim;i++)
    {
        TmaTrix[i].resize(dim,0);
        for(int j=0;j<dim;j++)
        {
            TmaTrix[i][j]=maTrix[j][i];
        }
        if(i<dim-1)
        {
            submaTrix[i].resize(dim-1,0);
        }
    }

    for(int x=0;x<dim;x++)
    {
        for(int y=0;y<dim;y++)
        {
            int subi=0;
            for(int i=0;i<dim;i++)
            {
                if(i==x)continue;
                int subj=0;
                for(int j=0;j<dim;j++)
                {
                    if(j==y)continue;
                    submaTrix[subi][subj]=TmaTrix[i][j];
                    subj++;
                }
                subi++;
            }
            double subdeT=0;
            calcDeterminant(subdeT,submaTrix,dim-1);
            invmaTrix[x][y]=std::pow(-1,x+y)*(subdeT/deTer);
        }
    }
}

void CalcWarpingIndex::initialEventualmTRE(double& ind0,
                      double& ind1,
                      const Field3D* v0,
                      const Field3D* vest,
                      const std::vector<std::vector<double> > tag_mr,
                      const std::vector<std::vector<double> > tag_us)
{
    ind0=0;ind1=0;
    int N=tag_mr.size();
    std::vector<double> sum1(N),sum2(N);
    std::vector<std::vector<double> > Ttag_us(N),uTtag_us(N);
    std::ofstream fileDir1("uTtag_us.txt");
    trilinearInterpolate(Ttag_us,v0,tag_us,false);
    trilinearInterpolate(uTtag_us,vest,Ttag_us,true);
    for(int i=0;i<N;i++)
    {
       for(int j=0;j<tag_mr[0].size();j++)
       {
            sum1[i]+=pow(tag_mr[i][j]-tag_us[i][j],2);
            sum2[i]+=pow(tag_mr[i][j]-uTtag_us[i][j],2);
            //fileDir1<<pow(uTtag_us[i][j]-tag_mr[i][j],2)<<" ";
       }
       ind0+=sqrt(sum1[i])/N;
       ind1+=sqrt(sum2[i])/N;       
       fileDir1<<sqrt(sum2[i])<<std::endl;
    }
    fileDir1.close();
}

void CalcWarpingIndex::initialEventualmTREV2(double& ind0,
                      double& ind1,
                      const Field3D* v0,
                      const Field3D* phiinv,
                      const std::vector<vector<double>> tag_mr,
                      const std::vector<vector<double>> tag_us)
{
    ind0=0;ind1=0;
    int N=tag_mr.size();
    std::vector<double> sum1(N),sum2(N);
    std::vector<std::vector<double> > Ttag_us(N),uTtag_us(N);
    std::ofstream fileDir1("uTtag_us.txt");
    trilinearInterpolate(Ttag_us,v0,tag_us,false);
    trilinearInterpolateV2(uTtag_us,phiinv,Ttag_us);
    for(int i=0;i<N;i++)
    {
       for(int j=0;j<tag_mr[0].size();j++)
       {
            sum1[i]+=pow(tag_mr[i][j]-tag_us[i][j],2);
            sum2[i]+=pow(tag_mr[i][j]-uTtag_us[i][j],2);
       }
       ind0+=sqrt(sum1[i])/N;
       ind1+=sqrt(sum2[i])/N;
       fileDir1<<sqrt(sum2[i])<<std::endl;
    }
    fileDir1.close();
}

void CalcWarpingIndex::initialEventualmTREV3(double& ind0,
                      double& ind1,
                      const Field3D* v0,
                      const Field3D* vest,
                      const std::vector<vector<double>> tag_mr,
                      const std::vector<vector<double>> tag_us)
{
    ind0=0;ind1=0;
    int N=tag_mr.size();
    std::vector<double> sum1(N),sum2(N);
    std::vector<std::vector<double> > Ttag_us(N),uTtag_us(N);
    std::ofstream fileDir1("uTtag_us.txt");
    trilinearInterpolate(Ttag_us,v0,tag_us,false);
    trilinearInterpolate(uTtag_us,vest,Ttag_us,true);
    for(int i=0;i<N;i++)
    {
       for(int j=0;j<tag_mr[0].size();j++)
       {
            sum1[i]+=pow(tag_mr[i][j]-tag_us[i][j],2);
            sum2[i]+=pow(tag_mr[i][j]-uTtag_us[i][j],2);
       }
       ind0+=sqrt(sum1[i])/N;
       ind1+=sqrt(sum2[i])/N;
       fileDir1<<sqrt(sum2[i])<<std::endl;
    }
    fileDir1.close();
}

void CalcWarpingIndex::calcDice_classic(std::vector<double>& diceScores,
                            const Image3D* I1seg,
                            const Image3D* I0seg,
                            const Image3D* deformI0seg)
{
    int initInt=0, finalInt=0, numI0=0, numI1=0, numdeformI0=0;
    int finalcsf=0, numI0csf=0, numI1csf=0, numdeformI0csf=0;
    int finalwm=0, numI0wm=0, numI1wm=0, numdeformI0wm=0;
    int finalgm=0, numI0gm=0, numI1gm=0, numdeformI0gm=0;
    double initDice=0, finalDice=0, csfDice=0, wmDice=0, gmDice=0;
    double rnI1=0, rnI0=0, rndI0=0;
    GridInfo griD = I1seg->grid();
    Vec3Di mSize = griD.size();
    int fsx = mSize.x;int fsy = mSize.y;int fsz = mSize.z;
    for (int zz = 0; zz < fsx; zz++)
    {
       for(int jj=0;jj<fsy;jj++)
       {
          for(int kk=0;kk<fsz;kk++)
          {
              rnI1=std::round(I1seg->get(zz,jj,kk));
              rnI0=std::round(I0seg->get(zz,jj,kk));
              rndI0=std::round(deformI0seg->get(zz,jj,kk));
              if((rnI1==rnI0)&&(rnI1!=0)&&(rnI0!=0))
              {
                  initInt++;
              }
              if((rnI1==rndI0)&&(rnI1!=0)&&(rndI0!=0))
              {
                  finalInt++;
                  switch((int)rnI1)
                  {
                    case 1:
                      finalcsf++;
                      break;
                    case 2:
                      finalgm++;
                      break;
                    case 3:
                      finalwm++;
                      break;
                    default:
                      finalwm++;
                      break;
                  }
              }
              if(rnI1!=0)
              {
                  numI1++;
                  switch((int)rnI1)
                  {
                    case 1:
                      numI1csf++;
                      break;
                    case 2:
                      numI1gm++;
                      break;
                    case 3:
                      numI1wm++;
                      break;
                    default:
                      numI1wm++;
                      break;
                  }
              }
              if(rnI0!=0)
              {
                  numI0++;
                  switch((int)rnI0)
                  {
                    case 1:
                      numI0csf++;
                      break;
                    case 2:
                      numI0gm++;
                      break;
                    case 3:
                      numI0wm++;
                      break;
                    default:
                      numI0wm++;
                      break;
                  }
              }
              if(rndI0!=0)
              {
                  numdeformI0++;
                  switch((int)rndI0)
                  {
                    case 1:
                      numdeformI0csf++;
                      break;
                    case 2:
                      numdeformI0gm++;
                      break;
                    case 3:
                      numdeformI0wm++;
                      break;
                    default:
                      numdeformI0wm++;
                      break;
                  }
              }
          }
       }
    }
    initDice=(2*(double)initInt)/((double)(numI0+numI1));
    finalDice=(2*(double)finalInt)/((double)(numdeformI0+numI1));
    csfDice=(2*(double)finalcsf)/((double)(numdeformI0csf+numI1csf));
    wmDice=(2*(double)finalwm)/((double)(numdeformI0wm+numI1wm));
    gmDice=(2*(double)finalgm)/((double)(numdeformI0gm+numI1gm));
    diceScores[0]=initDice;
    diceScores[1]=finalDice;
    diceScores[2]=csfDice;
    diceScores[3]=wmDice;
    diceScores[4]=gmDice;
}

void CalcWarpingIndex::calcDice_struct(std::vector<double>& diceScores,
                        const Image3D* I1segstruct,
                        const Image3D* I0segstruct,
                        const Image3D* deformI0segstruct,
                        const int Nstruct)
{
    std::vector<int> intersectCount (Nstruct, 0); //intersectCount[0]-> initial , intersectCount[Nstruct]-> final
    std::vector<int> numI1 (Nstruct);
    std::vector<int> numI0 (Nstruct);
    double rnI1=0, rnI0=0, rndI0=0;
    GridInfo griD = I1segstruct->grid();
    Vec3Di mSize = griD.size();
    int fsx = mSize.x;int fsy = mSize.y;int fsz = mSize.z;
    for (int zz = 0; zz < fsx; zz++)
    {
       for(int jj=0;jj<fsy;jj++)
       {
          for(int kk=0;kk<fsz;kk++)
          {
//              rnI1=std::round(I1segstruct->get(zz,jj,kk));
//              rnI0=std::round(I0segstruct->get(zz,jj,kk));
//              rndI0=std::round(deformI0segstruct->get(zz,jj,kk));

              rnI1=I1segstruct->get(zz,jj,kk);
              rnI0=I0segstruct->get(zz,jj,kk);
              rndI0=deformI0segstruct->get(zz,jj,kk);

              if((rnI1==rnI0)&&(rnI1!=0)&&(rnI0!=0))
              {
                  intersectCount[0]++;
              }
              if((rnI1==rndI0)&&(rnI1!=0)&&(rndI0!=0))
              {
                  //std::cout<<I1segstruct->get(zz,jj,kk)<<" "<<I0segstruct->get(zz,jj,kk)<<" "<<deformI0segstruct->get(zz,jj,kk)<<std::endl;
                  //std::cout<<rnI1<<" "<<rnI0<<" "<<rndI0<<" "<<zz<<" "<<jj<<" "<<kk<<std::endl;
                  intersectCount[Nstruct-1]++;
                  for(int ii=1;ii<Nstruct-1;ii++)
                  {
                      if(ii==(int)rnI1)
                      {
                          intersectCount[ii]++;
                      }
                  }
                  if((int)rnI1==Nstruct-1)
                  {
                      intersectCount[Nstruct-2]++;
                  }
              }
              if(rnI1!=0)
              {
                  numI1[0]++;
                  numI1[Nstruct-1]++;
                  for(int ii=1;ii<Nstruct-1;ii++)
                  {
                      if(ii==(int)rnI1)
                      {
                          numI1[ii]++;
                      }
                  }
                  if((int)rnI1==Nstruct-1)
                  {
                      numI1[Nstruct-2]++;
                  }

              }
              if(rnI0!=0)
              {
                  numI0[0]++;

              }
              if(rndI0!=0)
              {
                  numI0[Nstruct-1]++;
                  for(int ii=1;ii<Nstruct-1;ii++)
                  {
                      if(ii==(int)rndI0)
                      {
                          numI0[ii]++;
                      }
                  }
                  if((int)rndI0==Nstruct-1)
                  {
                      numI0[Nstruct-2]++;
                  }
              }
          }
       }
    }    
    for(int ii=0;ii<Nstruct;ii++)
    {
        //std::cout<<2*intersectCount[ii]<<" "<<numI0[ii]<<" "<<numI1[ii]<<std::endl;
        diceScores[ii]=
                (2*(double)intersectCount[ii])/((double)(numI0[ii]+numI1[ii]));
    }
}

void CalcWarpingIndex::calcDice_selective(std::vector<double>& diceScores,
                                          const Image3D* I1seg,
                                          const Image3D* I0seg,
                                          const Image3D* deformI0seg)
{
    // Dynamic Label Stat
    Vec2D<float> I0mm,I1mm,I0dmm;
    Opers::MaxMin(I0mm,*I0seg);
    Opers::MaxMin(I1mm,*I1seg);
    Opers::MaxMin(I0dmm,*deformI0seg);
    int maxVol=(int)(std::max(std::max(I0mm.x,I1mm.x),I0dmm.x));
    GridInfo griD = I0seg->grid();
    Vec3Di mSize = griD.size();
    int fsx = mSize.x;int fsy = mSize.y;int fsz = mSize.z;
//    int imageHist[(int)(maxVol+1)]={0};
//    int singularPoint=0;//, nonzeroLabel=0;

    int Nstruct =maxVol+2;
    std::vector<int> intersectCount (Nstruct, 0); //intersectCount[0]-> initial , intersectCount[Nstruct-1]-> final
    std::vector<int> numI1 (Nstruct, 0);
    std::vector<int> numI0 (Nstruct, 0);
    double rnI1=0, rnI0=0, rndI0=0;

    // Finding the number of pixels in each label
//    for(int i=0;i<fsx;i++)
//    {
//        for(int j=0;j<fsy;j++)
//        {
//            for(int k=0;k<fsz;k++)
//            {
//                if(isfinite(I0seg->get(i,j,k))||isfinite(I1seg->get(i,j,k))||isfinite(deformI0seg->get(i,j,k)))
//                {
//                  imageHist[(int)I0seg->get(i,j,k)]+=1;
//                }else
//                {
//                  singularPoint+=1;
//                }
//            }
//        }
//    }
//    for(int i=1;i<(int)(maxVol+1);i++)
//    {
//        if(imageHist[i]!=0)
//        {
//          nonzeroLabel++;
//        }
//    }
    //int Nstruct =nonzeroLabel+2;    

    // ////////////////////////////////////    
    for (int zz=0; zz<fsx; zz++)
    {
       for(int jj=0; jj<fsy; jj++)
       {
          for(int kk=0; kk<fsz; kk++)
          {
              rnI1=I1seg->get(zz,jj,kk);
              rnI0=I0seg->get(zz,jj,kk);
              rndI0=deformI0seg->get(zz,jj,kk);

              // Initial Dice nominator
              if((rnI1==rnI0)&&(rnI1!=0)&&(rnI0!=0))
              {
                  intersectCount[0]++;
              }

              // Final Dice nominator
              if((rnI1==rndI0)&&(rnI1!=0)&&(rndI0!=0))
              {
                  intersectCount[Nstruct-1]++;
                  for(int ii=1;ii<Nstruct-1;ii++)
                  {
                      if(ii==(int)rnI1)
                      {
                          intersectCount[ii]++;
                      }
                  }              
              }

              // I1 share for the denominator
              if(rnI1!=0)
              {
                  numI1[0]++;
                  numI1[Nstruct-1]++;
                  for(int ii=1;ii<Nstruct-1;ii++)
                  {
                      if(ii==(int)rnI1)
                      {
                          numI1[ii]++;
                      }
                  }           

              }

              // I0 share for the denominator
              if(rnI0!=0)
              {
                  numI0[0]++;

              }

              // I0d share for the denominator
              if(rndI0!=0)
              {
                  numI0[Nstruct-1]++;
                  for(int ii=1;ii<Nstruct-1;ii++)
                  {
                      if(ii==(int)rndI0)
                      {
                          numI0[ii]++;
                      }
                  }     
              }
          }
       }
    }

    for(int ii=0; ii<Nstruct; ii++)
    {
        if(intersectCount[ii]|numI0[ii]|numI1[ii])
        {
            diceScores.push_back((2*(double)intersectCount[ii])/((double)(numI0[ii]+numI1[ii])));
            //std::cout<<"\n"<<ii<<" "<<intersectCount[ii]<<" "<<numI0[ii]<<" "<<numI1[ii]<<std::endl;
        }
    }
}





















