# DiffeoRaptor
Diffeomorphic Multimodal Image Registration using RaPTOR

![alt text](https://github.com/nimamasoumi/DiffeoRaptor/blob/main/brain.png?raw=true)

This code is based on Fourier-Approximated Lie Algebras (FLASH). FLASH can be downloaded from https://bitbucket.org/FlashC/flashc/src/master/

The instructions to install DiffeoRaptor is similar to the one in the above link. You can find the same steps below:

# Installation for Linux users #

- Step 1:
  Install Armadillo
  
  The version of Armadillo >= 3.920
  http://arma.sourceforge.net/download.html

- Step 2:
  Install Pyca
  Pyca can be downloaded from https://bitbucket.org/scicompanat/pyca. 
  Follow the steps closely. Pay attention to software requirements and their version.
  You need specific version of your c++ compiler and other tools to successfully install Pyca.
  You can switch off the GPU functionality in Pyca since this is not required in our 
  code. 


- Step 2:
  After unzipping the DiffeoRaptor.zip, in your terminal:
  
  1) cd DiffeoRaptor
  
  2) ccmake .
  
     Specify Pyca source code, and Pyca binary path

  3) make 

# Run the program #

The DiffeoRaptor software package contains several functions for image processing as well as example script files to run the software.

1) In the folder DiffeoRaptor/Testing/runO3test you can find two example scripts. For running DiffeoRaptor, use "diffeoraptorO3_IXI.sh". For evaluating the Dice score use "processResults.sh". You need to modify both scripts based on your current directories and preffered parameter settings.

2) The scripts run the cxx files in the parent directory (DiffeoRaptor/Testing). In order, to modify the functionality, you will need to modify those cxx files.

Please cite our IJCARS paper and the FLASH paper in case you used the code in your project:

[1] N. Masoumi, H. Rivaz, M.O. Ahmad, and Y. Xiao, DiffeoRaptor: diffeomorphic inter-modal image registration using RaPTOR. Int J CARS (2022). https://doi.org/10.1007/s11548-022-02749-2

[2] M. Zhang and P. T. Fletcher, “Fast Diffeomorphic Image Registration via Fourier-Approximated Lie Algebras”, In International Journal of Computer Vision 2019, 127, 61–73. https://doi.org/10.1007/s11263-018-1099-x

[3] M. Zhang and P. T. Fletcher, “Finite-Dimensional Lie Algebras for Fast Diffeomorphic Image Registration”, In Information Processing in Medical Imaging. IPMI, 2015. Lecture Notes in Computer Science, vol 9123. Springer, Cham. https://doi.org/10.1007/978-3-319-19992-4_19

In case you have a question, don't hesitate to post in "Issuses" or contact the first author by:
n_masoum@encs.concordia.ca

