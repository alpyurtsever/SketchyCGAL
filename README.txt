---------------------------------------------------------------------------
                              SketchyCGAL
                   Scalable Semidefinite Programming
        A.Yurtsever - J.A. Tropp - O. Fercoq - M. Udell - V. Cevher
---------------------------------------------------------------------------

This toolbox includes the source code for the numerical experiments in 
[YTFUC2019].

Please follow <https://github.com/alpyurtsever/SketchyCGAL> for the 
updates. 

Please contact "alp.yurtsever@epfl.ch" or "alpy@mit.edu" for your questions 
and feedbacks.

---------------------------------------------------------------------------
---------------------------------------------------------------------------

INDEX

"./@memLog/memLog.m"
- implements a class handle to externally monitor memory usage of MATLAB.
  Works only on UNIX systems.

"./FilesMaxcut/data/G/"
"./FilesMaxcut/data/DIMACS10/"
- Empty folders. You should download the datasets from GSET and DIMACS10 
  benchmark to run MaxCut experiments.

"./FilesMaxcut/data/WALDSPURGER/C1.mat"
"./FilesMaxcut/data/WALDSPURGER/C2.mat"
"./FilesMaxcut/data/WALDSPURGER/C3.mat"
"./FilesMaxcut/data/WALDSPURGER/C4.mat"
"./FilesMaxcut/data/WALDSPURGER/C5.mat"
"./FilesMaxcut/data/WALDSPURGER/C6.mat"
"./FilesMaxcut/data/WALDSPURGER/C7.mat"
"./FilesMaxcut/data/WALDSPURGER/C8.mat"
"./FilesMaxcut/data/WALDSPURGER/C9.mat"
"./FilesMaxcut/data/WALDSPURGER/C10.mat"
- Datasets constructed by the script provided by Irene Waldspurger. These 
  datasets are the realizations of problems with spurious solutions for 
  Burer-Monteiro factorization methods. See [WW2019] for details.

"./FilesPhaseRetrieval/CreateData.m"
- This script generates data for the abstract phase retrieval experiments.

"./FilesPtychography/data/Cell160.tiff"
"./FilesPtychography/data/Cell320.tiff"
"./FilesPtychography/data/Cell640.tiff"
"./FilesPtychography/data/PIXNIO-39417-2650x1770"
"./FilesPtychography/data/PIXNIO-39417-2650x1770.txt"
- Cell images used in the Ptychography experiments. These images are 
  cropped from the image downloaded from "https://pixnio.com/science/
  microscopy-images/hemorrhagic-fever-marburg-virus/cytoarchitectural-
  histopathologic-changes-detected-in-a-kidney-sample-taken-from-a-marburg-
  patient" (Last accessed on December 05, 2019). We thank Dr. J. Lyle
  Conrad, USCDCP, for sharing this image with CC0 license.

"./FilesPtychography/data/CreateData"
- This script generates data for the synthetic transition matrices for the
  ptychography experiments. These matrices are different from the ones we 
  used in the actual experiments. Therefore the results can be different.

"./FilesPtychography/FilesQAP/munkres/munkres.m"
"./FilesPtychography/FilesQAP/munkres/license.txt"
- Munkres (Hungarian) algorithm used in rounding of QAP SDP solution. 
  We use the implementation by Yi Ciao. Reference:
  "Munkres' Assignment Algorithm, Modified for Rectangular Matrices"
  http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html

"./FilesPtychography/FilesQAP/PartialTrace/TrX.m"
"./FilesPtychography/FilesQAP/PartialTrace/license.txt"
- This function computes the partial trace from a state vector of density 
  matrix. We use the implementation by Tony Cubitt.

"./FilesPtychography/FilesQAP/qapread.m"
"./FilesPtychography/FilesQAP/tspread.m"
- Functions to load data from QAPLIB and TSPLIB.

"./FilesPtychography/FilesQAP/qapdata/"
"./FilesPtychography/FilesQAP/tspdata/"
- Empty folders. You need to download datasets from QAPLIB and TSPLIB to 
  These folders to run QAP experiments. 

"./FilesPtychography/solver/CGAL.m" 
"./FilesPtychography/solver/@NystromSketch/NystromSketch.m" 
"./FilesPtychography/solver/@NystromSketch/COPYRIGHT.txt" 
- Implementation of CGAL, ThinCGAL, and SketchyCGAL. NystromSketch class 
  is a minimal version of our SKETCH toolbox. You can find more details on  
  SKETCH toolbox at <https://github.com/alpyurtsever/SKETCH>.

"./FilesPtychography/utils/mexPrimitive2.c" 
"./FilesPtychography/utils/mexPrimitive3.c" 
"./FilesPtychography/utils/mexSparseMult.c"
- Generating sparse matrices is a time consuming process in MATLAB. To 
  avoid regenerating sparse matrices in the QAP experiment, we implement 
  these three MEX files, which directly computes the sparse matrix - dense 
  vector products from indexes. Primitive2 and Primitive3 are used to 
  implement the black-box oracles. mexSparseMult is used in nonnegativity 
  constraints.

"./plotter/Figure71_MaxCut.m"
"./plotter/Figure72_MaxCut.m"
"./plotter/Figure73_PhaseRetrieval.m"
"./plotter/Figure75_QAP.m"
- These scripts generate the Figures in the paper. Before running these 
  scripts, you should run the associated experiments. These scripts only 
  plot the figures from the saved results of these experiments. 

"./Test_MaxCut_CGAL.m" 
"./Test_MaxCut_Mosek.m" 
"./Test_MaxCut_SDPT3.m" 
"./Test_MaxCut_SDPNAL.m" 
"./Test_MaxCut_Sedumi.m" 
- These files implements the MaxCut experiments that we used to produce 
  Figure 7.1.

"./Test_MaxCut_CGAL_Figure72.m" 
- This file implements the MaxCut experiment to produce Figure 7.2.

"./Test_MaxCut_Waldspurger.m" 
- This file implements the numerical demonstration of spurious solutions
  for Burer-Monteiro factorization methods, present in the supplements. 
  See Figure SM6.2.

"./Test_PhaseRetrieval.m" 
- This file implements the abstract Phase Retrieval experiments to produce 
  Figure 7.3. You shuold run this file for all problem dimensions n, with 
  all methods, and from MC = 1 to 20, for 20 Monte-Carlo trials. 

"./Test_Ptychography.m" 
- This file implements the Ptchograph imaging experiments. You can run it 
  for different problem sizes n = {160,320,640}. 

"./Test_QAP.m" 
- This file implements the QAP experiments.

"./COPYRIGHT.txt" 
"./LICENSE.txt" 
- License information. 

"./README.txt" 
- This README file. 

---------------------------------------------------------------------------
---------------------------------------------------------------------------

NOTES

This code is tested on MATLAB R2018a. It might not work on other versions.

"Test_MaxCut_SDPNAL" requires SDPNAL+ (Tested with version 1.0).
Visit "https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/" to download 
SDPNAL+.

"Test_MaxCut_SDPT3" requires SDPT3. (Tested with version 4.0).
Visit "https://blog.nus.edu.sg/mattohkc/softwares/sdpt3/" to download 
SDPT3.

"Test_MaxCut_Sedumi" requires Sedumi. (Tested with version 1.3)
Visit "http://sedumi.ie.lehigh.edu/" to download Sedumi. 

"Test_MaxCut_manopt" and "Test_MaxCut_Waldspurger" require manopt. 
(Tested with version 5) 
Visit "https://www.manopt.org/" to download manopt.

"Test_MaxCut_Mosek" requires MoSeK. MoSeK requires a license, personal 
academic license for MoSeK is free. (Tested with version 8)
Visit "http://www.mosek.com/" to download MOSEK.

---------------------------------------------------------------------------
---------------------------------------------------------------------------

CITATION

If you find this toolbox useful in your research, please cite our work:

[YTFUC2019] A.Yurtsever, J.A. Tropp, O. Fercoq, M. Udell, V. Cevher
"Scalable Semidefinite Programming", 2019.

---------------------------------------------------------------------------
---------------------------------------------------------------------------

Last edit: Alp Yurtsever - December 06, 2019