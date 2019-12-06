// Implements the 3rd Primitive
// contact: Alp Yurtsever - alp.yurtsever@epfl.ch
// Last modified: 05 December 2019

#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize n, m, k, j;
    double *MprV, *Vpr, *Cpr;
    int *MprN, *MprI, *MprJ;
    
    n = mxGetM(prhs[0]);
    k = mxGetM(prhs[5]);
    m = *mxGetPr(prhs[4]);
    
    plhs[0] = mxCreateDoubleMatrix( m, 1, mxREAL );
    
    MprN = (int*)mxGetData(prhs[0]);
    MprI = (int*)mxGetData(prhs[1]);
    MprJ = (int*)mxGetData(prhs[2]);
    MprV = mxGetPr(prhs[3]);
    Vpr = mxGetPr(prhs[5]);
    Cpr = mxGetPr(plhs[0]);
    
    for( j=0; j<n; j++ ) {
        Cpr[*MprN++ - 1] += *MprV++ * Vpr[*MprI++ - 1] * Vpr[*MprJ++ - 1];
    }
    
}
