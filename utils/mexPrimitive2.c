// Implements the 2nd Primitive
// contact: Alp Yurtsever - alp.yurtsever@epfl.ch
// Last modified: 05 December 2019

#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize m, n, j;
    double *MprV, *Rpr, *Vpr, *Cpr;
    int *MprN, *MprI, *MprJ;
    
    n = mxGetM(prhs[0]);
    m = mxGetM(prhs[5]);
    
    plhs[0] = mxCreateDoubleMatrix( m, 1, mxREAL );
    
    MprN = (int*)mxGetData(prhs[0]);
    MprI = (int*)mxGetData(prhs[1]);
    MprJ = (int*)mxGetData(prhs[2]);
    MprV = mxGetPr(prhs[3]);
    Rpr = mxGetPr(prhs[4]);
    Vpr = mxGetPr(prhs[5]);
    Cpr = mxGetPr(plhs[0]);
    
    for( j=0; j<n; j++ ) {
        Cpr[*MprJ - 1] += 0.5 * (Rpr[*MprN - 1] * *MprV * Vpr[*MprI - 1]);  // Accumulate contribution of Vpr[j]
        Cpr[*MprI++ - 1] += 0.5 * (Rpr[*MprN++ - 1] * *MprV++ * Vpr[*MprJ++ - 1]);  // Accumulate contribution of Vpr[j]
    }
    
}
