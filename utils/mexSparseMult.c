// Computes the sparse matrix - vector multiplication 
// contact: Alp Yurtsever - alp.yurtsever@epfl.ch     
// Last modified: 26 November 2019                         

#include "mex.h"
#include <time.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize m, n, nnz, j, indi, indj;
    double *Mindi, *Mindj, *Mvals, *Xvals, *Ovals;
    
    m = mxGetScalar(prhs[3]);
    n = mxGetM(prhs[4]);
    nnz = mxGetM(prhs[0]);

    plhs[0] = mxCreateDoubleMatrix( m, 1, mxREAL );

    Ovals = mxGetPr(plhs[0]);
    Mindi = mxGetPr(prhs[0]);
    Mindj = mxGetPr(prhs[1]);
    Mvals = mxGetPr(prhs[2]);
    Xvals = mxGetPr(prhs[4]);

    for( j=0; j<nnz; j++ ) {
        indi = Mindi[j]-1;
        indj = Mindj[j]-1;
        if( indi >= m ) {
            mexErrMsgTxt("Index out of bound (i > m)");
        }
        if( indj >= m ) {
            mexErrMsgTxt("Index out of bound (j > n)");
        }
        Ovals[indi] += Mvals[j] * Xvals[indj];
    }
}
