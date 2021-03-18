// Computes the sparse matrix - vector multiplication 
// contact: Alp Yurtsever - alp.yurtsever@epfl.ch     
// Last modified: 26 November 2019                         

#include "mex.h"
#include <time.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize m, n, d, nnz, i, j, indi, indj;
    double *Mindi, *Mindj, *Mvals, *Xvals, *Ovals;
    
    m = mxGetScalar(prhs[3]);
    n = mxGetM(prhs[4]);
    d = mxGetN(prhs[4]);
    nnz = mxGetM(prhs[0]);

    plhs[0] = mxCreateDoubleMatrix( m, d, mxREAL );

    Ovals = mxGetPr(plhs[0]);
    Mindi = mxGetPr(prhs[0]);
    Mindj = mxGetPr(prhs[1]);
    Mvals = mxGetPr(prhs[2]);
    Xvals = mxGetPr(prhs[4]);
            
    for( i=0; i<d; i++){
        for( j=0; j<nnz; j++ ) {
            indi = i*m + Mindi[j]-1;
            indj = i*m + Mindj[j]-1;
            if( indi >= (i+1)*m ) {
                mexErrMsgTxt("Index out of bound (i > m)");
            }
            if( indj >= (i+1)*m ) {
                mexErrMsgTxt("Index out of bound (j > n)");
            }
            Ovals[indi] += Mvals[j] * Xvals[indj];
        }
    }
}
