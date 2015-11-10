/* 
    Qt = getQt(Kmt,M,pt);

 */
#include "mex.h"

/* The computational routine */
void getQt(double *Qt, mwSize K, mwSize M, const double *pt)
{
    /* index */
    mwSize i,j,s;    
    double hKij, hKijs;
    double tmp;
    
    /* set Q */
    for (j = 0; j < M && j <= K; j++){
        /* hKij = hyge(K,j,0,0) = 1; */
        hKij = 1.0;
        for (i = j; i >= 0; i--){
            /* hKijs = hyge(K,i,j-i,0) */
            hKijs = hKij;
            tmp = 0.0;
            for (s = j-i; s <= j; s++) {
                tmp += pt[s] * hKijs;
                /* hKijs = hyge(K,i,s+1,i+s+1-j) = hyge(K,i,s,i+s-j) * (j-s) * (s+1) / (K-s) / (i+s-j+1) */
                if (s < j) {
                    hKijs *= (double)(j-s) * (double)(s+1) / (double)(K-s) / (double)(i+s-j+1);
                }
            }
            Qt[j*(K+1)+i] = tmp;
            /* hKij = hyge(K,i-1,j-i+1,0) = hyge(K,i,j-i,0) * (K - i + 1) / (K - j + i) */
            if (i > 0) {
                hKij *= (double)(K - i + 1) / (double)(K - j + i);
            }
        }
    }
    
    for (j = M; j <= K; j++){
        /* hKij = hyge(K,j,0,0) = 1; */
        hKij = 1.0;
        for (i = j; i >= (j-M); i--){
            /* hKijs = hyge(K,i,j-i,0) */
            hKijs = hKij;
            tmp = 0.0;
            for (s = j-i; s <= M; s++) {
                tmp += pt[s] * hKijs;
                /* hKijs = hyge(K,i,s+1,i+s+1-j) = hyge(K,i,s,i+s-j) * (j-s) * (s+1) / (K-s) / (i+s-j+1) */
                if (s < M) {
                    hKijs *= (double)(j-s) * (double)(s+1) / (double)(K-s) / (double)(i+s-j+1);
                }
            }
            Qt[j*(K+1)+i] = tmp;
            /* hKij = hyge(K,i-1,j-i+1,0) = hyge(K,i,j-i,0) * (K - i + 1) / (K - j + i) */
            if (i > (j - M)) {
                hKij *= (double)(K - i + 1) / (double)(K - j + i);
            }
        }
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSize Kmt;
    mwSize M;
    double *pt;                 /* 1xN matrix */
    double *Qt;              /* output matrix */

    /* check for proper number of arguments */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Three inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    /* make sure the first input argument is scalar */
    if( mxGetNumberOfElements(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Input Kmt must be a scalar.");
    }
    
    /* make sure the second input argument is scalar */
    if( mxGetNumberOfElements(prhs[1])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input M must be a scalar.");
    }
    
    /* check that number of rows in third input argument is 1 */
    if(mxGetM(prhs[2])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input pt must be a row vector.");
    }
    
    /* get the value of Kmt and M  */
    Kmt = mxGetScalar(prhs[0]);
    M = mxGetScalar(prhs[1]);

    /* create a pointer to the real data in the input matrix  */
    pt = mxGetPr(prhs[2]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize)(Kmt+1),(mwSize)(Kmt+1),mxREAL);

    /* get a pointer to the real data in the output matrix */
    Qt = mxGetPr(plhs[0]);

    /* call the computational routine */
    getQt(Qt,Kmt,M,pt);
}
