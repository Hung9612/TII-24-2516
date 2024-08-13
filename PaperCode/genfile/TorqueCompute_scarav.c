/*
% TorqueCompute_scarav.c 
*/
#include "mex.h"
#include <math.h>
#include "TorqueCompute_scarav_lib.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
	double *Q;
	double *Qd;
	double *Qdd;
	double *Theta;
	double *tau;
	int N;


	if(nrhs!=4) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Wrong number of inputs.");
    }
    if(nlhs> 1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Too many outputs!");
    }
    int k;
    for (k=0;k<4;k++){
	    if( (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ) ||
	        (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ) ||
	        (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ) ||
	        (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) )) {
	        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrices must be type double.");
	    }
	}
    if(mxGetM(prhs[0])!=4 || mxGetM(prhs[1])!=4 || mxGetM(prhs[2])!=4) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Wrong number of rows in input matrices.");
    }
    if(mxGetM(prhs[3])!=8 || mxGetN(prhs[3])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Theta is the wrong size.");
    }

	Q = mxGetDoubles(prhs[0]);
	Qd = mxGetDoubles(prhs[1]);
	Qdd = mxGetDoubles(prhs[2]);
	Theta = mxGetDoubles(prhs[3]);
	N = mxGetN(prhs[0]);

 	plhs[0] = mxCreateDoubleMatrix(4, N, mxREAL);
	tau = mxGetPr(plhs[0]);

	TorqueCompute_scarav_lib(Q, Qd, Qdd, Theta, N, tau);
}