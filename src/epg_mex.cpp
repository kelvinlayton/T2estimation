#include "mex.h"
#include <math.h>
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    // Process input arguments
	//
	int nEchoes = (int)mxGetScalar(prhs[0]);
	double tau = (double)mxGetScalar(prhs[1]);
	double R1 = (double)mxGetScalar(prhs[2]);
	double *R2 = (double*)mxGetPr(prhs[3]);
	double alpha = (double)mxGetScalar(prhs[4]);	
	int nRates = (int)mxGetM(prhs[3]);

	// Define RF mixing matrix
	//
	double *T = (double*)mxCalloc(9,sizeof(double)); 
	double cosa = cos(alpha);
	double sina = sin(alpha);

	T[0] = 0.5*(1+cosa); T[1] = 0.5*(1-cosa); T[2] = sina;
	T[3] = T[1];         T[4] = T[0];         T[5] = -sina;
	T[6] = -0.5*sina;    T[7] = 0.5*sina;     T[8] = cosa;
	
	// T1 relaxation
	//
	double E1 = exp(-tau/2*R1);

	// Allocate state arrays
	//
	int nFStates = 2*nEchoes+1;
	double *Z = (double*)mxCalloc(nEchoes,sizeof(double));  // Z1 Z2 Z3 ...
	double *F = (double*)mxCalloc(nFStates,sizeof(double));	// F0  F1 F-1 F2 F-2 ...

	// Allocate output array
	//
	plhs[0] = mxCreateDoubleMatrix(nEchoes, nRates, mxREAL);
	double * output = mxGetPr(plhs[0]);

	// Iterate over different relaxation rates
	//
	for (int iRate=0; iRate<nRates; iRate++) 
	{
		// T2 relaxation
		//
		double E2 = exp(-tau/2*R2[iRate]);
		
		// Reset state arrays
		//
		memset(F,0,nFStates*sizeof(double));
		memset(Z,0,nEchoes*sizeof(double));
		F[nEchoes] = 1;

		// Iterate over echoes
		//
		for (int iEcho=0; iEcho<nEchoes; iEcho++) {

			// Precession and relaxation
			for(int i=0; i<2*nEchoes; i++) {
				F[i]=F[i+1]*E2;
			}
			for(int i=0; i<nEchoes; i++) {
				Z[i]=Z[i]*E1;
			}

			// RF mixing
			for (int i=1; i<nEchoes; i++) {
				double Fpos = F[nEchoes-i];
				double Fneg = F[nEchoes+i];
				double Zpos = Z[nEchoes-i];

				F[nEchoes-i] = T[0]*Fpos + T[1]*Fneg + T[2]*Zpos;
				F[nEchoes+i] = T[3]*Fpos + T[4]*Fneg + T[5]*Zpos;
				Z[nEchoes-i] = T[6]*Fpos + T[7]*Fneg + T[8]*Zpos;
			}

			// Precession and relaxation
			for(int i=0; i<2*nEchoes; i++) {
				F[i]=F[i+1]*E2;
			}
			for(int i=0; i<nEchoes; i++) {
				Z[i]=Z[i]*E1;
			}

			// Save echo amplitudes
			//
			output[iEcho + iRate*nEchoes] = F[nEchoes];
		}
	}

	// Free memory
	//
	mxFree(Z);
	mxFree(F);
	mxFree(T);

} 
