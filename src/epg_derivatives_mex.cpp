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

	// Define RF mixing matrix and it's derivative with respect to alpha
	//
	double *T = (double*)mxCalloc(9,sizeof(double)); 
	double *Td = (double*)mxCalloc(9,sizeof(double)); 

	double cosa = cos(alpha);
	double sina = sin(alpha);
	
	T[0] = 0.5*(1+cosa); T[1] = 0.5*(1-cosa); T[2] = sina;
	T[3] = T[1];         T[4] = T[0];         T[5] = -sina;
	T[6] = -0.5*sina;    T[7] = 0.5*sina;     T[8] = cosa;
	
	Td[0] = -0.5*sina;    Td[1] = 0.5*sina;     Td[2] = cosa;
	Td[3] = Td[1];        Td[4] = Td[0];        Td[5] = -cosa;
	Td[6] = -0.5*cosa;    Td[7] = 0.5*cosa;     Td[8] = -sina;

	// T1 relaxation and derivative with respect to R2
	//
	double E1 = exp(-tau/2*R1);
	double E1d = 0;

	// Allocate state arrays including intermediate states
	//
	int nFStates = 2*nEchoes+1;
	double *F = (double*)mxCalloc(nFStates,sizeof(double));	// F0  F1 F-1 F2 F-2 ...
	double *Z = (double*)mxCalloc(nEchoes,sizeof(double)); // Z1 Z2 Z3 ...
	
	double *FdR = (double*)mxCalloc(nFStates,sizeof(double));
	double *ZdR = (double*)mxCalloc(nEchoes,sizeof(double));
	double *FdRa = (double*)mxCalloc(nFStates,sizeof(double));
	double *ZdRa = (double*)mxCalloc(nEchoes,sizeof(double));

	double *FdA = (double*)mxCalloc(nFStates,sizeof(double));
	double *ZdA = (double*)mxCalloc(nEchoes,sizeof(double));
	double *FdAa = (double*)mxCalloc(nFStates,sizeof(double));
	double *ZdAa = (double*)mxCalloc(nEchoes,sizeof(double));

	// Allocate output arrays
	//
	plhs[0] = mxCreateDoubleMatrix(nEchoes, nRates, mxREAL);
	double * derivative_wrt_R = mxGetPr(plhs[0]);
	double * derivative_wrt_A = 0; 
	if (nlhs>1) {
		plhs[1] = mxCreateDoubleMatrix(nEchoes, nRates, mxREAL);
		derivative_wrt_A = mxGetPr(plhs[1]);
	}
	
	// Iterate over different relaxation rates
	//
	for (int iRate=0; iRate<nRates; iRate++) 
	{
		// T2 relaxation and derivative
		//
		double E2 = exp(-tau/2*R2[iRate]);
		double E2d = -tau/2*E2;

		// Reset all state arrays
		//
		memset(F,0,nFStates*sizeof(double));
		memset(Z,0,nEchoes*sizeof(double));
		memset(FdR,0,nFStates*sizeof(double));
		memset(ZdR,0,nEchoes*sizeof(double));
		memset(FdRa,0,nFStates*sizeof(double));
		memset(ZdRa,0,nEchoes*sizeof(double));
		memset(FdA,0,nFStates*sizeof(double));
		memset(ZdA,0,nEchoes*sizeof(double));
		memset(FdAa,0,nFStates*sizeof(double));
		memset(ZdAa,0,nEchoes*sizeof(double));

		F[nEchoes] = 1;
		FdRa[nEchoes] = 1;
		FdAa[nEchoes] = 1;

		// Iterate over echoes
		//
		for (int iEcho=0; iEcho<nEchoes; iEcho++) {

			// Precession and relaxation of all intermediate states
			//
			for(int i=0; i<2*nEchoes; i++) {
				F[i]=F[i+1]*E2;
				FdRa[i]=FdRa[i+1]*E2d;
				FdR[i]=FdR[i+1]*E2;
				FdAa[i]=FdAa[i+1]*E2;
				FdA[i]=FdA[i+1]*E2;
			}
			for(int i=0; i<nEchoes; i++) {
				Z[i]=Z[i]*E1;
				ZdRa[i]=ZdRa[i]*E1d;
				ZdR[i]=ZdR[i]*E1;
				ZdAa[i]=ZdAa[i]*E1;
				ZdA[i]=ZdA[i]*E1;
			}

			// RF mixing 
			//
			for (int i=1; i<nEchoes; i++) {
				double Fpos = F[nEchoes-i];
				double Fneg = F[nEchoes+i];
				double Zpos = Z[nEchoes-i];

				// Signal path
				F[nEchoes-i] = T[0]*Fpos + T[1]*Fneg + T[2]*Zpos;
				F[nEchoes+i] = T[3]*Fpos + T[4]*Fneg + T[5]*Zpos;
				Z[nEchoes-i] = T[6]*Fpos + T[7]*Fneg + T[8]*Zpos;

				// First term w.r.t R
				Fpos = FdRa[nEchoes-i];
				Fneg = FdRa[nEchoes+i];
				Zpos = ZdRa[nEchoes-i];
				FdRa[nEchoes-i] = T[0]*Fpos + T[1]*Fneg + T[2]*Zpos;
				FdRa[nEchoes+i] = T[3]*Fpos + T[4]*Fneg + T[5]*Zpos;
				ZdRa[nEchoes-i] = T[6]*Fpos + T[7]*Fneg + T[8]*Zpos;

				// Second term w.r.t R
				Fpos = FdR[nEchoes-i];
				Fneg = FdR[nEchoes+i];
				Zpos = ZdR[nEchoes-i];
				FdR[nEchoes-i] = T[0]*Fpos + T[1]*Fneg + T[2]*Zpos;
				FdR[nEchoes+i] = T[3]*Fpos + T[4]*Fneg + T[5]*Zpos;
				ZdR[nEchoes-i] = T[6]*Fpos + T[7]*Fneg + T[8]*Zpos;

				// First term w.r.t A
				Fpos = FdAa[nEchoes-i];
				Fneg = FdAa[nEchoes+i];
				Zpos = ZdAa[nEchoes-i];
				FdAa[nEchoes-i] = Td[0]*Fpos + Td[1]*Fneg + Td[2]*Zpos;
				FdAa[nEchoes+i] = Td[3]*Fpos + Td[4]*Fneg + Td[5]*Zpos;
				ZdAa[nEchoes-i] = Td[6]*Fpos + Td[7]*Fneg + Td[8]*Zpos;

				// Second term w.r.t A
				Fpos = FdA[nEchoes-i];
				Fneg = FdA[nEchoes+i];
				Zpos = ZdA[nEchoes-i];
				FdA[nEchoes-i] = T[0]*Fpos + T[1]*Fneg + T[2]*Zpos;
				FdA[nEchoes+i] = T[3]*Fpos + T[4]*Fneg + T[5]*Zpos;
				ZdA[nEchoes-i] = T[6]*Fpos + T[7]*Fneg + T[8]*Zpos;
			}

			// Precession and relaxation
			//
			for(int i=0; i<2*nEchoes; i++) {
				F[i]=F[i+1]*E2;
				FdRa[i]=FdRa[i+1]*E2;
				FdR[i]=FdR[i+1]*E2;
				FdAa[i]=FdAa[i+1]*E2;
				FdA[i]=FdA[i+1]*E2;
			}
			for(int i=0; i<nEchoes; i++) {
				Z[i]=Z[i]*E1;
				ZdRa[i]=ZdRa[i]*E1;
				ZdR[i]=ZdR[i]*E1;
				ZdAa[i]=ZdAa[i]*E1;
				ZdA[i]=ZdA[i]*E1;
			}

			// Combine intermediate states to get derivatives
			//
			for(int i=0; i<2*nEchoes; i++) {
				FdR[i] = 2*FdRa[i] + FdR[i];
				FdA[i] = FdAa[i] + FdA[i];
			}
			for(int i=0; i<nEchoes; i++) {
				ZdR[i] = 2*ZdRa[i] + ZdR[i];
				ZdA[i] = ZdAa[i] + ZdA[i];
			}

			// Copy current state vector to intermediate variables for next iteration
			//
			memcpy(FdRa,F,nFStates*sizeof(double));
			memcpy(ZdRa,Z,nEchoes*sizeof(double));
			memcpy(FdAa,F,nFStates*sizeof(double));
			memcpy(ZdAa,Z,nEchoes*sizeof(double));

			// Save derivatives of echo amplitudes with respect to R2 and alpha
			//
			derivative_wrt_R[iEcho + iRate*nEchoes] = FdR[nEchoes];
			if (derivative_wrt_A) {
				derivative_wrt_A[iEcho + iRate*nEchoes] = FdA[nEchoes];
			}
		}	// loop echoes
	}	// loop rates

	// Free memory
	//
	mxFree(Z);
	mxFree(F);
	mxFree(ZdR);
	mxFree(FdR);
	mxFree(ZdRa);
	mxFree(FdRa);
	mxFree(ZdA);
	mxFree(FdA);
	mxFree(ZdAa);
	mxFree(FdAa);
	mxFree(T);
	mxFree(Td);
} 
