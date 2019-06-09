#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <string.h>
//#include <cmath>
//#define _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <math.h>
#include "mex.h"
#include "engine.h"
#include <omp.h>
//static const int lat_len=100;//lattice length
//static const int isRing=1,isMagnet=1,fluxDiv=1;//isRing indicates whether we are working on a ring,isMagnet
//indicates whether we have a magnetic flux, fluxDiv is the number of flux point we want to find
//static const double lambda=1.0,t_amp=1.0;//potential amplitude,frequency,hopping amplitude
#define lat_len 50//lattice length
#define isRing 0//1 if the model is a ring, 0 otherwise
#define isMagnet 0//1 if there is a magnetic field, 0 otherwise
#define fluxDiv 10000// number of flux resolution
#define lambda 0.0//potential amplitude
#define t_amp -1.0//hopping amplitude
#define BSZ 1024//buffer length
#define isHartr 0//1 if there is a Hartree term,0 otherwise
#define isFock 0//1 if there is a Fock term,0 otherwise
//#define interest_len 99// the length that's interesting
#define isHubbard 0//1 if there is a Hubbard term,0 otherwise
#define converg_lim 1e-12//the limit in which we stop looking for energies
#define U_Hubb 1.0
#define U_HF 1.0
#define phi 0//2.199114857//cosine phase
#define isFlip 1
#define t_flip -0.1
#define p_alpha 0
#define q_alpha 5
#define f1 0.5//first row flux coefficient
#define f2 -0.5//second row flux coefficient
#define isNext 0
#define t_next 0.50
#define isSOC 0
#define isConstSOC 1
#define SOCDiv 200
#define SOCMax 2.2
//#define SOC_amp 1.0
#define FluxCase 5//8
#define IterLim 3000
#define isCurCalc 1
#define isCurPrint 1
#define isCurSumCalc 0
#define EE_part 0.5
#define isChargeDensityCalc 1
#define isEEcalc 1
#define isSymmetric 1.0//relevant only for impurity(FluxCase==5). defines whether the impurity is symmetric or not
#define isCondCalc 1//calculate conductivity according to Kubo formula
#define min_fill 2//the lowest value of filling. can't be less than 1
#define max_fill 20//the maximal value of filling, cant be more than 2*lat_len
#define ImpurityNum 50
#define max_impurity 2.0
#define Full2Pi 0//#define isMixed 0//take mixed states as trial wave functions

//mxArray *mexstateR = NULL, *mexstateI = NULL, *mexE = NULL, *mexSumState = NULL, *mexSumStateSqUp = NULL, *mexSumStateSqDown = NULL;
//mxArray *mexstate = NULL, *mexE = NULL;
//double *E, *SumStateSqUp, *SumStateSqDown, *SumStateSqUp0, *SumStateSqDown0, *stateR, *stateI, *SumStateR, *stateR0, *stateI0, *E0;
//double *E=NULL, *state=NULL, *E0=NULL,*state_origR=NULL,*state_origI=NULL;
//double *E, *state, *E0, *state_origR, *state_origI;
double fluxVal, alpha, SOC_amp,ImpureAmp;
int isMixed = 0;//take mixed states as trial wave functions, I do it as a global variable because i want to allow changes in the program 
char FileStrE[20] = "Energy.txt";//_impurity01_symetry
char FileStrCurLon[20] = "Cur_lon.txt";
char FileStrCurTrn[20] = "Cur_trn.txt";
char FileStrCurSum[20] = "Cur_sum.txt";
char FileStrCur[20] = "Cur_flip.txt";
char FileStrState[25] = "EE_half_system.txt"; 
char PathStr[150] = "C:\\Users\\Arik\\Desktop\\dummy\\";
//complex double *state, *SumState;

int sign(int x) {
	return (x > 0) - (x	 < 0);
}

int Factorial(int x) {
	int i, n = 1;
	for (i = 2; i <= x; i++) {
		n *= i;
	}
	return n; 
}

double ComplexAbs(double x, double y) {//x is the real part and y is the imaginary
	return sqrt(pow(x, 2) + pow(y, 2));
}

double ImagMul(double x1, double y1, double x2, double y2) {
	return x2 * y1 + y2 * x1;
}

double RealMul(double x1, double y1, double x2, double y2) {
	return x2 * x1 - y1 * y2;
}

double MatMulReal(double* xR, double* yR, double* xI, double* yI, int ind1, int ind2, int size1, char MulCase) {
	int i, j, temp1 = 0.0, temp2 = 0.0;
	double mul = 0.0;
	switch (MulCase) {
	case 'N':
		for (i = 0; i < size1; i++) {
			temp1 = i * size1 + ind1;
			temp2 = ind2 * size1 + i;
			mul += RealMul(xR[temp1], xI[temp1], yR[temp2], yI[temp2]);
		}
		break;
	case 'T':
		for (i = 0; i < size1; i++) {
			temp1 = i * size1 + ind1;
			temp2 = i * size1 + ind1;
			mul += RealMul(xR[temp1], xI[temp1], yR[temp1], yI[temp1]);
		}
		break;
	}

	return mul;
}

int isHermitian(double *HR[], double *HI[],int size) {
	int i,dif;
	if (FluxCase != 8 && FluxCase == 9)
		dif = 2;
	else
		dif = 1;
	for (i = dif; i < size-dif; i=i+2*dif) {
		if (HR[i-dif][i] != HR[i][i-dif] || HI[i-dif][i] != -HI[i][i-dif] || HR[i + dif][i] != HR[i][i + dif] || HI[i + dif][i] != -HI[i][i + dif]) {
			return 0;
		}
	}
	return 1;
}

double MatMulImag(double* xR, double* yR, double* xI, double* yI, int ind1, int ind2, int size1) {
	int i, j, temp1, temp2;
	double mul = 0.0;
	for (i = 0; i < size1; i++) {
		temp1 = i * 2 * lat_len + ind1;
		temp2 = ind2 * 2 * lat_len + i;
		mul += ImagMul(xR[temp1], xI[temp1], yR[temp2], yI[temp2]);
	}
	return mul;
}
//void init_array(double a[][lat_len])
void init_array(double* a[], int size1, int size2)//the function initializes the array given
{
	int i, j;
	for (i = 0; i<size1; i++) {
		for (j = 0; j<size2; j++) {
			a[i][j] = 0.0;
		}
	}
}

void set_array(double* ToSet[], double* Original[], int size1, int size2)//the function sets the 1st array given to the values of the 2nd array
{
	int i, j;
	for (i = 0; i<size1; i++) {
		for (j = 0; j<size2; j++) {
			ToSet[i][j] = Original[i][j];
		}
	}
}

double factorial(double f)
{
	if (f == 0)
		return 1;
	return(f * factorial(f - 1));
}

void EE_calc(double* stateI, double* stateR, int interest_len, Engine *ep, int dif) {
	int i, j, ind1, ind2, fill;
	int MatSize = (int) 2.0*EE_part*lat_len;
	double **corrMatR = (double**)malloc(MatSize * sizeof(double*));
	double **corrMatI = (double**)malloc(MatSize * sizeof(double*));
	for (i = 0; i < MatSize; i++) {
		corrMatR[i] = (double*)malloc(MatSize * sizeof(double));
		corrMatI[i] = (double*)malloc(MatSize * sizeof(double));
	}
	init_array(corrMatR, MatSize, MatSize);
	init_array(corrMatI, MatSize, MatSize);

	for (i = 0; i < MatSize; i++) {
		for (j = 0; j < MatSize; j++) {
			for (fill = 1; fill <= interest_len; fill++) {
				ind1 = 2 * lat_len*(fill - 1) + i * dif;
				ind2 = 2 * lat_len*(fill - 1) + j * dif;
				corrMatR[i][j] += RealMul(stateR[ind1], -stateI[ind1], stateR[ind2], stateI[ind2]);//RealMul(stateR[ind1],stateI[ind1],stateR[ind2],stateI[ind2]);
				corrMatI[i][j] += ImagMul(stateR[ind1], -stateI[ind1], stateR[ind2], stateI[ind2]);
			}
		}
	}
	mxArray *mxCorrI = NULL, *mxCorrR = NULL;
	mxCorrI = mxCreateNumericMatrix(MatSize, MatSize, mxDOUBLE_CLASS, mxREAL);
	mxCorrR = mxCreateNumericMatrix(MatSize, MatSize, mxDOUBLE_CLASS, mxREAL);
	for (i = 0; i < MatSize; i++) {
		for (j = 0; j < MatSize; j++) {
			((double*)mxGetPr(mxCorrR))[j + MatSize * i] = corrMatR[j][i];
			((double*)mxGetPr(mxCorrI))[j + MatSize * i] = corrMatI[j][i];
		}
	}
	engPutVariable(ep, "CorrR", mxCorrR);
	engPutVariable(ep, "CorrI", mxCorrI);
	engEvalString(ep, "corr=complex(CorrR,CorrI);");
	engEvalString(ep, "[~,cd]=eig(corr);");
	engEvalString(ep, "cd=diag(cd);");
	engEvalString(ep, "cd=cd((cd>1e-12) & (cd<1-1e-12));");
	if (dif == 1)//Vertical cut off(left part of the system)
		engEvalString(ep, "EEV(MyLen,fluxVal+1)= sum(-cd.*log(cd) - (1 - cd).*log(1 - cd))");
	else if (dif == 2)//Horizontal cut-off(top leg)
		engEvalString(ep, "EEH(MyLen,fluxVal+1)= sum(-cd.*log(cd) - (1 - cd).*log(1 - cd))");
	//engEvalString(ep, "S(MyLen,fluxVal)=EE_func(corr);");
	mxDestroyArray(mxCorrR); mxDestroyArray(mxCorrI);
	free(corrMatR); free(corrMatI);
}

double Hubbard_Exchange(int j, int ind1, int dif, int caseC, double *stateR, double *stateI) {
	double term1I, term1R, term2R, term2I;
	int temp;
	temp = (j << 1)*lat_len + ind1;
	term1I = ImagMul(stateR[temp], -stateI[temp], stateR[temp + dif], stateI[temp + dif]);
	term1R = RealMul(stateR[temp], -stateI[temp], stateR[temp + dif], stateI[temp + dif]);
	/*temp = (k << 1)*lat_len + ind1;
	term2I = ImagMul(term1R, term1I, stateR[temp + dif], -stateI[temp + dif]);
	term2R = RealMul(term1R, term1I, stateR[temp + dif], -stateI[temp + dif]);
	term1I = ImagMul(term2R, term2I, stateR[temp], stateI[temp]);
	term1R = RealMul(term2R, term2I, stateR[temp], stateI[temp]);*/
	if (caseC == 1)//real
		return term1R;
	else if (caseC == 2)
		return term1I;
}

void add_Hubbard_inter_aid(double* HR[], double* HI[], double* stateR, double* stateI, int ind1, int interest_len) {//,int dif) {
	int j, k, temp, dif;
	double norm1, norm2, term1R, term2R, term1I, term2I, NormFactor = 1.0;//Factorial(interest_len); (double)interest_len*(interest_len - 1) / 2.0;
	double HubbSum1R = 0.0, HubbSum2R = 0.0, HubbSum1I = 0.0, HubbSum2I = 0.0, HubbSum3 = 0.0, HubbSum4 = 0.0, HubbSum5 = 0.0, HubbSum6 = 0.0;
	/*for (j = 0; j < interest_len - 1; j++) {
	for (k = j + 1; k < interest_len; k++) {*/
	if (ind1 % 2 == 0)
		dif = 1;
	else
		dif = -1;
	//for (j = 0; j < lat_len; j++) {
	for (k = 0; k < interest_len; k++) {
		//if (j != k) {
		term1R = Hubbard_Exchange(k, ind1, dif, 1, stateR, stateI);
		HubbSum1R += term1R;//real part
		term1I = Hubbard_Exchange(k, ind1, dif, 2, stateR, stateI);
		HubbSum1I += term1I;//imaginary part

		norm1 = pow(stateR[2 * k*lat_len + ind1 + dif], 2) + pow(stateI[2 * k*lat_len + ind1 + dif], 2);
		//norm2 = pow(stateR[2 * k*lat_len + ind1 + 1], 2) + pow(stateI[2 * k*lat_len + ind1 + 1], 2);
		HubbSum3 += norm1;// *norm2;
						  //}
	}
	HR[ind1][ind1] += U_Hubb * HubbSum3 / NormFactor;
	//HR[ind1 + 1][ind1 + 1] += U_Hubb*HubbSum3/NormFactor;
	//HI[ind1][ind1] += U_Hubb*HubbSum5 / NormFactor;
	//HI[ind1 + 1][ind1 + 1] += U_Hubb*HubbSum6 / NormFactor;
	if (isFlip) {
		HR[ind1][ind1 + dif] -= U_Hubb * HubbSum1R / NormFactor;
		HI[ind1][ind1 + dif] -= U_Hubb * HubbSum1I / NormFactor;
	}
}

void add_Hartr_inter_aid(double* HR[], double* HI[], int ind1, int interest_len, double* stateR, double* stateI, int LengthFactor, int dif) {
	double norm1, norm2, NormFactor = 2.0, HartrSum1 = 0.0;
	int j, k, m;
	if (ind1 < LengthFactor - dif) {
		//for (j = 0; j < interest_len; j=j++) {
		for (k = 0; k < interest_len; k = k++) {
			norm1 = pow(ComplexAbs(stateR[LengthFactor * k + ind1 + dif], stateI[LengthFactor * k + ind1 + dif]), 2);
			norm2 = 1.0;// pow(ComplexAbs(stateR[LengthFactor * k + ind1 + dif], stateI[LengthFactor * k + ind1 + dif]), 2);
			HartrSum1 += norm1 * norm2;
		}
	}

	if (ind1 >= dif) {
		//for (j = 0; j < interest_len; j = j++) {
		for (k = 0; k < interest_len; k = k++) {
			norm1 = pow(ComplexAbs(stateR[LengthFactor * k + ind1 - dif], stateI[LengthFactor * k + ind1 - dif]), 2);
			HartrSum1 += norm1;// *norm2;
		}
	}
	HR[ind1][ind1] += U_HF * HartrSum1 / NormFactor;
}

void add_Fock_inter_aid(double* HR[], double* HI[], int ind1, int interest_len, double* stateR, double* stateI, int LengthFactor, int dif) {
	double term1R, term1I, term2R, term2I, FockSumR = 0.0, FockSumI = 0.0, NormFactor = 2.0;
	int i, j, temp;
	//for (i = 0; i < interest_len; i++) {
	//for (i = 0; i < lat_len; i++) {
	for (j = 0; j < interest_len; j++) {
		//			if (i != j) {
		temp = j * LengthFactor + ind1;
		term1I = ImagMul(stateR[temp], -stateI[temp], stateR[temp + dif], stateI[temp + dif]);
		term1R = RealMul(stateR[temp], -stateI[temp], stateR[temp + dif], stateI[temp + dif]);
		/*temp = i*LengthFactor + ind1;
		term2I = ImagMul(term1R, term1I, stateR[temp + dif], -stateI[temp + dif]);
		term2R = RealMul(term1R, term1I, stateR[temp + dif], -stateI[temp + dif]);
		term1I = ImagMul(term2R, term2I, stateR[temp], stateI[temp]);
		term1R = RealMul(term2R, term2I, stateR[temp], stateI[temp]);*/
		FockSumR += term1R;
		FockSumI += term1I;
		//			}
	}
	//}
	HR[ind1][ind1 + dif] -= U_HF * FockSumR / NormFactor;
	//HI[ind1][ind1 + dif] += sign(dif)*U_HF*FockSumI / NormFactor;
	//HR[ind1 + dif][ind1] -= U_HF*FockSumR / NormFactor;
	HI[ind1 + dif][ind1] -= U_HF*FockSumI / NormFactor;
}

void add_inter(double* HR[], double* HI[], double* HR0[], double* HI0[], int interest_len, double* stateR, double* stateI, int LengthFactor) {
	int i, j, k, dif;
	set_array(HR, HR0, LengthFactor, LengthFactor);
	set_array(HI, HI0, LengthFactor, LengthFactor);
	if (isHubbard) {
		//for (i = 0; i <= 2 * lat_len - 2; i = i + 2) {
		for (i = 0; i <= 2 * lat_len - 1; i = i++) {
			add_Hubbard_inter_aid(HR, HI, stateR, stateI, i, interest_len);
		}
	}
	/*if (isFlip)
	dif = 1;
	else
	dif = 2;*/
	//dif = 1;
	if (LengthFactor == 2 * lat_len)
		dif = 2;
	else
		dif = 1;

	if (isHartr) {
		for (i = 0; i <= LengthFactor - 1; i++) {
			add_Hartr_inter_aid(HR, HI, i, interest_len, stateR, stateI, LengthFactor, dif);
		}
	}
	if (isFock) {
		add_Fock_inter_aid(HR, HI, 0, interest_len, stateR, stateI, LengthFactor, dif);
		if (LengthFactor == 2 * lat_len)
			add_Fock_inter_aid(HR, HI, 1, interest_len, stateR, stateI, LengthFactor, dif);
		for (i = dif; i < LengthFactor - dif; i = i++) {
			add_Fock_inter_aid(HR, HI, i, interest_len, stateR, stateI, LengthFactor, dif);
			add_Fock_inter_aid(HR, HI, i, interest_len, stateR, stateI, LengthFactor, -dif);
		}
		if (LengthFactor == 2 * lat_len)
			add_Fock_inter_aid(HR, HI, LengthFactor - 2, interest_len, stateR, stateI, LengthFactor, -dif);
		add_Fock_inter_aid(HR, HI, LengthFactor - 1, interest_len, stateR, stateI, LengthFactor, -dif);
	}
}

void CurCalcV3(int fluxInd, double* CurLong[], double* CurTran[], int interest_len, Engine *ep) {//another version of calculating the current
	int i, j, temp, MatSize = 2 * lat_len;
	double curUP, curDOWN, finalCur, arg;
	engEvalString(ep, "temp=state0;");
	if (FluxCase == 3) {
		engEvalString(ep, "fluxPh=fluxPh/(2*lat_len);");
		engEvalString(ep, "CurTrnUp(1,fluxVal+1,MyLen-low_lim+1)=sum(-temp(2,1:MyLen).*conj(temp(1,1:MyLen))*exp(1i*fluxPh)+conj(temp(2,1:MyLen)).*temp(1,1:MyLen)*exp(-1i*fluxPh),2);");
		engEvalString(ep, "CurTrnUp(length(state0)/2,fluxVal+1,MyLen-low_lim+1)=sum(-temp(end-1,1:MyLen).*conj(temp(end,1:MyLen))*exp(1i*fluxPh)+conj(temp(end-1,1:MyLen)).*temp(end,1:MyLen)*exp(-1i*fluxPh),2);");
	}
	else if (FluxCase == 4) {
		engEvalString(ep, "f1=f1/lat_len;");
		engEvalString(ep, "f2=f2/lat_len;");
	}
	else if (FluxCase == 6) {
		engEvalString(ep, "CurTrnUp(:,fluxVal+1,MyLen-low_lim+1)=sum(-temp(1:2:end-1,1:MyLen).*conj(temp(2:2:end,1:MyLen))*exp(1i*f1*fluxPh)+conj(temp(1:2:end-1,1:MyLen)).*temp(2:2:end,1:MyLen),2)*exp(-1i*f1*fluxPh);");
		engEvalString(ep, "CurLonUp(:,fluxVal+1,MyLen-low_lim+1)=sum(-temp(1:2:end-3,1:MyLen).*conj(temp(3:2:end-1,1:MyLen))+conj(temp(1:2:end-3,1:MyLen)).*temp(3:2:end-1,1:MyLen),2);");
		engEvalString(ep, "CurLonDn(:,fluxVal+1,MyLen-low_lim+1)=sum(-temp(2:2:end-2,1:MyLen).*conj(temp(4:2:end,1:MyLen))+conj(temp(2:2:end-2,1:MyLen)).*temp(4:2:end,1:MyLen),2);");
	}
	else
		engEvalString(ep, "CurTrnUp(:,fluxVal+1,MyLen-low_lim+1)=sum(-temp(1:2:end-1,1:MyLen).*conj(temp(2:2:end,1:MyLen))+conj(temp(1:2:end-1,1:MyLen)).*temp(2:2:end,1:MyLen),2);");
	//engEvalString(ep, "fluxPh=fluxPh/lat_len;");
	if (FluxCase != 6) {
		engEvalString(ep, "CurLonUp(:,fluxVal+1,MyLen-low_lim+1)=sum(-temp(1:2:end-3,1:MyLen).*conj(temp(3:2:end-1,1:MyLen))*exp(-1i*f1*fluxPh)+conj(temp(1:2:end-3,1:MyLen)).*temp(3:2:end-1,1:MyLen)*exp(1i*f1*fluxPh),2);");
		engEvalString(ep, "CurLonDn(:,fluxVal+1,MyLen-low_lim+1)=sum(-temp(2:2:end-2,1:MyLen).*conj(temp(4:2:end,1:MyLen))*exp(-1i*f2*fluxPh)+conj(temp(2:2:end-2,1:MyLen)).*temp(4:2:end,1:MyLen)*exp(1i*f2*fluxPh),2);");
		if (isRing) {
			engEvalString(ep, "CurLonUp(:,fluxVal+1,MyLen-low_lim+1)= CurLonUp(:,fluxVal+1) + sum(temp(1,1:MyLen).*conj(temp(end-1,1:MyLen))*exp(-1i*fluxPh)-conj(temp(1,1:MyLen)).*temp(end-1,1:MyLen)*exp(1i*fluxPh),2);");
			engEvalString(ep, "CurLonDn(:,fluxVal+1,MyLen-low_lim+1)= CurLonDn(:,fluxVal+1) + sum(temp(2,1:MyLen).*conj(temp(end,1:MyLen))*exp(-1i*fluxPh)-conj(temp(2,1:MyLen)).*temp(end,1:MyLen)*exp(1i*fluxPh),2);");
		}
	}

	//engEvalString(ep, "temp=state0(1:2:end,:)+state0(2:2:end,:);");
	engEvalString(ep, "CurLonUp(:,fluxVal+1,MyLen-low_lim+1)=imag(CurLonUp(:,fluxVal+1,MyLen-low_lim+1));");
	engEvalString(ep, "CurLonDn(:,fluxVal+1,MyLen-low_lim+1)=imag(CurLonDn(:,fluxVal+1,MyLen-low_lim+1));");
	engEvalString(ep, "CurTrnUp(:,fluxVal+1,MyLen-low_lim+1)=imag(CurTrnUp(:,fluxVal+1,MyLen-low_lim+1));");
	//if (interest_len%2==0)
	//engEvalString(ep, "Cur_OprToPer(fluxVal+1,MyLen-low_lim+1)=f1.*sum(CurLonUp(1:MyLen/2,fluxVal+1,MyLen))+f2.*sum(CurLonDn(1:MyLen/2,fluxVal+1,MyLen));");
	engEvalString(ep, "Cur_OprToPer(fluxVal+1,MyLen-low_lim+1)=f1.*sum(CurLonUp(:,fluxVal+1,MyLen-low_lim+1))+f2.*sum(CurLonDn(:,fluxVal+1,MyLen-low_lim+1));");
	/*else
	engEvalString(ep, "Cur_OprToPer(fluxVal+1,MyLen-low_lim+1)=f1.*sum(CurLonUp(1:(MyLen+1)/2,fluxVal+1,MyLen-low_lim+1))+f2.*sum(CurLonDn(1:(MyLen-1)/2,fluxVal+1,MyLen-low_lim+1));");*/
}

void CurCalcV4(int fluxInd, double* CurLong[], double* CurTran[], int interest_len, Engine *ep, int CaseInd) {//another version of calculating the current
	int i, j, temp, MatSize = 2 * lat_len;
	double curUP, curDOWN, finalCur, arg, fluxPh;
	engEvalString(ep, "temp=state;");

	if (CaseInd == 1)//current calculation
		engEvalString(ep, "SumFactor=-1;");
	else if (CaseInd == 2)//charge density calculation
		engEvalString(ep, "SumFactor=1;");

	if (FluxCase == 3) {
		engEvalString(ep, "fluxPh=fluxPh/(2*lat_len);");
		engEvalString(ep, "temp2=temp(2,1:MyLen).*conj(temp(1,1:MyLen))*exp(1i*fluxPh);");
		engEvalString(ep, "tempTrnUp(1)=sum(SumFactor*temp2+conj(temp2),2);");

		engEvalString(ep, "temp2=temp(end-1,1:MyLen).*conj(temp(end,1:MyLen))*exp(1i*fluxPh);");
		engEvalString(ep, "tempTrnUp(length(state0)/2)=sum(SumFactor*temp2+conj(temp2),2);");
	}
	else if (FluxCase == 4) {
		engEvalString(ep, "f1=f1/lat_len;");
		engEvalString(ep, "f2=f2/lat_len;");
	}
	else if (FluxCase == 6) {
		engEvalString(ep, "tempTrnUp=sum(-temp(1:2:end-1,1:MyLen).*conj(temp(2:2:end,1:MyLen))*exp(1i*f1*fluxPh)+conj(temp(1:2:end-1,1:MyLen)).*temp(2:2:end,1:MyLen),2)*exp(-1i*f1*fluxPh);");
		engEvalString(ep, "tempLonUp=sum(-temp(1:2:end-3,1:MyLen).*conj(temp(3:2:end-1,1:MyLen))+conj(temp(1:2:end-3,1:MyLen)).*temp(3:2:end-1,1:MyLen),2);");
		engEvalString(ep, "tempLonDn=sum(-temp(2:2:end-2,1:MyLen).*conj(temp(4:2:end,1:MyLen))+conj(temp(2:2:end-2,1:MyLen)).*temp(4:2:end,1:MyLen),2);");
	}
	else {
		engEvalString(ep, "temp2=temp(1:2:end-1,1:MyLen).*conj(temp(2:2:end,1:MyLen));");
		engEvalString(ep, "tempTrnUp=sum(SumFactor*temp2+conj(temp2),2);");
	}
	//engEvalString(ep, "fluxPh=fluxPh/lat_len;");
	if (FluxCase != 6) {
		engEvalString(ep, "temp2=temp(1:2:end-3,1:MyLen).*conj(temp(3:2:end-1,1:MyLen))*exp(-1i*f1*fluxPh);");
		engEvalString(ep, "tempLonUp=sum(SumFactor*temp2+conj(temp2),2);");
		engEvalString(ep, "temp2=temp(2:2:end-2,1:MyLen).*conj(temp(4:2:end,1:MyLen))*exp(-1i*f2*fluxPh);");
		engEvalString(ep, "tempLonDn=sum(SumFactor*temp2+conj(temp2),2);");
		if (isRing) {
			engEvalString(ep, "temp2=temp(1,1:MyLen).*conj(temp(end-1,1:MyLen))*exp(-1i*f1*fluxPh);");
			engEvalString(ep, "tempLonUp= tempLonUp(:,fluxVal+1) + sum(temp2+SumFactor*conj(temp2),2);");
			engEvalString(ep, "temp2=temp(2,1:MyLen).*conj(temp(end,1:MyLen))*exp(-1i*f2*fluxPh);");
			engEvalString(ep, "tempLonDn= tempLonDn(:,fluxVal+1) + sum(temp2+SumFactor*conj(temp2);");
		}
	}

	//engEvalString(ep, "temp=state0(1:2:end,:)+state0(2:2:end,:);");
	if (CaseInd == 1) {
		engEvalString(ep, "CurLonUp(:,fluxVal+1,MyLen-low_lim+1)=imag(tempLonUp);");
		engEvalString(ep, "CurLonDn(:,fluxVal+1,MyLen-low_lim+1)=imag(tempLonDn);");
		engEvalString(ep, "CurTrnUp(:,fluxVal+1,MyLen-low_lim+1)=imag(tempTrnUp);");
	}
	else if (CaseInd == 2) {
		engEvalString(ep, "ChargeLonUp(:,fluxVal+1,MyLen-low_lim+1)=real(tempLonUp);");
		engEvalString(ep, "ChargeLonDn(:,fluxVal+1,MyLen-low_lim+1)=real(tempLonDn);");
		engEvalString(ep, "ChargeTrnUp(:,fluxVal+1,MyLen-low_lim+1)=real(tempTrnUp);");
	}

	//if (interest_len%2==0)
	//engEvalString(ep, "Cur_OprToPer(fluxVal+1,MyLen-low_lim+1)=f1.*sum(CurLonUp(1:MyLen/2,fluxVal+1,MyLen))+f2.*sum(CurLonDn(1:MyLen/2,fluxVal+1,MyLen));");
	engEvalString(ep, "Cur_OprToPer(fluxVal+1,MyLen-low_lim+1)=f1.*sum(CurLonUp(:,fluxVal+1,MyLen-low_lim+1))+f2.*sum(CurLonDn(:,fluxVal+1,MyLen-low_lim+1));");
	engEvalString(ep, "Charge_OprToPer(fluxVal+1,MyLen-low_lim+1)=f1.*sum(CurLonUp(:,fluxVal+1,MyLen-low_lim+1))+f2.*sum(CurLonDn(:,fluxVal+1,MyLen-low_lim+1));");
	/*else
	engEvalString(ep, "Cur_OprToPer(fluxVal+1,MyLen-low_lim+1)=f1.*sum(CurLonUp(1:(MyLen+1)/2,fluxVal+1,MyLen-low_lim+1))+f2.*sum(CurLonDn(1:(MyLen-1)/2,fluxVal+1,MyLen-low_lim+1));");*/
}

void ConductivityCalc(Engine *ep) {//another version of calculating the current
	engEvalString(ep, "CalcState=state;");
	engEvalString(ep, "GS=CalcState(:,1)';");
	engEvalString(ep, "term1=GS(3:2:end-1)*conj(CalcState(1:2:end-3,2:MyLen))*exp(-1i*f1*fluxPh);");
	engEvalString(ep, "term2=GS(1:2:end-3)*conj(CalcState(3:2:end-1,2:MyLen))*exp(1i*f1*fluxPh);");
	engEvalString(ep, "LonUp=term1-term2;");
	engEvalString(ep, "term1=GS(4:2:end)*conj(CalcState(2:2:end-2,2:MyLen))*exp(-1i*f2*fluxPh);");
	engEvalString(ep, "term2=GS(2:2:end-2)*conj(CalcState(4:2:end,2:MyLen))*exp(1i*f2*fluxPh);");
	engEvalString(ep, "LonDn=term1-term2;");
	engEvalString(ep, "term1=conj(GS(1:2:end-1))*CalcState(2:2:end,2:MyLen);");
	engEvalString(ep, "term2=conj(GS(2:2:end))*CalcState(1:2:end-1,2:MyLen);");
	engEvalString(ep, "Trn=term1-term2;");
	if (isRing) {
		engEvalString(ep, "term1=conj(GS(end-1))*CalcState(1,2:MyLen)*exp(-1i*f1*fluxPh);");
		engEvalString(ep, "term2=conj(GS(1))*CalcState(end-1,2:MyLen)*exp(1i*f1*fluxPh);");
		engEvalString(ep, "tempLonUp= tempLonUp(:,fluxVal+1) + sum(term1-term2,2);");
		engEvalString(ep, "term1=conj(GS(end))*CalcState(2,2:MyLen)*exp(-1i*f2*fluxPh);");
		engEvalString(ep, "term2=conj(GS(2))*CalcState(end,2:MyLen)*exp(1i*f2*fluxPh);");
		engEvalString(ep, "tempLonDn= tempLonDn(:,fluxVal+1) + sum(term1-term2,2);");
	}
	engEvalString(ep, "condUp(fluxVal+1,MyLen-1)=2*imag(dot(LonUp./((E(2:MyLen)'-E(1)).^2),Trn));");
	engEvalString(ep, "condDn(fluxVal+1,MyLen-1)=2*imag(dot(LonDn./((E(2:MyLen)'-E(1)).^2),Trn));");
}

void CurUpdate(double* CurLong[], double* CurTran[], Engine *ep) {//another version of calculating the current
	int i, j, temp, MatSize = 2 * lat_len, len;
	double curUP, curDOWN, finalCur, arg;
	int FluxLim;
	//double *tempCurLonUp, *tempCurLonDn, *tempCurTrnUp;
	if (isSOC == 1 && isConstSOC == 0)
		len = SOCDiv;
	else
		len = fluxDiv;
	if (isMagnet)
		FluxLim = fluxDiv;
	else
		FluxLim = 1;
	double *tempCurLonUp = malloc(MatSize * len * sizeof(double));
	double *tempCurLonDn = malloc(MatSize * len * sizeof(double));
	//double *tempCurTrnUp = malloc(MatSize * fluxDiv * sizeof(double));
	double *tempCurTrnUp = malloc(lat_len * len * sizeof(double));
	mxArray *MexCurTrnUp = engGetVariable(ep, "CurTrnUp");
	mxArray *MexCurLonDn = engGetVariable(ep, "CurLonDn");
	mxArray *MexCurLonUp = engGetVariable(ep, "CurLonUp");
	tempCurTrnUp = (double*)mxGetPr(MexCurTrnUp);
	tempCurLonUp = (double*)mxGetPr(MexCurLonUp);
	tempCurLonDn = (double*)mxGetPr(MexCurLonDn);
	for (i = 0; i < lat_len; i++) {
		for (j = 0; j < FluxLim; j++) {
			temp = i + j * lat_len;
			CurTran[i][j] = tempCurTrnUp[temp];
		}
	}
	for (i = 0; i < 2 * lat_len - 2; i = i + 2) {
		for (j = 0; j < FluxLim; j++) {
			temp = i / 2 + j * (lat_len - 1);
			CurLong[i][j] = tempCurLonUp[temp];
		}
	}
	for (i = 1; i < 2 * lat_len - 2; i = i + 2) {
		for (j = 0; j < FluxLim; j++) {
			temp = (i - 1) / 2 + j * (lat_len - 1);
			CurLong[i][j] = tempCurLonDn[temp];
		}
	}
	free(tempCurLonDn); free(tempCurLonUp); free(tempCurTrnUp);
	mxDestroyArray(MexCurTrnUp); mxDestroyArray(MexCurLonUp); mxDestroyArray(MexCurLonDn);
}

void NNHopping(double* HR[], double* HI[]) {//adds spin orbit coupling too the matrix
	int i;
	HR[0][4] = t_next;
	HR[1][5] = t_next;
	HR[2][6] = t_next;
	HR[3][7] = t_next;
	for (i = 4; i < 2 * lat_len - 4; i++) {
		HR[i][i - 4] = t_next;
		HR[i][i + 4] = t_next;
	}
	HR[2 * lat_len - 4][2 * lat_len - 8] = t_next;
	HR[2 * lat_len - 3][2 * lat_len - 7] = t_next;
	HR[2 * lat_len - 2][2 * lat_len - 6] = t_next;
	HR[2 * lat_len - 1][2 * lat_len - 5] = t_next;
	if (isRing) {
		HR[2 * lat_len - 2][2] = t_next;
		HR[2 * lat_len - 1][3] = t_next;
		HR[0][2 * lat_len - 4] = t_next;
		HR[1][2 * lat_len - 3] = t_next;
	}
}

void SOC_V1(double* HR[], double* HI[]) {//adds spin orbit coupling too the matrix
	int i;
	HR[0][3] = SOC_amp;
	HR[1][2] = -SOC_amp;
	for (i = 2; i < 2 * lat_len - 2; i = i + 2) {
		HR[i][i - 1] = -SOC_amp;
		HR[i][i + 3] = SOC_amp;
		HR[i + 1][i - 2] = SOC_amp;
		HR[i + 1][i + 2] = -SOC_amp;
	}
	HR[2 * lat_len - 2][2 * lat_len - 3] = -SOC_amp;
	HR[2 * lat_len - 1][2 * lat_len - 4] = SOC_amp;
	if (isRing) {
		HR[0][2 * lat_len - 1] = -SOC_amp;
		HR[1][2 * lat_len - 2] = SOC_amp;
		HR[2 * lat_len - 1][0] = SOC_amp;
		HR[2 * lat_len - 2][1] = -SOC_amp;
	}
}

void SOC_V2(double* HR[], double* HI[]) {//adds spin orbit coupling too the matrix
	int i;
	double phase;
	phase = (double)M_PI / lat_len;
	HR[0][3] = SOC_amp * sin(phase);
	HI[0][3] = -SOC_amp * cos(phase);
	HR[1][2] = -SOC_amp * sin(phase);
	HI[1][2] = -SOC_amp * cos(phase);
	for (i = 2; i < 2 * lat_len - 2; i = i + 2) {
		phase = (double)M_PI*(i + 1) / lat_len;
		HR[i][i - 1] = -SOC_amp * sin(phase);
		HI[i][i - 1] = -SOC_amp * cos(phase);
		HR[i][i + 3] = SOC_amp * sin(phase);
		HI[i][i + 3] = -SOC_amp * cos(phase);
		HR[i + 1][i - 2] = SOC_amp * sin(phase);
		HI[i + 1][i - 2] = -SOC_amp * cos(phase);
		HR[i + 1][i + 2] = -SOC_amp * sin(phase);
		HI[i + 1][i + 2] = -SOC_amp * cos(phase);
	}
	phase = (double)M_PI*(4 * lat_len - 3) / lat_len;
	HR[2 * lat_len - 2][2 * lat_len - 3] = -SOC_amp * sin(phase);
	HI[2 * lat_len - 2][2 * lat_len - 3] = -SOC_amp * cos(phase);
	HR[2 * lat_len - 1][2 * lat_len - 4] = SOC_amp * sin(phase);
	HI[2 * lat_len - 1][2 * lat_len - 4] = -SOC_amp * cos(phase);
	if (isRing) {
		HR[0][2 * lat_len - 1] = -SOC_amp * sin(phase);
		HI[0][2 * lat_len - 1] = -SOC_amp * cos(phase);
		HR[1][2 * lat_len - 2] = SOC_amp * sin(phase);
		HI[1][2 * lat_len - 2] = -SOC_amp * cos(phase);
		HR[2 * lat_len - 1][0] = SOC_amp * sin(phase);
		HI[2 * lat_len - 1][0] = SOC_amp * cos(phase);
		HR[2 * lat_len - 2][1] = -SOC_amp * sin(phase);
		HI[2 * lat_len - 2][1] = -SOC_amp * sin(phase);
	}
}

void PutValInMat(double* HR[], double* HI[], int ind, double LocFluxVal, double sigma, double f) {
	double temp;
	if (isRing == 1 && (isFlip == 0 || t_flip == 0 || FluxCase==8 || FluxCase == 9))
		LocFluxVal = LocFluxVal / lat_len;
	if (ind % 2 == 0) {
		temp = (double)ind / 2;
		HR[ind][ind] = lambda * cos(2.0*M_PI*sigma*(ind / 2 + 1) + phi);
		if (FluxCase >= 6) {
			HR[ind][ind + 1] = t_flip * cos((temp + 1)*(LocFluxVal + alpha * (temp - 1) * 2 * M_PI));
			HI[ind][ind + 1] = -t_flip * sin((temp + 1)*(LocFluxVal + alpha * (temp - 1) * 2 * M_PI));
			HR[ind + 1][ind] = t_flip * cos((temp + 1)*(LocFluxVal + alpha * (temp - 1) * 2 * M_PI));
			HI[ind + 1][ind] = t_flip * sin((temp + 1)*(LocFluxVal + alpha * (temp - 1) * 2 * M_PI));
		}
		else if (isFlip) {
			HR[ind][ind + 1] = t_flip;
			HR[ind + 1][ind] = t_flip;
		}
	}
	else {
		temp = (double)(ind - 1) / 2;
		HR[ind][ind] = lambda * cos(2.0*M_PI*sigma*((ind + 1) / 2) + phi);
	}
	if (ind < 2 * lat_len - 2) {//hopping to next site
		if (FluxCase < 6) {
			HR[ind][ind + 2] = t_amp * cos(f*(LocFluxVal + alpha * temp * 2 * M_PI));
			HI[ind][ind + 2] = t_amp * sin(f*(LocFluxVal + alpha * temp * 2 * M_PI));
		}
		else
			HR[ind][ind + 2] = -t_amp;
	}
	if (ind > 1) {
		if (FluxCase < 6) {
			HR[ind][ind - 2] = t_amp * cos(f*(LocFluxVal + alpha * (temp - 1) * 2 * M_PI));
			HI[ind][ind - 2] = -t_amp * sin(f*(LocFluxVal + alpha * (temp - 1) * 2 * M_PI));
		}
		else {
			HR[ind][ind - 2] = -t_amp;
		}
	}
}

void createMH_reg(double sigma, double* HR0[], double* HI0[], int isConst) {
	int i, ind1, dif, low_lim;
	int MatSize = 2 * lat_len;
	double Loc_f1, Loc_f2;//in case 4, i divide the value of f1 and f2 by lat_len
						  //double DoubleFlux=2*fluxVal;
	if (FluxCase == 4) {
		Loc_f1 = f1 / lat_len;
		Loc_f2 = f2 / lat_len;
	}
	else {
		Loc_f1 = f1;
		Loc_f2 = f2;
	}

	i = 0;
	PutValInMat(HR0, HI0, i, fluxVal, sigma, Loc_f1);
	i = 1;
	PutValInMat(HR0, HI0, i, fluxVal, sigma, Loc_f2);
	if (isConst) {
		dif = 2;
		low_lim = 2;
	}
	else {
		dif = 4;
		low_lim = 4;
	}
	for (i = low_lim; i < 2 * lat_len - 4; i = i + dif)//create the matrix
	{
		PutValInMat(HR0, HI0, i, fluxVal, sigma, f1);
		PutValInMat(HR0, HI0, i + 1, fluxVal, sigma, f2);
	}
	if ((lat_len / 2) % 2 == 1 && isConst == 0) {
		ind1 = MatSize - 2;
		PutValInMat(HR0, HI0, ind1, fluxVal, sigma, Loc_f1);
		ind1 = MatSize - 1;
		PutValInMat(HR0, HI0, ind1, fluxVal, sigma, Loc_f2);
	}
	else if (isConst == 1) {
		ind1 = MatSize - 2;
		PutValInMat(HR0, HI0, ind1, fluxVal, sigma, Loc_f1);
		ind1 = MatSize - 1;
		PutValInMat(HR0, HI0, ind1, fluxVal, sigma, Loc_f2);
		ind1 = MatSize - 4;
		PutValInMat(HR0, HI0, ind1, fluxVal, sigma, Loc_f1);
		ind1 = MatSize - 3;
		PutValInMat(HR0, HI0, ind1, fluxVal, sigma, Loc_f2);
	}
	else {
		ind1 = MatSize - 4;
		PutValInMat(HR0, HI0, ind1, fluxVal, sigma, Loc_f1);
		ind1 = MatSize - 3;
		PutValInMat(HR0, HI0, ind1, fluxVal, sigma, Loc_f2);
	}
	//HR0[10][11] = 0.01;
	//HR0[11][10] = 0.01;
}

void createMH_loop(double sigma, double* HR0[], double* HI0[]) {
	int i, ind1, dif, low_lim;
	int MatSize = 2 * lat_len;
	//double DoubleFlux=2*fluxVal
	fluxVal = fluxVal / (2 * lat_len);
	HR0[0][2] = t_amp * cos(f1*fluxVal);
	HI0[0][2] = -t_amp * sin(f1*fluxVal);
	HR0[0][1] = t_amp * cos(f1*fluxVal);
	HI0[0][1] = t_amp * sin(f1*fluxVal);
	HR0[1][0] = t_amp * cos(f1*fluxVal);
	HI0[1][0] = -t_amp * sin(f1*fluxVal);
	HR0[1][3] = t_amp * cos(f1*fluxVal);
	HI0[1][3] = t_amp * sin(f1*fluxVal);
	for (i = 2; i < 2 * lat_len - 2; i = i + 2)//create the matrix
	{
		HR0[i][i + 2] = t_amp * cos(f1*fluxVal);
		HI0[i][i + 2] = -t_amp * sin(f1*fluxVal);
		HR0[i][i - 2] = t_amp * cos(f1*fluxVal);
		HI0[i][i - 2] = t_amp * sin(f1*fluxVal);
	}
	for (i = 3; i < 2 * lat_len - 2; i = i + 2)//create the matrix
	{
		HR0[i][i + 2] = t_amp * cos(f1*fluxVal);
		HI0[i][i + 2] = t_amp * sin(f1*fluxVal);
		HR0[i][i - 2] = t_amp * cos(f1*fluxVal);
		HI0[i][i - 2] = -t_amp * sin(f1*fluxVal);
	}
	HR0[2 * lat_len - 1][2 * lat_len - 2] = t_amp * cos(f1*fluxVal);
	HI0[2 * lat_len - 1][2 * lat_len - 2] = t_amp * sin(f1*fluxVal);
	HR0[2 * lat_len - 1][2 * lat_len - 3] = t_amp * cos(f1*fluxVal);
	HI0[2 * lat_len - 1][2 * lat_len - 3] = -t_amp * sin(f1*fluxVal);
	HR0[2 * lat_len - 2][2 * lat_len - 4] = t_amp * cos(f1*fluxVal);
	HI0[2 * lat_len - 2][2 * lat_len - 4] = t_amp * sin(f1*fluxVal);
	HR0[2 * lat_len - 2][2 * lat_len - 1] = t_amp * cos(f1*fluxVal);
	HI0[2 * lat_len - 2][2 * lat_len - 1] = -t_amp * sin(f1*fluxVal);
}

void createMH_reversed(double sigma, double* HR0[], double* HI0[]) {
	int i, ind1;
	int MatSize = 2 * lat_len;
	//double DoubleFlux = 2 * fluxVal;
	for (i = 2; i < 2 * lat_len - 2; i = i + 4)//create the matrix
	{
		PutValInMat(HR0, HI0, i, fluxVal, sigma, f1);
		PutValInMat(HR0, HI0, i + 1, fluxVal, sigma, f2);
	}
}

void createMH_non(double sigma, double* HR0[], double* HI0[]) {// no flux in the system
	int i, ind1;
	int MatSize = 2 * lat_len;
	for (i = 2; i < 2 * lat_len - 2; i = i + 4)//create the matrix
	{
		PutValInMat(HR0, HI0, i, 0.0, sigma, f1);
	}
	if ((lat_len / 2) % 2 == 1) {
		ind1 = MatSize - 4;
		PutValInMat(HR0, HI0, ind1, 0.0, sigma, f1);
		ind1 = MatSize - 3;
		PutValInMat(HR0, HI0, ind1, 0.0, sigma, f2);
	}
	else {
		ind1 = MatSize - 2;
		PutValInMat(HR0, HI0, ind1, 0.0, sigma, f1);
		ind1 = MatSize - 1;
		PutValInMat(HR0, HI0, ind1, 0.0, sigma, f2);
	}
}

void RingCase(double* HR0[], double* HI0[], int LengthFactor) {
	//if (isFlip == 0 || t_flip==0.0)
	//fluxVal = fluxVal / lat_len;
	if (FluxCase != 8 && FluxCase !=9) {
		HR0[0][LengthFactor - 2] = t_amp * cos(f1*(fluxVal + alpha * 2 * M_PI));
		HI0[0][LengthFactor - 2] = -t_amp * sin(f1*(fluxVal + alpha * 2 * M_PI));
		HR0[LengthFactor - 2][0] = t_amp * cos(f1*(fluxVal + alpha * 2 * M_PI));
		HI0[LengthFactor - 2][0] = t_amp * sin(f1*(fluxVal + alpha * 2 * M_PI));
		HR0[1][LengthFactor - 1] = t_amp * cos(f2*(fluxVal + alpha * 2 * M_PI));// +2 * M_PI / fluxDiv);
		HI0[1][LengthFactor - 1] = -t_amp * sin(f2*(fluxVal + alpha * 2 * M_PI));// +2 * M_PI / fluxDiv);
		HR0[LengthFactor - 1][1] = t_amp * cos(f2*(fluxVal + alpha * 2 * M_PI));// + 2 * M_PI / fluxDiv);
		HI0[LengthFactor - 1][1] = t_amp * sin(f2*(fluxVal + alpha * 2 * M_PI));// +2 * M_PI / fluxDiv);
	}
	else
	{
		double LocFluxVal = (double)fluxVal / lat_len;
		HR0[0][LengthFactor - 1] = t_amp * cos(LocFluxVal + alpha * 2 * M_PI);// +2 * M_PI / fluxDiv);
		HI0[0][LengthFactor - 1] = -t_amp * sin(LocFluxVal + alpha * 2 * M_PI);// +2 * M_PI / fluxDiv);
		HR0[LengthFactor - 1][0] = t_amp * cos(LocFluxVal + alpha * 2 * M_PI);// + 2 * M_PI / fluxDiv);
		HI0[LengthFactor - 1][0] = t_amp * sin(LocFluxVal + alpha * 2 * M_PI);// +2 * M_PI / fluxDiv);
	}
}

void CreateSingleChain(double sigma, double* HR0[],double* HI0[]) {//creates the hamiltonian of a single chain
	int i, ind1, dif, low_lim;
	int MatSize = 2 * lat_len;
	double LocFluxVal = fluxVal / lat_len;
	HR0[0][0] = lambda * cos(2.0*M_PI*sigma + phi);
	HR0[0][1] = t_amp * cos(LocFluxVal + alpha * 2 * M_PI);
	HI0[0][1] = t_amp * sin(LocFluxVal + alpha * 2 * M_PI);
	for (i = 1; i < lat_len - 1; i++)//create the matrix
	{
		HR0[i][i] = lambda * cos(2.0*M_PI*sigma*(i+1) + phi);
		HR0[i][i + 1] = t_amp * cos(LocFluxVal + alpha * 2 * M_PI);
		HR0[i][i - 1] = t_amp * cos(LocFluxVal + alpha * 2 * M_PI);
		HI0[i][i + 1] = t_amp * sin(LocFluxVal + alpha * 2 * M_PI);
		HI0[i][i - 1] = -t_amp * sin(LocFluxVal + alpha * 2 * M_PI);
	}
	HR0[lat_len - 1][lat_len - 1] = lambda * cos(2.0*M_PI*sigma*(lat_len) + phi);
	HR0[lat_len - 1][lat_len - 2] = t_amp * cos(LocFluxVal + alpha * 2 * M_PI);
	HI0[lat_len - 1][lat_len - 2] = -t_amp * sin(LocFluxVal + alpha * 2 * M_PI);
}

void createMH(double sigma, double* HR0[], double* HI0[], int LengthFactor)
{
	int i;
	int MatSize = 2 * lat_len;
	//init_array(HR0,MatSize,MatSize);
	//init_array(HI0,MatSize,MatSize);
	switch (FluxCase) {
	case 0://flux identical for each square
		createMH_reg(sigma, HR0, HI0, 1);
		break;
	case 1://flux reversed in each consqutive square
		createMH_reg(sigma, HR0, HI0, 0);
		createMH_reversed(sigma, HR0, HI0);
		break;
	case 2://flux exists only in even squares
		createMH_reg(sigma, HR0, HI0, 0);
		createMH_non(sigma, HR0, HI0);
		break;
	case 3://one big loop, cut all the flips in the middle
		createMH_loop(sigma, HR0, HI0);
		break;
	case 4:
		createMH_reg(sigma, HR0, HI0, 1);
		break;
	case 5://One plaquette with an impurity
		createMH_reg(sigma, HR0, HI0, 1);
		//HR0[lat_len - 1][lat_len - 1] = 1.0;
		//HR0[lat_len - 2][lat_len - 2] = 1.0;
		if (lat_len % 2 == 0){
			HR0[lat_len][lat_len] = (2 * isSymmetric - 1)*ImpureAmp;
			HR0[lat_len + 1][lat_len + 1] = ImpureAmp;
		}
		else
		{
			HR0[lat_len][lat_len] = (2 * isSymmetric - 1)*ImpureAmp;
			HR0[lat_len - 1][lat_len - 1] = ImpureAmp;
		}
		//HR0[lat_len+1][lat_len + 1] = -ImpureAmp;
		//HR0[lat_len + 2][lat_len + 2] = -ImpureAmp;
		break;
	case 6://flux is on the vertical connection, not longitudinal
		createMH_reg(sigma, HR0, HI0, 1);
		break;
	case 7://flux is on the vertical connection, not longitudinal, with impurity
		createMH_reg(sigma, HR0, HI0, 1);
		HR0[lat_len][lat_len] = 0.01;
		HR0[lat_len - 1][lat_len - 1] = 0.01;
		HR0[lat_len + 1][lat_len + 1] = 0.01;
		HR0[lat_len + 2][lat_len + 2] = 0.01;
		break;
	case 8: //single chain
		CreateSingleChain(sigma, HR0,HI0);
		break;
	case 9: //single chain with impurity
		CreateSingleChain(sigma, HR0, HI0);
		HR0[0][0] = ImpureAmp;
		break;
	}
	if (isNext) {//add next nearest neighbour hopping terms
		NNHopping(HR0, HI0);
	}
	if (isSOC) {//add SOC terms
		if (isRing)
			SOC_V2(HR0, HI0);
		else
			SOC_V1(HR0, HI0);
	}
	if (isRing)
	{
		RingCase(HR0, HI0,LengthFactor);
	}
}

void FindInitStateMag(double sigma, double* r_H[], double* i_H[], double *E0[], Engine *ep, int LengthFactor)//,int interest_len)
{
	int i, j;
	mxArray  *mexR = NULL, *mexI = NULL;// , *mexE = NULL;//,*mexstate,*mexE;
	createMH(sigma, r_H, i_H,LengthFactor);
	mxArray *mexE = mxCreateNumericMatrix(LengthFactor, 1, mxDOUBLE_CLASS, mxREAL);
	mexR = mxCreateNumericMatrix(LengthFactor, LengthFactor, mxDOUBLE_CLASS, mxREAL);
	mexI = mxCreateNumericMatrix(LengthFactor, LengthFactor, mxDOUBLE_CLASS, mxREAL);
	for (i = 0; i < LengthFactor; i++) {
		for (j = 0; j < LengthFactor; j++) {
			((double*)mxGetPr(mexR))[j + LengthFactor * i] = r_H[j][i];
			((double*)mxGetPr(mexI))[j + LengthFactor * i] = i_H[j][i];
		}
	}
	engPutVariable(ep, "HR", mexR);
	engPutVariable(ep, "HI", mexI);
	//engPutVariable(ep, "interest_len", len);
	engEvalString(ep, "H=complex(HR,HI);");
	engEvalString(ep, "[state0,	E]=eig(H);");
	//engEvalString(ep, "[E,IX]=sort(diag(real(E)));");
	//engEvalString(ep, "state0=state0(:,IX);");
	mexE = engGetVariable(ep, "E");
	//for(i=0;i<2*lat_len;i++)
	//	*E0[i] = (mxGetPr(mexE))[i];
	for (i = 0; i < LengthFactor; i++) {
		E0[0][i] = ((double*)mxGetPr(mexE))[i];
	}
	mxDestroyArray(mexR); mxDestroyArray(mexI);// mxDestroyArray(mexE);
											   //MakeInitState (interest_len,ep);
}

//void MakeInitState(int interest_len,double **state,double **state_R,double **state_I Engine *ep) {
void MakeInitState(double **state_R, double **state_I, Engine *ep, int LengthFactor) {
	int i, j, ind;
	//engEvalString(ep, "state=abs(state0(:,1:MyLen));");
	//engEvalString(ep, "stateR=real(state0(:,1:MyLen));");
	//engEvalString(ep, "stateI=imag(state0(:,1:MyLen));");
	engEvalString(ep, "stateR=real(state0);");
	engEvalString(ep, "stateI=imag(state0);");
	mxArray *mexstateR = engGetVariable(ep, "stateR");
	mxArray *mexstateI = engGetVariable(ep, "stateI");
	for (i = 0; i < LengthFactor; i++) {
		for (j = 0; j < LengthFactor; j++) {
			state_R[0][i*LengthFactor + j] = ((double*)mxGetPr(mexstateR))[i*LengthFactor + j];
			state_I[0][i*LengthFactor + j] = ((double*)mxGetPr(mexstateI))[i * LengthFactor + j];
		}
	}
	mxDestroyArray(mexstateR); mxDestroyArray(mexstateI);
}

void mixedState(Engine *ep, int interest_len) {
	if (interest_len % 2 == 0)
		engEvalString(ep, "state=abs(state(:,1:MyLen));");
	else {
		engEvalString(ep, "state_rest=abs(state(:,MyLen));");
		engEvalString(ep, "state=abs(state(:,1:MyLen-1));");
	}

	engEvalString(ep, "temp=reshape(state,[],2);");
	engEvalString(ep, "state1=temp(:,2)+temp(:,1);");//the sum part
	engEvalString(ep, "state1=state1./sqrt(sum(abs(state1).^2));");//normalization
	engEvalString(ep, "state1=reshape(state1,size(state,1),[]);");

	engEvalString(ep, "state2=temp(:,2)-temp(:,1);");//the sum part
	engEvalString(ep, "state2=state2./sqrt(sum(abs(state2).^2));");//normalization
	engEvalString(ep, "state2=reshape(state2,size(state,1),[]);");

	engEvalString(ep, "state=[state1,state2];");
	if (interest_len % 2 != 0) {
		engEvalString(ep, "state=[state1,state2,state_rest];");
		engEvalString(ep, "clear state_rest;");
	}
	engEvalString(ep, "clear state1 state2 temp;");
}

void converg_func(double* HR[], double* HI[], double* HR0[], double* HI0[], int interest_len, double **E, double **stateR, double **stateI, Engine *ep, int LengthFactor) {
	int i, j;
	mwSize MatSize = (mwSize)LengthFactor;
	mxArray *mxHR = mxCreateNumericMatrix(MatSize, MatSize, mxDOUBLE_CLASS, mxREAL);//mxCreateNumericMatrix(2 * lat_len, 2 * lat_len, mxDOUBLE_CLASS, mxREAL);
	mxArray *mxHI = mxCreateNumericMatrix(MatSize, MatSize, mxDOUBLE_CLASS, mxREAL);
	add_inter(HR, HI, HR0, HI0, interest_len, *stateR, *stateI, LengthFactor);
	if (mxHR == NULL || mxHI == NULL) {
		printf("error in mx");
		return;
	}
	for (i = 0; i < LengthFactor; i++)
		for (j = 0; j < LengthFactor; j++) {
			((double*)mxGetPr(mxHR))[j + i * LengthFactor] = HR[j][i];
			((double*)mxGetPr(mxHI))[j + i * LengthFactor] = HI[j][i];
		}
	engPutVariable(ep, "HR", mxHR);
	engPutVariable(ep, "HI", mxHI);
	engEvalString(ep, "H=complex(HR,HI);");
	engEvalString(ep, "[state,E]=eig(H);");
	//engEvalString(ep, "[E,IX]=sort(diag(real(E)));");
	//engEvalString(ep, "state=state(:,IX);");
	//engEvalString(ep, "StateMat(:,:,fluxVal+1,MyLen-low_lim+1)=state;");//save the state data into a final 4D matrix
	if (isMixed)
		mixedState(ep, interest_len);
	/*else {
	engEvalString(ep, "stateR=real(state(:,1:MyLen));");
	engEvalString(ep, "stateI=imag(state(:,1:MyLen));");
	engEvalString(ep, "clear state;");
	}*/
	engEvalString(ep, "stateR=real(state);");
	engEvalString(ep, "stateI=imag(state);");
	//engEvalString(ep, "clear state;");

	mxArray *mexstateR = engGetVariable(ep, "stateR");
	mxArray *mexstateI = engGetVariable(ep, "stateI");
	mxArray *mexE = engGetVariable(ep, "E");
	for (i = 0; i < LengthFactor; i++) {
		E[0][i] = ((double*)mxGetPr(mexE))[i];
	}
	for (i = 0; i < LengthFactor; i++) {
		for (j = 0; j < LengthFactor; j++) {
			stateR[0][i*LengthFactor + j] = ((double*)mxGetPr(mexstateR))[i*LengthFactor + j];
			stateI[0][i *LengthFactor + j] = ((double*)mxGetPr(mexstateI))[i*LengthFactor + j];
		}
	}
	//*E = (double*)mxGetPr(mexE);
	mxDestroyArray(mexE);
	mxDestroyArray(mxHR);
	mxDestroyArray(mxHI);
	//mxDestroyArray(mexstate_origR);
	//mxDestroyArray(mexstate_origI);
	mxDestroyArray(mexstateR); mxDestroyArray(mexstateI);
}

void changing_flux(double* HR[], double* HI[], double* HR0[], double* HI0[], double** finalE, double* CurLon[], double* CurTrn[], int MaxFlux, int cmpnum, double sigma, Engine *ep, FILE* fpLog) {
	int i, count, CurLen, fluxInd, MatSize = 2 * lat_len, interest_len, count_dif, len_dif, LengthFactor;
	double dif, temp, difLast, matmul;
	mxArray *MexFlux = NULL, *mexlen = NULL, *mexPh = NULL, *mexLatLen = NULL, *mexf1 = NULL, *mexf2 = NULL, *mexHigh = NULL, *mexLow = NULL;
	double *E= (double*)calloc(MatSize, sizeof(double)) , *E0 = (double*)calloc(MatSize, sizeof(double)), *stateR= (double*)calloc(MatSize*MatSize, sizeof(double)) , *stateI= (double*)calloc(MatSize*MatSize, sizeof(double));//*state_origR,*state_origI
									 //#pragma omp parallel for
									 //for (fluxInd = 1; fluxInd < 2; fluxInd++) {

									 /*if (isFlip == 0 && isHubbard == 0)
									 len_dif = 2;//if there is no connection at all between the legs, odd fillings mean nothing
									 else
									 len_dif = 1;*/
	if (FluxCase != 8 && FluxCase != 9) {
		LengthFactor = 2 * lat_len;
	}
	else {
		LengthFactor = lat_len;
		if (LengthFactor < max_fill) {
			printf("filling value is not comaptible with case of two seperated chains - change maximal filling");
			return 0;
		}
	}

	//for (fluxInd = 0; fluxInd < MaxFlux; fluxInd++) {
	for (fluxInd = 1; fluxInd < 2; fluxInd++) {
		//E0 = (double*)calloc(MatSize, sizeof(double));
		//stateR = (double*)calloc(MatSize*MatSize, sizeof(double));
		//stateI = (double*)calloc(MatSize*MatSize, sizeof(double));
		memset(stateR, 0, MatSize*MatSize*sizeof(double));
		memset(stateI, 0, MatSize*MatSize*sizeof(double));
		memset(E0, 0, MatSize*sizeof(double));
		//state_origR = (double*)calloc(MatSize*MatSize, sizeof(double));
		//state_origI = (double*)calloc(MatSize*MatSize, sizeof(double));
		mexLatLen = mxCreateDoubleScalar((double)lat_len);
		engPutVariable(ep, "lat_len", mexLatLen);
		MexFlux = mxCreateDoubleScalar((double)fluxInd);
		engPutVariable(ep, "fluxVal", MexFlux);// this is the flux index
		fluxVal = (double)M_PI*(fluxInd << 1) / fluxDiv;
		mexPh = mxCreateDoubleScalar((double)fluxVal);
		engPutVariable(ep, "fluxPh", mexPh);// this is the flux phase
		mexHigh = mxCreateDoubleScalar((double)max_fill);
		engPutVariable(ep, "high_lim", mexHigh);// this is the highest filling
		mexLow = mxCreateDoubleScalar((double)min_fill);
		engPutVariable(ep, "low_lim", mexLow);// this is the lowest filling
		mexf2 = mxCreateDoubleScalar((double)f2);
		mexf1 = mxCreateDoubleScalar((double)f1);
		engPutVariable(ep, "f1", mexf1);
		engPutVariable(ep, "f2", mexf2);
		init_array(HR0, MatSize, MatSize);
		init_array(HI0, MatSize, MatSize);
		FindInitStateMag(sigma, HR0, HI0, &E0, ep, LengthFactor);//, interest_len);
		MakeInitState(&stateR, &stateI, ep, LengthFactor);
		mxDestroyArray(MexFlux); mxDestroyArray(mexf1); mxDestroyArray(mexf2);
		for (interest_len = min_fill; interest_len <= max_fill; interest_len++)//=len_dif)
		{
			cmpnum = interest_len;
			if (fluxInd % 10 == 0 && (interest_len-min_fill)%10==0) {
				system("cls");
				printf("impurity=%f,filling=%d,flux index=%d\n",ImpureAmp, interest_len, fluxInd);
			}
			mexlen = mxCreateDoubleScalar((double)interest_len);
			engPutVariable(ep, "MyLen", mexlen);
			dif = 1.0;
			count = 0;
			//MakeInitState(interest_len, &state, ep);
			temp = E0[cmpnum - 1];
			if (isHubbard || isFock || isHartr) {
				init_array(HR, MatSize, MatSize);
				init_array(HI, MatSize, MatSize);
				//E = (double*)calloc(MatSize, sizeof(double));
				memset(E, 0, MatSize*sizeof(double));
				count_dif = 0;
				difLast = 0;
				while ((dif > converg_lim && count < IterLim) && count_dif < 25)
				{
					converg_func(HR, HI, HR0, HI0, interest_len, &E, &stateR, &stateI, ep, LengthFactor);
					dif = fabs(E[cmpnum - 1] - temp);
					temp = E[cmpnum - 1];
					count++;
					if (count % 50 == 0) {
						printf("%d,%d,%.9f,%d\n", interest_len, count, dif, fluxInd);
						if (isHermitian(HR, HI, LengthFactor) == 0) {
							printf("matrix is not hermitian");
							return 0;
						}
						/*if (difLast - dif < dif / 100 && dif >= pow(10.0, -7.0) && isMixed == 0) {
						count_dif++;,
						}
						difLast = dif;*/
					}
				}
				if (count >= IterLim - 1)
					fprintf(fpLog, "flux:%d,interest len:%d did not converg. last difference is:%.9f\n", fluxInd, interest_len, dif);
				engEvalString(ep, "Emat(:,fluxVal+1,MyLen)=E;");
				for (i = 0; i < (lat_len << 1); i++)
					finalE[i][fluxInd] = E[i];
				//free(E); E = NULL;
			}
			else {
				engEvalString(ep, "Emat(:,fluxVal+1,MyLen)=E;");
				engEvalString(ep, "state=state0;");
				//engEvalString(ep, "[E,IX]=sort(E));");
				//engEvalString(ep, "state=state0(:,IX);");
				//engEvalString(ep, "StateMat(:,:,fluxVal+1,MyLen-low_lim+1)=state0;");//save the state data into a final 4D matrix
				for (i = 0; i < (lat_len << 1); i++)
					finalE[i][fluxInd] = E0[i];
			}
			if (isCurCalc) {//Calculates the current			
				CurLen = interest_len;
				CurCalcV4(fluxInd, CurLon, CurTrn, interest_len, ep, 1);
				//CurCalcV3(fluxInd, CurLon, CurTrn, interest_len, ep);
				//CurCalcV1(fluxInd, CurLon, CurTrn, CurLen, interest_len, state_origR, state_origI);
			}
			if (isChargeDensityCalc)
				CurCalcV4(fluxInd, CurLon, CurTrn, interest_len, ep, 2);
			//mxDestroyArray(mexle n); mexlen = NULL;
			if (isEEcalc) {
				EE_calc(stateI, stateR, interest_len, ep, 1);//Vertical EE
				EE_calc(stateI, stateR, interest_len, ep, 2);//Horizontal EE
			}
			if (isCondCalc && interest_len>1) {
				ConductivityCalc(ep);
			}
		}
		//free(E0); E0 = NULL;
		//free(stateR); stateR = NULL; free(stateI); stateI = NULL;
	}

	engEvalString(ep, "save(StrE,'Emat');");
	if (isCondCalc==0 || max_fill==1)
		engEvalString(ep, "save(StrCur,'CurLonUp','CurLonDn','CurTrnUp','Cur_OprToPer','ChargeLonUp','ChargeLonDn','ChargeTrnUp','Charge_OprToPer');");
	else
		engEvalString(ep, "save(StrCur,'CurLonUp','CurLonDn','CurTrnUp','Cur_OprToPer','ChargeLonUp','ChargeLonDn','ChargeTrnUp','Charge_OprToPer','condUp','condDn');");
	engEvalString(ep, "save(StrState,'EEH','EEV');");

	//if (isCurCalc) {//Calculates the current		
	//	CurUpdate(CurLon, CurTrn,ep);
	//}
	fclose(fpLog);
}

void changing_SOC(double* HR[], double* HI[], double* HR0[], double* HI0[], double* finalE[], double* CurLon[], double* CurTrn[], int cmpnum, double sigma, Engine *ep, FILE* fpLog) {
	int i, j, count, CurLen, SOCInd, interest_len, MatSize = 2 * lat_len, LengthFactor = 2 * lat_len;
	double dif, temp, matmul;
	mxArray *MexSOC = NULL, *mexlen = NULL, *MexFlux = NULL, *mexLatLen = NULL;
	double *E, *E0, *stateR, *stateI;//*state_origR, *state_origI,
	for (SOCInd = 0; SOCInd < SOCDiv; SOCInd++) {
		system("cls");
		E = (double*)calloc(MatSize, sizeof(double));
		E0 = (double*)calloc(MatSize, sizeof(double));
		stateR = (double*)calloc(MatSize*MatSize, sizeof(double));
		stateI = (double*)calloc(MatSize*MatSize, sizeof(double));
		//state_origR = (double*)malloc(MatSize*MatSize * sizeof(double));
		//state_origI = (double*)malloc(MatSize*MatSize * sizeof(double));
		SOC_amp = (double)SOCInd*SOCMax / SOCDiv;
		init_array(HR0, MatSize, MatSize);
		init_array(HI0, MatSize, MatSize);
		FindInitStateMag(sigma, HR0, HI0, &E0, ep, LengthFactor);//, interest_len);
		printf("%d\n", SOCInd);
		mexLatLen = mxCreateDoubleScalar((double)lat_len);
		engPutVariable(ep, "lat_len", mexLatLen);
		MexSOC = mxCreateDoubleScalar((double)SOCInd);
		engPutVariable(ep, "SOCVal", MexSOC);
		MexFlux = mxCreateDoubleScalar((double)fluxVal);
		engPutVariable(ep, "fluxVal", MexFlux);// this is the flux index
		MakeInitState(&stateR, &stateI, ep, LengthFactor);
		//if (isFock || isHartr) {
		for (interest_len = min_fill; interest_len <= max_fill; interest_len++)
		{
			if (SOCInd % 10 == 0) {
				printf("%d,%d\n", interest_len, SOCInd);
			}
			temp = (double)interest_len;
			mexlen = mxCreateDoubleScalar(temp);
			//mexlen= mxCreateDoubleMatrix(1, 1, mxREAL);
			//*mxGetPr(mexlen) = temp;
			engPutVariable(ep, "MyLen", mexlen);
			//mxDestroyArray(mexlen);
			dif = 1.0;
			count = 0;
			//MakeInitState(interest_len, ep);
			//MakeInitState(interest_len, &state_origR, &state_origI, &state, ep);
			if (isHubbard || isFock || isHartr) {
				init_array(HR, MatSize, MatSize);
				init_array(HI, MatSize, MatSize);
				temp = E0[cmpnum - 1];
				while (dif > converg_lim && count < IterLim)
				{
					converg_func(HR, HI, HR0, HI0, interest_len, &E, &stateR, &stateI, ep, LengthFactor);
					dif = fabs(E[cmpnum - 1] - temp);
					temp = E[cmpnum - 1];
					count++;
					if (count % 10 == 0)
						printf("%d,%d,%f,%d\n", interest_len, count, dif, SOCInd);
				}
				if (count >= IterLim - 1)
					fprintf(fpLog, "flux:%d,interest len:%d did not converg. last difference is:%.9f\n", SOCInd, interest_len, dif);
				engEvalString(ep, "Emat(:,SOCVal+1,MyLen)=E;");
				for (i = 0; i < 2 * lat_len; i++) {
					finalE[i][SOCInd] = E[i];
				}
			}
			else {
				engEvalString(ep, "Emat(:,SOCVal+1,MyLen)=E;");
				//engEvalString(ep, "Emat(:,SOCVal+1)=E;");
				for (i = 0; i < 2 * lat_len; i++) {
					finalE[i][SOCInd] = E0[i];
				}
			}
			if (isCurCalc) {//Calculates the current			
				CurLen = interest_len;
				CurCalcV4(SOCInd, CurLon, CurTrn, interest_len, ep, 1);
				//CurCalcV3(fluxInd, CurLon, CurTrn, interest_len, ep);
				//CurCalcV1(fluxInd, CurLon, CurTrn, CurLen, interest_len, state_origR, state_origI);
			}
			if (isChargeDensityCalc)
				CurCalcV4(SOCInd, CurLon, CurTrn, interest_len, ep, 2);
			//mxDestroyArray(mexle n); mexlen = NULL;
			if (isEEcalc) {
				EE_calc(stateI, stateR, interest_len, ep, 1);//Vertical EE
				EE_calc(stateI, stateR, interest_len, ep, 2);//Horizontal EE
			}
		}
		free(E0); E0 = NULL;
		free(stateR); stateR = NULL; free(stateI); stateI = NULL;
	}
	engEvalString(ep, "save(StrE,'Emat');");
	engEvalString(ep, "save(StrCur,'CurLonUp','CurLonDn','CurTrnUp','Cur_OprToPer','ChargeLonUp','ChargeLonDn','ChargeTrnUp','Charge_OprToPer');");
	engEvalString(ep, "save(StrState,'EEH','EEV');");
}

void printCurSum(Engine *ep, FILE *fpCurSum) {
	int len = max_fill - min_fill + 1;
	double *CurSum = malloc(fluxDiv * len * sizeof(double));
	int i, j, temp;
	//FILE *fpCurSum= _fsopen(FileStrCurSum, "wt", OF_SHARE_DENY_NONE);
	mxArray *mexCurSum = engGetVariable(ep, "Cur_OprToPer");
	CurSum = mxGetPr(mexCurSum);
	for (i = 0; i < len; i++) {
		for (j = 0; j < fluxDiv; j++) {
			temp = i * fluxDiv + j;
			fprintf(fpCurSum, "%.9f       ", CurSum[temp]);

		}
		fprintf(fpCurSum, "\n");
	}
	//fclose(fpCurSum);
}

void MyPrints(int rel_max, double* finalE[], double* CurLon[], double* CurTrn[]) {
	int i, j;
	FILE *fpEn;
	fpEn = _fsopen(FileStrE, "wt", OF_SHARE_DENY_NONE);
	for (i = 0; i<2 * lat_len; i++) {
		for (j = 0; j<rel_max; j++) {
			fprintf(fpEn, "%.9f      ", finalE[i][j]);
		}
		fprintf(fpEn, "\n");
	}

	if (isCurPrint) {//printf current values Calculates the currents
		FILE *fpCurLon = _fsopen(FileStrCurLon, "wt", OF_SHARE_DENY_NONE);
		for (i = 0; i<2 * lat_len - 2; i++) {
			for (j = 0; j<rel_max; j++) {
				fprintf(fpCurLon, "%.9f      ", CurLon[i][j]);
			}
			fprintf(fpCurLon, "\n");
		}
		fclose(fpCurLon);

		FILE *fpCurTrn = _fsopen(FileStrCurTrn, "wt", OF_SHARE_DENY_NONE);
		for (i = 0; i<lat_len; i++) {
			for (j = 0; j<rel_max; j++) {
				fprintf(fpCurTrn, "%.9f      ", CurTrn[i][j]);
			}
			fprintf(fpCurTrn, "\n");
		}
		fclose(fpCurTrn);
	}
	fclose(fpEn);
}

void PrintLog(FILE** fpLog) {
	char str1[35] = "Log_n";
	char StrFmt[50] = ".txt";
	char strN[50];
	char strFlux[50] = "_flux";
	char StrInter[50];
	char FinalStr[170];
	sprintf_s(strN, 20, "%d", lat_len);
	strcat_s(str1, 20, strN);
	if (isMagnet)
		strcat_s(str1, 20, strFlux);
	if (isHubbard) {
		sprintf_s(StrInter, 20, "%s", "WithHubb");
		strcat_s(str1, 35, StrInter);
	}
	else {
		sprintf_s(StrInter, 20, "%s", "noHubb");
		strcat_s(str1, 35, StrInter);
	}
	strcat_s(str1, 35, StrFmt);
	sprintf_s(FinalStr, 170, PathStr);
	strcat_s(FinalStr, 170, str1);
	*fpLog = _fsopen(FinalStr, "wt", OF_SHARE_DENY_NONE);

	fprintf(*fpLog, "no. of sites on each leg:%d\n", lat_len);
	fprintf(*fpLog, "filling values:%d to %d\n", min_fill, max_fill);
	fprintf(*fpLog, "flux:%d\n", isMagnet);
	if (isRing)
		fprintf(*fpLog, "this is a ring model\n", ImpureAmp);
	if (FluxCase == 5) {
		fprintf(*fpLog, "impurity:%f\n", ImpureAmp);
		//printf("%f\n", isSymmetric);
		fprintf(*fpLog, "Impurity Symmetric:%f\n", isSymmetric);
	}
	else
		fprintf(*fpLog, "no impurity\n");
	fprintf(*fpLog, "hopping:%f\n", t_amp);
	if (isFlip)
		fprintf(*fpLog, "flip:%f\n", t_flip);
	else
		fprintf(*fpLog, "no flip\n");
	if (isHubbard)
		fprintf(*fpLog, "Hubbard interaction:%f\n", U_Hubb);
	else
		fprintf(*fpLog, "no Hubbard interaction\n");
	if (isHartr && isFock)
		fprintf(*fpLog, "Longitudinal interaction(Hartree & Fock):%f\n", U_HF);
	else if (isHartr)
		fprintf(*fpLog, "Longitudinal interaction(only Hartree term):%f\n", U_HF);
	else if (isFock)
		fprintf(*fpLog, "Longitudinal interaction(only Fock term):%f\n", U_HF);
	else
		fprintf(*fpLog, "no Longitudinal interaction\n");
	fprintf(*fpLog, "number of flux points:%d\n", fluxDiv);
	fprintf(*fpLog, "flux value on first leg:%f\n flux value on second leg:%f\n", f1, f2);
	fprintf(*fpLog, "value of alpha is:%f\n", alpha);
	fprintf(*fpLog, "on site amplitude is:%f\n", lambda);
}

void CreateIter(Engine *ep, double sigma) {
	FILE *fpLog = NULL;
	//init_array(HI0, MatSize, MatSize);
	int i, MatSize = 2 * lat_len, rel_max, MaxFlux, cmpnum = lat_len;
	double **HR0 = (double**)malloc((lat_len << 1) * sizeof(double*));//the interaction-free part of the hamiltonian
	double **HR = (double**)malloc((lat_len << 1) * sizeof(double*));//the full hamiltonian
	double **HI0 = (double**)malloc((lat_len << 1) * sizeof(double*));//the interation free part of the hamiltonian
	double **HI = (double**)malloc((lat_len << 1) * sizeof(double*));//the full hamiltonian
	double **finalE = (double**)malloc((lat_len << 1) * sizeof(double*));
	double **CurLon = (double**)malloc((lat_len << 1) * sizeof(double*));
	double **CurTrn = (double**)malloc(lat_len * sizeof(double*));
	if (isConstSOC == 1 || isSOC == 0)
		rel_max = fluxDiv;
	else
		rel_max = SOCDiv;

	PrintLog(&fpLog);

	for (i = 0; i<(lat_len << 1); i++)
	{
		HR0[i] = (double*)malloc((lat_len << 1) * sizeof(double));
		HR[i] = (double*)malloc((lat_len << 1) * sizeof(double));
		HI0[i] = (double*)malloc((lat_len << 1) * sizeof(double));
		HI[i] = (double*)malloc((lat_len << 1) * sizeof(double));
		finalE[i] = (double*)malloc(rel_max * sizeof(double));
		CurLon[i] = (double*)malloc(rel_max * sizeof(double));
	}
	for (i = 0; i<lat_len; i++)
		CurTrn[i] = (double*)malloc(rel_max * sizeof(double));

	init_array(finalE, MatSize, rel_max);
	init_array(CurTrn, lat_len, rel_max);
	init_array(CurLon, MatSize, rel_max);
	if (isMagnet)
		if (Full2Pi)
			MaxFlux = (int)fluxDiv;
		else
			MaxFlux = (int)fluxDiv / 2;
	else {
		MaxFlux = 1;
	}

	if (isConstSOC == 1 || isSOC == 0) {
		changing_flux(HR, HI, HR0, HI0, finalE, CurLon, CurTrn, MaxFlux, cmpnum, sigma, ep, fpLog);
	}
	else {
		fluxVal = 0.0;
		changing_SOC(HR, HI, HR0, HI0, finalE, CurLon, CurTrn, cmpnum, sigma, ep, fpLog);
	}


	for (i = 0; i<2 * lat_len; i++)
	{
		free(HR0[i]);
		free(HR[i]);
		free(HI0[i]);
		free(HI[i]);
		free(finalE[i]);
		free(CurLon[i]);
		//free(CurTrn[i]);
	}
	for (i = 0; i<lat_len; i++)
		free(CurTrn[i]);
	free(HR);
	free(HR0);
	free(HI);
	free(HI0);
	free(CurLon);	
	free(finalE);
}

void ImpurityIter(Engine *ep, double sigma) {
	int i,MaxInd;
	char FldStr[170], ImpurityStr[20] = "impurity_", SlashStr[3] = "\\", ValStr[20];;
	char FinalStrE[170], FinalStrCur[170], FinalStrState[170];
	sprintf_s(FldStr, 170, PathStr);
	if (FluxCase == 5 || FluxCase==9)
		MaxInd = ImpurityNum + 1;
	else
		MaxInd = 1;
	for (i = 0; i < MaxInd; i++) {
		ImpureAmp = i * max_impurity / ImpurityNum;
		sprintf_s(ValStr, 20, "%.3f", ImpureAmp);
		strcat_s(ImpurityStr, 20, ValStr);
		//strcat_s(PathStr, 170, SlashStr);
		if (MaxInd > 1) {
			strcat_s(PathStr, 170, ImpurityStr);
			strcat_s(PathStr, 170, SlashStr);
			_mkdir(PathStr);
		}
		sprintf_s(FinalStrE, 170, PathStr);
		strcat_s(FinalStrE, 170, FileStrE);
		mxArray *mexStrE = mxCreateString(FinalStrE);
		engPutVariable(ep, "StrE", mexStrE);
		engEvalString(ep, "lastDot=find(StrE=='.',1,'last');StrE=strcat(StrE(1:lastDot-1),'.mat')");

		sprintf_s(FinalStrCur, 170, PathStr);
		strcat_s(FinalStrCur, 170, FileStrCur);
		mxArray *mexStrCur = mxCreateString(FinalStrCur);
		engPutVariable(ep, "StrCur", mexStrCur);
		engEvalString(ep, "lastDot=find(StrCur=='.',1,'last');StrCur=strcat(StrCur(1:lastDot-1),'.mat')");

		sprintf_s(FinalStrState, 170, PathStr);
		strcat_s(FinalStrState, 170, FileStrState);
		mxArray *mexStrState = mxCreateString(FinalStrState);
		engPutVariable(ep, "StrState", mexStrState);
		engEvalString(ep, "lastDot=find(StrState=='.',1,'last');StrState=strcat(StrState(1:lastDot-1),'.mat')");
		CreateIter(ep,sigma);

		sprintf_s(PathStr, 170, FldStr);
		sprintf_s(ImpurityStr, 20, "impurity_");
		engEvalString(ep, "clear all;");
	}
}

int main()
{
	double sigma = sqrt(90.5);//modualtion;
	//char buffer[BSZ  + 1];

	Engine *ep;//Open matlab engine session
	if (!(ep = engOpen(NULL)))
	{
		printf("cant open engine\n");
		exit(-1);
	}
	engEvalString(ep, "clear all;");

	alpha = (double)p_alpha / q_alpha;

	ImpurityIter(ep, sigma);

	printf("end\n");
	//	engOutputBuffer(ep,buffer,BSZ);
	//	engEvalString(ep, "whos");
	//	MessageBox ((HWND)NULL, (LPSTR)buffer, (LPSTR) "MATLAB - whos", MB_OK);
	//	printf("%s\n",buffer);
	//engClose(ep);
	return 0;
}
