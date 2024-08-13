/*
% _FUNCTIONNAME_.c
*/

#include <math.h>

void _FUNCTIONNAME_(const double *Q, const double *Qd, const double *Qdd, const double *Theta, int N, double *tau, double *Y)
{
	int pos;
	char k;
	int startInd, startIndY;
	double gi[_5DOF_];
	/*NAME_DEF*/

	/*SETUP1_CODE*/

	for (pos = 0; pos < N; pos++){

		/* Define start index */
		startInd = pos * _1DOF_;
		startIndY = pos * _1DOF_ * _ell_;

		for (k = 0; k < _1DOF_; k++){
			gi[k] = Q[startInd+k];
			gi[k+_1DOF_] = sin(gi[k]);
			gi[k+_2DOF_] = cos(gi[k]);
			gi[k+_3DOF_] = Qd[startInd+k];
			gi[k+_4DOF_] = Qdd[startInd+k];
		}

/*SETUP2_CODE*/

		/*TAU_CODE*/

	}
}
