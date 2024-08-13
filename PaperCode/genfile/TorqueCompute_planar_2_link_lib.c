/*
% TorqueCompute_planar_2_link_lib.c comput torque
*/

#include <math.h>

void TorqueCompute_planar_2_link_lib(const double *Q, const double *Qd, const double *Qdd, const double *Theta, int N, double *tau)
{
	int pos;
	char k;
	int startInd;
	double gi[10];
	double Y_2_d1_p1, Y_2_d2_p1, Y_2_d2_p2, Y_2_d2_p3, Y_2_d3_p1, Y_2_d4_p1, Y_2_d4_p2, Y_2_d5_p1, Y_2_d6_p1, Y_2_d7_p1, Y_6_d5_p1, Y_p1, Y_p2, Y_p3, Y_p4, Y_p5, Y_p7, gi2_p1, gi2_p2;


	

	for (pos = 0; pos < N; pos++){

		/* Define start index */
		startInd = pos * 2;

		for (k = 0; k < 2; k++){
			gi[k] = Q[startInd+k];
			gi[k+2] = sin(gi[k]);
			gi[k+4] = cos(gi[k]);
			gi[k+6] = Qd[startInd+k];
			gi[k+8] = Qdd[startInd+k];
		}

		gi2_p1 = gi[6]*gi[6];
		gi2_p2 = gi[7]*gi[7];

		Y_p1 = gi[9];
		Y_2_d1_p1 = 384.9444*gi[2];
		Y_2_d2_p1 = -2*gi[3];
		Y_2_d2_p2 = -gi[3];
		Y_2_d2_p3 = Y_2_d1_p1*gi[3];
		Y_2_d3_p1 = -384.9444*gi[4];
		Y_2_d4_p1 = 2*gi[5];
		Y_2_d4_p2 = Y_2_d3_p1*gi[5]+Y_2_d2_p3;
		Y_2_d5_p1 = Y_2_d2_p1*gi[6];
		Y_2_d6_p1 = Y_2_d5_p1*gi[7]+Y_2_d2_p2*gi2_p2+Y_2_d4_p2;
		Y_2_d7_p1 = Y_2_d4_p1*gi[8]+Y_2_d6_p1;
		Y_p3 = gi[5]*gi[9]+Y_2_d7_p1;
		Y_p5 = gi[8];
		Y_p7 = 9.81*gi[4];

		Y_p2 = gi[9]+gi[8];

		Y_6_d5_p1 = gi[3]*gi2_p1+Y_2_d4_p2;
		Y_p4 = gi[5]*gi[8]+Y_6_d5_p1;
		tau[startInd+0] = Y_p1*Theta[0]+Y_p3*Theta[1]+Y_p5*Theta[2]+Y_p7*Theta[3];
		tau[startInd+1] = Y_p2*Theta[0]+Y_p4*Theta[1];

	}
}