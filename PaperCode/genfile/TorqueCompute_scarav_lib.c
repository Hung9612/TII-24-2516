/*
% TorqueCompute_scarav_lib.c comput torque
*/

#include <math.h>

void TorqueCompute_scarav_lib(const double *Q, const double *Qd, const double *Qdd, const double *Theta, int N, double *tau)
{
	int pos;
	char k;
	int startInd;
	double gi[20];
	double Y_10_d3_p1, Y_10_d4_p1, Y_10_d4_p2, Y_10_d5_p1, Y_10_d5_p2, Y_10_d6_p1, Y_10_d7_p1, Y_10_d8_p1, Y_10_d9_p1, Y_13_d3_p1, Y_15_d3_p1, Y_15_d4_p1, Y_15_d5_p1, Y_15_d6_p1, Y_15_d7_p1, Y_16_d5_p1, Y_16_d6_p1, Y_16_d7_p1, Y_1_d1_p1, Y_1_d2_p1, Y_1_d3_p1, Y_1_d4_p1, Y_1_d5_p1, Y_2_d2_p1, Y_2_d2_p2, Y_2_d2_p3, Y_2_d2_p4, Y_2_d2_p5, Y_2_d3_p1, Y_2_d3_p2, Y_2_d3_p3, Y_2_d3_p4, Y_2_d3_p5, Y_2_d3_p6, Y_2_d3_p7, Y_2_d4_p1, Y_2_d4_p2, Y_2_d4_p3, Y_2_d4_p4, Y_2_d4_p5, Y_2_d4_p6, Y_2_d4_p7, Y_2_d5_p1, Y_2_d5_p2, Y_2_d6_p1, Y_2_d6_p2, Y_2_d7_p1, Y_2_d8_p1, Y_2_d9_p1, Y_3_d1_p1, Y_3_d2_p1, Y_3_d3_p1, Y_3_d3_p2, Y_3_d3_p3, Y_3_d4_p1, Y_3_d4_p2, Y_3_d4_p3, Y_3_d4_p4, Y_3_d4_p5, Y_3_d5_p1, Y_3_d5_p2, Y_3_d6_p1, Y_3_d6_p2, Y_3_d7_p1, Y_3_d8_p1, Y_3_d9_p1, Y_6_d1_p1, Y_6_d2_p1, Y_6_d3_p1, Y_6_d4_p1, Y_6_d5_p1, Y_8_d3_p1, Y_9_d3_p1, Y_9_d4_p1, Y_9_d4_p2, Y_9_d4_p3, Y_9_d4_p4, Y_9_d5_p1, Y_9_d5_p2, Y_9_d6_p1, Y_9_d7_p1, Y_9_d8_p1, Y_9_d9_p1, Y_p1, Y_p10, Y_p12, Y_p13, Y_p14, Y_p16, Y_p17, Y_p18, Y_p2, Y_p21, Y_p22, Y_p25, Y_p31, Y_p5, Y_p6, Y_p8, Y_p9, gi2_p1, gi2_p2, gi2_p3;


	

	for (pos = 0; pos < N; pos++){

		/* Define start index */
		startInd = pos * 4;

		for (k = 0; k < 4; k++){
			gi[k] = Q[startInd+k];
			gi[k+4] = sin(gi[k]);
			gi[k+8] = cos(gi[k]);
			gi[k+12] = Qd[startInd+k];
			gi[k+16] = Qdd[startInd+k];
		}

		gi2_p1 = gi[12]*gi[12];
		gi2_p2 = gi[13]*gi[13];
		gi2_p3 = gi[15]*gi[15];

				Y_1_d1_p1 = 0.5*gi[5];
		Y_1_d2_p1 = 0.5*gi[9];
		Y_1_d3_p1 = gi[9]*gi[12];
		Y_1_d4_p1 = Y_1_d3_p1*gi[13]+Y_1_d2_p1*gi2_p2;
		Y_1_d5_p1 = gi[5]*gi[16]+Y_1_d4_p1;
		Y_p1 = Y_1_d5_p1+Y_1_d1_p1*gi[17];

		Y_2_d2_p1 = gi[5]*gi[7];
		Y_2_d2_p2 = Y_1_d1_p1*gi[7];
		Y_2_d2_p3 = 0.577777981718733*gi[7];
		Y_2_d2_p4 = 0.5*gi[7];
		Y_2_d2_p5 = 0.288888990880427*gi[7];
		Y_2_d3_p1 = -0.288888990880427-0.5*gi[9];
		Y_2_d3_p2 = -0.5*gi[9]-0.577777981718733;
		Y_2_d3_p3 = -0.577777981718733-gi[9];
		Y_2_d3_p4 = gi[7]*gi[9]+Y_2_d2_p3;
		Y_2_d3_p5 = gi[7]*gi[9];
		Y_2_d3_p6 = Y_2_d2_p4*gi[9]+Y_2_d2_p5;
		Y_2_d3_p7 = Y_2_d2_p4*gi[9];
		Y_2_d4_p1 = Y_2_d2_p1+Y_2_d3_p3*gi[11];
		Y_2_d4_p2 = Y_2_d2_p2+Y_2_d3_p1*gi[11];
		Y_2_d4_p3 = Y_2_d2_p2+Y_2_d3_p2*gi[11];
		Y_2_d4_p4 = Y_2_d3_p4+gi[5]*gi[11];
		Y_2_d4_p5 = Y_2_d3_p5+gi[5]*gi[11];
		Y_2_d4_p6 = Y_2_d3_p6+Y_1_d1_p1*gi[11];
		Y_2_d4_p7 = Y_2_d3_p7+Y_1_d1_p1*gi[11];
		Y_2_d5_p1 = Y_2_d4_p4*gi[12];
		Y_2_d5_p2 = Y_2_d4_p5*gi[12];
		Y_2_d6_p1 = Y_2_d4_p4*gi[13]+Y_2_d5_p1;
		Y_2_d6_p2 = Y_2_d5_p2*gi[13]+Y_2_d4_p7*gi2_p2;
		Y_2_d7_p1 = Y_2_d6_p1*gi[15]+Y_2_d6_p2+Y_2_d4_p6*gi2_p3;
		Y_2_d8_p1 = Y_2_d4_p1*gi[16]+Y_2_d7_p1;
		Y_2_d9_p1 = Y_2_d8_p1+Y_2_d4_p3*gi[17];
		Y_p5 = Y_2_d9_p1+Y_2_d4_p2*gi[19];
		Y_3_d1_p1 = -gi[5];
		Y_3_d2_p1 = Y_3_d1_p1*gi[7];
		Y_3_d3_p1 = Y_2_d2_p3+Y_2_d2_p4*gi[9];
		Y_3_d3_p2 = gi[9]+0.577777981718733;
		Y_3_d3_p3 = 0.5*gi[9]+0.288888990880427;
		Y_3_d4_p1 = Y_3_d3_p1+Y_1_d1_p1*gi[11];
		Y_3_d4_p2 = Y_3_d3_p2*gi[11]+Y_3_d2_p1;
		Y_3_d4_p3 = gi[9]*gi[11]+Y_3_d2_p1;
		Y_3_d4_p4 = Y_3_d3_p3*gi[11]-Y_2_d2_p2;
		Y_3_d4_p5 = Y_1_d2_p1*gi[11]-Y_2_d2_p2;
		Y_3_d5_p1 = Y_3_d4_p2*gi[12];
		Y_3_d5_p2 = Y_3_d4_p3*gi[12];
		Y_3_d6_p1 = Y_3_d4_p2*gi[13]+Y_3_d5_p1;
		Y_3_d6_p2 = Y_3_d5_p2*gi[13]+Y_3_d4_p5*gi2_p2;
		Y_3_d7_p1 = Y_3_d6_p1*gi[15]+Y_3_d6_p2+Y_3_d4_p4*gi2_p3;
		Y_3_d8_p1 = Y_2_d4_p4*gi[16]+Y_3_d7_p1;
		Y_3_d9_p1 = Y_3_d8_p1+Y_3_d4_p1*gi[17];
		Y_p9 = Y_3_d9_p1+Y_2_d4_p6*gi[19];
		Y_p13 = gi[19];
		Y_p17 = gi[17];
		Y_6_d1_p1 = -2*gi[5];
		Y_6_d2_p1 = 2*gi[9];
		Y_6_d3_p1 = Y_6_d1_p1*gi[12];
		Y_6_d4_p1 = Y_6_d3_p1*gi[13]+Y_3_d1_p1*gi2_p2;
		Y_6_d5_p1 = Y_6_d2_p1*gi[16]+Y_6_d4_p1;
		Y_p21 = gi[9]*gi[17]+Y_6_d5_p1;
		Y_p25 = gi[16];

		Y_8_d3_p1 = -Y_1_d2_p1*gi2_p1;
		Y_p2 = Y_1_d1_p1*gi[16]+Y_8_d3_p1;

		Y_9_d3_p1 = -0.577777981718733-0.5*gi[9];
		Y_9_d4_p1 = Y_2_d2_p2+Y_9_d3_p1*gi[11];
		Y_9_d4_p2 = -0.288888990880427*gi[11];
		Y_9_d4_p3 = -0.577777981718733*gi[11];
		Y_9_d4_p4 = -Y_2_d3_p7-Y_1_d1_p1*gi[11];
		Y_9_d5_p1 = Y_2_d2_p3*gi[12];
		Y_9_d5_p2 = Y_9_d4_p4*gi2_p1;
		Y_9_d6_p1 = Y_2_d2_p3*gi[13]+Y_9_d5_p1;
		Y_9_d7_p1 = Y_9_d6_p1*gi[15]+Y_2_d2_p5*gi2_p3+Y_9_d5_p2;
		Y_9_d8_p1 = Y_9_d4_p1*gi[16]+Y_9_d7_p1;
		Y_9_d9_p1 = Y_9_d8_p1+Y_9_d4_p3*gi[17];
		Y_p6 = Y_9_d9_p1+Y_9_d4_p2*gi[19];

		Y_10_d3_p1 = Y_2_d2_p4*gi[9]+Y_2_d2_p3;
		Y_10_d4_p1 = Y_10_d3_p1+Y_1_d1_p1*gi[11];
		Y_10_d4_p2 = Y_2_d2_p2-Y_1_d2_p1*gi[11];
		Y_10_d5_p1 = -Y_9_d4_p3*gi[12];
		Y_10_d5_p2 = Y_10_d4_p2*gi2_p1;
		Y_10_d6_p1 = -Y_9_d4_p3*gi[13]+Y_10_d5_p1;
		Y_10_d7_p1 = Y_10_d6_p1*gi[15]-Y_9_d4_p2*gi2_p3+Y_10_d5_p2;
		Y_10_d8_p1 = Y_10_d4_p1*gi[16]+Y_10_d7_p1;
		Y_10_d9_p1 = Y_10_d8_p1+Y_2_d2_p3*gi[17];
		Y_p10 = Y_10_d9_p1+Y_2_d2_p5*gi[19];
		Y_p14 = Y_p13;

		Y_p18 = gi[17]+gi[16];

		Y_13_d3_p1 = gi[5]*gi2_p1;
		Y_p22 = gi[9]*gi[16]+Y_13_d3_p1;
		Y_p31 = gi[18]+96.2361041164153;

		Y_15_d3_p1 = -Y_2_d2_p5-Y_2_d2_p4*gi[9];
		Y_15_d4_p1 = Y_15_d3_p1-Y_1_d1_p1*gi[11];
		Y_15_d5_p1 = Y_15_d4_p1*gi2_p1;
		Y_15_d6_p1 = -Y_2_d2_p5*gi2_p2+Y_15_d5_p1-Y_9_d5_p1*gi[13];
		Y_15_d7_p1 = Y_2_d4_p2*gi[16]+Y_15_d6_p1;
		Y_p8 = Y_15_d7_p1+Y_9_d4_p2*gi[17];

		Y_16_d5_p1 = Y_2_d4_p2*gi2_p1;
		Y_16_d6_p1 = Y_16_d5_p1+Y_9_d4_p2*gi2_p2-Y_10_d5_p1*gi[13];
		Y_16_d7_p1 = Y_2_d4_p6*gi[16]+Y_16_d6_p1;
		Y_p12 = Y_16_d7_p1+Y_2_d2_p5*gi[17];

		Y_p16 = gi[19]+Y_p18;
		tau[startInd+0] = Y_p1*Theta[0]+Y_p5*Theta[1]+Y_p9*Theta[2]+Y_p13*Theta[3]+Y_p17*Theta[4]+Y_p21*Theta[5]+Y_p25*Theta[6];
		tau[startInd+1] = Y_p2*Theta[0]+Y_p6*Theta[1]+Y_p10*Theta[2]+Y_p14*Theta[3]+Y_p18*Theta[4]+Y_p22*Theta[5];
		tau[startInd+2] = Y_p31*Theta[7];
		tau[startInd+3] = Y_p8*Theta[1]+Y_p12*Theta[2]+Y_p16*Theta[3];

	}
}