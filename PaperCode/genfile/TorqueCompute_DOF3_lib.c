/*
% TorqueCompute_DOF3_lib.c comput torque
*/

#include <math.h>

void TorqueCompute_DOF3_lib(const double *Q, const double *Qd, const double *Qdd, const double *Theta, int N, double *tau)
{
	int pos;
	char k;
	int startInd;
	double gi[15];
	double Y_11_d10_p1, Y_11_d11_p1, Y_11_d3_p1, Y_11_d4_p1, Y_11_d5_p1, Y_11_d6_p1, Y_11_d6_p2, Y_11_d6_p3, Y_11_d7_p1, Y_11_d8_p1, Y_11_d8_p2, Y_11_d9_p1, Y_12_d3_p1, Y_15_d1_p1, Y_19_d2_p1, Y_19_d3_p1, Y_1_d1_p1, Y_1_d3_p1, Y_20_d2_p1, Y_20_d3_p1, Y_21_d3_p1, Y_21_d4_p1, Y_21_d4_p2, Y_21_d5_p1, Y_21_d5_p2, Y_21_d6_p1, Y_21_d7_p1, Y_22_d1_p1, Y_22_d2_p1, Y_22_d2_p2, Y_22_d3_p1, Y_22_d3_p2, Y_22_d4_p1, Y_22_d5_p1, Y_22_d5_p2, Y_22_d6_p1, Y_22_d7_p1, Y_23_d2_p1, Y_23_d3_p1, Y_23_d4_p1, Y_23_d5_p1, Y_23_d5_p2, Y_23_d6_p1, Y_23_d7_p1, Y_23_d8_p1, Y_24_d4_p1, Y_24_d4_p2, Y_24_d4_p3, Y_24_d5_p1, Y_24_d5_p2, Y_24_d6_p1, Y_24_d7_p1, Y_24_d8_p1, Y_25_d2_p1, Y_25_d3_p1, Y_25_d4_p1, Y_25_d4_p2, Y_25_d4_p3, Y_25_d5_p1, Y_25_d5_p2, Y_25_d6_p1, Y_26_d2_p1, Y_27_d3_p1, Y_28_d3_p1, Y_2_d10_p1, Y_2_d11_p1, Y_2_d1_p1, Y_2_d2_p1, Y_2_d3_p1, Y_2_d3_p2, Y_2_d3_p3, Y_2_d4_p1, Y_2_d5_p1, Y_2_d5_p2, Y_2_d5_p3, Y_2_d6_p1, Y_2_d6_p2, Y_2_d6_p3, Y_2_d6_p4, Y_2_d7_p1, Y_2_d8_p1, Y_2_d8_p2, Y_2_d9_p1, Y_32_d4_p1, Y_32_d5_p1, Y_32_d6_p1, Y_33_d3_p1, Y_33_d5_p1, Y_33_d6_p1, Y_33_d7_p1, Y_34_d3_p1, Y_34_d4_p1, Y_34_d5_p1, Y_34_d5_p2, Y_34_d6_p1, Y_34_d7_p1, Y_35_d5_p1, Y_35_d6_p1, Y_35_d7_p1, Y_36_d3_p1, Y_36_d3_p2, Y_36_d4_p1, Y_36_d4_p2, Y_36_d5_p1, Y_37_d4_p1, Y_37_d5_p1, Y_37_d6_p1, Y_37_d7_p1, Y_38_d4_p1, Y_38_d5_p1, Y_38_d5_p2, Y_38_d6_p1, Y_3_d2_p1, Y_3_d3_p1, Y_3_d4_p1, Y_4_d1_p1, Y_4_d2_p1, Y_4_d2_p2, Y_4_d3_p1, Y_4_d3_p2, Y_4_d4_p1, Y_4_d4_p2, Y_4_d4_p3, Y_4_d5_p1, Y_4_d5_p2, Y_4_d6_p1, Y_4_d7_p1, Y_4_d8_p1, Y_4_d9_p1, Y_5_d1_p1, Y_5_d1_p2, Y_5_d2_p1, Y_5_d2_p2, Y_5_d3_p1, Y_5_d3_p2, Y_5_d3_p3, Y_5_d3_p4, Y_5_d4_p1, Y_5_d4_p2, Y_5_d4_p3, Y_5_d5_p1, Y_5_d5_p2, Y_5_d6_p1, Y_5_d7_p1, Y_5_d8_p1, Y_5_d9_p1, Y_6_d2_p1, Y_6_d2_p2, Y_6_d3_p1, Y_6_d3_p2, Y_6_d3_p3, Y_6_d4_p1, Y_6_d4_p2, Y_6_d4_p3, Y_6_d4_p4, Y_6_d4_p5, Y_6_d5_p1, Y_6_d5_p2, Y_6_d6_p1, Y_6_d6_p2, Y_6_d7_p1, Y_6_d8_p1, Y_6_d9_p1, Y_7_d2_p1, Y_7_d2_p2, Y_7_d2_p3, Y_7_d3_p1, Y_7_d3_p2, Y_7_d3_p3, Y_7_d3_p4, Y_7_d4_p1, Y_7_d4_p2, Y_7_d4_p3, Y_7_d4_p4, Y_7_d4_p5, Y_7_d4_p6, Y_7_d5_p1, Y_7_d5_p2, Y_7_d6_p1, Y_7_d6_p2, Y_7_d7_p1, Y_7_d8_p1, Y_8_d3_p1, Y_8_d3_p2, Y_8_d3_p3, Y_8_d4_p1, Y_8_d4_p2, Y_8_d4_p3, Y_8_d5_p1, Y_8_d6_p1, Y_8_d6_p2, Y_8_d7_p1, Y_8_d8_p1, Y_9_d2_p1, Y_9_d3_p1, Y_9_d4_p1, Y_p1, Y_p10, Y_p11, Y_p12, Y_p13, Y_p14, Y_p15, Y_p16, Y_p17, Y_p18, Y_p19, Y_p2, Y_p20, Y_p21, Y_p22, Y_p23, Y_p24, Y_p25, Y_p26, Y_p28, Y_p29, Y_p31, Y_p32, Y_p33, Y_p34, Y_p35, Y_p37, Y_p4, Y_p40, Y_p43, Y_p44, Y_p46, Y_p47, Y_p49, Y_p5, Y_p6, Y_p7, Y_p8, gi2_p1, gi2_p2, gi2_p3, gi2_p4, gi2_p5, gi2_p6, gi2_p7;


	

	for (pos = 0; pos < N; pos++){

		/* Define start index */
		startInd = pos * 3;

		for (k = 0; k < 3; k++){
			gi[k] = Q[startInd+k];
			gi[k+3] = sin(gi[k]);
			gi[k+6] = cos(gi[k]);
			gi[k+9] = Qd[startInd+k];
			gi[k+12] = Qdd[startInd+k];
		}

		gi2_p1 = gi[4]*gi[4];
		gi2_p2 = gi[5]*gi[5];
		gi2_p3 = gi[7]*gi[7];
		gi2_p4 = gi[8]*gi[8];
		gi2_p5 = gi[9]*gi[9];
		gi2_p6 = gi[10]*gi[10];
		gi2_p7 = gi[11]*gi[11];

				Y_1_d1_p1 = -gi[4];

		Y_1_d3_p1 = Y_1_d1_p1*gi2_p6;
		Y_p1 = gi[7]*gi[13]+Y_1_d3_p1;
		Y_2_d1_p1 = -98.7036946115134*gi[3];
		Y_2_d2_p1 = 2*gi[4];
		Y_2_d3_p1 = gi[4]*gi[5];
		Y_2_d3_p2 = 2*gi[5];
		Y_2_d3_p3 = -98.7036946115134*gi[5];
		Y_2_d4_p1 = Y_2_d3_p3*gi[6];
		Y_2_d5_p1 = -gi[7];
		Y_2_d5_p2 = gi[5]*gi[7];
		Y_2_d5_p3 = Y_2_d4_p1*gi[7];
		Y_2_d6_p1 = Y_2_d5_p1*gi[8];
		Y_2_d6_p2 = -2*gi[8];
		Y_2_d6_p3 = Y_2_d2_p1*gi[8];
		Y_2_d6_p4 = Y_2_d5_p3+Y_2_d1_p1*gi[8];
		Y_2_d7_p1 = Y_2_d3_p2*gi[9];
		Y_2_d8_p1 = Y_2_d7_p1+Y_2_d6_p3*gi[10];
		Y_2_d8_p2 = Y_2_d5_p2*gi2_p6+Y_2_d6_p4;
		Y_2_d9_p1 = Y_2_d8_p1*gi[11]+Y_2_d5_p2*gi2_p7+Y_2_d8_p2;
		Y_2_d10_p1 = Y_2_d6_p2*gi[12]+Y_2_d9_p1;
		Y_2_d11_p1 = Y_2_d3_p1*gi[13]+Y_2_d10_p1;
		Y_p4 = Y_2_d11_p1+Y_2_d6_p1*gi[14];

		Y_3_d2_p1 = Y_2_d2_p1*gi[7];
		Y_3_d3_p1 = Y_3_d2_p1*gi[9];
		Y_3_d4_p1 = Y_3_d3_p1*gi[10];
		Y_p7 = gi2_p1*gi[12]+Y_3_d4_p1;
		Y_4_d1_p1 = 2*gi2_p1;
		Y_4_d2_p1 = Y_2_d2_p1*gi[5];
		Y_4_d2_p2 = Y_1_d1_p1*gi[5];
		Y_4_d3_p1 = Y_4_d2_p1*gi[7];
		Y_4_d3_p2 = Y_4_d1_p1-2*gi2_p3;
		Y_4_d4_p1 = Y_1_d1_p1*gi[8];
		Y_4_d4_p2 = -Y_3_d2_p1*gi[8];
		Y_4_d4_p3 = Y_4_d3_p2*gi[8];
		Y_4_d5_p1 = Y_4_d3_p1*gi[9];
		Y_4_d5_p2 = Y_4_d4_p3*gi[9];
		Y_4_d6_p1 = Y_4_d5_p2*gi[10]+Y_4_d2_p2*gi2_p6;
		Y_4_d7_p1 = Y_4_d5_p1*gi[11]+Y_4_d6_p1+Y_2_d3_p1*gi2_p7;
		Y_4_d8_p1 = Y_4_d4_p2*gi[12]+Y_4_d7_p1;
		Y_4_d9_p1 = Y_2_d5_p2*gi[13]+Y_4_d8_p1;
		Y_p10 = Y_4_d9_p1+Y_4_d4_p1*gi[14];
		Y_5_d1_p1 = 0.5*gi[4];
		Y_5_d1_p2 = -gi2_p1;
		Y_5_d2_p1 = Y_5_d1_p1*gi[5];
		Y_5_d2_p2 = Y_5_d1_p2*gi[5];
		Y_5_d3_p1 = Y_2_d3_p1*gi[7];
		Y_5_d3_p2 = 0.5*gi[7];
		Y_5_d3_p3 = gi[5]*gi2_p3+Y_5_d2_p2;
		Y_5_d3_p4 = gi[4]*gi[7];
		Y_5_d4_p1 = Y_5_d3_p2*gi[8];
		Y_5_d4_p2 = Y_5_d3_p4*gi[8];
		Y_5_d4_p3 = Y_5_d1_p1*gi[8];
		Y_5_d5_p1 = Y_5_d3_p3*gi[9];
		Y_5_d5_p2 = Y_5_d4_p2*gi[9];
		Y_5_d6_p1 = Y_5_d5_p1*gi[10]-Y_5_d4_p3*gi2_p6;
		Y_5_d7_p1 = Y_5_d6_p1+Y_5_d5_p2*gi[11]+Y_5_d4_p3*gi2_p7;
		Y_5_d8_p1 = Y_5_d3_p1*gi[12]+Y_5_d7_p1;
		Y_5_d9_p1 = Y_5_d8_p1+Y_5_d4_p1*gi[13];
		Y_p13 = Y_5_d9_p1+Y_5_d2_p1*gi[14];

		Y_6_d2_p1 = -gi2_p2;
		Y_6_d2_p2 = Y_2_d2_p1*gi2_p2;
		Y_6_d3_p1 = Y_6_d2_p1*gi2_p3;
		Y_6_d3_p2 = Y_6_d2_p2*gi[7];
		Y_6_d3_p3 = Y_2_d3_p2-Y_2_d3_p2*gi2_p3;
		Y_6_d4_p1 = Y_2_d3_p1*gi[8];
		Y_6_d4_p2 = -gi2_p4+Y_6_d3_p1;
		Y_6_d4_p3 = Y_2_d2_p1*gi2_p4;
		Y_6_d4_p4 = Y_6_d3_p3*gi[8];
		Y_6_d4_p5 = Y_2_d5_p2*gi[8];
		Y_6_d5_p1 = Y_6_d3_p2*gi[9];
		Y_6_d5_p2 = Y_6_d4_p4*gi[9];
		Y_6_d6_p1 = Y_6_d4_p3*gi[10]+Y_6_d5_p2;
		Y_6_d6_p2 = Y_6_d5_p1*gi[10]+Y_6_d4_p5*gi2_p6;
		Y_6_d7_p1 = Y_6_d6_p1*gi[11]+Y_6_d6_p2;
		Y_6_d8_p1 = Y_6_d4_p2*gi[12]+Y_6_d7_p1;
		Y_6_d9_p1 = Y_6_d4_p1*gi[13]+Y_6_d8_p1;
		Y_p16 = Y_6_d9_p1+Y_2_d5_p1*gi[14];

		Y_7_d2_p1 = Y_5_d1_p1*gi2_p2;
		Y_7_d2_p2 = -gi[5];
		Y_7_d2_p3 = 0.5*gi2_p2;
		Y_7_d3_p1 = gi[5]*gi2_p3+Y_7_d2_p2;
		Y_7_d3_p2 = gi2_p2+Y_6_d2_p1*gi2_p3;
		Y_7_d3_p3 = gi2_p3-1;
		Y_7_d3_p4 = Y_7_d2_p3*gi[7];
		Y_7_d4_p1 = Y_7_d3_p1*gi[8];
		Y_7_d4_p2 = Y_7_d2_p1-Y_5_d1_p1*gi2_p4;
		Y_7_d4_p3 = Y_4_d2_p1*gi[8];
		Y_7_d4_p4 = Y_7_d3_p2+Y_7_d3_p3*gi2_p4;
		Y_7_d4_p5 = Y_7_d3_p4-Y_5_d3_p2*gi2_p4;
		Y_7_d4_p6 = -Y_4_d3_p1*gi[8];
		Y_7_d5_p1 = Y_7_d4_p4*gi[9];
		Y_7_d5_p2 = Y_7_d4_p6*gi[9];
		Y_7_d6_p1 = Y_7_d4_p3*gi[10]+Y_7_d5_p1;
		Y_7_d6_p2 = Y_7_d4_p5*gi2_p6+Y_7_d5_p2*gi[10];
		Y_7_d7_p1 = Y_7_d6_p1*gi[11]+Y_7_d6_p2;
		Y_7_d8_p1 = Y_7_d4_p1*gi[12]+Y_7_d7_p1;
		Y_p19 = Y_7_d8_p1+Y_7_d4_p2*gi[13];

		Y_8_d3_p1 = gi2_p2+gi2_p2*gi2_p3;
		Y_8_d3_p2 = gi2_p3+1;
		Y_8_d3_p3 = 2*gi[7];
		Y_8_d4_p1 = Y_8_d3_p1+Y_8_d3_p2*gi2_p4;
		Y_8_d4_p2 = -Y_2_d2_p1*gi2_p4-Y_6_d2_p2;
		Y_8_d4_p3 = -Y_6_d3_p2-Y_3_d2_p1*gi2_p4;
		Y_8_d5_p1 = Y_8_d4_p3*gi[9];
		Y_8_d6_p1 = Y_8_d4_p2*gi[10];
		Y_8_d6_p2 = Y_8_d5_p1*gi[10];
		Y_8_d7_p1 = Y_8_d6_p1*gi[11]+Y_8_d6_p2;
		Y_8_d8_p1 = Y_8_d4_p1*gi[12]+Y_8_d7_p1;
		Y_p22 = Y_8_d8_p1+Y_8_d3_p3*gi[14];

		Y_9_d2_p1 = Y_5_d1_p2+gi2_p3;
		Y_9_d3_p1 = Y_9_d2_p1*gi[9];
		Y_9_d4_p1 = Y_9_d3_p1*gi[10];
		Y_p25 = Y_5_d3_p4*gi[12]+Y_9_d4_p1;

		Y_p28 = gi2_p3*gi[12]-Y_3_d4_p1;

		Y_11_d3_p1 = -Y_2_d1_p1*gi[5];
		Y_11_d4_p1 = -98.7036946115134*gi[6];
		Y_11_d5_p1 = Y_11_d4_p1*gi[7];
		Y_11_d6_p1 = gi[4]*gi[8];
		Y_11_d6_p2 = gi[7]*gi[8];
		Y_11_d6_p3 = Y_11_d5_p1*gi[8]+Y_11_d3_p1;
		Y_11_d7_p1 = -Y_2_d6_p2*gi[9];
		Y_11_d8_p1 = -Y_4_d2_p1*gi[10]+Y_11_d7_p1;
		Y_11_d8_p2 = Y_11_d6_p2*gi2_p6+Y_11_d6_p3;
		Y_11_d9_p1 = Y_11_d8_p1*gi[11]+Y_11_d6_p2*gi2_p7+Y_11_d8_p2;
		Y_11_d10_p1 = Y_2_d3_p2*gi[12]+Y_11_d9_p1;
		Y_11_d11_p1 = Y_11_d6_p1*gi[13]+Y_11_d10_p1;
		Y_p31 = Y_2_d5_p2*gi[14]+Y_11_d11_p1;

		Y_12_d3_p1 = gi[7]*gi2_p6;
		Y_p34 = gi[4]*gi[13]+Y_12_d3_p1;
		Y_p37 = gi[12];
		Y_p40 = 9.81*gi[6];
		Y_15_d1_p1 = 9.81*gi[4];
		Y_p43 = Y_15_d1_p1*gi[6];

		Y_p46 = Y_p40*gi[7];
		Y_p49 = 9.81*gi[3];

		Y_p2 = gi[7]*gi[12];

		Y_19_d2_p1 = -Y_2_d1_p1*gi[4];
		Y_19_d3_p1 = Y_19_d2_p1*gi[5];
		Y_p5 = Y_2_d3_p1*gi[12]+Y_19_d3_p1;

		Y_20_d2_p1 = Y_1_d1_p1*gi[7];
		Y_20_d3_p1 = Y_20_d2_p1*gi2_p5;
		Y_p8 = gi[13]+Y_20_d3_p1;

		Y_21_d3_p1 = gi2_p3+Y_5_d1_p2;
		Y_21_d4_p1 = Y_8_d3_p3*gi[8];
		Y_21_d4_p2 = Y_21_d3_p1*gi[8];
		Y_21_d5_p1 = Y_21_d4_p1*gi[9];
		Y_21_d5_p2 = Y_21_d4_p2*gi2_p5;
		Y_21_d6_p1 = Y_21_d5_p1*gi[11]+gi[8]*gi2_p7+Y_21_d5_p2;
		Y_21_d7_p1 = Y_2_d5_p2*gi[12]+Y_21_d6_p1;
		Y_p11 = gi[5]*gi[14]+Y_21_d7_p1;
		Y_22_d1_p1 = 0.5*gi2_p1;
		Y_22_d2_p1 = Y_22_d1_p1*gi[5];
		Y_22_d2_p2 = -0.5*gi[5];
		Y_22_d3_p1 = Y_22_d2_p1+Y_22_d2_p2*gi2_p3;
		Y_22_d3_p2 = Y_7_d2_p2*gi[7];
		Y_22_d4_p1 = 0.5*gi[8];
		Y_22_d5_p1 = Y_22_d3_p1*gi2_p5;
		Y_22_d5_p2 = Y_22_d3_p2*gi[9];
		Y_22_d6_p1 = Y_22_d5_p1+Y_22_d2_p2*gi2_p7+Y_22_d5_p2*gi[11];
		Y_22_d7_p1 = Y_5_d4_p1*gi[12]+Y_22_d6_p1;
		Y_p14 = Y_22_d4_p1*gi[14]+Y_22_d7_p1;

		Y_23_d2_p1 = Y_1_d1_p1*gi2_p2;
		Y_23_d3_p1 = Y_23_d2_p1*gi[7];
		Y_23_d4_p1 = -Y_2_d3_p2*gi[8];
		Y_23_d5_p1 = Y_23_d3_p1*gi2_p5;
		Y_23_d5_p2 = -Y_6_d2_p2*gi[9];
		Y_23_d6_p1 = Y_23_d4_p1*gi[10]+Y_23_d5_p2;
		Y_23_d7_p1 = Y_23_d6_p1*gi[11]+Y_23_d5_p1;
		Y_23_d8_p1 = Y_6_d4_p1*gi[12]+Y_23_d7_p1;
		Y_p17 = Y_23_d8_p1+Y_6_d2_p1*gi[13];

		Y_24_d4_p1 = gi[5]*gi[8];
		Y_24_d4_p2 = gi2_p4+Y_6_d2_p1;
		Y_24_d4_p3 = Y_5_d3_p1*gi[8];
		Y_24_d5_p1 = Y_7_d4_p3*gi[9];
		Y_24_d5_p2 = Y_24_d4_p3*gi2_p5;
		Y_24_d6_p1 = Y_24_d5_p1+Y_24_d4_p2*gi[10];
		Y_24_d7_p1 = Y_24_d6_p1*gi[11]+Y_24_d5_p2;
		Y_24_d8_p1 = Y_7_d4_p2*gi[12]+Y_24_d7_p1;
		Y_p20 = Y_24_d4_p1*gi[13]+Y_24_d8_p1;

		Y_25_d2_p1 = gi[4]*gi2_p2;
		Y_25_d3_p1 = Y_25_d2_p1*gi[7];
		Y_25_d4_p1 = gi2_p4+gi2_p2;
		Y_25_d4_p2 = Y_5_d3_p4*gi2_p4+Y_25_d3_p1;
		Y_25_d4_p3 = Y_2_d2_p1*gi2_p4+Y_6_d2_p2;
		Y_25_d5_p1 = Y_25_d4_p2*gi2_p5;
		Y_25_d5_p2 = Y_25_d4_p3*gi[9];
		Y_25_d6_p1 = Y_25_d5_p1+Y_25_d5_p2*gi[11];
		Y_p23 = Y_25_d4_p1*gi[13]+Y_25_d6_p1;

		Y_26_d2_p1 = -0.5*gi2_p3+Y_22_d1_p1;
		Y_p26 = Y_26_d2_p1*gi2_p5;

		Y_27_d3_p1 = Y_5_d3_p4*gi2_p5;
		Y_p29 = gi[13]+Y_27_d3_p1;

		Y_28_d3_p1 = Y_19_d2_p1*gi[8];
		Y_p32 = Y_11_d6_p1*gi[12]+Y_28_d3_p1;

		Y_p35 = gi[4]*gi[12];

		Y_p44 = Y_p49*gi[7];

		Y_p47 = -Y_p49*gi[4];

		Y_32_d4_p1 = Y_2_d1_p1*gi[7];
		Y_32_d5_p1 = Y_2_d4_p1+Y_32_d4_p1*gi[8];
		Y_32_d6_p1 = Y_7_d2_p2*gi2_p5+Y_32_d5_p1;
		Y_p6 = Y_2_d6_p1*gi[12]+Y_32_d6_p1;

		Y_33_d3_p1 = Y_4_d2_p2*gi[7];

		Y_33_d5_p1 = Y_33_d3_p1*gi2_p5;
		Y_33_d6_p1 = Y_33_d5_p1-Y_21_d5_p1*gi[10];
		Y_33_d7_p1 = Y_4_d4_p1*gi[12]+Y_33_d6_p1;
		Y_p12 = gi[5]*gi[13]+Y_33_d7_p1;

		Y_34_d3_p1 = -Y_5_d1_p1*gi[7];
		Y_34_d4_p1 = Y_34_d3_p1*gi[8];
		Y_34_d5_p1 = Y_2_d5_p2*gi[9];
		Y_34_d5_p2 = Y_34_d4_p1*gi2_p5;
		Y_34_d6_p1 = Y_34_d5_p1*gi[10]+Y_34_d5_p2;
		Y_34_d7_p1 = Y_5_d2_p1*gi[12]+Y_34_d6_p1;
		Y_p15 = Y_34_d7_p1+Y_22_d4_p1*gi[13];

		Y_35_d5_p1 = Y_7_d4_p1*gi2_p5;
		Y_35_d6_p1 = -Y_23_d5_p2*gi[10]+Y_24_d4_p1*gi2_p6+Y_35_d5_p1;
		Y_35_d7_p1 = Y_2_d5_p1*gi[12]+Y_35_d6_p1;
		Y_p18 = Y_35_d7_p1-gi[14];

		Y_36_d3_p1 = 0.5-0.5*gi2_p3;
		Y_36_d3_p2 = Y_7_d2_p3*gi2_p3-Y_7_d2_p3;
		Y_36_d4_p1 = Y_36_d3_p1*gi2_p4+Y_36_d3_p2;
		Y_36_d4_p2 = Y_7_d2_p3-0.5*gi2_p4;
		Y_36_d5_p1 = Y_36_d4_p1*gi2_p5;
		Y_p21 = Y_36_d5_p1+Y_36_d4_p2*gi2_p6-Y_24_d5_p1*gi[10];

		Y_37_d4_p1 = -Y_6_d2_p2-Y_2_d2_p1*gi2_p4;
		Y_37_d5_p1 = Y_37_d4_p1*gi[9];
		Y_37_d6_p1 = Y_37_d5_p1*gi[10];
		Y_37_d7_p1 = Y_8_d3_p3*gi[12]+Y_37_d6_p1;
		Y_p24 = Y_37_d7_p1+2*gi[14];

		Y_38_d4_p1 = Y_11_d3_p1*gi[7];
		Y_38_d5_p1 = -gi[8];
		Y_38_d5_p2 = Y_11_d4_p1*gi[8]+Y_38_d4_p1;
		Y_38_d6_p1 = Y_38_d5_p1*gi2_p5+Y_38_d5_p2;
		Y_p33 = Y_2_d5_p2*gi[12]+Y_38_d6_p1;
		tau[startInd+0] = Y_p1*Theta[0]+Y_p4*Theta[1]+Y_p7*Theta[2]+Y_p10*Theta[3]+Y_p13*Theta[4]+Y_p16*Theta[5]+Y_p19*Theta[6]+Y_p22*Theta[7]+Y_p25*Theta[8]+Y_p28*Theta[9]+Y_p31*Theta[10]+Y_p34*Theta[11]+Y_p37*Theta[12]+Y_p40*Theta[13]+Y_p43*Theta[14]+Y_p46*Theta[15]+Y_p49*Theta[16];
		tau[startInd+1] = Y_p2*Theta[0]+Y_p5*Theta[1]+Y_p8*Theta[2]+Y_p11*Theta[3]+Y_p14*Theta[4]+Y_p17*Theta[5]+Y_p20*Theta[6]+Y_p23*Theta[7]+Y_p26*Theta[8]+Y_p29*Theta[9]+Y_p32*Theta[10]+Y_p35*Theta[11]+Y_p44*Theta[14]+Y_p47*Theta[15];
		tau[startInd+2] = Y_p6*Theta[1]+Y_p12*Theta[3]+Y_p15*Theta[4]+Y_p18*Theta[5]+Y_p21*Theta[6]+Y_p24*Theta[7]+Y_p33*Theta[10];

	}
}