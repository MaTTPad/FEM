#pragma once
#include "GRID.h"
class Matrix
{
	double etaTab[4];
	double ksiTab[4];
	double Ntab[4][4];
	double dnPoDksi[4][4];
	double dnPoDeta[4][4];
	double jacobian[4][4];
	double detJ[4];
	double jacobianInverse[4][4];
	double dnPoDx[4][4];
	double dnPoDy[4][4];

	double dnPoDxDnPoDxTpkt1[4][4];
	double dnPoDxDnPoDxTpkt2[4][4];
	double dnPoDxDnPoDxTpkt3[4][4];
	double dnPoDxDnPoDxTpkt4[4][4];
	double dnPoDyDnPoDyTpkt1[4][4];
	double dnPoDyDnPoDyTpkt2[4][4];
	double dnPoDyDnPoDyTpkt3[4][4];
	double dnPoDyDnPoDyTpkt4[4][4];

	double dnPoDxDnPoDxTpkt1razydetJ[4][4];
	double dnPoDxDnPoDxTpkt2razydetJ[4][4];
	double dnPoDxDnPoDxTpkt3razydetJ[4][4];
	double dnPoDxDnPoDxTpkt4razydetJ[4][4];
	double dnPoDyDnPoDyTpkt1razydetJ[4][4];
	double dnPoDyDnPoDyTpkt2razydetJ[4][4];
	double dnPoDyDnPoDyTpkt3razydetJ[4][4];
	double dnPoDyDnPoDyTpkt4razydetJ[4][4];

	double kMacierzDetjpkt1[4][4];
	double kMacierzDetjpkt2[4][4];
	double kMacierzDetjpkt3[4][4];
	double kMacierzDetjpkt4[4][4];
	double matrixH[4][4];
public:
	Matrix();
	Matrix(GRID *a);
	~Matrix();

	double calculateN1(double eta, double ksi);
	double calculateN2(double eta, double ksi);
	double calculateN3(double eta, double ksi);
	double calculateN4(double eta, double ksi);
};

