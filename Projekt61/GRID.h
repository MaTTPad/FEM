#pragma once
#include "Element.h"
#include "Node.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

class GRID
{
public:
	double H, L, nH, nL;
	Node *nodes;
	Element *elements;

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
	GRID();
	~GRID();
	void showNodes();
	void showElements();
	void showAll();
	void showElementsNodes(int);
	void showNodesCords(int);

	void showMatrixH();
	void showMatrixH2();
	void showMatrixC();

};

