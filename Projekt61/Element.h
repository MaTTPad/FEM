#pragma once
#include "Edge.h"
class Element
{
public:
	int ID[4];
	double k, ro, c, alfa, tempotocz;
	double jacobian[4][4];
	double H[4][4], H2[4][4], matrixH[4][4];
	double matrixC[4][4], matrixP[4];
	Edge edges[4];

public:
	Element();
	~Element();
};

