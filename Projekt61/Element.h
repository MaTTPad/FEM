#pragma once
#include "Edge.h"
class Element
{
public:
	double ID[4], k, ro, c, alfa;
	double jacobian[4][4];
	double H[4][4], H2[4][4];
	double matrixC[4][4];
	Edge edges[4];

public:
	Element();
	~Element();
};

