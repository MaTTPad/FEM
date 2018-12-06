#pragma once
#include "Element.h"
#include "Node.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <nlohmann/json.hpp>

using json = nlohmann::json;
using namespace std;

class GRID
{
public:
	double H, L, nH, nL, initialTemp, simulationTime, simulationStepTime, ambientTemp, alfa, specificHeat, conductivity, density;
	Node *nodes;
	Element *elements;

	double etaTab[4];
	double ksiTab[4];

	double **globalMatrixH;
	double **globalMatrixC;
	double **MatrixHzDaszkiem;
	double *globalMatrixP;
	double *MatrixPzDaszkiem;



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
	void showMatrixHlokal();
	void showMatrixC();
	void calculateGlobalMatrixH();
	void calculateGlobalMatrixC();
	void showGlobalMatrixH();
	void showGlobalMatrixC();
	void calculateMatrixHzDaszkiem();
	void showMatrixHzDaszkiem();
	void showMatrixPzDaszkiem();

};

