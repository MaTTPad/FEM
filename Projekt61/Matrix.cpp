#include "Matrix.h"
#include "GRID.H"

const double eps = 1e-12; // sta³a przybli¿enia zera
bool gauss(int n, double ** AB, double * X);

Matrix::Matrix(GRID *grid)
{
	myGrid = grid;
	double elementsCount = (grid->nH - 1)*(grid->nL - 1);

	ksiTab[0] = -1 / sqrt(3);
	etaTab[0] = ksiTab[0];
	ksiTab[1] = -ksiTab[0];
	etaTab[1] = etaTab[0];
	ksiTab[2] = ksiTab[1];
	etaTab[2] = -etaTab[1];
	ksiTab[3] = -ksiTab[2];
	etaTab[3] = etaTab[2];

	double matrixH2[4][4];
	double matrixP[4];



	////////////////////////////////////////////////////////////////////// Uzupe³nianie Ntab

	for (int i = 0;i < 4;i++)
	{
		Ntab[i][0] = calculateN1(ksiTab[i], etaTab[i]);
		Ntab[i][1] = calculateN2(ksiTab[i], etaTab[i]);
		Ntab[i][2] = calculateN3(ksiTab[i], etaTab[i]);
		Ntab[i][3] = calculateN4(ksiTab[i], etaTab[i]);

	}

	/*for (int i = 0;i < 4;i++)
	{
	cout <<"N1="<< Ntab[i][0] << endl;
	cout << "N2="<<Ntab[i][1] << endl;   ///wyœwietlanie Ntab, poprawne
	cout << "N3="<<Ntab[i][2] << endl;
	cout << "N4="<<Ntab[i][3] << endl;
	cout << "=======================================" << endl;
	}
	*/

	///////////////////////////////////////////////////////////////////////////////////


	////////////////////////////////////////////////////////Wyliczanie dN/dKsi i dN/dEta


	for (int i = 0;i < 4;i++)
	{
		dnPoDksi[i][0] = -0.25*(1 - etaTab[i]);
		dnPoDksi[i][1] = 0.25*(1 - etaTab[i]);
		dnPoDksi[i][2] = 0.25*(1 + etaTab[i]);
		dnPoDksi[i][3] = -0.25*(1 + etaTab[i]);

		dnPoDeta[i][0] = -0.25*(1 - ksiTab[i]);
		dnPoDeta[i][1] = -0.25*(1 + ksiTab[i]);
		dnPoDeta[i][2] = 0.25*(1 + ksiTab[i]);
		dnPoDeta[i][3] = 0.25*(1 - ksiTab[i]);

	}
	/*
	for (int i = 0;i < 4;i++)
	{
	cout << "dN/dKsi "<<i<<"="<<dnPoDksi[i][0] << endl;
	cout << "dN/dKsi " << i << "=" << dnPoDksi[i][1] << endl;
	cout << "dN/dKsi " << i << "=" << dnPoDksi[i][2] << endl;
	cout << "dN/dKsi " << i << "=" << dnPoDksi[i][3] << endl;
	cout << endl;
	cout << "dN/dEta " << i << "=" << dnPoDeta[i][0] << endl;
	cout << "dN/dEta " << i << "=" << dnPoDeta[i][1] << endl;
	cout << "dN/dEta " << i << "=" << dnPoDeta[i][2] << endl;
	cout << "dN/dEta " << i << "=" << dnPoDeta[i][3] << endl;

	cout << "==================================" << endl;
	}
	*/

	///////////////////////////////////////////////////////////////////////////////////


	////////////////////////////////////////////////////////////////////// Wyliczanie Jacobianów

	for (int i = 0;i < elementsCount;i++)
	{
		for (int j = 0;j < 4;j++)
		{

			grid->elements[i].jacobian[j][0] = dnPoDksi[j][0] * grid->nodes[(int)grid->elements[j].ID[0]].x + dnPoDksi[j][1] * grid->nodes[(int)grid->elements[j].ID[1]].x + dnPoDksi[j][2] * grid->nodes[(int)grid->elements[j].ID[2]].x + dnPoDksi[j][3] * grid->nodes[(int)grid->elements[j].ID[3]].x;
			grid->elements[i].jacobian[j][2] = dnPoDeta[j][0] * grid->nodes[(int)grid->elements[j].ID[0]].x + dnPoDeta[j][1] * grid->nodes[(int)grid->elements[j].ID[1]].x + dnPoDeta[j][2] * grid->nodes[(int)grid->elements[j].ID[2]].x + dnPoDeta[j][3] * grid->nodes[(int)grid->elements[j].ID[3]].x;
			grid->elements[i].jacobian[j][1] = dnPoDksi[j][0] * grid->nodes[(int)grid->elements[j].ID[0]].y + dnPoDksi[j][1] * grid->nodes[(int)grid->elements[j].ID[1]].y + dnPoDksi[j][2] * grid->nodes[(int)grid->elements[j].ID[2]].y + dnPoDksi[j][3] * grid->nodes[(int)grid->elements[j].ID[3]].y;
			grid->elements[i].jacobian[j][3] = dnPoDeta[j][0] * grid->nodes[(int)grid->elements[j].ID[0]].y + dnPoDeta[j][1] * grid->nodes[(int)grid->elements[j].ID[1]].y + dnPoDeta[j][2] * grid->nodes[(int)grid->elements[j].ID[2]].y + dnPoDeta[j][3] * grid->nodes[(int)grid->elements[j].ID[3]].y;
		}



		for (int j = 0;j < 4;j++)
		{
			jacobian[j][0] = grid->elements[i].jacobian[j][0];
			jacobian[j][1] = grid->elements[i].jacobian[j][1];
			jacobian[j][2] = grid->elements[i].jacobian[j][2];
			jacobian[j][3] = grid->elements[i].jacobian[j][3];
		}

		//////////////////////////////////////////////////////////////////////////////////////////



		////////////////////////////////////////////////////////////////////////////// Wyznacznik detJ

		for (int j = 0;j < 4;j++)
			detJ[j] = grid->elements[0].jacobian[j][0] * grid->elements[0].jacobian[j][3] - grid->elements[0].jacobian[j][1] * grid->elements[0].jacobian[j][2];

		////////////////////////////////////////////////////////////////////////////////////////// J^-1

		for (int j = 0;j < 4;j++)
		{

			jacobianInverse[j][0] = jacobian[j][3] / detJ[j];
			jacobianInverse[j][1] = -jacobian[j][1] / detJ[j];
			jacobianInverse[j][2] = -jacobian[j][2] / detJ[j]; //TODO upewniæ sie, czy w excelu b³¹d
			jacobianInverse[j][3] = jacobian[j][0] / detJ[j];

		}

		/*
		for (int i = 0;i < 4;i++)
		{

		cout << jacobianInverse[i][0] << endl;
		cout << jacobianInverse[i][1] << endl;
		cout << jacobianInverse[i][2] << endl;
		cout << jacobianInverse[i][3] << endl;
		cout << "==============================" << endl;

		}
		*/


		/////////////////////////////////////////////////////////////////    Uzupe³nianie dn/dx i dn/dy

		for (int j = 0;j < 4;j++)
		{
			//dnPoDx[j][0] = jacobianInverse[j][0] * dnPoDksi[j][0] + jacobianInverse[j][1] * dnPoDeta[j][0];
			//dnPoDx[j][1] = jacobianInverse[j][0] * dnPoDksi[j][1] + jacobianInverse[j][1] * dnPoDeta[j][1];
			//dnPoDx[j][2] = jacobianInverse[j][0] * dnPoDksi[j][2] + jacobianInverse[j][1] * dnPoDeta[j][2];
			//dnPoDx[j][3] = jacobianInverse[j][0] * dnPoDksi[j][3] + jacobianInverse[j][1] * dnPoDeta[j][3];
			for (int k = 0;k < 4;k++)
			{
				dnPoDx[j][k] = jacobianInverse[j][0] * dnPoDksi[j][k] + jacobianInverse[j][1] * dnPoDeta[j][k];
				dnPoDy[j][k] = jacobianInverse[j][2] * dnPoDksi[j][k] + jacobianInverse[j][3] * dnPoDeta[j][k];
			}
		}

		/*
		for (int i = 0;i < 4;i++)
		{
		cout << dnPoDy[i][0] << endl;
		cout << dnPoDy[i][1] << endl;
		cout << dnPoDy[i][2] << endl;
		cout << dnPoDy[i][3] << endl;
		cout << "+=============================+"<<endl;

		}
		*/

		////////////////////////////////////////////////////////////////////////////////////////////////

		///////////////////////////////////////////////////////////////////////////////Uzupelnianie {dN/dx}{dN/dx}T i {dN/dy}{dN/dy}T

		for (int j = 0;j < 4;j++)
		{

			for (int k = 0;k < 4;k++)
			{
				dnPoDxDnPoDxTpkt1[j][k] = dnPoDx[0][j] * dnPoDx[0][k];
				dnPoDxDnPoDxTpkt2[j][k] = dnPoDx[1][j] * dnPoDx[1][k];
				dnPoDxDnPoDxTpkt3[j][k] = dnPoDx[2][j] * dnPoDx[2][k];
				dnPoDxDnPoDxTpkt4[j][k] = dnPoDx[3][j] * dnPoDx[3][k];


				dnPoDyDnPoDyTpkt1[j][k] = dnPoDy[0][j] * dnPoDy[0][k];
				dnPoDyDnPoDyTpkt2[j][k] = dnPoDy[1][j] * dnPoDy[1][k];
				dnPoDyDnPoDyTpkt3[j][k] = dnPoDy[2][j] * dnPoDy[2][k];
				dnPoDyDnPoDyTpkt4[j][k] = dnPoDy[3][j] * dnPoDy[3][k];

			}
		}

		/*
		for (int i = 0;i<4;i++)
		{

		for (int j = 0;j < 4;j++)
		{
		cout << dnPoDyDnPoDyTpkt1[i][j] << " ";
		}
		cout << endl;
		}
		cout << "========================" << endl;
		system("PAUSE");

		*/
		////////////////////////////////////////////////////////////////////////////////////////



		///////////////////////////////////////////////////////////////////{dN/dx}{dN/dx}T*DetJ i {dN/dy}{dN/dy}T*DetJ



		for (int j = 0;j < 4;j++)
		{
			for (int k = 0;k < 4;k++)
			{
				dnPoDxDnPoDxTpkt1razydetJ[j][k] = dnPoDxDnPoDxTpkt1[j][k] * detJ[j];
				dnPoDxDnPoDxTpkt2razydetJ[j][k] = dnPoDxDnPoDxTpkt2[j][k] * detJ[j];
				dnPoDxDnPoDxTpkt3razydetJ[j][k] = dnPoDxDnPoDxTpkt3[j][k] * detJ[j];
				dnPoDxDnPoDxTpkt4razydetJ[j][k] = dnPoDxDnPoDxTpkt4[j][k] * detJ[j];

				dnPoDyDnPoDyTpkt1razydetJ[j][k] = dnPoDyDnPoDyTpkt1[j][k] * detJ[j];
				dnPoDyDnPoDyTpkt2razydetJ[j][k] = dnPoDyDnPoDyTpkt2[j][k] * detJ[j];
				dnPoDyDnPoDyTpkt3razydetJ[j][k] = dnPoDyDnPoDyTpkt3[j][k] * detJ[j];
				dnPoDyDnPoDyTpkt4razydetJ[j][k] = dnPoDyDnPoDyTpkt4[j][k] * detJ[j];

			}
		}
		/*
		for (int i = 0;i<4;i++)
		{

		for (int j = 0;j < 4;j++)
		{
		cout << dnPoDyDnPoDyTpkt1razydetJ[i][j] << " ";
		}
		cout << endl;
		}
		cout << "========================" << endl;
		system("PAUSE");
		*/
		////////////////////////////////////////////////////////////////////



		//////////////////////////////////////////////////////////////////////// K*(     {dN/dx}{dN/dx}T  +  {dN/dy}{dN/dy}T)*DetJ

		for (int k = 0;k < 4;k++)
		{
			for (int j = 0;j < 4;j++)
			{
				kMacierzDetjpkt1[k][j] = grid->elements[i].k*(dnPoDxDnPoDxTpkt1razydetJ[k][j] + dnPoDyDnPoDyTpkt1razydetJ[k][j]);
				kMacierzDetjpkt2[k][j] = grid->elements[i].k*(dnPoDxDnPoDxTpkt2razydetJ[k][j] + dnPoDyDnPoDyTpkt2razydetJ[k][j]);
				kMacierzDetjpkt3[k][j] = grid->elements[i].k*(dnPoDxDnPoDxTpkt3razydetJ[k][j] + dnPoDyDnPoDyTpkt3razydetJ[k][j]);
				kMacierzDetjpkt4[k][j] = grid->elements[i].k*(dnPoDxDnPoDxTpkt4razydetJ[k][j] + dnPoDyDnPoDyTpkt4razydetJ[k][j]);

			}
		}



		//////////////////////////////////////////////////////////////////////////////////////// Macierz H

		//	cout << "Macierz H dla elementu numer: " << i << endl;
		for (int k = 0;k < 4;k++)
		{
			for (int j = 0;j < 4;j++)
			{
				matrixH[k][j] = kMacierzDetjpkt1[k][j] + kMacierzDetjpkt2[k][j] + kMacierzDetjpkt3[k][j] + kMacierzDetjpkt4[k][j];
			}
		}

		for (int k = 0;k < 4;k++) {
			for (int j = 0;j < 4;j++)
			{
				grid->elements[i].H[k][j] = matrixH[k][j];

			}
		}

		/*
		for (int k = 0;k < 4;k++) {
		for (int j = 0;j < 4;j++)
		{
		cout << elements[i].H[k][j] << " ";
		}
		cout << endl;
		}
		*/		//	cout << endl;



		/*
		for (int k = 0;k < 4;k++) {
		for (int j = 0;j < 4;j++)
		{
		cout << grid->elements[0].H[k][j] << " ";
		}
		cout << endl;
		}*/



		////////////////////////////////////////////// Macierz C //////////////////////////////////////////////////////////

		double nnt1cRoDetj[4][4], nnt2cRoDetj[4][4], nnt3cRoDetj[4][4], nnt4cRoDetj[4][4];

		for (int j = 0;j < 4;j++)
		{
			for (int k = 0;k < 4;k++)
			{
				nnt1cRoDetj[j][k] = Ntab[0][j] * Ntab[0][k] * detJ[0] * grid->elements[i].c*grid->elements[i].ro;
				nnt2cRoDetj[j][k] = Ntab[1][j] * Ntab[1][k] * detJ[1] * grid->elements[i].c*grid->elements[i].ro;
				nnt3cRoDetj[j][k] = Ntab[2][j] * Ntab[2][k] * detJ[2] * grid->elements[i].c*grid->elements[i].ro;
				nnt4cRoDetj[j][k] = Ntab[3][j] * Ntab[3][k] * detJ[3] * grid->elements[i].c*grid->elements[i].ro;
			}
		}

		for (int j = 0;j < 4;j++)
		{
			for (int k = 0;k < 4;k++)
			{
				grid->elements[i].matrixC[j][k] += (nnt1cRoDetj[j][k] + nnt2cRoDetj[j][k] + nnt3cRoDetj[j][k] + nnt4cRoDetj[j][k]);
			}
		}

		/*
		cout << "Macierz C" << endl;
		for (int j = 0;j < 4;j++)
		{
		for (int k = 0;k < 4;k++)
		{
		cout<<grid->elements[i].matrixC[j][k]<< " ";
		}
		cout << endl;
		}
		*/

		/////////////////////////////////////////////////////////// druga czêœæ macierzy H
		double ksiPow[8], etaPow[8];
		ksiPow[0] = -1 / sqrt(3);
		ksiPow[1] = -ksiPow[0];
		ksiPow[2] = 1;
		ksiPow[3] = 1;
		ksiPow[4] = 1 / sqrt(3);
		ksiPow[5] = -ksiPow[4];
		ksiPow[6] = -1;
		ksiPow[7] = -1;

		etaPow[0] = -1;
		etaPow[1] = -1;
		etaPow[2] = -1 / sqrt(3);
		etaPow[3] = -etaPow[2];
		etaPow[4] = 1;
		etaPow[5] = 1;
		etaPow[6] = 1 / sqrt(3);
		etaPow[7] = -etaPow[6];




		double Nshape[2][4];

		double sideLength[4];
		double detJ[4];
		sideLength[0] = sqrt(pow(grid->nodes[(int)grid->elements[i].ID[1]].x - grid->nodes[(int)grid->elements[i].ID[0]].x, 2) + pow(grid->nodes[(int)grid->elements[i].ID[1]].y - grid->nodes[(int)grid->elements[i].ID[0]].y, 2));
		sideLength[1] = sqrt(pow(grid->nodes[(int)grid->elements[i].ID[2]].x - grid->nodes[(int)grid->elements[i].ID[1]].x, 2) + pow(grid->nodes[(int)grid->elements[i].ID[2]].y - grid->nodes[(int)grid->elements[i].ID[1]].y, 2));
		sideLength[2] = sqrt(pow(grid->nodes[(int)grid->elements[i].ID[3]].x - grid->nodes[(int)grid->elements[i].ID[2]].x, 2) + pow(grid->nodes[(int)grid->elements[i].ID[3]].y - grid->nodes[(int)grid->elements[i].ID[2]].y, 2));
		sideLength[3] = sqrt(pow(grid->nodes[(int)grid->elements[i].ID[3]].x - grid->nodes[(int)grid->elements[i].ID[0]].x, 2) + pow(grid->nodes[(int)grid->elements[i].ID[3]].y - grid->nodes[(int)grid->elements[i].ID[0]].y, 2));

		detJ[0] = sideLength[0] / 2;
		detJ[1] = sideLength[1] / 2;
		detJ[2] = sideLength[2] / 2;
		detJ[3] = sideLength[3] / 2;

		double pc1[4][4], pc2[4][4], sum[4][4];

		for (int j = 0;j < 4;j++)
		{
			for (int k = 0;k < 4;k++)
			{
				matrixH2[j][k] = 0;
			}
		}
		for (int e = 0;e < 4;e++) //4 powierzchnie
		{
			Nshape[0][0] = calculateN1(ksiPow[2 * e], etaPow[2 * e]);
			Nshape[0][1] = calculateN2(ksiPow[2 * e], etaPow[2 * e]);
			Nshape[0][2] = calculateN3(ksiPow[2 * e], etaPow[2 * e]);
			Nshape[0][3] = calculateN4(ksiPow[2 * e], etaPow[2 * e]);

			Nshape[1][0] = calculateN1(ksiPow[2 * e + 1], etaPow[2 * e + 1]);
			Nshape[1][1] = calculateN2(ksiPow[2 * e + 1], etaPow[2 * e + 1]);
			Nshape[1][2] = calculateN3(ksiPow[2 * e + 1], etaPow[2 * e + 1]);
			Nshape[1][3] = calculateN4(ksiPow[2 * e + 1], etaPow[2 * e + 1]);

			for (int j = 0;j < 4;j++)
			{
				for (int k = 0;k < 4;k++)
				{
					pc1[j][k] = Nshape[0][j] * Nshape[0][k] * grid->elements[i].alfa;
					pc2[j][k] = Nshape[1][j] * Nshape[1][k] * grid->elements[i].alfa;
					sum[j][k] = (pc1[j][k] + pc2[j][k])*detJ[e];
					//cout << sum[j][k] << " ";
				}
				//cout << endl;
			}
			//cout << endl;

			if (grid->elements[i].edges[e].isBoundary())
			{
				for (int j = 0;j < 4;j++)
				{
					for (int k = 0;k < 4;k++)
					{
						matrixH2[j][k] += sum[j][k];
					}
				}
			}

		}
		for (int j = 0;j < 4;j++)
		{
			for (int k = 0;k < 4;k++)
			{
				grid->elements[i].H2[j][k] = matrixH2[j][k];
				//	cout << matrixH2[j][k] << " ";
			}
			//cout << endl;
		}
		//cout << endl;


		/////////////////////////////////////////////////////////////////////////////////////////koniec macierzy H2

		//////////////////////////////////////////////////////////////////////// £¹czenie macierzy H w ca³oœæ 



		for (int j = 0;j < 4;j++)
		{
			for (int k = 0;k < 4;k++)
			{
				grid->elements[i].matrixH[j][k] = grid->elements[i].H[j][k] + grid->elements[i].H2[j][k];
				//	cout << matrixH2[j][k] << " ";
			}
			//cout << endl;
		}
		//cout << endl;



		///////////////////////////////////////////////////////////////////////////////////////////



		//////////////////////////////////////////////////////////////////////////////////// Macierz P




		sideLength[0] = sqrt(pow(grid->nodes[(int)grid->elements[i].ID[1]].x - grid->nodes[(int)grid->elements[i].ID[0]].x, 2) + pow(grid->nodes[(int)grid->elements[i].ID[1]].y - grid->nodes[(int)grid->elements[i].ID[0]].y, 2));
		sideLength[1] = sqrt(pow(grid->nodes[(int)grid->elements[i].ID[2]].x - grid->nodes[(int)grid->elements[i].ID[1]].x, 2) + pow(grid->nodes[(int)grid->elements[i].ID[2]].y - grid->nodes[(int)grid->elements[i].ID[1]].y, 2));
		sideLength[2] = sqrt(pow(grid->nodes[(int)grid->elements[i].ID[3]].x - grid->nodes[(int)grid->elements[i].ID[2]].x, 2) + pow(grid->nodes[(int)grid->elements[i].ID[3]].y - grid->nodes[(int)grid->elements[i].ID[2]].y, 2));
		sideLength[3] = sqrt(pow(grid->nodes[(int)grid->elements[i].ID[3]].x - grid->nodes[(int)grid->elements[i].ID[0]].x, 2) + pow(grid->nodes[(int)grid->elements[i].ID[3]].y - grid->nodes[(int)grid->elements[i].ID[0]].y, 2));

		detJ[0] = sideLength[0] / 2;
		detJ[1] = sideLength[1] / 2;
		detJ[2] = sideLength[2] / 2;
		detJ[3] = sideLength[3] / 2;

		double ppc1[4], ppc2[4], psum[4];

		for (int j = 0;j < 4;j++)
		{

			matrixP[j] = 0;
		}


		for (int e = 0;e < 4;e++) //4 powierzchnie
		{
			Nshape[0][0] = calculateN1(ksiPow[2 * e], etaPow[2 * e]);
			Nshape[0][1] = calculateN2(ksiPow[2 * e], etaPow[2 * e]);
			Nshape[0][2] = calculateN3(ksiPow[2 * e], etaPow[2 * e]);
			Nshape[0][3] = calculateN4(ksiPow[2 * e], etaPow[2 * e]);

			Nshape[1][0] = calculateN1(ksiPow[2 * e + 1], etaPow[2 * e + 1]);
			Nshape[1][1] = calculateN2(ksiPow[2 * e + 1], etaPow[2 * e + 1]);
			Nshape[1][2] = calculateN3(ksiPow[2 * e + 1], etaPow[2 * e + 1]);
			Nshape[1][3] = calculateN4(ksiPow[2 * e + 1], etaPow[2 * e + 1]);

			for (int j = 0;j < 4;j++)
			{

				ppc1[j] = Nshape[0][j] * grid->elements[i].alfa*grid->elements[i].tempotocz;
				ppc2[j] = Nshape[1][j] * grid->elements[i].alfa*grid->elements[i].tempotocz;
				psum[j] = (ppc1[j] + ppc2[j])*detJ[e];

			}

			if (grid->elements[i].edges[e].isBoundary())
			{
				for (int j = 0;j < 4;j++)
				{
					matrixP[j] += psum[j];
				}
			}

		}
		for (int j = 0;j < 4;j++)
		{

			grid->elements[i].matrixP[j] = -matrixP[j];
			cout << -matrixP[j] << " ";
			cout << endl;
		}
		cout << endl;



		/////////////////////////////////////////////////////////////////////////////////////// Koniec macierzy P


		////////////////////////////////////////////////////////////////////////////////////// £¹czenie macierzy w macierz globaln¹ H

		/////	for (int j = 0;j < 4;j++)
		//		for (int k = 0;k < 4;k++)
		///			grid->globalMatrixH[(int)grid->elements[i].ID[j], (int)grid->elements[i].ID[k]]= grid->elements[i].matrixH[j][k];



		////////////////////////////////////////////////////////////////////////////////////////////////





	} //koniec pêtli po wszystkich elementach

	cout << "Simulation time" << endl;
	int simulationTime = 500;
	//cin >> simulationTime;
	int detau = 50;

	for (int lkrokow = 0;lkrokow<(simulationTime / detau);lkrokow++)
	{
		cout << "Iteracja nr " << lkrokow << endl;
		zerujMacierze();

		calculateGlobalMatrixC();
		calculateGlobalMatrixH();
		calculateMatrixHzDaszkiem();
		calculateGlobalMatrixP();
		calculateMatrixPzDaszkiem();


		double **AB, *X;
		int      n, i, j;

		// odczytujemy liczbê niewiadomych

		n = myGrid->nH*myGrid->nL;

		// tworzymy macierze AB i X

		AB = new double *[n];
		X = new double[n];

		for (i = 0; i < n; i++) AB[i] = new double[n + 1];

		// odczytujemy dane dla macierzy AB

		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++)
			{
				AB[i][j] = myGrid->MatrixHzDaszkiem[i][j];
			}
			AB[i][j] = myGrid->MatrixPzDaszkiem[i];
		}
		if (gauss(n, AB, X))
		{
			for (i = 0; i < n; i++)
				cout << "t" << i + 1 << " = " << X[i]
				<< endl;
		}
		else
			cout << "DZIELNIK ZERO\n";

		for (int i = 0;i < myGrid->nH*myGrid->nL;i++)
		{
			myGrid->nodes[i].t = X[i];
		}

		/*
		cout << "Przypisane wartosci temp: " << endl;
		for (int i = 0;i < myGrid->nH*myGrid->nL;i++)
		{
			cout << myGrid->nodes[i].t << endl;
		}*/


		// usuwamy macierze z pamiêci

		for (i = 0; i < n; i++) delete[] AB[i];
		delete[] AB;
		delete[] X;

	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}



Matrix::~Matrix()
{
}



double Matrix::calculateN1(double ksi, double eta)
{
	double N1 = 0.25*(1 - ksi)*(1 - eta);
	return N1;
}

double Matrix::calculateN2(double ksi, double eta)
{
	double N2 = 0.25*(1 + ksi)*(1 - eta);
	return N2;
}

double Matrix::calculateN3(double ksi, double eta)
{
	double N3 = 0.25*(1 + ksi)*(1 + eta);
	return N3;
}

double Matrix::calculateN4(double ksi, double eta)
{
	double N4 = 0.25*(1 - ksi)*(1 + eta);
	return N4;
}


Matrix::Matrix()
{
}




void Matrix::calculateGlobalMatrixH()
{
	cout << "Klasa matrix, obliczam macierz H globalna" << endl;
	for (int i = 0;i < (myGrid->nH - 1)*(myGrid->nL - 1);i++)
	{
		for (int j = 0;j < 4;j++)
			for (int k = 0;k < 4;k++) {
				int index1 = myGrid->elements[i].ID[j];
				int index2 = myGrid->elements[i].ID[k];

				myGrid->globalMatrixH[index1][index2] += myGrid->elements[i].matrixH[j][k];
			}
	}
}

void Matrix::calculateGlobalMatrixC()
{
	cout << "Klasa matrix, obliczam macierz C globalna" << endl;

	for (int i = 0;i < (myGrid->nH - 1)*(myGrid->nL - 1);i++)
	{
		for (int j = 0;j < 4;j++)
			for (int k = 0;k < 4;k++) {
				int index1 = myGrid->elements[i].ID[j];
				int index2 = myGrid->elements[i].ID[k];

				myGrid->globalMatrixC[index1][index2] += myGrid->elements[i].matrixC[j][k];
			}
	}
}

void Matrix::calculateMatrixHzDaszkiem()
{
	cout << "Klasa matrix, obliczam macierz H z daszkiem" << endl;

	//double dtau = 50;
	double **noweC = new double*[myGrid->nH*myGrid->nL];
	for (int i = 0;i < myGrid->nH*myGrid->nL;i++)
		noweC[i] = new double[myGrid->nH*myGrid->nL];

	for (int i = 0;i<myGrid->nH*myGrid->nL;i++)
		for (int j = 0;j < myGrid->nH*myGrid->nL;j++)
		{
			noweC[i][j] = myGrid->globalMatrixC[i][j] / myGrid->elements[0].dtau;
		}

	for (int i = 0;i < (myGrid->nH)*(myGrid->nL);i++)
	{
		for (int j = 0;j < (myGrid->nH)*(myGrid->nL);j++)
		{
			myGrid->MatrixHzDaszkiem[i][j] = myGrid->globalMatrixH[i][j] + noweC[i][j];
		}
	}
}

void Matrix::calculateGlobalMatrixP()
{
	cout << "Klasa matrix, obliczam macierz P globalna" << endl; //chyba coœ z globaln¹ zjebane DX

	for (int i = 0;i < (myGrid->nH - 1)*(myGrid->nL - 1);i++)
	{
		for (int k = 0;k < 4;k++) {
			int index1 = myGrid->elements[i].ID[k];
			myGrid->globalMatrixP[index1] += myGrid->elements[i].matrixP[k];
		}
	}
}

void Matrix::calculateMatrixPzDaszkiem()
{
	cout << "Klasa matrix, obliczam macierz P z daszkiem" << endl;

	double **noweC = new double*[myGrid->nH*myGrid->nL];
	double *noweC2 = new double[myGrid->nH*myGrid->nL];

	for (int i = 0;i < myGrid->nH*myGrid->nL;i++)
	{
		noweC2[i] = 0;
		noweC[i] = new double[myGrid->nH*myGrid->nL];
	}
	for (int i = 0;i<myGrid->nH*myGrid->nL;i++)
		for (int j = 0;j < myGrid->nH*myGrid->nL;j++)
		{
			noweC[i][j] = myGrid->globalMatrixC[i][j] / myGrid->elements[0].dtau;
		}
	/*
	cout << "Nowa macierz C:" << endl; //ok
	for (int i = 0;i < myGrid->nH*myGrid->nL;i++) {
	for (int j = 0;j < myGrid->nH*myGrid->nL;j++)
	{
	cout << noweC[i][j] << " ";
	}
	cout << endl;
	}*/
	///teraz pomno¿yæ * {t0}

	///zrobiæ wektor z temp w ka¿dym wêŸle, bêdzie proœciej ;)

	double *nodesTemp = new double[myGrid->nH*myGrid->nL];

	int counter = 0;
	for (int i = 0;i < (myGrid->nH)*(myGrid->nL);i++)
	{
		/*(for (int j = 0;j < 4;j++)
		{
			//nodesTemp[counter] = myGrid->nodes[(myGrid->elements[i].ID[j])].t;
			counter++;
		}*/
		nodesTemp[i] = myGrid->nodes[i].t;
		
	}

	/*cout << "Temperatury w wezlach:" << endl;
	for (int i = 0;i < myGrid->nH*myGrid->nL;i++)
	{
	cout << "Temp w wezle nr " <<i<<": "<< nodesTemp[i] << endl; //SZTOS
	}
	*/
	for (int i = 0;i < myGrid->nH*myGrid->nL;i++)
	{
		for (int j = 0;j < myGrid->nH*myGrid->nL;j++)
		{
			noweC2[i] += (noweC[i][j] * nodesTemp[j]);
		}
	}
	/*
	cout << "NOWE C2:" << endl;
	for (int i = 0;i < myGrid->nH*myGrid->nL;i++)
	{
	cout << noweC2[i] << " ";
	}
	cout << endl;
	///////////////////////////////////////////////////////

	*/
	//cout << "Globalna p" << endl;

	for (int i = 0;i < (myGrid->nH)*(myGrid->nL);i++)
	{
		//cout << myGrid->globalMatrixP[i] << " ";
		myGrid->MatrixPzDaszkiem[i] = noweC2[i] - myGrid->globalMatrixP[i]; //albo +=

	}
	//cout << endl;








}


bool gauss(int n, double ** AB, double * X)
{

	int i, j, k;
	double m, s;

	// eliminacja wspó³czynników

	for (i = 0; i < n - 1; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			if (fabs(AB[i][i]) < eps) return false;
			m = -AB[j][i] / AB[i][i];
			for (k = i + 1; k <= n; k++)
				AB[j][k] += m * AB[i][k];
		}
	}

	// wyliczanie niewiadomych

	for (i = n - 1; i >= 0; i--)
	{
		s = AB[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= AB[i][j] * X[j];
		if (fabs(AB[i][i]) < eps) return false;
		X[i] = s / AB[i][i];
	}
	return true;
}


void Matrix::zerujMacierze()
{
	for (int i = 0;i < myGrid->nH*myGrid->nL;i++)
	{
		for (int j = 0;j < myGrid->nH*myGrid->nL;j++)
		{
			myGrid->globalMatrixH[i][j] = 0;
			myGrid->globalMatrixC[i][j] = 0;

		}
		myGrid->globalMatrixP[i] = 0;
	}
}
