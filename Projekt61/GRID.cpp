#include "GRID.h"


GRID::GRID()
{
	//wczytanie z pliku
	fstream plik;
	plik.open("MES.txt", ios::in);

	if (plik.good())
	{
			plik >> H; // wysokoœæ
			plik >> L; //d³ugoœæ
			plik >> nH; // liczba wêz³ów po wysokoœci
			plik >> nL; // liczba wêz³ów po d³ugoœci
		

		int nodesCount = nH*nL;
		int elementsCount = (nH - 1)*(nL - 1);

		globalMatrixH = new double*[nodesCount];
		globalMatrixC = new double*[nodesCount];
		MatrixHzDaszkiem = new double*[nodesCount];
		globalMatrixP = new double[nodesCount];
		MatrixPzDaszkiem = new double[nodesCount];
		for (int i = 0;i < nodesCount;i++)
		{
			globalMatrixH[i] = new double[nodesCount];
			globalMatrixC[i] = new double[nodesCount];
			MatrixHzDaszkiem[i] = new double[nodesCount];

			globalMatrixP[i] = 0;
			MatrixPzDaszkiem[i] = 0;
		}

		for (int i = 0;i < nodesCount;i++)
			for (int j = 0;j < nodesCount;j++)
			{
				globalMatrixH[i][j] = 0;
				globalMatrixC[i][j] = 0;
			}



		ksiTab[0] = -1 / sqrt(3);
		etaTab[0] = ksiTab[0];
		ksiTab[1] = -ksiTab[0];
		etaTab[1] = etaTab[0];
		ksiTab[2] = ksiTab[1];
		etaTab[2] = -etaTab[1];
		ksiTab[3] = -ksiTab[2];
		etaTab[3] = etaTab[2];

		cout << "Siatka stworzona!" << endl;
		cout << "H=" << H << endl;
		cout << "L=" << L << endl;
		cout << "nH=" << nH << endl;
		cout << "nL=" << nL << endl;

		for (int i = 0;i < 4;i++)
		{
			cout << "Ksi" << i << "=" << ksiTab[i] << endl;
			cout << "Eta" << i << "=" << etaTab[i] << endl;
		}

		cout << "_______________________________________________________" << endl;

		nodes = new Node[nodesCount]; //wêz³y numerowane od 0
		elements = new Element[elementsCount];
		double x = 0;
		double y = 0;


		/////////////////////////////////////////////////////////////uzupe³nianie wspó³rzêdnych wêz³ów
		int licznik = 0;
		for (int k = 0;k < nodesCount;k++)
		{
			nodes[k].x = x;
			nodes[k].y = y;
			if (x == 0 || x == L || y == 0 || y == H) {
				nodes[k].boundary=true;
			}
			else
				nodes[k].boundary = false;
			y += (H / (nH - 1));
			licznik++;

			if (licznik >= nH)
			{
				y = 0;
				x += (L / (nL - 1));
				licznik = 0;
			}

		}



		//for (int i = 0;i < nodesCount;i++)
		//	{
		//	cout << "Wezel numer " << i << ": (" << nodes[i].x << "," << nodes[i].y << ")" << endl;
		//}
		///////////////////////////////////////////////////////////////////////////// wêz³y nale¿¹ce do elementu

		double val = 0;
		licznik = 0;
		for (int k = 0;k < elementsCount;k++)
		{
			elements[k].ID[0] = val;
			elements[k].ID[1] = val + nH;
			elements[k].ID[2] = val + nH + 1;
			elements[k].ID[3] = val + 1;

			elements[k].edges[0].nodes[0] = nodes[(int)elements[k].ID[0]];
			elements[k].edges[0].nodes[1] = nodes[(int)elements[k].ID[1]];
			elements[k].edges[0].setBoundary();
			elements[k].edges[1].nodes[0] = nodes[(int)elements[k].ID[1]];
			elements[k].edges[1].nodes[1] = nodes[(int)elements[k].ID[2]];
			elements[k].edges[1].setBoundary();

			elements[k].edges[2].nodes[0] = nodes[(int)elements[k].ID[2]];
			elements[k].edges[2].nodes[1] = nodes[(int)elements[k].ID[3]];
			elements[k].edges[2].setBoundary();

			elements[k].edges[3].nodes[0] = nodes[(int)elements[k].ID[3]];
			elements[k].edges[3].nodes[1] = nodes[(int)elements[k].ID[0]];
			elements[k].edges[3].setBoundary();
/*
			cout << "Edges element: " <<k<< endl;
			cout << elements[k].edges[0].isBoundary() << " ";
			cout << elements[k].edges[1].isBoundary() << " ";
			cout << elements[k].edges[2].isBoundary() << " ";
			cout << elements[k].edges[3].isBoundary() << " ";

			cout << endl << endl;
	*/		
			licznik++;
			val++;
			if (licznik == (nH - 1))
			{
				val++;
				licznik = 0;
			}
		}


	}



	else
	{
		cout << "Blad otwarcia pliku." << endl;
	}
}



GRID::~GRID()
{
}


void GRID::showElements()
{
	double elementsCount = (nH - 1)*(nL - 1);

	for (int i = 0;i < elementsCount;i++)
	{
		cout << "Element numer " << i << ": " << elements[i].ID[0] << " " << elements[i].ID[1] << " " << elements[i].ID[2] << " " << elements[i].ID[3] << endl;
	}
}

void GRID::showElementsNodes(int id)
{
	double elementsCount = (nH - 1)*(nL - 1);
	if (id < elementsCount)
	{
		cout << "Wezly nalezace do elementu numer " << id << ": " << elements[id].ID[0] << " " << elements[id].ID[1] << " " << elements[id].ID[2] << " " << elements[id].ID[3] << endl;
		cout << "Ich wspolrzedne:" << endl;
		cout << "Wezel nr " << elements[id].ID[0] << "; x=" << nodes[(int)elements[id].ID[0]].x << "; y=" << nodes[(int)elements[id].ID[0]].y << endl;
		cout << "Wezel nr " << elements[id].ID[1] << "; x=" << nodes[(int)elements[id].ID[1]].x << "; y=" << nodes[(int)elements[id].ID[1]].y << endl;
		cout << "Wezel nr " << elements[id].ID[2] << "; x=" << nodes[(int)elements[id].ID[2]].x << "; y=" << nodes[(int)elements[id].ID[2]].y << endl;
		cout << "Wezel nr " << elements[id].ID[3] << "; x=" << nodes[(int)elements[id].ID[3]].x << "; y=" << nodes[(int)elements[id].ID[3]].y << endl;

	}
	else
		cout << "Nie ma elementu o takim numerze." << endl;
}


void GRID::showNodesCords(int id)
{
	double elementsCount = (nH)*(nL);
	if (id < elementsCount)
		cout << "Wspolrzedne wezla numer " << id << ": (" << nodes[id].x << "," << nodes[id].y << ")" << endl;
	else
		cout << "Nie ma wezla o takim numerze." << endl;
}

void GRID::showNodes()
{
	double nodesCount = nH*nL;

	for (int i = 0;i < nodesCount;i++)
	{
		cout << "Wezel numer " << i << ": (" << nodes[i].x << "," << nodes[i].y << ")" << endl;
	}
}



void GRID::showAll()
{
	double nodesCount = nH*nL;
	double elementsCount = (nH - 1)*(nL - 1);

	for (int i = 0;i < nodesCount;i++)
	{
		cout << "Wezel numer " << i << ": " << nodes[i].x << " " << nodes[i].y << endl;
	}

	for (int i = 0;i < elementsCount;i++)
	{
		cout << "Element numer " << i << ": " << elements[i].ID[0] << " " << elements[i].ID[1] << " " << elements[i].ID[2] << " " << elements[i].ID[3] << endl;
	}
}

void GRID::showMatrixH()
{
	for (int i = 0;i < ((nH - 1)*(nL - 1));i++) {
		cout << "Macierz H elementu numer: " << i << endl;
		for (int k = 0;k < 4;k++) {
			for (int j = 0;j < 4;j++)
			{
				cout << elements[i].H[k][j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
}

void GRID::showMatrixH2()
{
	for (int i = 0;i < ((nH - 1)*(nL - 1));i++) {
		cout << "Macierz H2 elementu numer: " << i << endl;
		for (int k = 0;k < 4;k++) {
			for (int j = 0;j < 4;j++)
			{
				cout << elements[i].H2[k][j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
}

void GRID::showMatrixC()
{
	for (int i = 0;i < ((nH - 1)*(nL - 1));i++) {
		cout << "Macierz C elementu numer: " << i << endl;
		for (int k = 0;k < 4;k++) {
			for (int j = 0;j < 4;j++)
			{
				cout << elements[i].matrixC[k][j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
}

void GRID::showMatrixHlokal()
{
	for (int i = 0;i < ((nH - 1)*(nL - 1));i++) {
		cout << "Macierz lokalna H elementu numer: " << i << endl;
		for (int k = 0;k < 4;k++) {
			for (int j = 0;j < 4;j++)
			{
				cout << elements[i].matrixH[k][j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
}


void GRID::calculateGlobalMatrixH()
{

	for (int i = 0;i < (nH-1)*(nL-1) ;i++)
	{
		for (int j = 0;j < 4;j++)
			for (int k = 0;k < 4;k++) {
				int index1 = elements[i].ID[j];
				int index2 = elements[i].ID[k];

				globalMatrixH[index1] [index2] +=elements[i].matrixH[j][k];
			}
	}
}

void GRID::calculateGlobalMatrixC()
{

	for (int i = 0;i < (nH - 1)*(nL - 1);i++)
	{
		for (int j = 0;j < 4;j++)
			for (int k = 0;k < 4;k++) {
				int index1 = elements[i].ID[j];
				int index2 = elements[i].ID[k];

				globalMatrixC[index1][index2] += elements[i].matrixC[j][k];
			}
	}
}

void GRID::calculateMatrixHzDaszkiem()
{
	double dtau = 50;
	double **noweC = new double*[nH*nL];
	for (int i = 0;i < nH*nL;i++)
		noweC[i] = new double[nH*nL];

	for(int i=0;i<nH*nL;i++)
		for (int j = 0;j < nH*nL;j++)
		{
			noweC[i][j] = globalMatrixC[i][j] / dtau;
		}

	for (int i = 0;i < (nH)*(nL);i++)
	{
		for (int j = 0;j < (nH)*(nL);j++)
		{
			MatrixHzDaszkiem[i][j] = globalMatrixH[i][j] + noweC[i][j];
		}
	}
}


void GRID:: showGlobalMatrixH()
{

		//calculateGlobalMatrixH();
		cout << "Macierz globalna H:" << endl;
		for (int k = 0;k < nH*nL;k++) {
			for (int j = 0;j < nH*nL;j++)
			{
				cout << globalMatrixH[k][j] << " ";
			}
			cout << endl;
		}
}

void GRID::showGlobalMatrixC()
{

	//calculateGlobalMatrixC();
	cout << "Macierz globalna C:" << endl;
	for (int k = 0;k < nH*nL;k++) {
		for (int j = 0;j < nH*nL;j++)
		{
			cout << globalMatrixC[k][j] << " ";
		}
		cout << endl;
	}
}

void GRID::showMatrixHzDaszkiem()
{

	//calculateMatrixHzDaszkiem();
	cout << "Macierz h z daszkiem:" << endl;
	for (int k = 0;k < nH*nL;k++) {
		for (int j = 0;j < nH*nL;j++)
		{
			cout << MatrixHzDaszkiem[k][j] << " ";
		}
		cout << endl;
	}
}

void GRID::showMatrixPzDaszkiem()
{

	cout << "Macierz p z daszkiem:" << endl;
	for (int k = 0;k < nH*nL;k++) {
		cout << MatrixPzDaszkiem[k] << " ";
		cout << endl;
	}
}