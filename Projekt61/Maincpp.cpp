#include <iostream>
#include "Element.h"
#include "Node.h"
#include "GRID.h"
#include "Matrix.h"

using namespace std;

int main()

{
	GRID grid;
//	grid.showMatrixH();
	Matrix xD(&grid);

//	grid.showMatrixH();
	int choice;
	bool loop = true;
	while (loop)
	{
		cout << "_______________________________________________________" << endl;
		cout << "1. Wyswietl wszystkie wspolrzedne wezlow nalezacych do siatki." << endl;
		cout << "2. Wyswietl wszystkie numery wezlow elementow nalezacych do siatki." << endl;
		cout << "3. Wyswietl numery wezlow nalezacych  do danego elementu i ich wspolrzedne." << endl;
		cout << "4. Wyswietl wspolrzedne danego wezla." << endl;
		cout << "5. Wyswietl wszystko." << endl;
		cout << "6. Wyswietl macierz H." << endl;
		cout << "7. Wyswietl macierz C." << endl;
		cout << "8. Wyswietl macierz H2." << endl;
		cout << "9. Wyswietl sumarna lokalna macierz H." << endl;
		cout << "10. Wyswietl globalna macierz H." << endl;
		cout << "11. Wyswietl globalna macierz C." << endl;
		cout << "12. Wyswietl macierz h z daszkiem." << endl;
		cout << "13. Wyswietl macierz p z daszkiem." << endl;


		cout << "14. Wyjdz." << endl;

		cout << "Wybierz: ";
		cin >> choice;
		switch (choice)
		{
		case 1:
		{
			grid.showNodes();
			break;
		}

		case 2:
		{
			grid.showElements();
			break;
		}

		case 3:
		{
			cout << "Podaj numer elementu: ";
			int nr;
			cin >> nr;
			grid.showElementsNodes(nr);
			break;
		}

		case 4:
		{
			cout << "Podaj numer wezla: ";
			int nr;
			cin >> nr;
			grid.showNodesCords(nr);
			break;
		}

		case 5:
		{
			grid.showAll();
			break;
		}


		case 6:
		{
			grid.showMatrixH();
			break;
		}

		case 7:
		{
			grid.showMatrixC();
			break;
		}

		case 8:
		{
			grid.showMatrixH2();
			break;
		}

		case 9:
		{
			grid.showMatrixHlokal();
			break;
		}

		case 10:
		{
			grid.showGlobalMatrixH();
			break;
		}

		case 11:
		{
			grid.showGlobalMatrixC();
			break;
		}

		case 12:
		{
			grid.showMatrixHzDaszkiem();
			break;
		}

		case 13:
		{
			grid.showMatrixPzDaszkiem();
			break;
		}

		case 14:
		{
			loop = false;
			break;
		}

		default:
		{
			cout << "Zly wybor." << endl;
		}

		}
	}

}