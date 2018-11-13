#include "Element.h"



Element::Element()
{
	k = 30;
	ro = 7800;
	c = 700;
	alfa = 25;
	for (int j = 0;j < 4;j++)
	{
		for (int k = 0;k < 4;k++)
		{
			matrixC[j][k]=0;
			H[j][k] = 0;
		}
	}
}


Element::~Element()
{
}
