#include "Node.h"



Node::Node()
{
}

Node::Node(int initialTemp)
{
	t = initialTemp;
}


Node::~Node()
{
}

bool Node::isBoundary() {
	return boundary;
}
