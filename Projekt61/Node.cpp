#include "Node.h"



Node::Node()
{
	t = 100;
}


Node::~Node()
{
}

bool Node::isBoundary() {
	return boundary;
}
