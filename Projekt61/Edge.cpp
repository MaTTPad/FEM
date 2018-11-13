#include "Edge.h"



Edge::Edge()
{

}

Edge::Edge(Node node1, Node node2)
{
	nodes[0] = node1;
	nodes[1] = node2;
	setBoundary();
}

void Edge::setBoundary()
{
	if (nodes[0].isBoundary() && nodes[1].isBoundary()) {
		boundary = true;
	}
	else {
		boundary = false;
	}
}

bool Edge::isBoundary() {
	return boundary;
}


Edge::~Edge()
{
}
