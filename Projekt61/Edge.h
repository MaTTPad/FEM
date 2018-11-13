#pragma once
#include "Node.h"
class Edge
{
public:
	Node nodes[2];
	bool boundary;

	Edge();
	Edge(Node node1, Node node2);
	void setBoundary();
	bool isBoundary();
	~Edge();
};

