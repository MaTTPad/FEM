#pragma once
class Node
{
public:
	double x, y, t;
	bool boundary;

public:
	bool isBoundary();
	Node();
	Node(int);
	~Node();
};

