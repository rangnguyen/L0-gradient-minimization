/*
Author:		Rang M. H. Nguyen
Reference:	Rang M. H. Nguyen, Michael S. Brown
			Fast and Effective L0 Gradient Minimization by Region Fusion
			ICCV 2015
Date:		Dec 1st, 2015
*/
//-------------------------------

#pragma once
#include<vector>

using namespace std;
class Face
{
public:
	float* normal;
	int* vertices;
	float* center;

	vector<int> neighbours;

	Face(void);
	Face(float normal[3], int vertices[3]);
	void setNormal(float n1, float n2, float n3);
	void setVertex(int* ver);
	void setCenter(float* v1, float* v2, float* v3);
	void insertNeighbour(int idx);
	~Face(void);
};

