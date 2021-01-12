/*
Author:		Rang M. H. Nguyen
Reference:	Rang M. H. Nguyen, Michael S. Brown
			Fast and Effective L0 Gradient Minimization by Region Fusion
			ICCV 2015
Date:		Dec 1st, 2015
*/
//-------------------------------

#include "Face.h"


Face::Face(void)
{
	normal = new float[3]; 
	vertices = new int[3];
	center = new float[3];
}

Face::Face(float norm[3], int ver[3])
{
	for(int i = 0; i < 3; i++)
	{
		this->normal[i] = norm[i];
		this->vertices[i] = ver[i];
	}
}

void Face::setNormal(float n1, float n2, float n3)
{
	normal[0] = n1; normal[1] = n2; normal[2] = n3;
}

void Face::setVertex(int* ver)
{
	for(int i = 0; i < 3; i++)
	{
		vertices[i] = ver[i];
	}
}

void Face::insertNeighbour(int idx)
{
	neighbours.push_back(idx);
}


void Face::setCenter(float* v1, float* v2, float* v3)
{
	this->center = new float[3];
	for(int k = 0; k < 3; k++)
	{
		center[k] = (v1[k] + v2[k] + v3[k])/3;
	}
}

Face::~Face(void)
{
	delete[] normal;
	delete[] vertices;
	delete[] center;
}
