/*
Author:		Rang M. H. Nguyen
Reference:	Rang M. H. Nguyen, Michael S. Brown
			Fast and Effective L0 Gradient Minimization by Region Fusion
			ICCV 2015
Date:		Dec 1st, 2015
*/
//-------------------------------

#pragma once
#include <windows.h>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <string>
#include "Face.h"
#include <iostream>
#include <GL/gl.h>
#include <GL/glu.h>
#include "LinkedList.h"
#include <time.h>
 
using namespace std;
class Model_PLY 
{
public:
    int Model_PLY::Load(char *filename);
    void Model_PLY::Save(char*);
    float* Model_PLY::calculateNormal( float *coord1, float *coord2, float *coord3 );
	void Draw();
	void updateVertices(int );
	void addNoise(float noise);
	bool deNoise(float t, int maxSize, int maxLoop);
	int LoadObj(char *filename);
	void ScalingBox(void);

	void updateNormals();

    Model_PLY();
	~Model_PLY(void);

 
    float* Faces_Triangles;
    float* Faces_Quads;
    float* Vertex_Buffer;
    float* Normals;
	Face* faces; 
	int** hashFaces; 
	int* nHashFaces;
	int maxSizeHashFace;
    
	int TotalConnectedTriangles;    
    int TotalConnectedQuads;    
    int TotalConnectedPoints;
    int TotalFaces;
	float f3ObjCentre[3];
	float fObjScale; 
};

