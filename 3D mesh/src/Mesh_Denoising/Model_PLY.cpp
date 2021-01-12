/*
Author:		Rang M. H. Nguyen
Reference:	Rang M. H. Nguyen, Michael S. Brown
			Fast and Effective L0 Gradient Minimization by Region Fusion
			ICCV 2015
Date:		Dec 1st, 2015
*/
//-------------------------------

#include "Model_PLY.h"
#include <fstream>


Model_PLY::Model_PLY()
{
	maxSizeHashFace = 100;
}
 
 
float* Model_PLY::calculateNormal( float *coord1, float *coord2, float *coord3 )
{
   /* calculate Vector1 and Vector2 */
   float va[3], vb[3], vr[3], val;
   va[0] = coord1[0] - coord2[0];
   va[1] = coord1[1] - coord2[1];
   va[2] = coord1[2] - coord2[2];
 
   vb[0] = coord1[0] - coord3[0];
   vb[1] = coord1[1] - coord3[1];
   vb[2] = coord1[2] - coord3[2];
 
   /* cross product */
   vr[0] = va[1] * vb[2] - vb[1] * va[2];
   vr[1] = vb[0] * va[2] - va[0] * vb[2];
   vr[2] = va[0] * vb[1] - vb[0] * va[1];
 
   /* normalization factor */
   val = sqrt( vr[0]*vr[0] + vr[1]*vr[1] + vr[2]*vr[2] );
 
    float norm[3];
    norm[0] = vr[0]/val;
    norm[1] = vr[1]/val;
    norm[2] = vr[2]/val;
 
 
    return norm;
}
 
 
 
int Model_PLY::Load(char* filename)
{
    this->TotalConnectedTriangles = 0; 
    this->TotalConnectedQuads = 0;
    this->TotalConnectedPoints = 0;
 
    char* pch = strstr(filename,".ply");
 
    if (pch != NULL)
    {
       FILE* file = fopen(filename,"r");
	   fseek(file,0,SEEK_SET);        
 
       if (file)
       {
            int i = 0;   
            int temp = 0;
            int quads_index = 0;
            int triangle_index = 0;
            int normal_index = 0;
            char buffer[1000];
 
 
            fgets(buffer,300,file);            // ply
 
 
            // READ HEADER
            // -----------------
 
            // Find number of vertexes
            while (  strncmp( "element vertex", buffer,strlen("element vertex")) != 0  )
            {
                fgets(buffer,300,file);            // format
            }
            strcpy(buffer, buffer+strlen("element vertex"));
            sscanf(buffer,"%i", &this->TotalConnectedPoints);
 
 
            // Find number of vertexes
            fseek(file,0,SEEK_SET);
            while (  strncmp( "element face", buffer,strlen("element face")) != 0  )
            {
                fgets(buffer,300,file);            // format
            }
            strcpy(buffer, buffer+strlen("element face"));
            sscanf(buffer,"%i", &this->TotalFaces);
 
 
            // go to end_header
            while (  strncmp( "end_header", buffer,strlen("end_header")) != 0  )
            {
                fgets(buffer,300,file);            // format
            }
 
            //----------------------
			Vertex_Buffer = new float[3*TotalConnectedPoints];        
 
			Faces_Triangles = new float[9*TotalFaces];
			Normals  = new float[9*TotalFaces];
 
            // read verteces
            i =0;
            for (int iterator = 0; iterator < this->TotalConnectedPoints; iterator++)
            {
                fgets(buffer,300,file);
 
                sscanf(buffer,"%f %f %f", &Vertex_Buffer[i], &Vertex_Buffer[i+1], &Vertex_Buffer[i+2]);
                i += 3;
            }
 
			// Rang
			faces = new Face[this->TotalFaces];			
			hashFaces = new int*[this->TotalConnectedPoints];
			nHashFaces = new int[this->TotalConnectedPoints];
			for(int j = 0; j < this->TotalConnectedPoints; j++)
			{
				hashFaces[j] = new int[maxSizeHashFace];
				nHashFaces[j] = 0;
			}
			
            // read faces
            i =0;
            for (int iterator = 0; iterator < this->TotalFaces; iterator++)
            {
                fgets(buffer,300,file);
 
                    if (buffer[0] == '3')
                    {
 
                        int vertex1 = 0, vertex2 = 0, vertex3 = 0;
                        //sscanf(buffer,"%i%i%i\n", vertex1,vertex2,vertex3 );
                        buffer[0] = ' ';
                        sscanf(buffer,"%i%i%i", &vertex1,&vertex2,&vertex3 );
                        

                        //  vertex == punt van vertex lijst
                        // vertex_buffer -> xyz xyz xyz xyz
                        // printf("%f %f %f ", Vertex_Buffer[3*vertex1], Vertex_Buffer[3*vertex1+1], Vertex_Buffer[3*vertex1+2]);
 
              
 
                        float coord1[3] = { Vertex_Buffer[3*vertex1], Vertex_Buffer[3*vertex1+1],Vertex_Buffer[3*vertex1+2]};
                        float coord2[3] = {Vertex_Buffer[3*vertex2],Vertex_Buffer[3*vertex2+1],Vertex_Buffer[3*vertex2+2]};
                        float coord3[3] = {Vertex_Buffer[3*vertex3],Vertex_Buffer[3*vertex3+1],Vertex_Buffer[3*vertex3+2]};
                        float *norm = this->calculateNormal(coord1, coord2, coord3);
 

                        TotalConnectedTriangles += 3;

						// MR from here						
						// Search neighbour 						

						int vertices[3] = {vertex1, vertex2, vertex3};	
						faces[iterator].setNormal(norm[0], norm[1], norm[2]);
						faces[iterator].setVertex(vertices);	
						faces[iterator].setCenter(coord1, coord2, coord3);
						for(int k = 0; k < 3; k++)
						{	
							for(int iv = 0; iv < nHashFaces[vertices[k]]; iv++)
							{
								faces[iterator].insertNeighbour(hashFaces[vertices[k]][iv]);
								faces[hashFaces[vertices[k]][iv]].insertNeighbour(iterator);
							}
						}						
						
						for(int k = 0; k < 3; k++)
						{
							hashFaces[vertices[k]][nHashFaces[vertices[k]]] = iterator;
							nHashFaces[vertices[k]] += 1;
						}					

						// to here

                    }
 
 
                    i += 3;
            } 
			
			ScalingBox();
            fclose(file);
        }
 
      else { printf("File can't be opened\n"); }
    } else {
      printf("File does not have a .PLY extension. ");    
    }   
    return 0;
}

void Model_PLY::addNoise(float noise)
{
	float normal[3];
	for(int i = 0; i < TotalConnectedPoints; i++)
	{
		int ci = 3*i;
		float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		for(int k = 0; k < 3; k++)
			normal[k] = 0;
		for(int j = 0; j < nHashFaces[i]; j++)
		{
			for(int k = 0; k < 3; k++)
				normal[k] += faces[hashFaces[i][j]].normal[k];
		}
		r = noise*(r-0.5);
		if(nHashFaces[i] > 1)
		{
			for(int k = 0; k < 3; k++)
				Vertex_Buffer[ci+k] += r*normal[k]/nHashFaces[i];
		}
	}

	//for(int i = 0; i < 3*TotalConnectedPoints; i++)
	//{
	//	float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	//	r = noise*(r-0.5);
	//	Vertex_Buffer[i] += r*Vertex_Buffer[i];
	//}
	updateNormals();	
}

void Model_PLY::updateNormals()
{
	int triangle_index = 0;
	int normal_index = 0;
	for (int iterator = 0; iterator < this->TotalFaces; iterator++)
	{   
		int idx1 = 3 * faces[iterator].vertices[0];
		int idx2 = 3 * faces[iterator].vertices[1];
		int idx3 = 3 * faces[iterator].vertices[2];

		Faces_Triangles[triangle_index] = Vertex_Buffer[idx1];
		Faces_Triangles[triangle_index+1] = Vertex_Buffer[idx1+1];
		Faces_Triangles[triangle_index+2] = Vertex_Buffer[idx1+2];
		Faces_Triangles[triangle_index+3] = Vertex_Buffer[idx2];
		Faces_Triangles[triangle_index+4] = Vertex_Buffer[idx2+1];
		Faces_Triangles[triangle_index+5] = Vertex_Buffer[idx2+2];
		Faces_Triangles[triangle_index+6] = Vertex_Buffer[idx3];
		Faces_Triangles[triangle_index+7] = Vertex_Buffer[idx3+1];
		Faces_Triangles[triangle_index+8] = Vertex_Buffer[idx3+2];
 
		float coord1[3] = { Faces_Triangles[triangle_index], Faces_Triangles[triangle_index+1],Faces_Triangles[triangle_index+2]};
		float coord2[3] = {Faces_Triangles[triangle_index+3],Faces_Triangles[triangle_index+4],Faces_Triangles[triangle_index+5]};
		float coord3[3] = {Faces_Triangles[triangle_index+6],Faces_Triangles[triangle_index+7],Faces_Triangles[triangle_index+8]};
		float *norm = this->calculateNormal( coord1, coord2, coord3 );
 
		Normals[normal_index] = norm[0];
		Normals[normal_index+1] = norm[1];
		Normals[normal_index+2] = norm[2];
		Normals[normal_index+3] = norm[0];
		Normals[normal_index+4] = norm[1];
		Normals[normal_index+5] = norm[2];
		Normals[normal_index+6] = norm[0];
		Normals[normal_index+7] = norm[1];
		Normals[normal_index+8] = norm[2];
 
        normal_index += 9;
 
        triangle_index += 9;	
		faces[iterator].setNormal(norm[0], norm[1], norm[2]);
     }
}

void Model_PLY::updateVertices(int MAX_LOOP)
{

	float sumV3[3];
	for(int iter = 0; iter < MAX_LOOP; iter++)
	{
		for(int iv = 0; iv < TotalConnectedPoints; iv++)
		{
			int ci = 3*iv;
			for(int b = 0; b < 3; b++)
			{
				sumV3[b] = 0;
			}

			for(int k = 0; k < nHashFaces[iv]; k++)
			{
				float sum = 0;
				Face* faceTemp = &faces[hashFaces[iv][k]];
				float center[3];
				int va = 3*faceTemp->vertices[0];
				int vb = 3*faceTemp->vertices[1];
				int vc = 3*faceTemp->vertices[2];
				// compute the center
				for (int b=0; b<3; b++)
				{
					center[b] = (Vertex_Buffer[va+b] + Vertex_Buffer[vb+b] + Vertex_Buffer[vc+b])/3;
				}
				// compute the dot product
				for(int b =0; b < 3; b++)
				{					
					sum += faceTemp->normal[b]*(center[b] - Vertex_Buffer[ci+b]);
				}

				for(int b=0; b<3; b++)
				{
					sumV3[b] += faceTemp->normal[b]*sum;
				}
			}
			if(nHashFaces[iv] > 0)
			{
				for(int b = 0; b < 3; b++)
				{
					Vertex_Buffer[ci+b] += sumV3[b]/nHashFaces[iv]; 
				}
			}
			
		}
		cout<<"Iter: "<<iter<<endl;
	}
}

void Model_PLY::Save(char* filename)
{
	//FILE* file = fopen(filename,"w");
	ofstream outputFile;
	outputFile.open(filename, std::ofstream::out);
	if(outputFile.is_open())
	{
		outputFile<< "ply\n";
		outputFile<< "format ascii 1.0\n";
		outputFile<< "comment zipper output\n";
		outputFile<< "element vertex " <<TotalConnectedPoints<<endl;
		outputFile<< "property float x\n";
		outputFile<< "property float y\n";
		outputFile<< "property float z\n";
		outputFile<< "element face "<<TotalFaces<<endl;;
		outputFile<< "property list uchar int vertex_indices\n";
		outputFile<< "end_header\n";

		for (int iv = 0; iv < TotalConnectedPoints; iv++)
		{
			int ci = 3*iv;
			outputFile<<Vertex_Buffer[ci]<<" "<<Vertex_Buffer[ci+1]<<" "<<Vertex_Buffer[ci+2]<<endl;
		}

		for(int ifa = 0; ifa < TotalFaces; ifa++)
		{
			outputFile<<3<<" "<<faces[ifa].vertices[0]<<" "<<faces[ifa].vertices[1]<<" "<<faces[ifa].vertices[2]<<endl;
		}
		outputFile<<std::flush;
		outputFile.close();
		cout<<"File is already saved."<<endl;
	}
	else
	{
		cout<<"File cannot be opened for writing."<<endl;
	}
}


void Model_PLY::Draw()
{
	updateNormals();
    glEnableClientState(GL_VERTEX_ARRAY);    
    glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3,GL_FLOAT,    0,Faces_Triangles);    
    glNormalPointer(GL_FLOAT, 0, Normals);
    glDrawArrays(GL_TRIANGLES, 0, TotalConnectedTriangles);    
    glDisableClientState(GL_VERTEX_ARRAY);    
    glDisableClientState(GL_NORMAL_ARRAY);
}

Model_PLY::~Model_PLY(void)
{
	delete[] Faces_Triangles;
    delete[] Faces_Quads;
    delete[] Vertex_Buffer;
    delete[] Normals;
	delete[] faces; 
	delete[] nHashFaces;
	for(int j = 0; j < this->TotalConnectedPoints; j++)
	{
		delete[] hashFaces[j];
	}
	delete[] hashFaces;		
}

bool Model_PLY::deNoise(float t, int maxSize, int maxLoop)
{
	
	int M = TotalFaces;
	int bands = 3;
	double total = 0;
	double time;
	double step = (double)1/maxLoop;
	double ct = 0;
	double curThresh;
	clock_t start, end;

	float* Y = new float[bands*M];
	float* sumY = new float[bands*M];
	
	for(int i = 0; i<M; i++)
	{	
		int ci = i * bands;
		for(int j = 0; j < bands; j++)
		{
			sumY[ci+j] = Y[ci+j] = faces[i].normal[j];
		}
		
	}
	

	LinkedList* G = new LinkedList[M];
	int* NB = new int[M*maxSize];
	int* nNB = new int[M];
	int* maxNB = new int[M];
	int* W = new int[M];
	int* IDX = new int[M];
	int* lutIDX = new int[M];
	int* hashID = new int[M];
	bool* pagefault = new bool[M];
	int** startpage = new int*[M];

	// Init data
	for(int i = 0; i < M; i++)
	{
		G[i].insert(i);
		IDX[i] = i;
		lutIDX[i] = i;
		W[i] = 1;
		maxNB[i] = maxSize;
		hashID[i] = -1;
		pagefault[i] = false;
		startpage[i] = &NB[i*maxSize];
	}

	// Init neighbourship
	for(int i = 0; i < M; i++)
	{
		int curI = maxSize * i;		
		nNB[i] = 0;
		for(int j = 0; j < faces[i].neighbours.size(); j++)
		{
			int faceID = faces[i].neighbours[j];
			int found;
			for(found = 0; found < nNB[i]; found = found + 2)
			{
				if(NB[curI+found] == faceID) break;
			}
			if(found < nNB[i])
			{
				NB[curI+found + 1] += 3;
			}
			else
			{
				NB[curI+nNB[i]] = faceID;	
				NB[curI+nNB[i]+1] = 1;
				nNB[i] += 2;
			}
		}
	}

	// All other runs	
	int iter = 0, inc = -1;
	while(ct <= 1 && M > 1)
	{
		int maxNBnum = 0;
		inc = -1;		
		double curThresh = pow(ct,2.2)*t;
		//double curThresh = ct*t;
		ct +=step;
		for(int i = 0; i < M; i++)
		{
			int idx1 = IDX[i];				
			int* startI1 = startpage[idx1];

			int FIX_LOOP_TIMES = nNB[idx1];			
						
			// create hashIDX			
			for(int hi = 0; hi < nNB[idx1]; hi = hi+2)
			{
				hashID[*(startI1+hi)] = hi;
			}
			

			for(int j = 0; j < FIX_LOOP_TIMES; j=j+2)
			{				
				int idx2 = *(startI1+j);		
				int* startI2 = startpage[idx2];
				int rIdx1 = bands*idx1;
				int rIdx2 = bands*idx2;
				
				int len = *(startI1+j+1);					
				int sumW = W[idx1]+W[idx2];

				
				float dx = Y[rIdx1  ] - Y[rIdx2  ];
				float dy = Y[rIdx1+1] - Y[rIdx2+1];
				float dz = Y[rIdx1+2] - Y[rIdx2+2];
				double d = dx*dx + dy*dy + dz*dz;
				
				if(d*W[idx1]*W[idx2] < curThresh*len*sumW)
				{			
					
					// Join and erase mean set
					double tempSum = 0;
					double tempV3[3];
					//CHANGE normalize the normal (April 15)
					for(int b = 0; b < bands; b++)
					{						
						tempV3[b] = W[idx1]*Y[rIdx1+b]+W[idx2]*Y[rIdx2+b];
						tempSum += tempV3[b]*tempV3[b];
					}		

					//START HERE: Normalize the normal vector
					tempSum = sqrt(tempSum);
					for(int b = 0; b < bands; b++)
					{
						Y[rIdx1+b] = tempV3[b]/tempSum;
					}
					//END HERE

					// Join and erase weigh set
					W[idx1] = sumW;
									
					G[idx1].append(G[idx2]);	
					
					// Erase idx2 from idx1
					if(j != nNB[idx1]-2)
					{						
						swap(*(startI1+j), *(startI1+nNB[idx1]-2));
						swap(*(startI1+j+1), *(startI1+nNB[idx1]-1));
						hashID[*(startI1+j)] = j;
						j = j - 2;
					}

					nNB[idx1] -= 2;		
					if(nNB[idx1] < FIX_LOOP_TIMES)
						FIX_LOOP_TIMES -= 2;	

					////Update hashID
					hashID[idx2] = -1;				
					
					// INSERT					
					for(int t =0 ; t<nNB[idx2]; t = t+2)
					{
						int aa = *(startI2+t);
						int la = *(startI2+t+1);
						int* startIa = startpage[aa];
						if (aa == idx1)
							continue;
						
						// Check if aa is the common neighbor of idx1 and idx2
						int find = hashID[aa];						
						if (find > -1) // If YES
						{
							*(startI1+find+1) += la;							
							int k = 0;
							while(*(startIa+k) != idx1)
							{
								k += 2;
							}							
							// Update the new connection number
							*(startIa+k+1) = *(startI1+find+1);	

							// Erase idx2 from its neighbors
							k = 0;
							while(*(startIa+k) != idx2)
							{
								k += 2;
							}	
							if(k != nNB[aa]-2)
							{
								swap(*(startIa+k  ), *(startIa+nNB[aa]-2));
								swap(*(startIa+k+1), *(startIa+nNB[aa]-1));								
							}
							nNB[aa] -= 2;	
						}	
						else
						{
							if(nNB[idx1] >= maxNB[idx1])
							{
								if(pagefault[idx1] || maxNB[idx1+maxNB[idx1]/maxSize] != 0)
								{
									// PAGE FAULT									
									int*temp = new int [2*maxNB[idx1]];
									maxNB[idx1] = 2*maxNB[idx1];
									// Copy to new page
									for(int ii = 0; ii < nNB[idx1]; ii++)
									{
										temp[ii] = *(startI1+ii);
									}
									if(pagefault[idx1])
									{
										delete[] startpage[idx1];
									}
									startpage[idx1] = &temp[0];
									startI1 = startpage[idx1];
									pagefault[idx1] = true;	
									//cout<<"Page fault!\n";
								}
								else
								{
									maxNB[idx1] += maxSize;
								}
							}				
							
							*(startI1+nNB[idx1])   = aa;
							*(startI1+nNB[idx1]+1) = la;
							//Update hashID
							hashID[aa] = nNB[idx1];
							nNB[idx1] += 2;
							
							int k = 0;
							while(*(startIa+k) != idx2)
							{
								k += 2;
							}							
							
							*(startIa+k) = idx1;
						}
					}					

					// DELETE!	
					maxNB[idx2] = 0;
					nNB[idx2] = 0;
					if(pagefault[idx2])
						delete[] startpage[idx2];
					M = M - 1;					
					int p_idx2 = lutIDX[idx2];
					lutIDX[IDX[M]] = p_idx2;
					swap(IDX[p_idx2], IDX[M]);						
				}				
				
			}
			if (nNB[idx1] > maxNBnum)
				maxNBnum = nNB[idx1];			
			// clean hashIDX
			for(int hi = 0; hi < nNB[idx1]; hi = hi+2)
			{
				hashID[*(startI1+hi)] = -1;
			}
		}
		
		iter++;
		cout<<"Iteration: "<<iter<<", M = "<<M<<", MaxNbNum = "<<maxNBnum<<endl;
	}
	// restore image
	//start = clock();
	for(int i = 0; i < M; i++)
	{	
		int idx1 = IDX[i];
		if(pagefault[idx1])
		{
			delete[] startpage[idx1];
		}

		int ridx = bands*idx1;
		Node* temp = G[idx1].pHead;
		while(temp != NULL)
		{
			int lidx = temp->value;			
			for(int b = 0; b < bands; b++)
				faces[lidx].normal[b] = Y[ridx+b];
			temp = temp->next;
		}
	}
	// DELETE
	delete[] Y;
	delete [] NB;
	delete [] nNB;		
	delete[] W;
	delete[] sumY;
	delete[] IDX;
	delete[] lutIDX;
	delete[] hashID;
	delete[] G;
	delete[] pagefault;
	delete[] startpage;
	return true;
}


int Model_PLY::LoadObj(char* filename)
{
    this->TotalConnectedTriangles = 0; 
    this->TotalConnectedQuads = 0;
    this->TotalConnectedPoints = 0;
 
    char* pch = strstr(filename,".obj");
 
    if (pch != NULL)
    {
       FILE* file = fopen(filename,"r");
	}
	return 0;
}

void Model_PLY::ScalingBox(void)
{
    int i,j;
    float box[2][3];

    box[0][0] = box[0][1] = box[0][2] = FLT_MAX;
    box[1][0] = box[1][1] = box[1][2] = -FLT_MAX;
	for (i=0;i<TotalConnectedPoints;i++)
    {
		int ci = 3*i;
        for (j=0;j<3;j++)
        {
			if (box[0][j]>Vertex_Buffer[ci+j])
                box[0][j] = Vertex_Buffer[ci+j];
            if (box[1][j]< Vertex_Buffer[ci+j])
                box[1][j] = Vertex_Buffer[ci+j];
        }
    }
    f3ObjCentre[0] = (box[0][0]+box[1][0])/2.0;
    f3ObjCentre[1] = (box[0][1]+box[1][1])/2.0;
    f3ObjCentre[2] = (box[0][2]+box[1][2])/2.0;

    fObjScale = box[1][0]-box[0][0];	
	for (i = 1; i < 3; i ++)
	{
		if(fObjScale < box[1][i] - box[0][i])
			fObjScale = box[1][i] - box[0][i];
	}
    fObjScale /=2.0;
    for (i=0;i<TotalConnectedPoints;i++)
    {
		int ci = 3*i;
		for(j=0; j < 3; j++)
		{
			Vertex_Buffer[ci+j] = (Vertex_Buffer[ci+j] - f3ObjCentre[j])/fObjScale;
		}       
    }
}