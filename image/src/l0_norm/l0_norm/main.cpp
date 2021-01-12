/*
Author:		Rang M. H. Nguyen
Reference:	Rang M. H. Nguyen, Michael S. Brown
			Fast and Effective L0 Gradient Minimization by Region Fusion
			ICCV 2015
Date:		Dec 1st, 2015
*/
//-------------------------------

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <time.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include "LinkedList.h"
#include <string>
#include <fstream>

using namespace cv;
using namespace std;

Mat l0_norm(Mat img, double t, int maxSize, int maxLoop);
double objective_function(Mat& img1, Mat& img2, double t);

int main(int argc, char* argv[])
{
	if (argc < 5) 
	{ // Check the value of argc. If not enough parameters have been passed, inform user and exit.
        cout << "Usage is -i <infile> -t <lambda> -o <outfile>[Option]\n"; // Inform the user of how to use the program
        return 0;
    } 
	else 
	{ 
        string filein = "";
		string fileout = "";
		double t = 0;
        for (int i = 1; i < argc; i++) {
            if (i + 1 != argc) // Check that we haven't finished parsing already
                if (strcmp(argv[i],"-i") == 0) {
                    filein = string(argv[i + 1]);
					i++;
                } else if (strcmp(argv[i], "-t") == 0) {
                    t = stod(argv[i + 1]);
					i++;
                } else if (strcmp(argv[i], "-o") == 0) {
                    fileout = string(argv[i + 1]);
					i++;
                } else {
                    cout << "Not enough or invalid arguments, please try again.\n";
                    return 0;
            }
        }
		ifstream testfile(filein);
		if(!testfile.good())
		{
			cout <<  "Could not find or open the image\n";
			return -1;
		}
		Mat img = imread(filein, CV_LOAD_IMAGE_COLOR);   // Read the file
		Mat img1 = img.clone();

		double tc = t*255*255;
		if(!img.data )                              // Check for invalid input
		{
			cout <<  "Could not open or find the image\n";
			return -1;
		}
		clock_t start = clock();
		img = l0_norm(img, tc, 32, 100);
		clock_t end = clock();
		double time = (double) (end-start) / CLOCKS_PER_SEC;	
		cout<<"Elapsed time: "<<time<<endl; 
		cout<<"Objective function = "<<objective_function(img1, img, t)<<endl;

		if(fileout == "")
		{
			std::ostringstream ss;
			ss << t;
			int found=filein.find('.');
			fileout = filein.substr(0,found) + "_l0_" + ss.str() + ".png"; 
		}
		imwrite(fileout, img);
	
		waitKey(0);                                          // Wait for a keystroke in the window
		return 0;
	}
}


void createNeighbour(int rows, int cols, int maxSize, int*& NB, int*& nNB)
{
	int inc = 0;

	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < cols; j++)
		{			
			int curI = inc*maxSize-1;
			int count = 0;
			if(j+1 < cols)
			{
				NB[++curI] = i*cols+j+1;	
				NB[++curI] = 1;
				count += 2;
			}
			if(i+1 < rows)
			{
				NB[++curI] = (i+1)*cols+j;
				NB[++curI] = 1;
				count += 2;
			}
			if(j > 0)
			{
				NB[++curI] = i*cols+j-1;	
				NB[++curI] = 1;
				count += 2;
			}
			if(i > 0)
			{
				NB[++curI] = (i-1)*cols+j;	
				NB[++curI] = 1;
				count += 2;
			}
			nNB[inc] = count;
			inc++;
		}
	}
}


Mat l0_norm(Mat img, double t, int maxSize, int maxLoop)
{
	int rows = img.rows, cols = img.cols, bands = img.channels();
	int M = rows * cols;
	double total = 0;
	clock_t start, end;
	double time;
	double step = (double)1/maxLoop;
	double ct = 0;
	double curThresh;

	float* Y = new float[bands*M];
	float* sumY = new float[bands*M];
	for(int i = 0; i<bands*M; i++)
	{
		sumY[i] = Y[i] = img.data[i];
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
	
	int iter = 0, inc = -1;
	createNeighbour(rows, cols, maxSize, NB, nNB);


	// All other runs	
	while(iter <= maxLoop && M > 1)
	{
		int maxNBnum = 0;
		inc = -1;		
		//double curThresh = pow(1.5,iter-maxLoop)*t;
		//double curThresh = pow(ct,2.2)*t;
		double curThresh = ct*t;
		ct +=step;
		for(int i = 0; i < M; i++)
		{
			int idx1 = IDX[i];				
			int* startI1 = startpage[idx1];

			int FIX_LOOP_TIMES = nNB[idx1];		
			
			//			
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
				
				float dx = Y[rIdx1  ] - Y[rIdx2  ];
				float dy = Y[rIdx1+1] - Y[rIdx2+1];
				float dz = Y[rIdx1+2] - Y[rIdx2+2];
				double d = dx*dx + dy*dy + dz*dz;
				
				int sumW = W[idx1]+W[idx2];
				if(d*W[idx1]*W[idx2] <= curThresh*len*sumW) // Eq. 12 in the paper
				{			
					
					// Join and erase mean set
					for(int b = 0; b < bands; b++)
					{
						sumY[rIdx1+b] += sumY[rIdx2+b];
						Y[rIdx1+b] = sumY[rIdx1+b]/sumW;
					}					
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
					
					
					// Insert neighbour of idx2 into idx1					
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
						else // If NO
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

		if(iter == 0) cout<<"Original number of groups = "<<M<<endl;
		iter++;

	}
	cout<<"Final number of groups = "<<M<<endl;
	// restore image
	for(int i = 0; i < M; i++)
	{	
		int idx1 = IDX[i];
		if(pagefault[idx1])
		{
			delete[] startpage[idx1];
		}
		Node* temp = G[idx1].pHead;
		while(temp != NULL)
		{
			int lidx = bands*temp->value;
			int ridx = bands*idx1;
			for(int b = 0; b < bands; b++)
				img.data[lidx+b] = Y[ridx+b];
			temp = temp->next;
		}
	}
	// clear the memory
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


	return img;
}


double  objective_function(Mat& img1, Mat& img2, double t)
{
	double f = 0;
	int rows = img1.rows;
	int cols = img1.cols;
	int bands = img1.channels();
	float maxValue = 255*255;
	int curI = 0;
	for(int i = 0; i < rows; i ++)
	{
		for(int j = 0; j < cols; j++)
		{			
			float dx = img1.data[curI] - img2.data[curI];
			float dy = img1.data[curI+1] - img2.data[curI+1];
			float dz = img1.data[curI+2] - img2.data[curI+2];

			f += (dx*dx + dy*dy + dz*dz)/maxValue;
			if(j < cols-1)
			{
				int nextJ = curI+3;
				if(img2.data[curI] != img2.data[nextJ] || img2.data[curI+1] != img2.data[nextJ+1] || img2.data[curI+2] != img2.data[nextJ+2])
				{
					f += t;
				}
			}

			if(i < rows-1)
			{
				int nextI = curI+3*cols;
				if(img2.data[curI] != img2.data[nextI] || img2.data[curI+1] != img2.data[nextI+1] || img2.data[curI+2] != img2.data[nextI+2])
				{
					f += t;
				}
			}
			curI = curI + 3;
		}
	}

	return f;
}