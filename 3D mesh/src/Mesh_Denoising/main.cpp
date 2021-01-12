/*
Author:		Rang M. H. Nguyen
Reference:	Rang M. H. Nguyen, Michael S. Brown
			Fast and Effective L0 Gradient Minimization by Region Fusion
			ICCV 2015
Date:		Dec 1st, 2015
*/
//-------------------------------
#include "Model_PLY.h"
#include <windows.h>  // for MS Windows
#include <GL/glut.h>  // GLUT, include glu.h and gl.h
#include "Model_PLY.h"
#include <iostream>
#include <string>
#include <fstream>

using namespace std;
 
/* Global variables */
char title[] = "3D Shapes";

Model_PLY obj;

//light 0
GLfloat light0_ambient[] = {1, 1 , 1, 1.0};
GLfloat light0_position[] = {0, 2.0, 0.0, 1.0 };


// Gold
GLfloat silver_ambient[] = {0.24725, 0.1995, 0.0745, 1.0};
GLfloat silver_diffuse[] = {0.75164, 0.60648, 0.22648, 1.0};
GLfloat silver_specular[] = {0.628281	, 0.555802, 0.366065, 1.0};
GLfloat	silver_shiness = 0.6;

string filein = "";
string fileout = "";
double lambda = 0;
int id = 0;
/* Initialize OpenGL Graphics */
void initGL() {
   glClearColor(1.0f, 1.0f, 1.0f, 0.0f); 
   glClearDepth(1.0f);               
   glEnable(GL_DEPTH_TEST);  
   glDepthFunc(GL_LEQUAL);    
   //glShadeModel(GL_SMOOTH);   
   glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  


   // Insert
   	glEnable(GL_NORMALIZE);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);

	glEnable(GL_LIGHT0);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
	glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
	//glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0.2); 
	//glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.5);

	// Material 
	glMaterialfv(GL_FRONT, GL_AMBIENT, silver_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, silver_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, silver_specular);
	glMaterialf(GL_FRONT, GL_SHININESS, silver_shiness * 128.0);
	
	clock_t start, end;
	double time;
	obj.Load(const_cast<char*>(filein.c_str()));
	//obj.addNoise(0.03);
	//obj.Save("diamond_03.ply");
	start = clock();
	obj.deNoise(lambda, 128, 100);
	end = clock();
	time = (double) (end-start) / CLOCKS_PER_SEC;
	cout<<"Elapsed time denoise: "<<time<<endl; 
	obj.updateVertices(100);
	obj.Save(const_cast<char*>(fileout.c_str()));
}
 

void display() {
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   glMatrixMode(GL_MODELVIEW);     
   glLoadIdentity();                 
   // gluLookAt(2,2,-2,0,0,0,-1, -1, 0); for fandisk
   gluLookAt(1, 1.5, 2, 0, 0, 0, 0, 1, 0);
   obj.Draw();     
   glutSwapBuffers();  
}

void reshape(GLsizei width, GLsizei height) 
{  
	float rate = 2;
	if (height == 0) height = 1;                
	GLfloat aspect = (GLfloat)width / (GLfloat)height;
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);  
	glLoadIdentity();            
	if (width >= height) 
	{   
		glOrtho(-rate * aspect, rate * aspect, -rate, rate, 0.1, 100);
	} 
	else 
	{
		glOrtho(-rate, rate, -rate / aspect, rate / aspect, 0.1, 100);
	}
}

void keyPressed (unsigned char key, int x, int y) 
{  
	if(key == 13)
	{
		glutDestroyWindow (id);
		exit(0);
	}
	glutPostRedisplay();
}

/* Main function: GLUT runs as a console application starting at main() */
int main(int argc, char** argv) {
	if (argc < 7) 
	{ // Check the value of argc. If not enough parameters have been passed, inform user and exit.
        cout << "Usage is -i <infile> -t <lambda> -o <outfile>\n"; // Inform the user of how to use the program
        return 0;
    } 
	else 
	{ 
        for (int i = 1; i < argc; i++) {
            if (i + 1 != argc) // Check that we haven't finished parsing already
                if (strcmp(argv[i],"-i") == 0) {
                    filein = string(argv[i + 1]);
					i++;
                } else if (strcmp(argv[i], "-t") == 0) {
                    lambda = stod(argv[i + 1]);
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
			cout <<  "Could not find or open the mesh file\n";
			return -1;
		}
	}

   glutInit(&argc, argv);            
   glutInitDisplayMode(GLUT_DOUBLE); 
   glutInitWindowSize(640, 480);  
   glutInitWindowPosition(50, 50); 
   id = glutCreateWindow(title);          
   glutDisplayFunc(display);       
   glutReshapeFunc(reshape);    
   glutKeyboardFunc(keyPressed);
   initGL();                       
   glutMainLoop();                 
   return 0;
}