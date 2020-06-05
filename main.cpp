#include <cstdlib>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "PARTICLE_SYSTEM.h"

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
using namespace Eigen;
using namespace std;

// Initialize constants
int N = 2500;		// Number of particles
const float r = 0.02;	// Radius of each particle
const float m = 3.0;	// Mass of each particles
const float dt = 0.001;   // Timestep

bool animate = true;
bool spin = false;
bool single_cloud = true;

// Particle system variable
PARTICLE_SYSTEM particleSystem;

///////////////////////////////////////////////////////////////////////////////
// The drawing function
///////////////////////////////////////////////////////////////////////////////
void displayCallback()
{
	// clear away the previous frame
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	
	// draw the particle system
	particleSystem.draw();
	
	// swap the buffers
	glutSwapBuffers();
}

///////////////////////////////////////////////////////////////////////////////
// The projection function
///////////////////////////////////////////////////////////////////////////////
void reshapeCallback(int w, int h)
{
	// set the viewport resolution (w x h)
	glViewport(0, 0, (GLsizei) w, (GLsizei) h);
	
	// Make ensuing transforms affect the projection matrix
	glMatrixMode(GL_PROJECTION);
	
	// set the projection matrix to an orthographic view
	glLoadIdentity();
	glOrtho(-8, 8, -8, 8, -8, 8);
	
	// set the matric mode back to modelview
	glMatrixMode(GL_MODELVIEW);
	
	// set the lookat transform
	glLoadIdentity();
	gluLookAt(0, 0, 1, 0, 0, 0, 0, 1, 0);
}

///////////////////////////////////////////////////////////////////////////////
// Keyboard handling function
///////////////////////////////////////////////////////////////////////////////
void keyboardCallback(unsigned char key, int x, int y)
{
	/* this is the keyboard event handler
	   the x and y parameters are the mouse 
	   coordintes when the key was struck */
	switch (key)
	{
	case 'q':
	case 'Q':
		exit(0);
		break;
	case ' ':
		animate = !animate;
		break;
	}

	glutPostRedisplay();
}

///////////////////////////////////////////////////////////////////////////////
// Idle command processing function
///////////////////////////////////////////////////////////////////////////////
void idleCallback()
{
	if (!animate)
		return;

	particleSystem.construct_tree();
	particleSystem.compute_gravity_forces();
	particleSystem.compute_fluid_forces();
	particleSystem.leapfrog_step(dt, 0);

	glutPostRedisplay();
}

///////////////////////////////////////////////////////////////////////
// process the command line
///////////////////////////////////////////////////////////////////////
string toUpper(const string& input)
{
	string copy(input);
	for (unsigned int x = 0; x < input.length(); x++)
		copy[x] = std::toupper(copy[x]);
	return copy;
}

void readCommandLine(int argc, char** argv)
{
	if (argc == 2) {
		N = stoi(argv[1]);
	}
	else if (argc == 3) {
		string flag2 = toUpper(argv[2]);
		N = stoi(argv[1]);
		if (!flag2.compare("-SPIN")) {
			spin = true;
		}
		else if (!flag2.compare("-DOUBLE"))
			single_cloud = false;
	}
	else if (argc == 4) {
		N = stoi(argv[1]);
		spin = true;
		single_cloud = false;
	}
}

int main(int argc, char **argv)
{
	readCommandLine(argc, argv);

	// initialize GLUT
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE| GLUT_RGBA);
	glutInitWindowSize(1000, 1000); 
	glutInitWindowPosition(100, 100);
	glutCreateWindow(argv[0]);
	glShadeModel(GL_SMOOTH);
	
	// point GLUT to our callback functions
	glutDisplayFunc(displayCallback); 
	glutIdleFunc(idleCallback); 
	glutReshapeFunc(reshapeCallback);
	glutKeyboardFunc(keyboardCallback);
	
	// set background to black
	glClearColor(0.0, 0.0, 0.0, 0.0);

	// Initialize particle system
	particleSystem = PARTICLE_SYSTEM(N, r, m, spin, single_cloud);
	// Take first step
	particleSystem.construct_tree();
	particleSystem.compute_gravity_forces();
	particleSystem.compute_fluid_forces();
	particleSystem.leapfrog_step(dt, 1);

	glutMainLoop();

	return 0;
}