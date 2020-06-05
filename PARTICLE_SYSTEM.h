#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "BH_TREE.h"
#include "Eigen/Dense"
#include <vector>
using namespace std;
using namespace Eigen;

class PARTICLE_SYSTEM {
public:
	PARTICLE_SYSTEM();
	PARTICLE_SYSTEM(const int particle_num, const float radius, const float mass, bool spin, bool single_cloud);

	// Compute radial distance between two particles
	float distance(PARTICLE p1, PARTICLE p2);

	// Construct Barnes-Hut tree
	void construct_tree();

	// Compute gravitational force on each particle
	void compute_gravity_forces();

	// Compute fluid forces on each particle
	void compute_fluid_forces();

	// Leapfrog integration for each particle
	void leapfrog_step(float dt, int firstStep); 

	// Draw all particles
	void draw();

	// Draw density field
	void draw_density();

	// accessors
	vector<PARTICLE>& particles() { return _particles; };
	float& mass() 				  { return _mass; };
	float& radius()			      { return _radius; };
	int& size()					  { return N; };

private:
	vector<PARTICLE> _particles;
	BH_TREE _tree;
	QUAD _screen;
	int N;
	float _radius;
	float _mass;
};