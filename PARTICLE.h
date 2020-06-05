#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "Eigen/Dense"
using namespace Eigen;

const float G = 0.4f; // Gravity constant 0.05
const float H = 0.4;      // Smoothing length
const float A = 8;  // Gas constant  2
const float EPS = 0.2;   // Gravity epsilon smoothing
const float pi = 3.14159265358979f;

class PARTICLE {
public:
	PARTICLE();
	PARTICLE(const Vector3d& position, const Vector3d& velocity, const float radius, const float mass, bool spin);
	
	// draw to OGL
	void draw();

	// add gravity force due to particle p
	void add_gravity(PARTICLE* p);
	
	// accessors
	Vector3d& position() 	 { return _position; };
	Vector3d& velocity() 	 { return _velocity; };
	Vector3d& halfvelocity() {return _halfvelocity;} ;
	Vector3d& force()	 	 { return _force; };
	float& radius()   { return _radius; };
	float& mass()     { return _mass; };
	float& density()  { return _density; };
	float& pressure() { return _pressure; };

private:  
	Vector3d _position;
	Vector3d _velocity;
	Vector3d _halfvelocity;
	Vector3d _force;
	float _radius;
	float _mass;
	float _density;
	float _pressure;
	float _h;
};

PARTICLE* add(PARTICLE* a, PARTICLE* b);

// Compute cubic spline kernel
float spline_kernel(float r, const float h);

// Compute gradient of cubic spline kernel
Vector3d grad_spline_kernel(Vector3d x, float r, const float h);