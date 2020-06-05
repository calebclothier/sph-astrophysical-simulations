#include "PARTICLE.h"
#include <cmath>
#include <random>
#include <iostream>

using namespace Eigen;
using namespace std;

// Sign of float for velocity field initialization
float sgn(float v) {
	if (v < 0) return -1.f;
	if (v > 0) return 1.f;
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
PARTICLE::PARTICLE() 
{
	_position.setZero();
	_velocity.setZero();
	_halfvelocity.setZero();
	_force.setZero();
	_radius = 5;
	_mass = 1.0;
	_density = 0;
	_pressure = 0;
}

PARTICLE::PARTICLE(const Vector3d& position, const Vector3d& velocity, const float radius, const float mass, bool spin) :
	_position(position), _velocity(velocity), _radius(radius), _mass(mass)
{

	// Velocity spin test initialization
	if (spin == true) {
		_velocity.setZero();
		default_random_engine generator;
	    normal_distribution<double> distribution(0.10, 0.01);
	    float x = distribution(generator);
	    float y = distribution(generator);
		_velocity[0] = - x * sgn(_position[1]) * sqrt(abs(_position[1])) * exp(-abs(_position[1]) / 5);
		_velocity[1] = y * sgn(_position[0]) * sqrt(abs(_position[0])) * exp(-abs(_position[0]) / 5);
	}


	_halfvelocity.setZero();
	_force.setZero();
	_density = 0;
	_pressure = 0;
}

///////////////////////////////////////////////////////////////////////////////
// OGL drawing
///////////////////////////////////////////////////////////////////////////////
void PARTICLE::draw() 
{
	glPushMatrix();
	glTranslatef(_position[0], _position[1], _position[2]);
	float energy = 0.05 * log(_density + 1);
	if (energy > 1)
		glColor3f(1,1,1); 
	else
		glColor3f(0.8, energy, 0.1);
	if (_radius < 0.1)
		glutSolidSphere(_radius, 5, 5);
	else
		glutSolidSphere(_radius, 20, 20);
	glPopMatrix();
}

///////////////////////////////////////////////////////////////////////////////
// Add gravity to particle forces
///////////////////////////////////////////////////////////////////////////////
void PARTICLE::add_gravity(PARTICLE* p) 
{
	Vector3d separation = _position - p->position();
	_force -= (G * _mass * p->mass() / sqrt(pow(separation.squaredNorm() + EPS, 3))) * separation;
}

///////////////////////////////////////////////////////////////////////////////
// Create new COM particle for Barnes Hut tree
///////////////////////////////////////////////////////////////////////////////
PARTICLE* add(PARTICLE* a, PARTICLE* b)
{
	PARTICLE* com_particle = new PARTICLE();
	if ((a == NULL) && (b == NULL)) {
		return NULL;
	}
	else if (a == NULL) {
		com_particle->mass() = b->mass();
		com_particle->position() = b->position();
		return com_particle;
	}
	else if (b == NULL) {
		com_particle->mass() = a->mass();
		com_particle->position() = a->position();
		return com_particle;
	}
	else {
		com_particle->mass() = a->mass() + b->mass();
		com_particle->position() = (a->mass() * a->position() + b->mass() * b->position()) / com_particle->mass();
		return com_particle;
	}
}

///////////////////////////////////////////////////////////////////////////////
// Cubic spline kernel for interpolation
///////////////////////////////////////////////////////////////////////////////
float spline_kernel(float r, const float h) 
{
	float q = r / h;
	float s = 10.f / (7.f * pi * h * h);
	if (q > 2)
		return 0;
	else if (q > 1)
		return (s / 4.f) * pow(2-q, 3);
	else
		return s * (1 - 1.5 * q * q * (1 - 0.5*q));
}

Vector3d grad_spline_kernel(Vector3d x, float r, const float h)
{
	float q = r / h;
	float s = 10.f / (7.f * pi * h * h);
	if (r > 1e-12) {
		if (q > 2)
			return 0 * x;
		else if (q > 1)
			return (-0.75 * s * (2-q) * (2-q) / (r * h)) * x;
		else 
			return (-3.0 * s * q * (1-0.75*q) / (r * h)) * x;
	}
	else
		return 0 * x;
}



