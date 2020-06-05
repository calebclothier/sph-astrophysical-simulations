#include "PARTICLE_SYSTEM.h"
#include <random>
#include <cmath>
#include <iostream>

using namespace std;
using namespace Eigen;

float total_mass = 1500;

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
PARTICLE_SYSTEM::PARTICLE_SYSTEM() 
{
	_particles.push_back(PARTICLE(Vector3d(0,0,0), Vector3d(0,0,0), _radius, _mass, 0));
	_tree = BH_TREE();
	_screen = QUAD(Vector3d(-16, -16, 0), 32.f);
	N = 1;
	_mass = 1.0;
	_radius = 1.0;
}

PARTICLE_SYSTEM::PARTICLE_SYSTEM(const int particle_num, const float radius, const float mass, bool spin, bool single_cloud) :
	N(particle_num), _radius(radius), _mass(mass)
{

	_mass = total_mass / N;

	// Initialize empty BH tree and screen quadrant
	_tree = BH_TREE();
	_screen = QUAD(Vector3d(-16, -16, 0), 32.f);

	// Initialize distribution function
	default_random_engine generator;

	// Normal single cloud distribution
    normal_distribution<double> distribution(0, 1.5);

    // Plummer distribiution
    uniform_real_distribution<double> phiDist(0, 2*pi);
    uniform_real_distribution<double> thetaDist(-1, 1);
    uniform_real_distribution<double> rDist(0, 3);

    // Two cloud distribution
    normal_distribution<double> leftCloud(-2, 1);
    normal_distribution<double> rightCloud(2, 1);
    uniform_real_distribution<double> noise(-1, 1);

    float x;
    float y;
    float v_x = 0;
    float v_y = 0;
    float r;
    int i = 0;

    while (i < (N-1)) {

    	if (single_cloud == true) {

	    	// Normal  distribution
	    	// x = distribution(generator);
	    	// y = distribution(generator);

	    	// Plummer distribution
	    	float phi = phiDist(generator);
	    	float theta = acos(thetaDist(generator));
	    	/*
	    	r = 1 / sqrt(pow(rDist(generator), -2.f/3.f) - 1);
	    	if (i < 5 * N / 6)
	    		r *= 3;
	    	else
	    		r *= 0.8;
	    	*/
	    	r = rDist(generator);
	    	x = r * sin(theta) * cos(phi);
	    	y = r * sin(theta) * sin(phi);

	    	if ((r > .8) && (r < 7)) {
    			_particles.push_back(PARTICLE(Vector3d(x,y,0), Vector3d(v_x,v_y,0), _radius, _mass, spin));
    			i++;
    		}
	    }

    	// Two cloud distribution
    	else {
	    	if ((i % 2) == 0) {

	    		float phi = phiDist(generator);
	    		float theta = acos(thetaDist(generator));
	    		r = rDist(generator);
	    		x = r * sin(theta) * cos(phi);
	    		y = r * sin(theta) * sin(phi);
	    		x += 2.3;
	    		y += 2.3;
	    		v_x = -5 + noise(generator);
	    		v_y = -2 + noise(generator);

	    	}
	    	else {

	    		float phi = phiDist(generator);
	    		float theta = acos(thetaDist(generator));
	    		r = rDist(generator);
	    		x = r * sin(theta) * cos(phi);
	    		y = r * sin(theta) * sin(phi);
	    		x -= 2.3;
	    		y -= 2.3;
	    		v_x = 5 + noise(generator);
	    		v_y = 2 + noise(generator);

	    	}
	    	_particles.push_back(PARTICLE(Vector3d(x,y,0), Vector3d(v_x,v_y,0), _radius, _mass, spin));
    		i++;
	    }
    }

    if (single_cloud == true)
    	_particles.push_back(PARTICLE(Vector3d(-0.01,-0.01,0), Vector3d(0,0,0), 0.15, 1000, spin));
    else
    	x = rightCloud(generator);
	    y = rightCloud(generator);
	    _particles.push_back(PARTICLE(Vector3d(x,y,0), Vector3d(0,0,0), _radius, _mass, spin));


    cout << "Particle number: " << _particles.size() << "\n";

}


///////////////////////////////////////////////////////////////////////////////
// Compute distance between two particles
///////////////////////////////////////////////////////////////////////////////
float PARTICLE_SYSTEM::distance(PARTICLE p1, PARTICLE p2) 
{
	return sqrt(pow(p1.position()[0] - p2.position()[0], 2) + pow(p1.position()[1] - p2.position()[1], 2) + pow(p1.position()[2] - p2.position()[2], 2));
}


///////////////////////////////////////////////////////////////////////////////
// Construct Barnes-Hut tree
///////////////////////////////////////////////////////////////////////////////
void PARTICLE_SYSTEM::construct_tree() {
	// Create new Barnes-Hut Tree
	_tree = BH_TREE(_screen);
	for (int i = 0; i < _particles.size(); i++) {
		if (_screen.contains(&_particles[i]))
			_tree.insert(&_particles[i]);
		_particles[i].force().setZero();
	}
}

///////////////////////////////////////////////////////////////////////////////
// Compute gravitational forces between particles
///////////////////////////////////////////////////////////////////////////////
void PARTICLE_SYSTEM::compute_gravity_forces() 
{
	// Add gravitational forces
	for (int i = 0; i < _particles.size(); i++) {
		if (_screen.contains(&_particles[i]))
			_tree.add_gravity_force(&_particles[i]);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Compute fluid forces on each particle
///////////////////////////////////////////////////////////////////////////////
void PARTICLE_SYSTEM::compute_fluid_forces() {

	// Compute densities and pressures
	for (int i = 0; i < _particles.size(); i++) {
		if (_screen.contains(&_particles[i])) {
			_particles[i].density() = 0;
			_tree.compute_density(&_particles[i]);
			_particles[i].pressure() = A * pow(_particles[i].density(), 4.f / 3.f);
		}
	}

	// Compute forces
	for (int i = 0; i < _particles.size(); i++) {
		if (_screen.contains(&_particles[i])) {
			_tree.add_hydro_force(&_particles[i]);
		}
	}

}


///////////////////////////////////////////////////////////////////////////////
// Update positions and velocities using leapfrog time integration
///////////////////////////////////////////////////////////////////////////////
void PARTICLE_SYSTEM::leapfrog_step(float dt, int firstStep) 
{

	if (firstStep) {
		for (int i = 0; i < _particles.size(); i++) {
			Vector3d a = _particles[i].force() / _particles[i].mass();
			_particles[i].halfvelocity() = _particles[i].velocity() + (a * (dt / 2.f));
			_particles[i].position() += _particles[i].halfvelocity() * dt;
		}
	}
	else {
		for (int i = 0; i < _particles.size(); i++) {

			/*
			if (i == _particles.size()) {
				break;
			}
			*/
			
			Vector3d a = _particles[i].force() / _particles[i].mass();
			Vector3d v_old = _particles[i].halfvelocity();
			_particles[i].halfvelocity() += (a * dt);
			_particles[i].position() += (_particles[i].halfvelocity() * dt);
			_particles[i].velocity() = (v_old + _particles[i].halfvelocity()) / 2.f;

			/*
			// Add to star if close and not moving
			if (((_particles[i].position() - _particles[_particles.size() - 1].position()).norm() < _particles[_particles.size() - 1].radius() / 1.3) && (_particles[i].velocity().norm() < 0.3)) {
				if (i != (_particles.size() - 1)) {

					Vector3d star_velocity = _particles[_particles.size() - 1].velocity();
					Vector3d star_position = _particles[_particles.size() - 1].position();
					Vector3d part_velocity = _particles[i].velocity();
					Vector3d part_position = _particles[i].position();
					float star_mass = _particles[_particles.size() - 1].mass();
					float part_mass = _particles[i].mass();

					_particles[_particles.size() - 1].radius() += 0.00005;
					_particles[_particles.size() - 1].mass() += _particles[i].mass();
					_particles[_particles.size() - 1].velocity() = (star_velocity * star_mass + part_velocity * part_mass) / (star_mass + part_mass);
					_particles[_particles.size() - 1].position() = (star_position * star_mass + part_position * part_mass) / (star_mass + part_mass);

					_particles.erase(_particles.begin() + i);
				}
			}
			*/
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
// Draw particles in OpenGL
///////////////////////////////////////////////////////////////////////////////
void PARTICLE_SYSTEM::draw() 
{
	for (int i = 0; i < _particles.size(); i++) {
    	_particles[i].draw();
    }
}

///////////////////////////////////////////////////////////////////////////////
// Draw density field in OpenGL
///////////////////////////////////////////////////////////////////////////////
void PARTICLE_SYSTEM::draw_density() {

	// First, an unclever brute force way to draw the density field
	int _res = 512;
	float h = 16.f / (_res - 1);
	int smoothing_distance = round(H / h) + 20;
	float m = 1 / h;
	
	MatrixXf density_field(_res, _res);
	density_field.setZero();

	float x;
	float y;
	int x_index;
	int y_index;
	float distance;
	float density;
	float max_density = 10;

	// Calculate density at each location in window by looping over particles
	for (int i = 0; i <= _particles.size(); i++) {

		// Get x and y coordinates of particle (x: -8 -> 8, y: -8 -> 8)
		x = _particles[i].position()[0];
		y = _particles[i].position()[1];

		x_index = round(m * (x - 8)) + (_res - 1);
		y_index = round(m * (y - 8)) + (_res - 1);

		// For all neighboring cells within the smoothing length of the particle, update density accordingly
		for (int i = -smoothing_distance; i <= smoothing_distance; i++) {
			for (int j = -smoothing_distance; j <= smoothing_distance; j++) {

				if (((x_index + i) > _res) || ((x_index + i) < 0)) {
					continue;
				}
				if (((y_index + j) > _res) || ((y_index + j) < 0)) {
					continue;
				}

				distance = sqrt(pow(i * h, 2) + pow(j * h, 2));
				density_field(x_index + i, y_index + j) += _particles[i].mass() * spline_kernel(distance, H);
			}
		}
	}
	
	// Draw density field

	glBegin(GL_POINTS);
	for (int i = 0; i < _res; i++) {
		for (int j = 0; j < _res; j++) {

			density = density_field(i, j) / max_density;
			glColor3f(density, density, density);

			x = (i - (_res - 1)) / m + 8.f;
			y = (j - (_res - 1)) / m + 8.f;
			glVertex2f(x, y);

		}
	}

	glEnd();

}










