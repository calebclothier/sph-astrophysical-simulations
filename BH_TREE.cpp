#include "BH_TREE.h"
#include <iostream>

using namespace std;

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
QUAD::QUAD() {
	BL_corner = Vector3d(0,0,0);
	side_length = 1;
}
QUAD::QUAD(Vector3d corner, float side) : BL_corner(corner), side_length(side) {}

///////////////////////////////////////////////////////////////////////////////
// Function to determine whether a particle falls within the quadrant
///////////////////////////////////////////////////////////////////////////////
bool QUAD::contains(PARTICLE* p)
{
	float x = p->position()[0];
	float y = p->position()[1];
	if (((BL_corner[0] <= x) && (x <= BL_corner[0] + side_length)) &&
	    ((BL_corner[1] <= y) && (y <= BL_corner[1] + side_length)))
	    return true;
	else
		return false;
}

///////////////////////////////////////////////////////////////////////////////
// Returns vector to center of quadrant
///////////////////////////////////////////////////////////////////////////////
Vector3d QUAD::center() {
	return Vector3d(BL_corner[0] + side_length / 2.f, BL_corner[1] + side_length / 2.f, 0);
}

///////////////////////////////////////////////////////////////////////////////
// Functions to create and return sub-quadrants
///////////////////////////////////////////////////////////////////////////////
QUAD QUAD::TL()
{
	Vector3d corner;
	corner[0] = BL_corner[0];
	corner[1] = BL_corner[1] + (side_length / 2.f);
	corner[2] = BL_corner[2];
	return QUAD(corner, side_length / 2.f);
}

QUAD QUAD::TR()
{
	Vector3d corner;
	corner[0] = BL_corner[0] + (side_length / 2.f);
	corner[1] = BL_corner[1] + (side_length / 2.f);
	corner[2] = BL_corner[2];
	return QUAD(corner, side_length / 2.f);
}

QUAD QUAD::BL()
{
	return QUAD(BL_corner, side_length / 2.f);
}

QUAD QUAD::BR()
{
	Vector3d corner;
	corner[0] = BL_corner[0] + (side_length / 2.f);
	corner[1] = BL_corner[1];
	corner[2] = BL_corner[2];
	return QUAD(corner, side_length / 2.f);
}

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
BH_TREE::BH_TREE() 
{
	quad = QUAD();
	particle = NULL;
	TL = NULL;
	TR = NULL;
	BL = NULL;
	BR = NULL;
}
BH_TREE::BH_TREE(QUAD q) : quad(q)
{
	// Assign particle and sub-quadrant pointers to null
	particle = NULL;
	TL = NULL;
	TR = NULL;
	BL = NULL;
	BR = NULL;
}

///////////////////////////////////////////////////////////////////////////////
// Returns true if particle p lies within the quadrant q associated with tree
///////////////////////////////////////////////////////////////////////////////
bool BH_TREE::contains(PARTICLE* p) 
{
	return quad.contains(p);
}

///////////////////////////////////////////////////////////////////////////////
// Insert particle into tree
///////////////////////////////////////////////////////////////////////////////
void BH_TREE::insert(PARTICLE* p)
{
	// If no particle in this node, add particle
	if (particle == NULL) {
		particle = p;
	}
	// If internal node, update its particle and insert new particle recursively
	else if (TL || TR || BL || BR) {
		particle = add(particle, p);
		if ((TL->quad).contains(p))
			TL->insert(p);
		else if ((TR->quad).contains(p))
			TR->insert(p);
		else if ((BL->quad).contains(p))
			BL->insert(p);
		else
			BR->insert(p);
	}
	// External node case
	else {
		// Create new sub-trees 
		QUAD TLq = quad.TL();
		QUAD TRq = quad.TR();
		QUAD BLq = quad.BL();
		QUAD BRq = quad.BR();

		TL = new BH_TREE(TLq);
		TR = new BH_TREE(TRq);
		BL = new BH_TREE(BLq);
		BR = new BH_TREE(BRq);

		// Recursively insert particles into sub-trees
		if (TL->contains(p))
			TL->insert(p);
		else if (TR->contains(p))
			TR->insert(p);
		else if (BL->contains(p))
			BL->insert(p);
		else if (BR->contains(p))
			BR->insert(p);

		if (TL->contains(particle))
			TL->insert(particle);
		else if (TR->contains(particle))
			TR->insert(particle);
		else if (BL->contains(particle))
			BL->insert(particle);
		else if (BR->contains(particle))
			BR->insert(particle);

		// Update particle at node
		particle = add(add(TL->particle, TR->particle), add(BL->particle, BR->particle));

	}
}

///////////////////////////////////////////////////////////////////////////////
// Add gravitational forces to particle
///////////////////////////////////////////////////////////////////////////////
void BH_TREE::add_gravity_force(PARTICLE* p)
{
	// If external node and not particle p
	if (TL == NULL) {
		if ((p != particle) && (particle != NULL)) {
			p->add_gravity(particle);
		}
	}
	// If internal node and s/d < theta
	else if (quad.length() / (particle->position() - p->position()).norm() < THETA) {
		p->add_gravity(particle);
	}
	// Else run recursively on sub-trees
	else {
		TL->add_gravity_force(p);
		TR->add_gravity_force(p);
		BL->add_gravity_force(p);
		BR->add_gravity_force(p);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Add hydrodynamical forces to particle
///////////////////////////////////////////////////////////////////////////////
void BH_TREE::add_hydro_force(PARTICLE* p)
{
	// If leaf node, add force due to particle
	if (TL == NULL) {
		
		if (particle != NULL) {
			
			Vector3d rij = p->position() - particle->position();
			float r2 = rij.squaredNorm();
			float r = sqrt(r2);

			if ((r > H) || p == particle) {
				;
			}
			else {

				// Physical pressure term
				float pressure_term = (particle->pressure() / pow(particle->density(), 2)) + (p->pressure() / pow(p->density(), 2));

				// Artificial viscosity
				float wij = rij.dot(p->velocity() - particle->velocity());
				float II;
				if (wij >= 0) {
					float mu = H * wij / (r2 + 0.01*H*H);
					II = (-ALPHA*mu + BETA*mu*mu) / (0.5*(p->density() + particle->density()));
				}
				else {
					II = 0;
				}

				if ((p->mass() < 10) && (particle->mass() < 10)) {
					Vector3d force = p->mass() * particle->mass() * (pressure_term + II) * grad_spline_kernel(rij, r, H);
					p->force() -= force;

					/*
					if (force.norm() < 5)
						p->force() -= force;
					else {
						p->force() -= 5 * force / force.norm();
					}
					*/


				}

			}
		}
	}
	else {

		if (((TL->quad.center() - p->position()).norm() - H) < TL->quad.length() / sqrt(2)) {
			TL->add_hydro_force(p);
		}
		if (((TR->quad.center() - p->position()).norm() - H) < TR->quad.length() / sqrt(2)) {
			TR->add_hydro_force(p);
		}
		if (((BL->quad.center() - p->position()).norm() - H) < BL->quad.length() / sqrt(2)) {
			BL->add_hydro_force(p);
		}
		if (((BR->quad.center() - p->position()).norm() - H) < BR->quad.length() / sqrt(2)) {
			BR->add_hydro_force(p);
		}

	}
}

///////////////////////////////////////////////////////////////////////////////
// Compute density and pressure at each particle
///////////////////////////////////////////////////////////////////////////////
void BH_TREE::compute_density(PARTICLE* p)
{

	// If leaf node, update density and pressure
	if (TL == NULL) {
		if (particle != NULL) {
			float radius = (particle->position() - p->position()).norm();
			p->density() += particle->mass() * spline_kernel(radius, H);
		}
	}
	else {
		if (((TL->quad.center() - p->position()).norm() - H) < TL->quad.length() / sqrt(2)) {
			TL->compute_density(p);
		}
		if (((TR->quad.center() - p->position()).norm() - H) < TR->quad.length() / sqrt(2)) {
			TR->compute_density(p);
		}
		if (((BL->quad.center() - p->position()).norm() - H) < BL->quad.length() / sqrt(2)) {
			BL->compute_density(p);
		}
		if (((BR->quad.center() - p->position()).norm() - H) < BR->quad.length() / sqrt(2)) {
			BR->compute_density(p);
		}
	}
}






