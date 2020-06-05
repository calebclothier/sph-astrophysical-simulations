#include "PARTICLE.h"

// Length of gravitational cutoff
const float THETA = 1.8;
// Artificial viscosity parameters
const float ALPHA = 0;
const float BETA = 2*ALPHA;

class QUAD {
public:
	QUAD();
	QUAD(Vector3d corner, float side);
	
	// Returns true if particle is within quadrant, false otherwise
	bool contains(PARTICLE* p);
	// Returns quadrant sidelength
	float length() { return side_length; };
	// Returns vector to center of quadrant
	Vector3d center();
	
	// Creates and returns sub-quadrants 
	QUAD TL();
	QUAD TR();
	QUAD BL();
	QUAD BR();

private:
	Vector3d BL_corner;
	float side_length;
};


class BH_TREE {
public:
	BH_TREE();
	BH_TREE(QUAD q);

	// Returns true if particle p lies within the quadrant q associated with tree
	bool contains(PARTICLE* p);

	// Inserts particle into Barnes-Hut tree
	void insert(PARTICLE* p);

	// Approximates the gravitational force acting on particle p
	void add_gravity_force(PARTICLE* p);

	// Add hydrodynamical force acting on particle p
	void add_hydro_force(PARTICLE* p);

	// Calculates and updates density at each particle
	void compute_density(PARTICLE* p);

private:
	PARTICLE* particle;
	QUAD quad;
	BH_TREE* TL;
	BH_TREE* TR;
	BH_TREE* BL;
	BH_TREE* BR;
};