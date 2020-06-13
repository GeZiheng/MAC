#pragma once
#include <iostream>
#include <array>
#include <vector>
#include <Eigen>

class ParticleSystem {					// sparse grid for pressure and velocity
public:
	int							particle_num;		// number of particles
	std::vector<double>			positions;			// storage for grid positions
	std::vector<double>			velocities;			// storage for grid velocities

public:
	ParticleSystem() {
		particle_num = 0;
	}
	~ParticleSystem() {
		positions.clear();
		velocities.clear();
	}
};