#pragma once
#include <string>
#include "Grid.h"
//#include "ParticleSystem.h"
//#include "PoissonDisk.h"

class Simulator {
private:
	Grid<double> 	p_grid, u_grid, v_grid;			// grid for pressure and velocity (p: grid centers, u,v: grid edges)
	//ParticleSystem	pts_sys;						// particle system, stores particles for advection
	std::vector<int> 	labels;						// label each cell as fluid (non-negative), solid (-1), or air (-2) 
	double			t;								// time
	int 			frame;							// current number of frame
	double			rho;							// density of fluid

public:
	int				res;							// resolution of grid
	double			dx;								// grid dx width
	double			dt;								// time step length
	double			frame_dt;						// frame dt
	double			gravity;						// gravity
	double			end_t;							// end time
	int				num_frames;						// total number of frames

private:
	void projectPressure();							// solve pressure in projection equation
	void updateVelocity();							// update velocity with pressure
	void advection();								// semi-Lagrangian advection
	void applyBC();									// apply boundary condition
	void loadData();								// load data from file
	void writeData();								// write data into file

public:
	Simulator();
	~Simulator();
	void init();									// initialize parameters
	void advanceOneTimeStep();						// execute a whole cycle in one time step
	double getTime();								// return current time
};