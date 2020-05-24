#pragma once
#include <string>
#include "Grid.h"

class Simulator {
private:
	int 			res;							// number of cells along each dimension
	Grid<double> 	p_grid, u_grid, v_grid;			// grid for pressure and velocity (p: grid centers, u,v: grid edges)
	std::vector<int> 	labels;						// label each cell as fluid (non-negative), solid (-1), or air (-2) 
	double			t;								// time
	double			dt;								// dt for each timestep
	double			frame_dt;						// dt for each frame
	double			dx;								// dx for grid spacing
	int 			num_frames;						// total number of frames
	int 			frame;							// current number of frame

private:
	void projectPressure();							// projection step
	void advectVelocity();							// advection step
	void applyBC();									// apply boundary condition
	void writeData();								// write particles data into file

public:
	Simulator();
	~Simulator();
	void init();									// initialize parameters
	void advanceOneTimeStep();						// execute a whole cycle in one time step
	double getTime();								// return current time
};