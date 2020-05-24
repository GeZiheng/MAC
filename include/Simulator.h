#pragma once
#include <string>
#include "Grid.h"

class Simulator {
private:
	Grid<double> 	p_grid, u_grid, v_grid;			// grid for pressure and velocity (p: grid centers, u,v: grid edges)
	std::vector<int> 	labels;						// label each cell as fluid (non-negative), solid (-1), or air (-2) 
	double			t;								// time
	int 			frame;							// current number of frame

public:
	int				res;							// resolution of grid
	double			dx;								// grid dx width
	double			dt;								// time step length
	double			frame_dt;						// frame dt
	double			gravity;						// gravity
	double			end_t;							// end time
	int				num_frames;						// total number of frames

private:
	void projectPressure();							// projection step
	void advectVelocity();							// advection step
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