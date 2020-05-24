#pragma once
#include <Eigen>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "Simulator.h"
#include "Parameters.h"

using namespace Eigen;

Simulator::Simulator() {
}

Simulator::~Simulator() {
}

void Simulator::projectPressure() {
	int N = p_grid.active_cells;
	VectorXd p_vec(N), rhs(N);
	SparseMatrix<double> A(N, N);
	Index grid_index;
	int i, j, q;
	/* Build A and rhs for linear system */
	for(int p = 0; p < N; p++) {
		/* Calculate A entry for current node */
		A(p,p) = 4 * dt / (rho * dx * dx);
		/* Find grid index for current node */
		grid_index = p_grid.id_list(p);
		i = grid_index[0];
		j = grid_index[1];
		/* Calculate A entries for each neighbor nodes */
		q = labels(p_grid.flatIndex({i+1, j}));
		A(p, q) = - dt / (rho * dx * dx);
		q = labels(p_grid.flatIndex({i-1, j}));
		A(p, q) = - dt / (rho * dx * dx);
		q = labels(p_grid.flatIndex({i, j+1}));
		A(p, q) = - dt / (rho * dx * dx);
		q = labels(p_grid.flatIndex({i, j-1}));
		A(p, q) = - dt / (rho * dx * dx);
		/* Calculate rhs for each neighbor nodes */
		q = labels(p_grid.flatIndex({i+1, j}));
		rhs(p) -= u_grid.data(q) / dx;
		q = labels(p_grid.flatIndex({i-1, j}));
		rhs(p) += u_grid.data(q) / dx;
		q = labels(p_grid.flatIndex({i, j+1}));
		rhs(p) += v_grid.data(q) / dx;
		q = labels(p_grid.flatIndex({i, j-1}));
		rhs(p) += v_grid.data(q) / dx;
	}
	/* Solve linear system using conjugate gradient method */
	ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
	cg.compute(A);
	p_vec = cg.solve(rhs);
	/* Load data into pressure grid */
	p_grid.loadFromVector(N, p_vec);
}

void Simulator::advectVelocity() {
	int N = p_grid.active_cells;					// number of fluid cells
	Index grid_index;
	int i, j, q;
	/* Update velocity (u,v) grids */
	for(int p = 0; p < N; p++) {
		/* Find grid index for current node */
		grid_index = p_grid.id_list(p);
		i = grid_index[0];
		j = grid_index[1];
		/* update u grid nodes */
		q = labels(p_grid.flatIndex({i+1, j}));
		u_grid.data[q] = u_grid.data[q] - dt / (rho * dx) * (p_grid.data[q] - p_grid.data[p]);
		/* update v grid nodes */
		q = labels(p_grid.flatIndex({i, j+1}));
		v_grid.data[q] = v_grid.data[q] - dt / (rho * dx) * (p_grid.data[q] - p_grid.data[p]);
	}
}

void Simulator::loadData() {
	std::ifstream myfile;
	std::string s, str;
	int i;
	myfile.open("../config.txt");

	// get res and dx from file
	str.clear();
	getline(myfile, s);
	for (i = 11; i < s.length(); i++)
		str.push_back(s[i]);
	res = stoi(str);
	dx = 1.0 / res;

	// get dt from file
	str.clear();
	getline(myfile, s);
	for (i = 8; i < s.length(); i++)
		str.push_back(s[i]);
	dt = stod(str);

	// get frame dt from file
	str.clear();
	getline(myfile, s);
	for (i = 11; i < s.length(); i++)
		str.push_back(s[i]);
	int frame_rate = stoi(str);
	frame_num = 0;
	frame_dt = 1.0 / double(frame_rate);

	// get gravity from file
	str.clear();
	getline(myfile, s);
	for (i = 8; i < s.length(); i++)
		str.push_back(s[i]);
	gravity = stod(str);

	// get end_t from file
	str.clear();
	getline(myfile, s);
	for (i = 6; i < s.length(); i++)
		str.push_back(s[i]);
	end_t = stod(str);

	// get num_frames from file
	str.clear();
	getline(myfile, s);
	for (i = 6; i < s.length(); i++)
		str.push_back(s[i]);
	end_t = stod(str);

	myfile.close();
}

void Simulator::writeData() {
	std::ofstream myfile;
	std::ostringstream filename;
	int N = p_grid.active_cells;					// number of fluid cells
	Index grid_index;
	int i, j, q;
	filename << "../output/frame" << std::setfill('0') << std::setw(4) << frame_num << ".dat";
	myfile.open(filename.str());
	myfile << "variables=\"x\",\"y\",\"p\",\"u\",\"v\"" << std:: endl;
	for (int p = 0; p < N; p++) {
		grid_index = p_grid.id_list(p);
		i = grid_index[0];
		j = grid_index[1];
		myfile << std::scientific << std::setprecision(16) << p_grid.origin[0] + i * p_grid.dx << ' ' << p_grid.origin[1] + j * p_grid.dx << ' ' << p_grid.data[p] << ' ' << u_grid.data[p] << ' ' << v_grid.data[p] << std::endl;
	}
	myfile.close();
}

void Simulator::init() {
	loadData();
	t = 0;
	frame = 1;
	dx = double(1) / double(res);
}

void Simulator::advanceOneTimeStep() {
	projectPressure();				// solve for pressure grid
	advectVelocity();				// update velocity grids using pressure grid
	t += dt;						// increase time by dt
	if(t + dt > frame_dt * frame) {
		dt = frame_dt * frame - t;	// shrink dt to match frame dt
		writeData();				// write data for each frame
		frame++;					// increase frame number
	}
}

void applyBC() {
}

double Simulator::getTime() {
	return t;
}