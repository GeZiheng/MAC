#pragma once
#include <Eigen>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "Simulator.h"

using namespace Eigen;

Simulator::Simulator() {
}

Simulator::~Simulator() {
}

void Simulator::projectPressure() {
	int N = p_grid.active_cells;
	VectorXd p_vec(N), rhs(N);
	SparseMatrix<double> A(N, N);
	Vec2i grid_index;
	int i, j, q;
	std::vector<Eigen::Triplet<double>> triplets;
	/* Build A and rhs for linear system */
	for(int p = 0; p < N; p++) {
		/* Calculate A entry for current node */
		triplets.push_back(Eigen::Triplet<double>(p, p, 4 * dt / (rho * dx * dx)));
		/* Find grid index for current node */
		grid_index = p_grid.id_list[p];
		i = grid_index[0];
		j = grid_index[1];
		/* Calculate A entries for each neighbor nodes */
		q = labels[p_grid.flatIndex({i+1, j})];
		triplets.push_back(Eigen::Triplet<double>(p, q, -dt / (rho * dx * dx)));
		q = labels[p_grid.flatIndex({i-1, j})];
		triplets.push_back(Eigen::Triplet<double>(p, q, -dt / (rho * dx * dx)));
		q = labels[p_grid.flatIndex({i, j+1})];
		triplets.push_back(Eigen::Triplet<double>(p, q, -dt / (rho * dx * dx)));
		q = labels[p_grid.flatIndex({i, j-1})];
		triplets.push_back(Eigen::Triplet<double>(p, q, -dt / (rho * dx * dx)));
		/* Calculate rhs for each neighbor nodes */
		q = labels[p_grid.flatIndex({i+1, j})];
		rhs(p) -= u_grid.data[q] / dx;
		q = labels[p_grid.flatIndex({i-1, j})];
		rhs(p) += u_grid.data[q] / dx;
		q = labels[p_grid.flatIndex({i, j+1})];
		rhs(p) += v_grid.data[q] / dx;
		q = labels[p_grid.flatIndex({i, j-1})];
		rhs(p) += v_grid.data[q] / dx;
	}
	A.setFromTriplets(triplets.begin(), triplets.end());
	/* Solve linear system using conjugate gradient method */
	ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
	cg.compute(A);
	p_vec = cg.solve(rhs);
	/* Load data into pressure grid */
	p_grid.loadFromVector(N, p_vec);
}

void Simulator::updateVelocity() {
	int N = p_grid.active_cells;					// number of fluid cells
	Vec2i grid_index;
	int i, j, q;
	/* Update velocity (u,v) grids */
	for(int p = 0; p < N; p++) {
		/* Find grid index for current node */
		grid_index = p_grid.id_list[p];
		i = grid_index[0];
		j = grid_index[1];
		/* update u grid nodes */
		q = labels[p_grid.flatIndex({i+1, j})];
		u_grid.data[q] = u_grid.data[q] - dt / (rho * dx) * (p_grid.data[q] - p_grid.data[p]);
		/* update v grid nodes */
		q = labels[p_grid.flatIndex({i, j+1})];
		v_grid.data[q] = v_grid.data[q] - dt / (rho * dx) * (p_grid.data[q] - p_grid.data[p]);
	}
}

void Simulator::advection() {
	std::vector<Vec2i> grid_id;
	std::vector<Vec2d> weights;
	/* Semi-Lagrangian advection for u-grid */
	for (int p = 0; p < u_grid.active_cells; p++) {
		double u = u_grid.data[p];
		int i = u_grid.id_list[p][0];
		int j = u_grid.id_list[p][1];
	}
	/* Semi-Lagrangian advection for v-grid */
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
	dx = double(1) / double(res);

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
	num_frames = 0;
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

	// get density from file
	str.clear();
	getline(myfile, s);
	for (i = 8; i < s.length(); i++)
		str.push_back(s[i]);
	rho = stod(str);

	// get num_frames from file
	str.clear();
	getline(myfile, s);
	for (i = 8; i < s.length(); i++)
		str.push_back(s[i]);
	num_frames = stoi(str);

	myfile.close();
}

void Simulator::writeData() {
	std::ofstream myfile;
	std::ostringstream filename;
	int N = p_grid.active_cells;					// number of fluid cells
	Vec2i grid_index;
	int i, j;
	filename << "../output/frame" << std::setfill('0') << std::setw(4) << num_frames << ".dat";
	myfile.open(filename.str());
	myfile << "variables=\"x\",\"y\",\"p\",\"u\",\"v\"" << std:: endl;
	for (int p = 0; p < N; p++) {
		grid_index = p_grid.id_list[p];
		i = grid_index[0];
		j = grid_index[1];
		myfile << std::scientific << std::setprecision(16) << p_grid.origin[0] + i * p_grid.dx << ' ' << p_grid.origin[1] + j * p_grid.dx << ' ' << p_grid.data[p] << ' ' << u_grid.data[p] << ' ' << v_grid.data[p] << std::endl;
	}
	myfile.close();
}

void Simulator::applyBC() {
}

void Simulator::init() {
	loadData();
	t = 0;
	frame = 1;
	p_grid = Grid<double>({ res + 1,res + 1 }, { 0,0 }, 0, dx);
	u_grid = Grid<double>({ res + 1,res + 1 }, { 0,0 }, 0, dx);
	v_grid = Grid<double>({ res + 1,res + 1 }, { 0,0 }, 0, dx);
}

void Simulator::advanceOneTimeStep() {
	projectPressure();				// solve for pressure grid
	updateVelocity();				// update velocity grids using pressure grid
	t += dt;						// increase time by dt
	if(t + dt > frame_dt * frame) {
		dt = frame_dt * frame - t;	// shrink dt to match frame dt
		writeData();				// write data for each frame
		frame++;					// increase frame number
	}
}

double Simulator::getTime() {
	return t;
}