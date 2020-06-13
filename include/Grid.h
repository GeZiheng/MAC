#pragma once
#include <iostream>
#include <array>
#include <vector>
#include <Eigen>

using Vec2i = std::array<int, 2>;
using Vec2d = std::array<double, 2>;

template<typename T>
class Grid {					// sparse grid for pressure and velocity
public:
	int							active_cells;		// number of active cells in the grid
	std::vector<T>				data;				// storage for grid data
	std::vector<Vec2i>			id_list;			// index list for grid nodes
	Vec2i 						gridN;				// size of grid by dimension
	Vec2d						origin;				// point of origin of grid
	double						dx;					// dx for grid spacing

public:
	Grid() {
	}
	Grid(Vec2i in_grid_N, Vec2d in_origin, int in_active_cells, double in_dx) {			// initialize grid with gridN and origin
		gridN[0] = in_grid_N[0];
		gridN[1] = in_grid_N[1];
		origin[0] = in_origin[0];
		origin[1] = in_origin[1];
		active_cells = in_active_cells;
		data.resize(active_cells);
		id_list.resize(active_cells);
		dx = in_dx;
	}
	~Grid() {
		data.clear();
		id_list.clear();
	}
	int flatIndex(Vec2i index) {			// calculate flat index from (i,j)
		return index[0] * gridN[1] + index[1];
	}
	Vec2i gridIndex(int i) {				// calculate (i,j) from flat index
		Vec2i index;
		index[0] = i / gridN[0];
		index[1] = i - index[0] * gridN[0];
		return index;
	}
	void loadFromVector(int N, Eigen::VectorXd data_vec) {				// load grid data from an Eigen Vector
		for (int i = 0; i < N; i++) {
			data[i] = data_vec[i];
		}
	}

	void interpolate(Vec2d pt, std::vector<Vec2i>& grid_id, std::vector<double>& weights) {		// find weights for (linear) interpolation on a particle
		grid_id.clear();
		weights.clear();
		int bn_i = floor(pt[0] / dx);
		int bn_j = floor(pt[1] / dx);
		double wx = pt[0] / dx - i;
		double wy = pt[1] / dx - j;
		// downleft node
		grid_id.push_back({ bn_i, bn_j });
		weights.push_back(wx * wy);
		// downright node
		grid_id.push_back({ bn_i, bn_j+1 });
		weights.push_back(wx * (1-wy));
		// upleft node
		grid_id.push_back({ bn_i+1, bn_j });
		weights.push_back((1-wx) * wy);
		// upright node
		grid_id.push_back({ bn_i+1, bn_j+1 });
		weights.push_back((1-wx) * (1-wy));
	}
};