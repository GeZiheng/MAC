#pragma once
#include <iostream>
#include <array>
#include <vector>
#include <Eigen>

typedef std::array<int, 2> Index;
typedef std::array<double, 2> Vec2d;

template<typename T>
class Grid {					// sparse grid for pressure and velocity
public:
	int							active_cells;		// number of active cells in the grid
	std::vector<T>				data;				// storage for grid data
	std::vector<Index>			id_list;			// index list for grid nodes
	Index 						gridN;				// size of grid by dimension
	Vec2d						origin;				// point of origin of grid
	double						dx;					// dx for grid spacing

public:
	Grid() {
	}
	Grid(Index in_grid_N, Vec2d in_origin, int in_active_cells, double in_dx) {			// initialize grid with gridN and origin
		gridN[0] = in_grid_N[0];
		gridN[1] = in_grid_N[1];
		origin[0] = in_origin[0];
		origin[1] = in_origin[1];
		active_cells = in_active_cells;
		data.resize(active_cells);
		index_list.resize(active_cells);
		dx = in_dx;
	}
	~Grid() {
		data.clear();
		id_list.clear();
	}
	int flatIndex(Index index) {			// calculate flat index from (i,j)
		return index[0] * gridN[1] + index[1];
	}
	Index gridIndex(int i) {				// calculate (i,j) from flat index
		Index index;
		index[0] = i / gridN[0];
		index[1] = i - index[0] * gridN[0];
		return index;
	}
	void loadFromVector(int N, Eigen::VectorXd data_vec) {				// load grid data from an Eigen Vector
		for (int i = 0; i < N; i++) {
			data[i] = data_vec[i];
		}
	}
};