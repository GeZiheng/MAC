#pragma once
#include <cmath>

const double Pi = 4.0 * atan(1);

extern int res;						// grid res
extern double dx;					// grid dx width
extern double dt;					// time step length
extern double frame_dt;				// frame dt
extern double rho;					// density (constant, not changing)
extern double gravity;				// gravity
extern int frame_num;				// number of frames