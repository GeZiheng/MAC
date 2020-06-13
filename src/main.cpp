#include <iostream>
#include <fstream>
#include <GL/glut.h>
#include "Simulator.h"

using namespace std;

int main(int argc, char **argv) {
	Simulator sim;
	sim.init();
	while (sim.getTime() < sim.end_t)
		sim.advanceOneTimeStep();
	return 0;
}