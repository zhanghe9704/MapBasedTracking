#ifndef TRACKING_HPP
#define TRACKING_HPP

#include <vector>
#include <iostream>
#include "DNewtonSolver.hpp"
#include "gf.hpp"
#include "map.hpp"
#include "../tpsa_lib/include/da.h"

using std::vector;

//The equations to solve.
void gf2Eqns (double * xf, void * p, double * f);

//Symplectic tracking
//dim - dimenstion of the problem, xi - initial coordinates, eqns - equations of GF2 to solve, map - the truncated map
//xf - final coordinates, nIter - max iteration time for Newton solver, delta - tolerance for Newton solver
//nTrk - trakcing times, outToFile - output to file or screen, filename - file to output data, nFreq - output every nFreq times
int sympTrack(const int dim, double * xi, vector<Map> &map, int type, const int nTrk, double * xf,  bool outToFile=false,
              char * filename=NULL, const int nFreq=1, const int nIter=10, const double delta=1e-16);
int mapTrack(const int dim, double * xi, vector<Map> &map, const int nTrk,  double * xf, bool outToFile=false, char * filename=NULL, const int nFreq=1);

void gfun(const int dim, vector<Map> &tr_map, int type, vector<Map> &eqns);

void trans_madx_to_cosy(double gamma, std::vector<DAVector> &trans);
void trans_cosy_to_madx(double gamma, std::vector<DAVector> &trans);
#endif
