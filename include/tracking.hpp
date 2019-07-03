#ifndef TRACKING_HPP
#define TRACKING_HPP

#include <vector>
#include <iostream>
#include "coordinates_transfer.hpp"
#include "DNewtonSolver.hpp"
#include "gf.hpp"
#include "map.hpp"
#include "../tpsa_lib/include/da.h"

using std::vector;

//The equations to solve.
void gf2Eqns (double * xf, void * p, double * f);

/** \brief Perform symplectic tracking.
 *
 * \param dim - Dimension of the dynamic system.
 * \param xi - Initial coordinates.
 * \param map - Truncated map
 * \param type - Type of the generating function, 1, 2, 3, or 4.
 * \param nTrk - Number of turns to track.
 * \param xf - Final coordinates.
 * \param outToFile - Whether to write the tracking result to a file.
 * \param filename - File to save the tracking result.
 * \param nFreq - Save the result every nFreq turns.
 * \param nIter - Max number of iterations in Newton solver.
 * \param delta - Tolerance of the Newton solver.
 * \return 0.
 *
 */
int sympTrack(const int dim, double * xi, vector<Map> &map, int type, const int nTrk, double * xf,  bool outToFile=false,
              char * filename=NULL, const int nFreq=1, const int nIter=10, const double delta=1e-16);

/** \brief Perform symplectic tracking.
 *
 * \param dim - Dimension of the dynamic system.
 * \param xi - Initial coordinates.
 * \param map - Truncated map
 * \param eqns - The differential equations derived from the generating function.
 * \param type - Type of the generating function, 1, 2, 3, or 4.
 * \param n_turn - Number of turns to track.
 * \param xf - Final coordiantes.
 * \param nIter - Max number of iterations in Newton solver.
 * \param delta - Tolerance of the Newton solver.
 * \return 0.
 *
 */
int sympTrack(const int dim, double * xi, vector<Map> &tr_map, vector<Map> &eqns, int type, const int n_turn, double * xf,
              const int n_iter=10, const double delta=1e-16);

/** \brief Perform tracking using the truncated map.
 *
 * \param dim - Dimension of the dynamic system.
 * \param xi - Initial coordinates.
 * \param map - Truncated map
 * \param nTrk - Number of turns to track.
 * \param xf - Final coordinates.
 * \param outToFile - Whether to write the tracking result to a file.
 * \param filename - File to save the tracking result.
 * \param nFreq - Save the result every nFreq turns.
 * \return 0.
 *
 */
int mapTrack(const int dim, double * xi, vector<Map> &map, const int nTrk,  double * xf, bool outToFile=false, char * filename=NULL, const int nFreq=1);

/** \brief Derive the differential equations for symplectic tracking.
 *
 * \param dim - Dimension of the dynamic system.
 * \param tr_map - Truncated map.
 * \param type - Type of the generating function, 1, 2, 3, or 4.
 * \param eqns - Equations for symplectic tracking.
 * \return None.
 *
 */
void gfun(const int dim, vector<Map> &tr_map, int type, vector<Map> &eqns);

/** \brief Convert a map that takes PTC coordinates into one that takes MAD-X coordinates.
 *
 * \param[in] ptc_map - map that takes PTC coordinates.
 * \param[out] madx_map - map that takes MADX-X coordinates.
 * \return None.
 *
 */
void map_ptc_to_madx(vector<Map>& ptc_map, vector<Map>& madx_map);

/** \brief Convert a map that takes COSY Infinity 9.x coordinates into one that takes MAD-X coordinates.
 *
 * \param[in] cosy_map - map that takes COSY Infinity 9.x coordinates
 * \param[in] gamma - Lorentz factor of the reference particle.
 * \param[out] madx_map - map that takes MAD-X coordinates.
 * \return
 *
 */
void map_cosy_to_madx(vector<Map>& cosy_map, double gamma, vector<Map>& madx_map);
#endif
