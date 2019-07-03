#ifndef COORDINATES_TRANSFER_HPP_INCLUDED
#define COORDINATES_TRANSFER_HPP_INCLUDED

#include <cmath>
#include <vector>
#include "../tpsa_lib/include/da.h"

//X_madx = trans(X_cosy)
void trans_madx_to_cosy(double gamma, std::vector<DAVector> &trans);
//X_cosy = trans(X_madx)
void trans_cosy_to_madx(double gamma, std::vector<DAVector> &trans);
//X_madx = trans(X_ptc)
void trans_madx_to_ptc(std::vector<DAVector> &trans);
//X_ptc = trans(X_madx)
void trans_ptc_to_madx(std::vector<DAVector> &trans);
//X_cosy = trans(X_ptc)
void trans_cosy_to_ptc(double gamma, std::vector<DAVector> &trans);
//X_ptc = trans(X_cosy)
void trans_ptc_to_cosy(double gamma, std::vector<DAVector> &trans);

void da_map_ptc_to_madx(std::vector<DAVector>& ptc_map, std::vector<DAVector>& madx_map);
void da_map_cosy_to_madx(std::vector<DAVector>& cosy_map, double gamma, std::vector<DAVector>& madx_map);

#endif // COORDINATES_TRANSFER_HPP_INCLUDED

