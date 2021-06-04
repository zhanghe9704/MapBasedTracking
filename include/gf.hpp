#ifndef GF_HPP_INCLUDED
#define GF_HPP_INCLUDED

#include "map.hpp"
#include "../tpsa_lib/include/da.h"

void map2da(Map& m, DAVector& d);
void da2map(DAVector& d, Map& m);
void generating_function(std::vector<DAVector>& ivecs, int dim, const int type, DAVector& gf) ;
void gf_eqns(std::vector<DAVector>& map, int dim, const int type, std::vector<DAVector>& eqns);
#endif // GF_HPP_INCLUDED
