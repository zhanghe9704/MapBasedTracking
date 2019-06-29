#ifndef GF_HPP_INCLUDED
#define GF_HPP_INCLUDED

#include "map.hpp"
#include "../tpsa_lib/include/da.h"

void map2da(Map& m, DAVector& d);
void da2map(DAVector& d, Map& m);
void generating_function(std::vector<DAVector>& ivecs, int dim, const int type, DAVector& gf) ;
void generating_function_1(std::vector<DAVector>& ivecs, int dim, DAVector& gf);
void generating_function_2(std::vector<DAVector>& ivecs, int dim, DAVector& gf);
void generating_function_3(std::vector<DAVector>& ivecs, int dim, DAVector& gf);
void generating_function_4(std::vector<DAVector>& ivecs, int dim, DAVector& gf);
void eqns_gf2(DAVector& gf2, int dim, std::vector<DAVector>& eqns);
void gf_eqns(std::vector<DAVector>& map, int dim, const int type, std::vector<DAVector>& eqns);
#endif // GF_HPP_INCLUDED
