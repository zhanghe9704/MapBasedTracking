#include "../include/coordinates_transfer.hpp"

//X_madx = trans(X_cosy)
void trans_madx_to_cosy(double gamma, std::vector<DAVector> &trans) {
    trans.resize(6);
    trans.at(0) = da[0];
    trans.at(1) = da[1];
    trans.at(2) = da[2];
    trans.at(3) = da[3];

    double tmp = sqrt((gamma-1)/(gamma+1));
    trans.at(4) = da[4]*tmp;
    trans.at(5) = da[5]/tmp;
}

//X_cosy = trans(X_madx)
void trans_cosy_to_madx(double gamma, std::vector<DAVector> &trans) {
    trans.resize(6);
    trans.at(0) = da[0];
    trans.at(1) = da[1];
    trans.at(2) = da[2];
    trans.at(3) = da[3];

    double tmp = sqrt((gamma-1)/(gamma+1));
    trans.at(4) = da[4]*tmp;
    trans.at(5) = da[5]/tmp;
}

//X_madx = trans(X_ptc)
void trans_madx_to_ptc(std::vector<DAVector> &trans) {
    trans.resize(6);
    trans.at(0) = da[0];
    trans.at(1) = da[1];
    trans.at(2) = da[2];
    trans.at(3) = da[3];
    trans.at(4) = da[5];
    trans.at(5) = -da[4];
}

//X_ptc = trans(X_madx)
void trans_ptc_to_madx(std::vector<DAVector> &trans) {
   trans_madx_to_ptc(trans);
}

//X_cosy = trans(X_ptc)
void trans_cosy_to_ptc(double gamma, std::vector<DAVector> &trans) {
    trans.resize(6);
    trans.at(0) = da[0];
    trans.at(1) = da[1];
    trans.at(2) = da[2];
    trans.at(3) = da[3];

    double tmp = sqrt((gamma-1)/(gamma+1));
    trans.at(4) = -da[5]*tmp;
    trans.at(5) = da[4]/tmp;
}

//X_ptc = trans(X_cosy)
void trans_ptc_to_cosy(double gamma, std::vector<DAVector> &trans) {
    trans.resize(6);
    trans.at(0) = da[0];
    trans.at(1) = da[1];
    trans.at(2) = da[2];
    trans.at(3) = da[3];

    double tmp = sqrt((gamma-1)/(gamma+1));
    trans.at(4) = -da[5]/tmp;
    trans.at(5) = da[4]*tmp;
}

void da_map_ptc_to_madx(std::vector<DAVector>& ptc_map, std::vector<DAVector>& madx_map) {
    std::vector<DAVector> trans;
    trans_ptc_to_madx(trans);
    da_composition(ptc_map, trans,  madx_map);
    std::swap(madx_map.at(4), madx_map.at(5));
    madx_map.at(4) *= -1;
}

void da_map_cosy_to_madx(std::vector<DAVector>& cosy_map, double gamma, std::vector<DAVector>& madx_map) {
    std::vector<DAVector> trans;
    trans_cosy_to_madx(gamma, trans);
    da_composition(cosy_map, trans, madx_map);
    double coef = sqrt((gamma+1)/(gamma-1));
    madx_map.at(4) *= coef;
    coef = 1/coef;
    madx_map.at(5) *= coef;
}
