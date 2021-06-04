#include "../include/gf.hpp"

#include <cassert>
#include <iostream>
#include <limits>
#include <vector>

void map2da(Map& m, DAVector& d) {
    int dim = m.getDim();
    int n = m.getTerm();
    d.reset_const(0);
    int *c = new int[2*dim];
    for(int i=0; i<n; ++i) {
        for(int j=0; j<2*dim; ++j) {
            c[j] = m.getOrders(i,j);
        }
        double coef = m.getCoef(i);
        d.set_element(c, coef);
    }
    delete[] c;
}

void da2map(DAVector& d, Map& m) {
    m.reset();
    unsigned int n = d.length();
    int dim = DAVector::dim();
    unsigned int *c = new unsigned int[2*dim];
    double elem = 0;
    std::vector<int> vec(dim, 0);
    int k = 0;
    int total_order = 0;
    for(unsigned int i=0; i<n; ++i) {
        d.element(i, c, elem);
        if (elem!=0) {
            m.addTerm();
            m.addCoef(elem);
            int order = 0;
            for(int j = 0; j<dim; ++j) {
                m.setOrders(k,j, c[j]);
                order += c[j];
            }
            if (order>total_order) total_order = order;
            ++k;
        }
    }
    m.setTotalOrder(total_order);
    delete[] c;
}

DAVector potential(std::vector<DAVector> &field, unsigned int dim) {
    dim *= 2;
    DAVector ptl;
    for(unsigned int i=0; i<dim; ++i) {
        DAVector ff = field.at(i);
        for(unsigned int j=i+1; j<dim; ++j) {
            DAVector s;
            da_substitute_const(ff, j, 0, s);
            ff = s;
        }
        ptl = ptl + da_int(ff, i);
    }
    return ptl;
}

void gf_equation_1(std::vector<DAVector>& ivecs, int dim, std::vector<DAVector>& ovecs) {
    dim *= 2;
    assert(dim==ivecs.size()&&dim<=DAVector::dim()&&"Wrong dimension of map in gf_equation_1!");
    for(auto& v : ivecs) assert(v.con()<std::numeric_limits<double>::min()&&"Constant part of the input map is NOT zero in gf_equation_1!");

    std::vector<DAVector> n2(dim);
    std::vector<DAVector> tmp(dim);

    ovecs.resize(dim);
    for(int i=0; i<dim; i+=2) {
        ovecs.at(i) = da[i];
        ovecs.at(i+1) = ivecs.at(i);
        n2.at(i) = da[i+1];
        n2.at(i+1) = -1*ivecs.at(i+1);
    }
    inv_map(ovecs, dim, tmp);
    for(int i=0; i<dim; ++i) ovecs.at(i) = tmp.at(i);

    da_composition(n2, ovecs, tmp);
    for(int i=0; i<dim; ++i) ovecs.at(i) = tmp.at(i);
}

void gf_equation_2(std::vector<DAVector>& ivecs, int dim, std::vector<DAVector>& ovecs) {
    dim *= 2;
    assert(dim==ivecs.size()&&dim<=DAVector::dim()&&"Wrong dimension of map in gf_equation_1!");
    for(auto& v : ivecs) assert(v.con()<std::numeric_limits<double>::min()&&"Constant part of the input map is NOT zero in gf_equation_1!");

    std::vector<DAVector> n2(dim);
    std::vector<DAVector> tmp(dim);

    ovecs.resize(dim);
    for(int i=0; i<dim; i+=2) {
        ovecs.at(i) = da[i];
        ovecs.at(i+1) = ivecs.at(i+1);
        n2.at(i) = da[i+1];
        n2.at(i+1) = ivecs.at(i);
    }
    inv_map(ovecs, dim, tmp);
    for(int i=0; i<dim; ++i) ovecs.at(i) = tmp.at(i);

    da_composition(n2, ovecs, tmp);
    for(int i=0; i<dim; ++i) ovecs.at(i) = tmp.at(i);
}

void gf_equation_3(std::vector<DAVector>& ivecs, int dim, std::vector<DAVector>& ovecs) {
    dim *= 2;
    assert(dim==ivecs.size()&&dim<=DAVector::dim()&&"Wrong dimension of map in gf_equation_1!");
    for(auto& v : ivecs) assert(v.con()<std::numeric_limits<double>::min()&&"Constant part of the input map is NOT zero in gf_equation_1!");

    std::vector<DAVector> n2(dim);
    std::vector<DAVector> tmp(dim);

    ovecs.resize(dim);
    for(int i=0; i<dim; i+=2) {
        ovecs.at(i) = da[i+1];
        ovecs.at(i+1) = ivecs.at(i);
        n2.at(i) = -1*da[i];
        n2.at(i+1) = -1*ivecs.at(i+1);
    }
    inv_map(ovecs, dim, tmp);
    for(int i=0; i<dim; ++i) ovecs.at(i) = tmp.at(i);

    da_composition(n2, ovecs, tmp);
    for(int i=0; i<dim; ++i) ovecs.at(i) = tmp.at(i);
}

void gf_equation_4(std::vector<DAVector>& ivecs, int dim, std::vector<DAVector>& ovecs) {
    dim *= 2;
    assert(dim==ivecs.size()&&dim<=DAVector::dim()&&"Wrong dimension of map in gf_equation_1!");
    for(auto& v : ivecs) assert(v.con()<std::numeric_limits<double>::min()&&"Constant part of the input map is NOT zero in gf_equation_1!");

    std::vector<DAVector> n2(dim);
    std::vector<DAVector> tmp(dim);

    ovecs.resize(dim);
    for(int i=0; i<dim; i+=2) {
        ovecs.at(i) = da[i+1];
        ovecs.at(i+1) = ivecs.at(i+1);
        n2.at(i) = -1*da[i];
        n2.at(i+1) = ivecs.at(i);
    }
    inv_map(ovecs, dim, tmp);
    for(int i=0; i<dim; ++i) ovecs.at(i) = tmp.at(i);

    da_composition(n2, ovecs, tmp);
    for(int i=0; i<dim; ++i) ovecs.at(i) = tmp.at(i);
}

void gf_eqns(std::vector<DAVector>& ivecs, int dim, const int type, std::vector<DAVector>& ovecs) {
    assert(2*dim<=DAVector::dim()&&"Wrong dimension of map in gf_eqns!");
    assert(ivecs.size()==ovecs.size()&&"Input map and output eqns have different dimension in gf_eqns!");
    switch(type) {
    case 1 : {
        gf_equation_1(ivecs, dim, ovecs);
        break;
    }
    case 2 : {
        gf_equation_2(ivecs, dim, ovecs);
        break;
    }
    case 3 : {
        gf_equation_3(ivecs, dim, ovecs);
        break;
    }
    case 4 : {
        gf_equation_4(ivecs, dim, ovecs);
        break;
    }
    default : {
        std::cout<< "Wrong selection of generating function!" <<std::endl;
    }
    }

    da_restore_order();
    DAVector gf = potential(ovecs, dim);

//    gf.print();

    switch(type) {
    case 1 : {
//        for(int i=0; i<2*dim; i+=2) ovecs.at(i) = da_der(gf, i+1);
//        for(int i=1; i<2*dim; i+=2) ovecs.at(i) = -1*da_der(gf, i-1);
        for(int i=0; i<2*dim; i+=2) ovecs.at(i) = da_der(gf, i);
        for(int i=1; i<2*dim; i+=2) ovecs.at(i) = -1*da_der(gf, i);
//        gf_equation_1(ivecs, dim, ovecs);
        break;
    }
    case 2 : {
        for(int i=0; i<2*dim; i+=2) ovecs.at(i) = da_der(gf, i+1);
        for(int i=1; i<2*dim; i+=2) ovecs.at(i) = da_der(gf, i-1);
        break;
    }
    case 3 : {
        for(int i=0; i<2*dim; i+=2) ovecs.at(i) = -1*da_der(gf, i+1);
        for(int i=1; i<2*dim; i+=2) ovecs.at(i) = -1*da_der(gf, i-1);
//        gf_equation_3(ivecs, dim, ovecs);
        break;
    }
    case 4 : {
        for(int i=0; i<2*dim; i+=2) ovecs.at(i) = -1*da_der(gf, i);
        for(int i=1; i<2*dim; i+=2) ovecs.at(i) = da_der(gf, i);
//        gf_equation_4(ivecs, dim, ovecs);
        break;
    }
    default : {
        std::cout<< "Wrong selection of generating function!" <<std::endl;
    }
    }
    da_change_order(DAVector::order()-1);
}

//
void generating_function(std::vector<DAVector>& ivecs, int dim, const int type, DAVector& gf) {
    da_change_order(DAVector::order()-1);
    std::vector<DAVector> ovecs(dim);
    switch(type) {
    case 1 : {
        gf_equation_1(ivecs, dim, ovecs);
        break;
    }
    case 2 : {
        gf_equation_2(ivecs, dim, ovecs);
        break;
    }
    case 3 : {
        gf_equation_3(ivecs, dim, ovecs);
        break;
    }
    case 4 : {
        gf_equation_4(ivecs, dim, ovecs);
        break;
    }
    default : {
        std::cout<< "Wrong selection of generating function!" <<std::endl;
    }
    }
    da_restore_order();
    gf = potential(ovecs, dim);
}

