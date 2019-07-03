#include <fstream>
#include <iostream>
#include "../include/map.hpp"
#include "../include/tracking.hpp"

/******
* In this example, we read the truncated map from an output of COSY Infinity 9.x
* and track the particle with a given initial coordinates using the map.
* Different with example_01, in this example, we separate the symplectic tracking
* into two steps: (1) create the differential equations and (2) perform tracking.
******/
int main() {
    int dim = 2;
    int da_order = 3;  //This should be 1 above the order of the map for gf calculation.
    int da_dim = 2*dim;
    int da_scratch = 1000;
    da_init(da_order, da_dim, da_scratch);  //Initialize the da environment.
    da_change_order(da_order-1);            //Reduce the order by one.

    //Read the truncated map from file.
    vector<Map> trunc_map(da_dim, Map(dim, da_full_length()));
    char filename[100] = "4D_2ndOrder_COSY.txt";
    readCOSYMap(dim, filename, trunc_map);

    //Print the map to screen.
    for(auto &v: trunc_map)
        std::cout<<v<<std::endl;

    //Step 1: Create the differential equations.
    vector<Map> eqns(da_dim, Map(dim, da_full_length()));
    int gf_type = 1; //Choose the type of the generating function.
    gfun(dim, trunc_map, gf_type, eqns);

    //Step 2: Symplectic tracking.
    char output_symp[100] = "exp03_track_symplectic_type_1_cosy.txt";
    std::ofstream output;
	output.open(output_symp);

    double xi[4] = {3.3000000000000000E-02,   0.000000000000000,  3.3000000000000000E-02,   0.000000000000000};
    double xf[4] = {0,0,0,0};

    int n_turns = 1000000;
    int n_freq = 1000;
    for(int i=0; i<n_turns; ++i) {
        if(i%n_freq==0) {
            output.precision(10);
            output<<std::noshowpos;
            output<<std::scientific;
            output<<i<<' ';
            output<<std::showpos;
            for(int j=0;j<2*dim;++j){
                output<<xi[j]<<' ';
            }
            output<<std::endl;
        }

        sympTrack(dim, xi, trunc_map, eqns, gf_type, 1, xf);
        for(int j=0; j<2*dim; ++j) xi[j] = xf[j];
    }
    output.precision(10);
    output<<std::noshowpos;
    output<<std::scientific;
    output<<n_turns<<' ';
    output<<std::showpos;
    for(int j=0;j<2*dim;++j){
        output<<xi[j]<<' ';
    }
    output<<std::endl;
    output.close();

    return 0;
}

