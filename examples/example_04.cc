#include <iostream>
#include "../include/map.hpp"
#include "../include/tracking.hpp"

/******
* In this example, we demonstrate how to convert a map that takes PTC coordinates
* into one that takes MAD-X coordinates and use it for tracking.
******/

int main() {
    int dim = 3;
    int da_order = 2;
    int da_dim = 2*dim;
    int da_scratch = 1000;
    da_init(da_order, da_dim, da_scratch);  //Initialize the da environment.

    //Read the truncated map from file.
    vector<Map> ptc_map(da_dim, Map(dim, da_full_length()));
    char filename[100] = "6D_2ndOrder_MADX.txt";
    readMadXMap(dim, filename, ptc_map);
    //Print the truncated map to screen.
    std::cout<<"PTC map: "<<std::endl;
    for(auto &v: ptc_map)
        std::cout<<v<<std::endl;

    //Tracking with PTC map
    char output[100] = "exp04_truncated_map_ptc_tracking.txt";
    double xi[6] = {0.3000000000000000E-04,   0.000000000000000,  0.3000000000000000E-03,   0.000000000000000, 0.3000000000000000E-03,   0.000000000000000};
    double xf[6] = {0,0,0,0,0,0};

    int n_turns = 10;
    bool output_to_file = true;
    int dn = 1;
    mapTrack(dim, xi, ptc_map, n_turns, xf, output_to_file, output, dn);

    //Convert into MAD-X map.
    vector<Map> madx_map(da_dim, Map(dim, da_full_length()));
    map_ptc_to_madx(ptc_map, madx_map);
    std::cout<<"MAD-X map: "<<std::endl;
    for(auto &v: madx_map)
        std::cout<<v<<std::endl;

    //Tracking with the MAD-X map
    char output2[100] = "exp04_truncated_map_madx_tracking.txt";
    double xi2[6] = {0.3000000000000000E-04,   0.000000000000000,  0.3000000000000000E-03,   0.000000000000000,   0.000000000000000, 0.3000000000000000E-03};

    mapTrack(dim, xi2, madx_map, n_turns, xf, output_to_file, output2, dn);



}
