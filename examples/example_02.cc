#include <iostream>
#include "../include/map.hpp"
#include "../include/tracking.hpp"

/******
* In this example, we read the truncated map from an output of MAD-X
* and track the particle with a given initial coordinates using the map.
******/

int main() {
    int dim = 3;
    int da_order = 3;  //This should be 1 above the order of the map for gf calculation.
    int da_dim = 2*dim;
    int da_scratch = 1000;
    da_init(da_order, da_dim, da_scratch);  //Initialize the da environment.
    da_change_order(da_order-1);            //Reduce the order by one.

    //Read the truncated map from file.
    vector<Map> trunc_map(da_dim, Map(dim, da_full_length()));
    char filename[100] = "6D_2ndOrder_MADX.txt";
    readMadXMap(dim, filename, trunc_map);
    //Print the truncated map to screen.
    for(auto &v: trunc_map)
        std::cout<<v<<std::endl;

    //Convert the map into a group of da vectors.
    vector<DAVector> da_map(da_dim);
    for(int i=0; i<da_dim; ++i) {
            map2da(trunc_map.at(i), da_map.at(i));
            da_map.at(i) = da_map.at(i) - da_map.at(i).con();   //Remove the constant part!
    }
    //Print the map to screen again as DA vectors.
    for(auto&v: da_map) v.print();

    char output[100] = "track_truncated_map_madx.txt";
    double xi[6] = {0.3000000000000000E-04,   0.000000000000000,  0.3000000000000000E-03,   0.000000000000000, 0.3000000000000000E-03,   0.000000000000000};
    double xf[6] = {0,0,0,0,0,0};

     int n_turns = 10000;
    bool output_to_file = true;
    mapTrack(dim, xi, trunc_map, n_turns, xf, output_to_file, output);

    int gf_type = 2; //Choose the type of the generating function.
    char output_symp[100] = "track_symplectic_type_2_madx.txt";
    sympTrack(dim, xi, trunc_map, gf_type, n_turns, xf, output_to_file, output_symp);

    return 0;
}
