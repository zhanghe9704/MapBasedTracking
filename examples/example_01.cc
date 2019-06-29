#include <iostream>
#include "../include/map.hpp"
#include "../include/tracking.hpp"

/******
* In this example, we read the truncated map from an output of COSY Infinity 9.x
* and track the particle with a given initial coordinates using the map.
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

    //Convert the map into a group of da vectors.
    vector<DAVector> da_map(da_dim);
    for(int i=0; i<da_dim; ++i) {
            map2da(trunc_map.at(i), da_map.at(i));
            da_map.at(i) = da_map.at(i) - da_map.at(i).con();   //Remove the constant part!
    }
    //Print the map to screen again as DA vectors.
    for(auto&v: da_map) v.print();

    char output[100] = "track_truncated_map_cosy.txt";

    double xi[4] = {3.3000000000000000E-02,   0.000000000000000,  3.3000000000000000E-02,   0.000000000000000};
    double xf[4] = {0,0,0,0};

    int n_turns = 1000000;
    bool output_to_file = true;
    int dn = 1000;
    mapTrack(dim, xi, trunc_map, n_turns, xf, output_to_file, output, dn);

    int gf_type = 1; //Choose the type of the generating function.
    char output_symp[100] = "track_symplectic_type_1_cosy.txt";
    sympTrack(dim, xi, trunc_map, gf_type, n_turns, xf, output_to_file, output_symp, dn);

    return 0;
}
