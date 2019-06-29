#include "../include/tracking.hpp"
#include <cmath>
#include <iomanip>
#include <fstream>


using std::cout;
using std::endl;


//Tracking by truncated map
int mapTrack(const int dim, double * xi, vector<Map> &map, const int nTrk, double * xf, bool outToFile, char * filename, const int nFreq){
	std::ofstream outfile;
	if(outToFile) outfile.open(filename);

	std::ostream & output = (outToFile?outfile:cout);

	output.precision(10);
    output<<std::noshowpos;
    output<<std::scientific;
    output<<0<<' ';
    output<<std::showpos;
    for(int j=0;j<2*dim;++j){
        output<<xi[j]<<' ';
    }
    output<<endl;

	for(int i=0; i<2*dim;++i) xf[i] = xi[i];

	for(int i=1; i<=nTrk; ++i){
		applyMap(dim, xf, map, xf);
		if (i%nFreq==0){
            output<<std::noshowpos;
			output<<i<<' ';
            output<<std::showpos;
			for(int j=0;j<2*dim;++j){
				output<<xf[j]<<' ';
			}
			output<<endl;
		}
	}
	if(outToFile) outfile.close();
	return 0;
}

//Parameters for the function gf2Eqns
struct gfParas{
	vector<Map> * eqns;
	double * xi;
	int dim;
	double * x;
};

//The equations to solve.
void gf2Eqns (double * xf, void * p, double * f){
	int dim = ((gfParas *) p)->dim;
	double * x = ((gfParas *) p)->x;
//	double * x = new double [2*dim];

	for(int i=0; i<dim; ++i){
		x[2*i] = ((gfParas *) p)->xi[2*i];
		x[2*i+1] = xf[2*i+1];
	}
	applyMap(dim, x, *(((gfParas *) p)->eqns), f);

	for(int i=0; i<dim; ++i){
		f[2*i] = xf[2*i]-f[2*i];
		f[2*i+1] = ((gfParas *) p)->xi[2*i+1]-f[2*i+1];
	}
//	delete[] x;
}

void gf1Eqns (double * xf, void * p, double * f){
	int dim = ((gfParas *) p)->dim;
	double * x = ((gfParas *) p)->x;
//	double * x = new double [2*dim];

	for(int i=0; i<dim; ++i){
		x[2*i] = ((gfParas *) p)->xi[2*i];
		x[2*i+1] = xf[2*i];
	}
	applyMap(dim, x, *(((gfParas *) p)->eqns), f);

	for(int i=0; i<dim; ++i){
		f[2*i] = ((gfParas *) p)->xi[2*i+1]-f[2*i] ;
		f[2*i+1] = xf[2*i+1]-f[2*i+1];
	}
//	delete[] x;
}

void gf3Eqns (double * xf, void * p, double * f){
	int dim = ((gfParas *) p)->dim;
	double * x = ((gfParas *) p)->x;
//	double * x = new double [2*dim];

	for(int i=0; i<dim; ++i){
		x[2*i] = ((gfParas *) p)->xi[2*i+1];
		x[2*i+1] = xf[2*i];
	}
	applyMap(dim, x, *(((gfParas *) p)->eqns), f);

	for(int i=0; i<dim; ++i){
		f[2*i] = xf[2*i+1]-f[2*i];
		f[2*i+1] = ((gfParas *) p)->xi[2*i]-f[2*i+1];
	}

//	delete[] x;
}

void gf4Eqns (double * xf, void * p, double * f){
	int dim = ((gfParas *) p)->dim;
	double * x = ((gfParas *) p)->x;
//	double * x = new double [2*dim];

	for(int i=0; i<dim; ++i){
		x[2*i] = ((gfParas *) p)->xi[2*i+1];
		x[2*i+1] = xf[2*i+1];
	}
	applyMap(dim, x, *(((gfParas *) p)->eqns), f);

	for(int i=0; i<dim; ++i){
		f[2*i] = ((gfParas *) p)->xi[2*i]-f[2*i];
		f[2*i+1] = xf[2*i]-f[2*i+1];
	}
//	delete[] x;
}

//Symplectic tracking
//dim - dimenstion of the problem, xi - initial coordinates, eqns - equations of GF2 to solve, map - the truncated map
//xf - final coordinates, nIter - max iteration time for Newton solver, delta - tolerance for Newton solver
//nTrk - trakcing times, outToFile - output to file or screen, filename - file to output data, nFreq - output every nFreq times
int sympTrack(const int dim, double * xi, vector<Map> &tr_map, int type, const int nTrk, double * xf, bool outToFile,
              char * filename, const int nFreq, const int nIter, const double delta){

	std::vector<DAVector> da_eqns(dim*2);
	std::vector<DAVector> da_map(dim*2);
	vector<Map> eqns(dim*2, Map(dim, da_full_length()));
	for(int i=0; i<dim*2; ++i) {
        map2da(tr_map.at(i), da_map.at(i));
        da_map.at(i) = da_map.at(i) - da_map.at(i).con();
	}
	gf_eqns(da_map, dim, type, da_eqns);

	for(int i=0; i<dim*2; ++i)
        da2map(da_eqns.at(i), eqns.at(i));

//    for(auto& v: da_eqns) v.print();

    double * xt = new double [2*dim];
	struct gfParas p = {&eqns, xi, dim, xt};

	std::ofstream outfile;
	if(outToFile) outfile.open(filename);

	std::ostream & output = (outToFile?outfile:cout);

	double * x = new double[2*dim];
	for(int i=0; i<2*dim;++i) x[i] = xi[i];
	void (*pFunc)(double *, void *, double *) = nullptr;
	switch(type) {
	case 1 : {
	    pFunc = &gf1Eqns;
        break;
	}
	case 2 : {
	    pFunc = &gf2Eqns;
        break;
	}
	case 3 : {
	    pFunc = &gf3Eqns;
        break;
	}
	case 4 : {
	    pFunc = &gf4Eqns;
        break;
	}
	default : {
	    std::cout<< "Wrong selection of generating function!" <<std::endl;
	}
	}



	output.precision(10);
    output<<std::noshowpos;
    output<<std::scientific;
    output<<0<<' ';
    output<<std::showpos;
    for(int j=0;j<2*dim;++j){
        output<<xi[j]<<' ';
    }
    output<<endl;


	for(int i=1; i<=nTrk; ++i){
		applyMap(dim, x, tr_map, xf);
		NewtonIter(pFunc, xf, &p, 2*dim, nIter, delta);
		for(int j=0; j<2*dim; ++j) x[j] = xf[j];
		p = {&eqns, x, dim, xt};

		if (i%nFreq==0){
            output<<std::noshowpos;
			output<<i<<' ';
            output<<std::showpos;
			for(int j=0;j<2*dim;++j){
//				output.setf(std::ios::scientific);
//				output<<std::setprecision(8);
				output<<xf[j]<<' ';
			}
			output<<endl;
		}
	}

	if(outToFile) outfile.close();
	delete[] x;
	delete[] xt;
	return 0;
}


int sympTrack(const int dim, double * xi, vector<Map> &tr_map, vector<Map> &eqns, int type, const int n_turn, double * xf,
              const int n_iter, const double delta){

    vector<double> vec_x(2*dim);
    vector<double> vec_xt(2*dim);
    double * xt = &vec_xt[0];
    double * x = &vec_x[0];

//    double * xt = new double [2*dim];
	struct gfParas p = {&eqns, xi, dim, xt};

//	double * x = new double[2*dim];
	for(int i=0; i<2*dim;++i) x[i] = xi[i];
	void (*pFunc)(double *, void *, double *) = nullptr;
	switch(type) {
	case 1 : {
	    pFunc = &gf1Eqns;
        break;
	}
	case 2 : {
	    pFunc = &gf2Eqns;
        break;
	}
	case 3 : {
	    pFunc = &gf3Eqns;
        break;
	}
	case 4 : {
	    pFunc = &gf4Eqns;
        break;
	}
	default : {
	    std::cout<< "Wrong selection of generating function!" <<std::endl;
	}
	}

	for(int i=1; i<=n_turn; ++i){
		applyMap(dim, x, tr_map, xf);
		NewtonIter(pFunc, xf, &p, 2*dim, n_iter, delta);
		for(int j=0; j<2*dim; ++j) x[j] = xf[j];
		p = {&eqns, x, dim, xt};
	}

//	delete[] x;
//	delete[] xt;
	return 0;
}

void gfun(const int dim, vector<Map> &tr_map, int type, vector<Map> &eqns){

	int da_dim = 2*dim;
	int da_order = 0;
	for(auto& v: tr_map) {
//        std::cout<<v;
//        std::cout<<"=================="<<std::endl;
        if(v.getTotalOrder()>da_order)
            da_order = v.getTotalOrder();
	}
	int da_scratch = 1000;

    da_order += 1;
	da_init(da_order, da_dim, da_scratch);
	da_change_order(da_order-1);

	std::vector<DAVector> da_eqns(da_dim);
	std::vector<DAVector> da_map(da_dim);

	for(int i=0; i<da_dim; ++i) {
        map2da(tr_map.at(i), da_map.at(i));
        da_map.at(i) = da_map.at(i) - da_map.at(i).con();
	}
	gf_eqns(da_map, dim, type, da_eqns);

	for(int i=0; i<da_dim; ++i)
        da2map(da_eqns.at(i), eqns.at(i));

    da_clear();

}

void trans_madx_to_cosy(double gamma, std::vector<DAVector> &trans) {
    trans.resize(6);
    trans.at(0) = da[0];
    trans.at(1) = da[1];
    trans.at(2) = da[2];
    trans.at(3) = da[3];

    double tmp = sqrt((gamma-1)/(gamma+1));
    trans.at(4) = da[4]*tmp;
    trans.at(5) = da[5]*gamma/tmp;
}

void trans_cosy_to_madx(double gamma, std::vector<DAVector> &trans) {
    trans.resize(6);
    trans.at(0) = da[0];
    trans.at(1) = da[1];
    trans.at(2) = da[2];
    trans.at(3) = da[3];

    double tmp = sqrt((gamma-1)/(gamma+1));
    trans.at(4) = da[4]/tmp;
    trans.at(5) = da[5]*tmp;
}

