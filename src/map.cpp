

#include "../include/map.hpp"
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>


using std::cout;
using std::endl;
using std::vector;

//constructor of the class Map
Map::Map(const int dim, const int maxTerm):dim(dim),orders(maxTerm,vector<int>(2*dim)){
	nTerm = 0;
	coef.reserve(maxTerm);
}

//output the map
void Map::outputMap(){
	for(int i=0;i<nTerm;++i){
		cout << i+1 <<' '<<coef[i];
		for (int j=0;j<2*dim;++j){
			cout<<' '<<orders[i][j];
		}
		cout << endl;
	}

}

void Map::reset() {
    nTerm = 0;
    order = 0;
    std::fill(coef.begin(), coef.end(), 0);
    for (auto& sub : orders) {
        std::fill(sub.begin(), sub.end(), 0);
    }
}

void readCOSYMap(int dim, char fileName[], vector<Map> &myMap){
	std::ifstream infile;
	infile.open(fileName);
	auto iter=myMap.begin();
	int total_order = 0;

	std::string line;
	std::string start = "COEFFICIENT";
	std::string ending = "--";
	while(std::getline(infile, line)) {
        if(line.empty()) continue;
        auto found = line.find(start);
        if(found!=std::string::npos) continue;
        found = line.find(ending);
        if(found!=std::string::npos) {
            iter->setTotalOrder(total_order);
            ++iter;
            total_order = 0;
            continue;
        }

        iter->addTerm();
        double coef = std::stod(line.substr(8,24));
        iter->addCoef(coef);
        int current_order = std::stoi(line.substr(32,4));
        if (current_order>total_order) total_order = current_order;
        for(int i=0; i<2*dim; ++i) {
            int order;
            if(i==0) order = std::stoi(line.substr(36,2));
            else if(i==1) order = std::stoi(line.substr(38,2));
            else if(i==2) order = std::stoi(line.substr(41,2));
            else if(i==3) order = std::stoi(line.substr(43,2));
            else if(i==4) order = std::stoi(line.substr(46,2));
            else if(i==5) order = std::stoi(line.substr(48,2));
            iter->setOrders(iter->getTerm()-1, i, order);
        }
	}
	if(iter!=myMap.end()) iter->setTotalOrder(total_order);
	infile.close();
}

void readMadXMap(int dim, char filename[], vector<Map> &my_map) {
    std::ifstream infile;
	infile.open(filename);
	auto iter=my_map.end();
	int total_order = 0;

	std::string line;
	std::string start = "COEFFICIENT";
	std::string pass_1 = "etall";
	std::string pass_2 = "*****";
	bool skip = false;

	while(std::getline(infile, line)) {
        if (skip) {
            skip = false;
            continue;
        }
        if(line.empty()) continue;
        if(line=="  0.000000000000000E+000") break;
        auto found = line.find(start);
        if(found!=std::string::npos) {
            if(iter==my_map.end()) {
                iter = my_map.begin();
            }
            else {
                iter->setTotalOrder(total_order);
                ++iter;
                total_order = 0;
            }
            continue;
        }

        found = line.find(pass_1);
        if(found!=std::string::npos) continue;
        found = line.find(pass_2);
        if(found!=std::string::npos) continue;
        if(line.find_first_not_of(' ') == std::string::npos) continue;

        iter->addTerm();
        double coef = std::stod(line.substr(8,21));
        iter->addCoef(coef);
        int current_order = std::stoi(line.substr(30,4));
        if (current_order>total_order) total_order = current_order;

        for(int i=0; i<2*dim; ++i) {
            int order;
            if(i==0) order = std::stoi(line.substr(37,2));
            else if(i==1) order = std::stoi(line.substr(39,2));
            else if(i==2) order = std::stoi(line.substr(42,2));
            else if(i==3) order = std::stoi(line.substr(44,2));
            else if(i==4) order = std::stoi(line.substr(47,2));
            else if(i==5) order = std::stoi(line.substr(49,2));
             iter->setOrders(iter->getTerm()-1, i, order);
        }
        skip = true;
	}
	if(iter!=my_map.end()) iter->setTotalOrder(total_order);
	infile.close();
}

void readMadX5DMap(int dim, char filename[], vector<Map> &my_map, bool rmv_const) {
    std::ifstream infile;
	infile.open(filename);
	auto iter=my_map.end();
	int total_order = 0;

	std::string line;
	std::string start = "COEFFICIENT";
	std::string pass_1 = "etall";
	std::string pass_2 = "*****";
	bool skip = false;
	bool check_const = rmv_const;

	while(std::getline(infile, line)) {
        if (skip) {
            skip = false;
            continue;
        }
        if(line.empty()) continue;
        if(line=="  0.000000000000000E+000") break;
        auto found = line.find(start);
        if(found!=std::string::npos) {
            if(iter==my_map.end()) {
                iter = my_map.begin();
            }
            else {
                iter->setTotalOrder(total_order);
                ++iter;
                total_order = 0;
                check_const = rmv_const;
            }
            continue;
        }

        found = line.find(pass_1);
        if(found!=std::string::npos) continue;
        found = line.find(pass_2);
        if(found!=std::string::npos) continue;
        if(line.find_first_not_of(' ') == std::string::npos) continue;

        int current_order = std::stoi(line.substr(30,4));
        if(check_const) {
            if(current_order==0) {
                check_const = false;
                skip = true;
                continue;
            }
        }

        iter->addTerm();
        double coef = std::stod(line.substr(8,21));
        iter->addCoef(coef);
//        int current_order = std::stoi(line.substr(30,4));
        if (current_order>total_order) total_order = current_order;

        for(int i=0; i<2*dim; ++i) {
            int order;
            if(i==0) order = std::stoi(line.substr(37,2));
            else if(i==1) order = std::stoi(line.substr(39,2));
            else if(i==2) order = std::stoi(line.substr(42,2));
            else if(i==3) order = std::stoi(line.substr(44,2));
            else if(i==4) order = std::stoi(line.substr(47,2));
//            else if(i==5) order = std::stoi(line.substr(49,2));
            else if(i==5) order = 0;
             iter->setOrders(iter->getTerm()-1, i, order);
        }
        skip = true;
	}
	if(iter!=my_map.end()) iter->setTotalOrder(total_order);
	infile.close();
}

//Output a map
std::ostream& operator<<(std::ostream& os, Map& map){
	for(int i=0;i<map.getTerm();++i){
		os << i+1 <<' '<<map.getCoef(i);
		for (int j=0;j<2*map.getDim();++j){
			os<<' '<<map.getOrders(i,j);
		}
		os << endl;
	}
	return os;
}

//Apply the map once
void applyMap(const int dim, double * xi, vector<Map> &map, double * xf){
	vector<double> xt(2*dim);
	int order = 0;
	for(auto& m: map)
        if (m.getTotalOrder()>order) order = m.getTotalOrder();
    vector<vector<double> > xn(2*dim, vector<double>(order+1));
    for(int i=0; i<2*dim; ++i) {
        xn[i][0] = 1;
        for(int j=1; j<order+1; ++j) {
            xn[i][j] = xi[i]*xn[i][j-1];
        }
    }

    for (int i=0; i<2*dim; ++i){
		for (int j=0; j<map[i].getTerm(); ++j){
			double term=map[i].getCoef(j);
			for (int k=0; k<2*dim; ++k){
                term *= xn[k][map[i].getOrders(j,k)];
			}
			xt[i]+=term;
		}
	}

	for(int i=0; i<2*dim; ++i){
		xf[i] = xt[i];
	}
}




