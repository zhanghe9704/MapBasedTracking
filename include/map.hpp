/**************************************************

***************************************************/

#ifndef MAP_HPP
#define MAP_HPP

#include <vector>
#include <iostream>

using std::vector;

class Map{
	const int dim;   //Dimension of the problem
	int nTerm;  //number of terms in the map
	int order;  //highest order
	std::vector<std::vector<int> > orders;	//orders of all variables for each term
	std::vector<double> coef;		//coef for each term

	public:
		int addTerm(){ nTerm += 1; return nTerm;}
		int getDim(){return dim;}
		int getTerm(){return nTerm;		}
		void setTotalOrder(int n) {order = n;}
		int getTotalOrder() {return order;}
		void addOrders(const int i, const int order){ orders[i].push_back(order);		}
		void setOrders(const int i, const int j, const int order){orders[i][j]=order;		}
		void setTerm(const int n){ nTerm = n;		}
		void addCoef(double value){coef.push_back(value);		}
		bool coefEmpty(){return coef.empty();		}	//Check whether the coef vector is empty
		double getCoef(const int i){return coef[i];		}
		int getOrders(const int i, const int j){ return orders[i][j];		}
		void outputMap();							//Output the map
		int coefSize(){return coef.size();		}   //Check the number of coef, output length of the coef vector
		void reset();
		Map(const int dim, const int maxTerm);

};

//Output a map
std::ostream& operator<<(std::ostream& os, Map& map);

//Apply the map once
void applyMap(const int dim, double * xi, vector<Map> &map, double * xf);

//Read the truncated map from the COSY output file.
void readCOSYMap(int dim, char fileName[], vector<Map> &myMap);

//Read the truncated map from MadX output file.
void readMadXMap(int dim, char filename[], vector<Map> &my_map);

//Read the 5D truncated map from MadX output file.
void readMadX5DMap(int dim, char filename[], vector<Map> &my_map, bool rmv_const);

#endif
