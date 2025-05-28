#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
#include <map>
#include <Eigen/Dense>
#include <cmath>

using namespace std;
using namespace Eigen;
namespace Utils{

	/*array<int,3> to_array(const Vector3d& v){
		double eps = 1e-3;
		return {static_cast<int>(v[0]/eps),static_cast<int>(v[1]/eps),static_cast<int>(v[2]/eps)};
	}*/
		
	vector<Vector3d> punti_triangolazione (Vector3d A, Vector3d B, Vector3d C, int b) {
		vector<Vector3d> points;
		for (int k = 0; k <= b; k++) {
			for (int j = 0; j <= b - k; j++) {
				double c_A = 1.0 - (double)k / b - (double)j / b;
				double c_B = (double)k / b;
				double c_C = (double)j / b;
				Vector3d P = c_A * A + c_B * B + c_C * C;
				points.push_back(P.normalized());
			}
		}
		return(points);
	}

	void file_vertici(const vector<Vector3d>& points, 
					  map<array<int,3> , int>& mappa_vertici, 
					  int& id_vertice, 
					  ofstream& s_g_Cell0Ds) {
						  
		for(int z = 0; z < points.size(); z++) {
			double eps = 1e-3;
			array<int,3> key = {static_cast<int>(points[z][0]/eps),static_cast<int>(points[z][1]/eps),static_cast<int>(points[z][2]/eps)};    
			if (mappa_vertici.find(key) == mappa_vertici.end()) {
				mappa_vertici[key] = id_vertice;
				double eps = 1e-3;
				s_g_Cell0Ds << id_vertice << " " << key[0]*eps << " " << key[1]*eps << " " << key[2]*eps << "\n";
				id_vertice++;							// assegna nuovo ID a vertice se non esiste						
			}
		}
	}
}