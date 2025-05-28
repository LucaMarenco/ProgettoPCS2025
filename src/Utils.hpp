#pragma once

#include <iostream>
#include <array>
#include <vector>
#include <map>
#include <Eigen/Dense>
#include <cmath>


using namespace std;
using namespace Eigen; 

namespace Utils{

//array<int,3> to_array(const Vector3d& v);

vector<Vector3d> punti_triangolazione (Vector3d A, Vector3d B, Vector3d C, int b);

void file_vertici(vector<Vector3d>& points, map<array<int,3> , int>& mappa_vertici, int& id_vertice, ofstream& s_g_Cell0Ds);
}