#pragma once

#include <iostream>
#include <array>
#include <vector>
#include <map>
#include <Eigen/Dense>
#include <cmath>
#include <optional>


using namespace std;
using namespace Eigen; 



array<int,3> to_array(const Vector3d& v);

vector<Vector3d> punti_triangolazione(Vector3d A, Vector3d B, Vector3d C, int b);

bool file_vertici(const vector<Vector3d>& points, 
				  map<array<int,3> , int>& mappa_vertici, 
				  int& id_vertice, 
				  ofstream& s_g_Cell0Ds);

bool file_lati(const vector<Vector3d>& points, 
			   map<pair<array<int,3>, array<int,3>>, int>& mappa_lati,
			   map<array<int,3> , int>& mappa_vertici,
			   int& id_lato,
			   int& b,
			   ofstream& s_g_Cell1Ds);
			   
bool file_facce(const vector<Vector3d>& points,
				map<int, pair<Vector3i, Vector3i>>& mappa_facce,
				map<pair<array<int,3>, array<int,3>>, int>& mappa_lati,
			    map<array<int,3> , int>& mappa_vertici,
				int& id_faccia,
				int& b,
				ofstream& s_g_Cell2Ds);

bool file_poliedro(int& F_s_g,
				   int& V_s_g,
				   int& L_s_g,
				   ofstream& s_g_Cell3Ds);

optional<Vector3d> calcola_intersezione(Vector3d A, Vector3d B, Vector3d C, Vector3d D);

void punti_lungo_il_lato(int b, Vector3d A, Vector3d B, vector<Vector3d>& punti_l_l);

vector<Vector3d> punti_triangolazione_II(Vector3d A, Vector3d B, Vector3d C, int b);

bool file_vertici_II(const vector<Vector3d>& points, 
				  map<array<int,3> , int>& mappa_vertici, 
				  int& id_vertice, 
				  ofstream& s_g_Cell0Ds);

void punti_lungo_il_lato_normalizzati(int b, Vector3d A, Vector3d B, vector<Vector3d>& punti_l_l);
