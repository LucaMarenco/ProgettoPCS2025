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

vector<Vector3d> punti_triangolazione(const Vector3d& A, const Vector3d& B, const Vector3d& C, int b);

bool file_vertici(const vector<Vector3d>& points, 
				  map<array<int,3> , int>& mappa_vertici, 
				  int& id_vertice, 
				  ostream& s_g_Cell0Ds, bool dual);

bool file_lati(const vector<Vector3d>& points, 
			   map<pair<array<int,3>, array<int,3>>, int>& mappa_lati,
			   map<array<int,3> , int>& mappa_vertici,
			   int& id_lato,
			   int b,
			   ostream& s_g_Cell1Ds, bool duale);
			   
bool file_facce(const vector<Vector3d>& points,
				map<int, pair<Vector3i, Vector3i>>& mappa_facce,
				map<pair<array<int,3>, array<int,3>>, int>& mappa_lati,
			    map<array<int,3> , int>& mappa_vertici,
				int& id_faccia,
				int b);

bool file_poliedro(int& F_s_g,
				   int& V_s_g,
				   int& L_s_g,
				   ostream& s_g_Cell3Ds);		   

MatrixXd Cell0DsConverter(int V_s_g, map<array<int,3> , int> mappa_vertici);

MatrixXi Cell1DsConverter(int L_s_g, map<array<int,3> , int> mappa_vertici , map<pair<array<int,3>, array<int,3>>, int> mappa_lati) ;

optional<Vector3d> calcola_intersezione(const Vector3d& A, const Vector3d& B, const Vector3d& C, const Vector3d& D);

vector<Vector3d> punti_lungo_i_lati(int b, const Vector3d& A, const Vector3d& B, const Vector3d& C);

vector<Vector3d> punti_triangolazione_II(const Vector3d& A, const Vector3d& B, const Vector3d& C, int b);

vector<Vector3d> trova_k_punti_vicini(const Vector3d& punto, const vector<Vector3d>& punti, size_t k);

bool file_lati_II(const vector<Vector3d>& punti_unici,
                  map<pair<array<int,3>, array<int,3>>, int>& mappa_lati,
                  map<array<int,3> , int>& mappa_vertici,
                  int& id_lato,
                  int& b,
                  ostream& s_g_Cell1Ds, bool duale);

vector<int> dijkstra(int n, vector<vector<int>>& adiac_nodi, vector<vector<double>>& adiac_pesi, int start, int end) ;

bool file_facce_II(const vector<Vector3d>& punti_unici,
				map<array<array<int, 3>, 3>, int>& mappa_facce_2,
				map<pair<array<int,3>, array<int,3>>, int>& mappa_lati,
			        map<array<int,3> , int>& mappa_vertici,
				int& id_faccia,
				int& b,
				ostream& s_g_Cell2Ds,
				map<int, pair<Vector3i, Vector3i>>& mappa_facce , bool duale);
				
pair<vector<Vector3d>, map< int, Vector3d>> file_vertici_duale(int F_s_g, map<int, pair<Vector3i, Vector3i>>& mappa_facce, map<array<int,3> , int>& mappa_vertici, map<array<int,3> , int>& mappa_vertici_duale, ofstream& s_g_Cell0Ds );
map<pair<array<int,3>, array<int,3>>, int> file_lati_duale(vector<Vector3d>& baricentri, map<array<int,3>, int>& mappa_vertici_duale, int& id_lato, ofstream& s_g_Cell1Ds) ;
map<int, pair<vector<int>, vector<int>>> file_facce_duale( map<int, pair<Vector3i, Vector3i>>& mappa_facce, map<pair<array<int,3>, array<int,3>>, int>& mappa_lati_duale, map<array<int,3> , int>& mappa_vertici, map< int, Vector3d>& mappa_baricentri, int& id_faccia_duale,ofstream& s_g_Cell2Ds);