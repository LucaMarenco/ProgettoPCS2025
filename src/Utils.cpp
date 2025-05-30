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

array<int,3> to_array(const Vector3d& v){
	double eps = 1e-3;
	return {static_cast<int>(v[0]/eps),static_cast<int>(v[1]/eps),static_cast<int>(v[2]/eps)};
}
	
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

bool file_vertici(const vector<Vector3d>& points, 
				  map<array<int,3> , int>& mappa_vertici, 
				  int& id_vertice, 
				  ofstream& s_g_Cell0Ds) {
					  
	for(int z = 0; z < points.size(); z++) {
		double eps = 1e-3;
		array<int,3> key = to_array(points[z]);    
		if (mappa_vertici.find(key) == mappa_vertici.end()) {
			mappa_vertici[key] = id_vertice;
			double eps = 1e-3;
			s_g_Cell0Ds << id_vertice << " " << key[0]*eps << " " << key[1]*eps << " " << key[2]*eps << "\n";
			id_vertice++;							// assegna nuovo ID a vertice se non esiste						
		}
	}
	return true;
}

bool file_lati(const vector<Vector3d>& points, 
			   map<pair<array<int,3>, array<int,3>>, int>& mappa_lati,
			   map<array<int,3> , int>& mappa_vertici,
			   int& id_lato,
			   int& b,
			   ofstream& s_g_Cell1Ds) {
				   
	int d = 0;
	for (int f = 0; f <= b; f++) {   
		for (int j = d; j < d + b - f; j++) {
			auto key_j = to_array(points[j]);
			auto key_j1 = to_array(points[j+1]);
			auto key_jb = to_array(points[j+b+1-f]);
			
			if (mappa_lati.find({min({key_j,key_jb}), max({key_j,key_jb})}) != mappa_lati.end()){}
			else{
				mappa_lati.insert({{min({key_j,key_jb}), max({key_j,key_jb})},id_lato});
				s_g_Cell1Ds << id_lato << " " <<mappa_vertici[min({key_j,key_jb})] << " "<<  mappa_vertici[max({key_j,key_jb})]  <<"\n";
				id_lato++;
			}
			
			if (mappa_lati.find({min({key_j,key_j1}), max({key_j,key_j1})}) != mappa_lati.end()){}
			else{
				mappa_lati.insert({{min({key_j,key_j1}), max({key_j,key_j1})} ,id_lato});
				s_g_Cell1Ds << id_lato << " " << mappa_vertici[min({key_j,key_j1})]  << " " << mappa_vertici[max({key_j,key_j1})] << "\n";
				id_lato++;
			}
			
			if (mappa_lati.find({min({key_jb,key_j1}), max({key_jb,key_j1})}) != mappa_lati.end()){}
			else{
				mappa_lati.insert({{min({key_jb,key_j1}), max({key_jb,key_j1})},id_lato});
				s_g_Cell1Ds << id_lato << " " <<mappa_vertici[min({key_jb,key_j1})] << " " << mappa_vertici[max({key_jb,key_j1})]<< "\n";
				id_lato++;
			}	
		}
		d = d + b + 1 - f;	
	}
	return true;
}			
			   
bool file_facce(const vector<Vector3d>& points,
				map<int, pair<Vector3i, Vector3i>>& mappa_facce,
				map<pair<array<int,3>, array<int,3>>, int>& mappa_lati,
			    map<array<int,3> , int>& mappa_vertici,
				int& id_faccia,
				int& b,
				ofstream& s_g_Cell2Ds) {
	
	
	Vector3i id_vertici_faccia;
	Vector3i id_lati_faccia;
	
	int d = 0;
	for (int i = 0; i < b; i++) {   
		for (int j = d; j < d + b - i; j++){
			auto key_j = to_array(points[j]);
			auto key_j1 = to_array(points[j+1]);
			auto key_jb = to_array(points[j+b+1-i]);
							
			if( i==0){
				id_vertici_faccia = Vector3i{mappa_vertici[key_j], mappa_vertici[key_j1], mappa_vertici[key_jb]}; 
				if (mappa_lati.find({min({key_j,key_jb}),max({key_j,key_jb})}) != mappa_lati.end()){
					id_lati_faccia= Vector3i{mappa_lati[{min({key_j,key_j1}),max({key_j,key_j1})}], mappa_lati[{min({key_j1,key_jb}),max({key_j1,key_jb})}], mappa_lati[{min({key_j,key_jb}),max({key_j,key_jb})}]};
					mappa_facce.insert({id_faccia, {id_vertici_faccia, id_lati_faccia}}); 
				}
				id_faccia++;
			}
			else {
				auto key_jb_m = to_array(points[j-b-1+i]);
				
				id_vertici_faccia = Vector3i{ mappa_vertici[key_j], mappa_vertici[key_j1], mappa_vertici[key_jb]};  
				id_lati_faccia = Vector3i{mappa_lati[{min({key_j,key_j1}),max({key_j,key_j1})}], mappa_lati[{min({key_j1,key_jb}),max({key_j1,key_jb})}], mappa_lati[{min({key_j,key_jb}),max({key_j,key_jb})}]};
				mappa_facce.insert({id_faccia, {id_vertici_faccia, id_lati_faccia}}); 
				id_faccia++;
				
				id_vertici_faccia = Vector3i{ mappa_vertici[key_j], mappa_vertici[key_j1], mappa_vertici[key_jb_m]}; 
				id_lati_faccia = Vector3i{ mappa_lati[{min({key_j,key_j1}),max({key_j,key_j1})}], mappa_lati[{min({key_j1,key_jb_m}),max({key_j1,key_jb_m})}], mappa_lati[{min({key_j,key_jb_m}),max({key_j,key_jb_m})}]};
				mappa_facce.insert({id_faccia, {id_vertici_faccia, id_lati_faccia}});
				id_faccia++;	
			}
		}
		d = d + b + 1 - i;
	}
	return true;
}

bool file_poliedro(int& F_s_g,
				   int& V_s_g,
				   int& L_s_g,
				   ofstream& s_g_Cell3Ds) {
	s_g_Cell3Ds << 0 << " " << V_s_g << " " << L_s_g << " " << F_s_g << " ";
	for(int i = 0; i < V_s_g; i++) {
		 s_g_Cell3Ds << i << " ";
	}
	for(int i = 0; i < L_s_g; i++) {
		s_g_Cell3Ds << i << " ";
	}
	for(int i = 0; i < F_s_g; i++) {
		s_g_Cell3Ds << i << " ";
	}
	return true;
}
					