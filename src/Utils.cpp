#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
#include <map>
#include <Eigen/Dense>
#include <cmath>
#include <optional>

using namespace std;
using namespace Eigen;

array<int,3> to_array(const Vector3d& v){
	double eps = 1e-3;
	return {static_cast<int>(v[0]/eps), static_cast<int>(v[1]/eps), static_cast<int>(v[2]/eps)};
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

optional<Vector3d> calcola_intersezione(Vector3d A, Vector3d B, Vector3d C, Vector3d D)  
{
    Vector3d dir1 = B - A;  // Vettore AB
    Vector3d dir2 = D - C;  // Vettore CD
    Vector3d AC = A - C;    // Vettore tra A e C

    double a = dir1.dot(dir1);
    double b = dir1.dot(dir2);
    double c = dir2.dot(dir2);
    double d = dir1.dot(AC);
    double e = dir2.dot(AC);

    double denominatore = a * c - b * b;
    double t = (b * e - c * d) / denominatore;
    double s = (a * e - b * d) / denominatore;

    // Controllo se l'intersezione cade dentro i segmenti [0,1]
    if (t < 0.0 || t > 1.0 || s < 0.0 || s > 1.0) {
        return nullopt; // nessuna intersezione sui segmenti (non restituisce niente)
    }

    Vector3d punto_intersezione = A + t * dir1;
    return punto_intersezione; // Intersezione tra AB e CD
}

void punti_lungo_il_lato(int b, Vector3d A, Vector3d B, vector<Vector3d>& punti_l_l) { // Per prendere i punti che si vanno a formare
	double lunghezza_lato = (B - A).norm();                                            // lungo i lati della faccia che ho triangolato
	punti_l_l.push_back(A);
	double suddivisione = lunghezza_lato / (2 * b);
	for (int i = 1; i < 2 * b; i++) {
		punti_l_l.push_back(((suddivisione * i) * (B - A) / lunghezza_lato) + A);
	}
} 

void punti_lungo_il_lato_normalizzati(int b, Vector3d A, Vector3d B, vector<Vector3d>& punti_l_l) {
	double lunghezza_lato = (B - A).norm();
	punti_l_l.push_back(A);
	double suddivisione = lunghezza_lato / (2 * b);
	for (int i = 1; i < 2 * b; i++) {
		punti_l_l.push_back(((suddivisione * i * (B - A) / lunghezza_lato) + A).normalized());
	}
} 

					
vector<Vector3d> punti_triangolazione_II(Vector3d A, Vector3d B, Vector3d C, int b) {
	vector<Vector3d> points;
	punti_lungo_il_lato_normalizzati(b, A, B, points);
	punti_lungo_il_lato_normalizzati(b, B, C, points); // Punti lungo i lati della faccia già proiettati sulla sfera
	punti_lungo_il_lato_normalizzati(b, C, A, points);
	vector<Vector3d> punti_lungo_lato;
	punti_lungo_il_lato(b, A, B, punti_lungo_lato);
	punti_lungo_il_lato(b, B, C, punti_lungo_lato); // Punti lungo i lati della faccia
	punti_lungo_il_lato(b, C, A, punti_lungo_lato);
	vector<pair<Vector3d, Vector3d>> verticali; // Coppie di estremi che formano i lati verticali
	int d_base = 1;
	for(int j = 2; j < 4 * b - 1; j++) {  // Prendo i lati verticali per poi calcolare l'intersezione di ogni lato obliquo con 
		if(j % 2 == 0) {                  // tutte le verticali prendendo così tutti i punti 
			verticali.push_back({punti_lungo_lato[j], punti_lungo_lato[6 * b - d_base]});
			d_base++;
		}	
	}
	int d_inclinato_1 = 0;
	int d_inclinato_2 = 1;
	for(int i = 0; i < 6 * b; i++) {
		if (i % 2 == 0 && i < 2 * b) { // Qua prendo i lati obliqui che hanno come uno degli estremi i punti sul lato "sinistro"
			for(int k = 0; k < verticali.size(); k++) {
				auto intersezione = calcola_intersezione(punti_lungo_lato[i], punti_lungo_lato[3 * b - d_inclinato_1], verticali[k].first, verticali[k].second);  		
				if (intersezione) {
					points.push_back(intersezione->normalized());
				}
			}
			d_inclinato_1++;
		}
		else if (i % 2 == 0 && i > 4 * b) {
			for(int k = 0; k < verticali.size(); k++) { // Qua prendo i lati obliqui che hanno come uno degli estremi i punti sulla base
				auto intersezione = calcola_intersezione(punti_lungo_lato[i], punti_lungo_lato[4 * b - d_inclinato_2], verticali[k].first, verticali[k].second);
				if (intersezione) {
					points.push_back(intersezione->normalized());
				}
			}
			d_inclinato_2++;
		}		
	}
	return points;
}

bool file_vertici_II(const vector<Vector3d>& points, 
				  map<array<int,3> , int>& mappa_vertici, 
				  int& id_vertice, 
				  ofstream& s_g_Cell0Ds) {
					  
	for(int z = 0; z < points.size(); z++) {
		double eps = 1e-3;
		array<int,3> key = to_array(points[z]);    
		if (mappa_vertici.find(key) == mappa_vertici.end()) {
			mappa_vertici[key] = id_vertice;
			s_g_Cell0Ds << id_vertice << " " << key[0]*eps << " " << key[1]*eps << " " << key[2]*eps << "\n";
			id_vertice++;							// assegna nuovo ID a vertice se non esiste						
		}
	}
	return true;
}

	

			   

/*
Alternativa per suddividere le facce della triangolazione II
std::vector<Vec3> subdivideFace(const Vec3& A, const Vec3& B, const Vec3& C, int freq) {
    std::vector<Vec3> points;
    for (int i = 0; i <= freq; ++i) {
        Vec3 AB = A * (1 - i / (double)freq) + B * (i / (double)freq);
        Vec3 AC = A * (1 - i / (double)freq) + C * (i / (double)freq);
        for (int j = 0; j <= i; ++j) {
            Vec3 P = AB * (1 - j / (double)i) + AC * (j / (double)i);
            P.normalize(); // proiezione sulla sfera
            points.push_back(P);
        }
    }
    return points;
}*/