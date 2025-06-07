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
#include <queue>
#include <climits>
#include <set>
#include <algorithm>

using namespace std;
using namespace Eigen;

array<int,3> to_array(const Vector3d& v){
	double eps = 1e-3;
	return {static_cast<int>(v[0]/eps), static_cast<int>(v[1]/eps), static_cast<int>(v[2]/eps)};
}
	
vector<Vector3d> punti_triangolazione(const Vector3d& A, const Vector3d& B, const Vector3d& C, int b) {
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
				  ostream& s_g_Cell0Ds) {		  
	for(size_t z = 0; z < points.size(); z++) {
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
			   ostream& s_g_Cell1Ds) {				   
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
				int& b) {		
	Vector3i id_vertici_faccia;
	Vector3i id_lati_faccia;	
	int d = 0;
	for (int i = 0; i < b; i++) {   
		for (int j = d; j < d + b - i; j++){
			auto key_j = to_array(points[j]);
			auto key_j1 = to_array(points[j+1]);
			auto key_jb = to_array(points[j+b+1-i]);
							
			if(i == 0){
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
				   ostream& s_g_Cell3Ds) {
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


MatrixXd Cell0DsConverter(int V_s_g, map<array<int,3> , int> mappa_vertici) {
	MatrixXd Cell0DsCoordinates(3, V_s_g);
	for (const auto& [coord,id] : mappa_vertici){
		for (int d = 0; d < 3; ++d) {
			Cell0DsCoordinates(d,id) = coord[d];
		}
	}		
	return Cell0DsCoordinates;
}


MatrixXi Cell1DsConverter(int L_s_g, map<array<int,3> , int> mappa_vertici, map<pair<array<int,3>, array<int,3>>, int> mappa_lati) {
	MatrixXi Cell1DsExtrema(2, L_s_g);
	for (const auto& [coord,id] : mappa_lati){
		Cell1DsExtrema(0,id) = mappa_vertici[coord.first];
		Cell1DsExtrema(1,id) = mappa_vertici[coord.second];
	}		
	return Cell1DsExtrema;
}



optional<Vector3d> calcola_intersezione(const Vector3d& A, const Vector3d& B, const Vector3d& C, const Vector3d& D)  
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
        return std::nullopt; // nessuna intersezione sui segmenti
    }

    Vector3d punto_intersezione = A + t * dir1;
    return punto_intersezione; // Intersezione tra AB e CD
}

vector<Vector3d> punti_lungo_i_lati(int b, const Vector3d& A, const Vector3d& B, const Vector3d& C) { // Per prendere i punti che si vanno a formare
	vector<Vector3d> points;                                                     // lungo i lati della faccia che ho triangolato                              
	points.push_back(A);
	double frequenza = 2 * b;
	for (int i = 1; i < 2 * b; i++) {
		points.push_back(A * (1 - i / (double)frequenza) + B * (i / (double)frequenza));
	}
	points.push_back(B);
	for (int i = 1; i < 2 * b; i++) {
		points.push_back(B * (1 - i / (double)frequenza) + C * (i / (double)frequenza));
	}
	points.push_back(C);
	for (int i = 1; i < 2 * b; i++) {
		points.push_back(C * (1 - i / (double)frequenza) + A * (i / (double)frequenza));
	}
	return points;
} 

vector<Vector3d> punti_lungo_i_lati_normalizzati(int b, const Vector3d& A, const Vector3d& B, const Vector3d& C) {
	vector<Vector3d> punti_l_l;
	punti_l_l.push_back(A);
	double frequenza = 2 * b;
	for (int i = 1; i < 2 * b; i++) {
		punti_l_l.push_back((A * (1 - i / (double)frequenza) + B * (i / (double)frequenza)).normalized());
	}
	punti_l_l.push_back(B);
	for (int i = 1; i < 2 * b; i++) {
		punti_l_l.push_back((B * (1 - i / (double)frequenza) + C * (i / (double)frequenza)).normalized());
	}
	punti_l_l.push_back(C);
	for (int i = 1; i < 2 * b; i++) {
		punti_l_l.push_back((C * (1 - i / (double)frequenza) + A * (i / (double)frequenza)).normalized());
	}
	return punti_l_l;
} 

					
vector<Vector3d> punti_triangolazione_II(const Vector3d& A, const Vector3d& B, const Vector3d& C, int b) {
	vector<Vector3d> points = punti_lungo_i_lati_normalizzati(b, A, B, C); // Punti lungo i lati della faccia già proiettati sulla sfera
	vector<Vector3d> punti_lungo_lato = punti_lungo_i_lati(b, A, B, C);
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
			for(size_t k = 0; k < verticali.size(); k++) {
				auto intersezione = calcola_intersezione(punti_lungo_lato[i], punti_lungo_lato[3 * b - d_inclinato_1], verticali[k].first, verticali[k].second);  		
				if (intersezione) {
					points.push_back(intersezione->normalized());
				}
			}
			d_inclinato_1++;
		}
		else if (i % 2 == 0 && i > 4 * b) {
			for(size_t k = 0; k < verticali.size(); k++) { // Qua prendo i lati obliqui che hanno come uno degli estremi i punti sulla base
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


// Questa triangolazione non è normalizzata e non contiene duplicati

vector<Vector3d> punti_triangolazione_II_n_n(const Vector3d& A, const Vector3d& B, const Vector3d& C, int b) {
	vector<Vector3d> points = punti_lungo_i_lati(b, A, B, C); // Punti lungo i lati della faccia già proiettati sulla sfera
	vector<Vector3d> punti_lungo_lato = punti_lungo_i_lati(b, A, B, C);
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
			for(size_t k = 0; k < verticali.size(); k++) {
				auto intersezione = calcola_intersezione(punti_lungo_lato[i], punti_lungo_lato[3 * b - d_inclinato_1], verticali[k].first, verticali[k].second);  		
				if (intersezione) {
					points.push_back(*intersezione);
				}
			}
			d_inclinato_1++;
		}
		else if (i % 2 == 0 && i > 4 * b) {
			for(size_t k = 0; k < verticali.size(); k++) { // Qua prendo i lati obliqui che hanno come uno degli estremi i punti sulla base
				auto intersezione = calcola_intersezione(punti_lungo_lato[i], punti_lungo_lato[4 * b - d_inclinato_2], verticali[k].first, verticali[k].second);
				if (intersezione) {
					points.push_back(*intersezione);
				}
			}
			d_inclinato_2++;
		}		
	}
	vector<Vector3d> punti_unici;
	map<array<int,3> , int> mappa_vertici;
	int id_vertice = 0;
	for(size_t z = 0; z < points.size(); z++) {
		array<int,3> key = to_array(points[z]);    
		if (mappa_vertici.find(key) == mappa_vertici.end()) {
			mappa_vertici[key] = id_vertice;
			punti_unici.push_back(points[z]);									
			id_vertice++;
		}
	}
	return punti_unici;
}

vector<Vector3d> trova_punti_vicini(const Vector3d& punto, const vector<Vector3d>& punti) {
    vector<pair<double, size_t>> distanze_con_indice;

    for (size_t i = 0; i < punti.size(); ++i) {
        if ((punti[i] - punto).norm() > 1e-12) {
            double distanza = (punti[i] - punto).norm();
            distanze_con_indice.emplace_back(distanza, i);
        }
    }

    sort(distanze_con_indice.begin(), distanze_con_indice.end());

    // Prendi indici dei 6 punti più vicini
    vector<size_t> vicini_indici;
    for (size_t i = 0; i < min(size_t(6), distanze_con_indice.size()); ++i) {
        vicini_indici.push_back(distanze_con_indice[i].second);
    }

    // Ordina i 6 punti in modo tale che ogni punto sia il più vicino possibile al precedente
    vector<Vector3d> ordinati;
    set<size_t> usati;

    // Parto dal primo
    size_t attuale = vicini_indici[0];
    ordinati.push_back(punti[attuale]);
    usati.insert(attuale);

    for (size_t step = 1; step < vicini_indici.size(); ++step) {
        double min_d = numeric_limits<double>::max();
        size_t prossimo = -1;

        for (size_t idx : vicini_indici) {
            if (usati.count(idx)) continue;
            double d = (punti[attuale] - punti[idx]).norm();
            if (d < min_d) {
                min_d = d;
                prossimo = idx;
            }
        }

        if (prossimo != size_t(-1)) {
            ordinati.push_back(punti[prossimo]);
            usati.insert(prossimo);
            attuale = prossimo;
        }
    }

    return ordinati;
}


bool file_lati_II(const vector<Vector3d>& punti_unici, 
                 map<pair<array<int,3>, array<int,3>>, int>& mappa_lati,
                 map<array<int,3>, int>& mappa_vertici,
                 int& id_lato,
                 int& b,
                 ostream& s_g_Cell1Ds) 
{
    // Lati sul bordo della faccia
    for (int i = 0; i < 6 * b - 1; i++) {
        auto key_j = to_array(punti_unici[i].normalized());
        auto key_j1 = to_array(punti_unici[i + 1].normalized());
        if (mappa_lati.find({key_j, key_j1}) == mappa_lati.end() &&
            mappa_lati.find({key_j1, key_j}) == mappa_lati.end()) 
        {
            mappa_lati.insert({{key_j, key_j1}, id_lato});
            mappa_lati.insert({{key_j1, key_j}, id_lato});
            s_g_Cell1Ds << id_lato << " " << mappa_vertici[key_j] << " " << mappa_vertici[key_j1] << "\n";
            id_lato++;
        }
    }

    // Chiudo il ciclo con il lato tra primo e ultimo punto sul bordo

	auto key_j = to_array(punti_unici[0].normalized());
	auto key_j1 = to_array(punti_unici[6 * b - 1].normalized());
	if (mappa_lati.find({key_j, key_j1}) == mappa_lati.end() &&
		mappa_lati.find({key_j1, key_j}) == mappa_lati.end()) 
	{
		mappa_lati.insert({{key_j, key_j1}, id_lato});
		mappa_lati.insert({{key_j1, key_j}, id_lato});
		s_g_Cell1Ds << id_lato << " " << mappa_vertici[key_j] << " " << mappa_vertici[key_j1] << "\n";
		id_lato++;
	}
	for(size_t j = 0; j < punti_unici.size() - 6 * b; j++) {
		vector<Vector3d> vicini = trova_punti_vicini(punti_unici[j + 6*b], punti_unici);
		auto key_j = to_array(punti_unici[6 * b + j].normalized());
		for(size_t k = 0; k < vicini.size(); k++){
			auto key_j1 = to_array(vicini[k].normalized());
			if (mappa_lati.find({key_j, key_j1}) == mappa_lati.end() &&
				mappa_lati.find({key_j1, key_j}) == mappa_lati.end()) 
			{
				mappa_lati.insert({{key_j1, key_j}, id_lato});
				mappa_lati.insert({{key_j, key_j1}, id_lato});
				s_g_Cell1Ds << id_lato << " " << mappa_vertici[key_j] << " " << mappa_vertici[key_j1] << "\n";
				id_lato++;
			}
		}
		
	}
    return true;
}

using P = pair<double, int>;

vector<int> dijkstra(int n, vector<vector<int>>& adiac_nodi, vector<vector<double>>& adiac_pesi, int start, int end) {
    vector<double> dist(n, numeric_limits<double>::infinity());
    vector<int> pred(n, -1);
    priority_queue<P, vector<P>, greater<>> pq;

    dist[start] = 0.0;
    pq.emplace(0.0, start);

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();

        if (d > dist[u]) continue;

        for (size_t i = 0; i < adiac_nodi[u].size(); ++i) {
            int v = adiac_nodi[u][i];
            double w = adiac_pesi[u][i];

            if (dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w;
                pred[v] = u;
                pq.emplace(dist[v], v);
            }
        }
    }

    // Ricostruzione cammino
    vector<int> path;
    if (dist[end] == numeric_limits<double>::infinity()) {
        return path; // Nessun cammino
    }

    for (int v = end; v != -1; v = pred[v]) {
        path.push_back(v);
    }

    reverse(path.begin(), path.end());
    return path;
}

bool file_facce_II(const vector<Vector3d>& punti_unici,
				map<array<array<int, 3>, 3>, int>& mappa_facce,
				map<pair<array<int,3>, array<int,3>>, int>& mappa_lati,
			    map<array<int,3> , int>& mappa_vertici,
				int& id_faccia,
				int& b,
				ostream& s_g_Cell2Ds) {	
	for(size_t i = 0; i < punti_unici.size() - 6*b; i++) {
		vector<Vector3d> vicini = trova_punti_vicini(punti_unici[i + 6*b], punti_unici);
		auto key_j = to_array(punti_unici[6*b + i].normalized());
		for(size_t j = 0; j < vicini.size() - 1; j++) {
			auto key_j1 = to_array(vicini[j].normalized());
			auto key_j2 = to_array(vicini[j + 1].normalized());
			if(mappa_facce.find({key_j, key_j1,key_j2}) == mappa_facce.end() && mappa_facce.find({key_j, key_j2,key_j1}) == mappa_facce.end() &&
				mappa_facce.find({key_j1, key_j,key_j2}) == mappa_facce.end() && mappa_facce.find({key_j1, key_j2,key_j}) == mappa_facce.end() &&
				mappa_facce.find({key_j2, key_j1,key_j}) == mappa_facce.end() && mappa_facce.find({key_j2, key_j,key_j1}) == mappa_facce.end())	{
				mappa_facce.insert({{key_j,key_j1,key_j2},id_faccia});
				mappa_facce.insert({{key_j,key_j2,key_j1},id_faccia});
				mappa_facce.insert({{key_j1,key_j2,key_j},id_faccia});
				mappa_facce.insert({{key_j1,key_j,key_j2},id_faccia});
				mappa_facce.insert({{key_j2,key_j1,key_j},id_faccia});
				mappa_facce.insert({{key_j2,key_j,key_j1},id_faccia});
				s_g_Cell2Ds << id_faccia << " 3 " << "3 " << mappa_vertici[key_j] << " "<<  mappa_vertici[key_j1] << " "<<
				mappa_vertici[key_j2]<<" "<< mappa_lati[{key_j,key_j1}]<<" "<<mappa_lati[{key_j1,key_j2}]<<" "<<mappa_lati[{key_j2,key_j}]<<"\n";
				id_faccia++;
			}
		}
		auto key_j1 = to_array(vicini[0].normalized());
		auto key_j2 = to_array(vicini[5].normalized());
		if(mappa_facce.find({key_j, key_j1,key_j2}) == mappa_facce.end() && mappa_facce.find({key_j, key_j2,key_j1}) == mappa_facce.end() &&
			mappa_facce.find({key_j1, key_j,key_j2}) == mappa_facce.end() && mappa_facce.find({key_j1, key_j2,key_j}) == mappa_facce.end() &&
			mappa_facce.find({key_j2, key_j1,key_j}) == mappa_facce.end() && mappa_facce.find({key_j2, key_j,key_j1}) == mappa_facce.end())	{
			mappa_facce.insert({{key_j,key_j1,key_j2},id_faccia});
			mappa_facce.insert({{key_j,key_j2,key_j1},id_faccia});
			mappa_facce.insert({{key_j1,key_j2,key_j},id_faccia});
			mappa_facce.insert({{key_j1,key_j,key_j2},id_faccia});
			mappa_facce.insert({{key_j2,key_j1,key_j},id_faccia});
			mappa_facce.insert({{key_j2,key_j,key_j1},id_faccia});
			s_g_Cell2Ds << id_faccia << " 3 " << "3 " << mappa_vertici[key_j] << " "<<  mappa_vertici[key_j1] << " "<<
			mappa_vertici[key_j2]<<" "<< mappa_lati[{key_j,key_j1}]<<" "<<mappa_lati[{key_j1,key_j2}]<<" "<<mappa_lati[{key_j2,key_j}]<<"\n";
			id_faccia++;
		}
	}

	return true;
}

/*
bool file_vertici_duale(int F_s_g, ) {
	vector<Vector3d> baricentri(F_s_g);
	double eps = 1e-3;
	for(int i = 0; i < F_s_g; i++) {
		Vector3d id_vertici_faccia = mappa_facce[i].first();
		for (const auto& [key, valore] : mappa_vertici) {
			if(key == id_vertici_faccia[0]) {
				Vector3d A(valore[0]*eps, valore[1]*eps, valore[2]*eps);
			}
			if(key == id_vertici_faccia[1]) {
				Vector3d B(valore[0]*eps, valore[1]*eps, valore[2]*eps);
			}
			if(key == id_vertici_faccia[2]) {
				Vector3d C(valore[0]*eps, valore[1]*eps, valore[2]*eps);
			}	
		}
		Vector3d baricentro = punti_triangolazione_II(A, B, C, 1).end();
		baricentri[i] = baricentro;
	}
	for(int j = 0; j < baricentri.size(); j++) {
		s_g_Cell0Ds << j << " " << baricentro[0] << " " << baricentro[1] << " " baricentro[2] << "\n";
	}
}
*/