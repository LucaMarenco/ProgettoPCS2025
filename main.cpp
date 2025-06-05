/*Consideriamo di partire che abbiamo già dei file con le celle che rappresentano i solidi platonici con p = 3, 
quindi tetraedro, ottaedro e icosaedro; e consideriamo anche di aver già aperto i file in lettura 
nel programma (presi in input).
Dati in input p,q,b,c bisogna innanzitutto, perché sia un solido geodetico, avere p = 3.
Faccio una funzione void (perché non ritorna niente) che, dati in input i 4 interi p,q,b e c,
crea 4 file txt che rappresentano il poliedro geodetico che si ottiene a partire da un certo solido platonico:
*/
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <array>
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <optional>
#include "Utils.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;



int main(int argc, char *argv[])
{
	if (argc == 5 || argc ==7 ){
		int b;
		int c;
		int p;
		int q;
		if (stoi(argv[3]) != 0 && stoi(argv[4]) == 0){
			b = stoi(argv[3]);
			c = stoi(argv[4]);
		}
		else if (stoi(argv[3]) == 0 && stoi(argv[4]) != 0){
			b = stoi(argv[4]);
			c = stoi(argv[3]);
		}
		else if (stoi(argv[3]) >= 1 && stoi(argv[3]) == stoi(argv[4])){
			b = stoi(argv[3]);
			c = stoi(argv[4]);
		}
		
		if (stoi(argv[1]) == 3){
			p = stoi(argv[1]);
			q = stoi(argv[2]);
		}
		else if (stoi(argv[2]) == 3){
			p = stoi(argv[2]);
			q = stoi(argv[1]);
		}
		
		
		
		vector<vector<int>> tCells = {
			{0, 3, 3, 0, 1, 2, 0, 3, 1},
			{1, 3, 3, 0, 1, 3, 0, 4, 2},   // Questo è il Cell2 del tetraedro di partenza
			{2, 3, 3, 0, 3, 2, 2, 5, 1},
			{3, 3, 3, 1, 2, 3, 3, 5, 4}
		};
		
		vector<vector<int>> oCells = {
			{0, 3, 3, 0, 1, 2, 0, 4, 1},
			{1, 3, 3, 0, 2, 3, 1, 5, 2},   // Questo è il Cell2 dell'ottaedro di partenza
			{2, 3, 3, 0, 3, 4, 2, 6, 3},
			{3, 3, 3, 0, 4, 1, 3, 7, 0},
			{4, 3, 3, 5, 2, 1, 9, 4, 8},
			{5, 3, 3, 5, 3, 2, 10, 5, 9},
			{6, 3, 3, 5, 4, 3, 11, 6, 10},
			{7, 3, 3, 5, 1, 4, 8, 7, 11}
		};
		
		vector<vector<int>> iCells = {
			{0, 3, 3, 0, 1, 8, 0, 7, 3},
			{1, 3, 3, 0, 8, 4, 3, 19, 1},   // Questo è il Cell2 dell'icosaedro di partenza
			{2, 3, 3, 0, 4, 5, 1, 18, 2},
			{3, 3, 3, 0, 5, 10, 2, 21, 4},
			{4, 3, 3, 0, 10, 1, 4, 8, 0},
			{5, 3, 3, 1, 10, 7, 8, 25, 6},
			{6, 3, 3, 1, 7, 6, 6, 22, 5},
			{7, 3, 3, 1, 6, 8, 5, 23, 7},
			{8, 3, 3, 2, 9, 3, 12, 16, 9},
			{9, 3, 3, 2, 4, 9, 10, 20, 12},
			{10, 3, 3, 2, 5, 4, 11, 18, 10},
			{11, 3, 3, 2, 11, 5, 13, 28, 11},
			{12, 3, 3, 2, 3, 11, 9, 17, 13},
			{13, 3, 3, 3, 6, 7, 14, 22, 15},
			{14, 3, 3, 3, 7, 11, 15, 26, 17},
			{15, 3, 3, 3, 9, 6, 16, 24, 14},
			{16, 3, 3, 7, 10, 11, 25, 29, 26},
			{17, 3, 3, 4, 8, 9, 19, 27, 20},
			{18, 3, 3, 5, 11, 10, 28, 29, 21},
			{19, 3, 3, 6, 9, 8, 24, 27, 23},
		};

		// Mappa che associa l'ID della faccia al vettore dei suoi vertici
		map<int, vector<int>> tCell2DsVertices;
		map<int, vector<int>> oCell2DsVertices;
		map<int, vector<int>> iCell2DsVertices;
		
		for (const auto& cell : tCells) {
			int face_id = cell[0]; // ID della faccia
			int num_vertices = cell[1]; // Numero di vertici

			// Estrai i vertici dalla riga
			vector<int> tvertices(cell.begin() + 3, cell.begin() + 3 + num_vertices);

			// Aggiungi alla mappa
			tCell2DsVertices[face_id] = tvertices;
		}
		
		for (const auto& cell : oCells) {
			int face_id = cell[0]; // ID della faccia
			int num_vertices = cell[1]; // Numero di vertici

			// Estrai i vertici dalla riga
			vector<int> overtices(cell.begin() + 3, cell.begin() + 3 + num_vertices);

			// Aggiungi alla mappa
			oCell2DsVertices[face_id] = overtices;
		}
		
		for (const auto& cell : iCells) {
			int face_id = cell[0]; // ID della faccia
			int num_vertices = cell[1]; // Numero di vertici

			// Estrai i vertici dalla riga
			vector<int> ivertices(cell.begin() + 3, cell.begin() + 3 + num_vertices);

			// Aggiungi alla mappa
			iCell2DsVertices[face_id] = ivertices;
		}
		
		vector<pair<int, Vector3d>> tVertices = {
			{0, { 0.57735027,  0.57735027,  0.57735027}},
			{1, {-0.57735027, -0.57735027,  0.57735027}},    // Questo è il Cell0 del tetraedro di partenza
			{2, {-0.57735027,  0.57735027, -0.57735027}},
			{3, { 0.57735027, -0.57735027, -0.57735027}}
		}; 
		
		vector<pair<int, Vector3d>> oVertices = {
			{0, {0,  0, 1}},
			{1, {0, -1, 0}},    // Questo è il Cell0 dell'ottaedro di partenza
			{2, { 1, 0, 0}},
			{3, { 0, 1, 0}},
			{4, {-1, 0, 0}},
			{5, {0, 0, -1}}
		};
		
		vector<pair<int, Vector3d>> iVertices = {
			{0, {0.52573111, 0.85065081, 0.00000000}},
			{1, {-0.52573111, 0.85065081, 0.00000000}},
			{2, {0.52573111, -0.85065081, 0.00000000}},
			{3, {-0.52573111, -0.85065081, 0.00000000}},		// Questo è il Cell0 dell'icosaedro di partenza
			{4, {0.85065081, 0.00000000, 0.52573111}},
			{5, {0.85065081, 0.00000000, -0.52573111}},
			{6, {-0.85065081, 0.00000000, 0.52573111}},
			{7, {-0.85065081, 0.00000000, -0.52573111}},
			{8, {0.00000000, 0.52573111, 0.85065081}},
			{9, {0.00000000, -0.52573111, 0.85065081}},
			{10, {0.00000000, 0.52573111, -0.85065081}},
			{11, {0.00000000, -0.52573111, -0.85065081}}
		}; 
		
		map<int, Vector3d> tCell0DsCoordinates; // Mappa che associa l'ID di un vertice alle sue coordinate
		map<int, Vector3d> oCell0DsCoordinates;
		map<int, Vector3d> iCell0DsCoordinates;
		
		for (const auto& [id, coords] : tVertices) {
			tCell0DsCoordinates[id] = coords;
		}
		
		for (const auto& [id, coords] : oVertices) {
			oCell0DsCoordinates[id] = coords;
		}
		
		for (const auto& [id, coords] : iVertices) {
			iCell0DsCoordinates[id] = coords;
		}
		
        int F_s_g = 0;
		int V_s_g = 0;
		int L_s_g = 0;
		MatrixXd Cell0DsCoordinates;
		MatrixXi Cell1DsExtrema;
		if (p == 3) {
			if (b >= 1 && c == 0 ){
				int T = b * b + b * c + c * c; 
				ofstream s_g_Cell0Ds("s_g_Cell0Ds.txt");
				ofstream s_g_Cell1Ds("s_g_Cell1Ds.txt");
				ofstream s_g_Cell2Ds("s_g_Cell2Ds.txt");
				ofstream s_g_Cell3Ds("s_g_Cell3Ds.txt");
				s_g_Cell0Ds << "Id " << "x " << "y " << "z" << "\n";
				s_g_Cell1Ds << "Id " << "start_vertex " << "end_vertex" << "\n"; 
				s_g_Cell2Ds << "Id " << "num_vertices " << "num_edges " << "vertices " << "edges " << "\n";
				s_g_Cell3Ds << "Id " << "num_vertices " << "num_edges " << "num_faces " << "vertices " << "edges " << "faces " << "\n";
				if (q == 3) {
					int F = 4;
					F_s_g = 4 * T;
					V_s_g = 2 * T + 2;
					L_s_g = 6 * T;
					int id_vertice = 0;
					int id_lato = 0;
					int id_faccia = 0;
					map<array<int,3> , int> mappa_vertici;    
					map<pair<array<int,3>, array<int,3>>, int> mappa_lati;
					map<int, pair<Vector3i, Vector3i>> mappa_facce;
					//Siamo nel caso tetraedro:
					//in questo caso il programma deve anche restituirmi un poliedro di Goldberg di classe I:
					for(int i = 0; i < F; i++) {
						int id_A = tCell2DsVertices[i][0];
						int id_B = tCell2DsVertices[i][1]; 
						int id_C = tCell2DsVertices[i][2];
						Vector3d A = tCell0DsCoordinates[id_A];
						Vector3d B = tCell0DsCoordinates[id_B];
						Vector3d C = tCell0DsCoordinates[id_C];	
						vector<Vector3d> points = punti_triangolazione(A,B,C,b);
						if(!file_vertici(points, mappa_vertici, id_vertice, s_g_Cell0Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						};
						
						
						if(!file_lati(points, mappa_lati, mappa_vertici, id_lato, b, s_g_Cell1Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						};	
						
						if(!file_facce(points, mappa_facce, mappa_lati, mappa_vertici, id_faccia, b, s_g_Cell2Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						};					
					}
					
					for (const auto& [key, value] : mappa_facce)
					{
						s_g_Cell2Ds << key << " " << "3 3 " << value.first.transpose()<< " " << value.second.transpose() << endl;
					}
					
					if(!file_poliedro(F_s_g, V_s_g, L_s_g, s_g_Cell3Ds))
					{
						cerr << "errore nella compilazione del file" << endl;
						return 1;
					};
					
					
					Cell0DsCoordinates = Cell0DsConverter(V_s_g, mappa_vertici);
					Cell1DsExtrema = Cell1DsConverter(L_s_g, mappa_vertici, mappa_lati);

				}
				
				else if(q == 4) {
					int F = 8;
					F_s_g = 8 * T;
					V_s_g = 4 * T + 2;
					L_s_g = 12 * T;
					int id_vertice = 0;
					int id_lato = 0;
					int id_faccia = 0;
					map<array<int,3> , int> mappa_vertici;    
					map<pair<array<int,3>, array<int,3>>, int> mappa_lati;
					map<int, pair<Vector3i, Vector3i>> mappa_facce;
					//Siamo nel caso ottaedro:
					//in questo caso il programma deve anche restituirmi un poliedro di Goldberg di classe I:
					for(int i = 0; i < F; i++) {
						int id_A = oCell2DsVertices[i][0];
						int id_B = oCell2DsVertices[i][1]; 
						int id_C = oCell2DsVertices[i][2];
						Vector3d A = oCell0DsCoordinates[id_A];
						Vector3d B = oCell0DsCoordinates[id_B];
						Vector3d C = oCell0DsCoordinates[id_C];	
						vector<Vector3d> points = punti_triangolazione(A,B,C,b);
						if(!file_vertici(points, mappa_vertici, id_vertice, s_g_Cell0Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						};
						
						
						if(!file_lati(points, mappa_lati, mappa_vertici, id_lato, b, s_g_Cell1Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						};	
						
						if(!file_facce(points, mappa_facce, mappa_lati, mappa_vertici, id_faccia, b, s_g_Cell2Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						};					
					}
					
					for (const auto& [key, value] : mappa_facce)
					{
						s_g_Cell2Ds << key << " " << "3 3 " << value.first.transpose()<< " " << value.second.transpose() << endl;
					}
					
					if(!file_poliedro(F_s_g, V_s_g, L_s_g, s_g_Cell3Ds))
					{
						cerr << "errore nella compilazione del file" << endl;
						return 1;
					};
					
					Cell0DsCoordinates = Cell0DsConverter(V_s_g, mappa_vertici);
					Cell1DsExtrema = Cell1DsConverter(L_s_g, mappa_vertici, mappa_lati);
					
				}
				
				else if(q == 5) {
					int F = 20;
					F_s_g = 20 * T;
					V_s_g = 10 * T + 2;
					L_s_g = 30 * T;
					int id_vertice = 0;
					int id_lato = 0;
					int id_faccia = 0;
					map<array<int,3> , int> mappa_vertici;    
					map<pair<array<int,3>, array<int,3>>, int> mappa_lati;
					map<int, pair<Vector3i, Vector3i>> mappa_facce;
					//Siamo nel caso icosaedro:
					//in questo caso il programma deve anche restituirmi un poliedro di Goldberg di classe I:
					for(int i = 0; i < F; i++) {
						int id_A = iCell2DsVertices[i][0];
						int id_B = iCell2DsVertices[i][1]; 
						int id_C = iCell2DsVertices[i][2];
						Vector3d A = iCell0DsCoordinates[id_A];
						Vector3d B = iCell0DsCoordinates[id_B];
						Vector3d C = iCell0DsCoordinates[id_C];	
						vector<Vector3d> points = punti_triangolazione(A,B,C,b);
						if(!file_vertici(points, mappa_vertici, id_vertice, s_g_Cell0Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						};
						
						
						if(!file_lati(points, mappa_lati, mappa_vertici, id_lato, b, s_g_Cell1Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						};	
						
						if(!file_facce(points, mappa_facce, mappa_lati, mappa_vertici, id_faccia, b, s_g_Cell2Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						};					
					}
					
					for (const auto& [key, value] : mappa_facce)
					{
						s_g_Cell2Ds << key << " " << "3 3 " << value.first.transpose()<< " " << value.second.transpose() << endl;
					}
					
					if(!file_poliedro(F_s_g, V_s_g, L_s_g, s_g_Cell3Ds))
					{
						cerr << "errore nella compilazione del file" << endl;
						return 1;
					};
					
					Cell0DsCoordinates = Cell0DsConverter(V_s_g, mappa_vertici);
					Cell1DsExtrema = Cell1DsConverter(L_s_g, mappa_vertici, mappa_lati);
				
				}							
			}
					// Triangolazione classe II
			else if(b == c && b >= 1) {
				int T = b * b + b * c + c * c; 
				ofstream s_g_Cell0Ds("s_g_Cell0Ds.txt");
				ofstream s_g_Cell1Ds("s_g_Cell1Ds.txt");
				ofstream s_g_Cell2Ds("s_g_Cell2Ds.txt");
				ofstream s_g_Cell3Ds("s_g_Cell3Ds.txt");
				s_g_Cell0Ds << "Id " << "x " << "y " << "z" << "\n";
				s_g_Cell1Ds << "Id " << "start_vertex " << "end_vertex" << "\n"; 
				s_g_Cell2Ds << "Id " << "num_vertices " << "num_edges " << "vertices " << "edges " << "\n";
				s_g_Cell3Ds << "Id " << "num_vertices " << "num_edges " << "num_faces " << "vertices " << "edges " << "faces " << "\n";
				if (q == 3) {
					int F = 4;
					F_s_g = 4 * (3 * b * b + 3 * b);
					V_s_g = 4 + 6 * (2 * b - 1) + 4 * (3 * b * b / 2.0 - 3 * b / 2.0 + 1);
					L_s_g = 6 * (2 * b) + 4 * (9 * b * b / 2.0 + 3 * b / 2.0);
					int id_vertice = 0;
					int id_lato = 0;
					int id_faccia = 0;
					map<array<int,3> , int> mappa_vertici;
					map<pair<array<int,3>, array<int,3>>, int> mappa_lati;
					map<array<array<int, 3>, 3>, int> mappa_facce;
					for(int i = 0; i < F; i++) {
						int id_A = tCell2DsVertices[i][0];
						int id_B = tCell2DsVertices[i][1]; 
						int id_C = tCell2DsVertices[i][2];
						Vector3d A = tCell0DsCoordinates[id_A];
						Vector3d B = tCell0DsCoordinates[id_B];
						Vector3d C = tCell0DsCoordinates[id_C];
						vector<Vector3d> points = punti_triangolazione_II(A, B, C, b);
						vector<Vector3d> punti_n_n = punti_triangolazione_II_n_n(A, B, C, b);
						if(!file_vertici_II(points, mappa_vertici, id_vertice, s_g_Cell0Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						};

						if (!file_lati_II(punti_n_n, mappa_lati, mappa_vertici, id_lato, b, s_g_Cell1Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						}; 
						if (!file_facce_II(punti_n_n, mappa_facce,mappa_lati, mappa_vertici, id_faccia, b, s_g_Cell2Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						}; 
						if(!file_poliedro(F_s_g, V_s_g, L_s_g, s_g_Cell3Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						};
						Cell0DsCoordinates = Cell0DsConverter(V_s_g, mappa_vertici);
						Cell1DsExtrema = Cell1DsConverter(L_s_g, mappa_vertici, mappa_lati);	
					} 
				}
				else if (q == 4) {
					int F = 8;
					F_s_g = 8 * (3 * b * b + 3 * b);
					V_s_g = 6 + 12 * (2 * b - 1) + 8 * (3 * b * b / 2.0 - 3 * b / 2.0 + 1);
					L_s_g = 12 * (2 * b) + 8 * (9 * b * b / 2.0 + 3 * b / 2.0);
					int id_vertice = 0;
					int id_lato = 0;
					int id_faccia = 0;
					map<array<int,3> , int> mappa_vertici;
					map<pair<array<int,3>, array<int,3>>, int> mappa_lati;
					map<array<array<int, 3>, 3>, int> mappa_facce;
					for(int i = 0; i < F; i++) {
						int id_A = oCell2DsVertices[i][0];
						int id_B = oCell2DsVertices[i][1]; 
						int id_C = oCell2DsVertices[i][2];
						Vector3d A = oCell0DsCoordinates[id_A];
						Vector3d B = oCell0DsCoordinates[id_B];
						Vector3d C = oCell0DsCoordinates[id_C];
						vector<Vector3d> points = punti_triangolazione_II(A, B, C, b);
						vector<Vector3d> punti_n_n = punti_triangolazione_II_n_n(A, B, C, b);
						if(!file_vertici_II(points, mappa_vertici, id_vertice, s_g_Cell0Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						};

						if (!file_lati_II(punti_n_n, mappa_lati, mappa_vertici, id_lato, b, s_g_Cell1Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						}; 
						if (!file_facce_II(punti_n_n, mappa_facce,mappa_lati, mappa_vertici, id_faccia, b, s_g_Cell2Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						}; 
						if(!file_poliedro(F_s_g, V_s_g, L_s_g, s_g_Cell3Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						};

						Cell0DsCoordinates = Cell0DsConverter(V_s_g, mappa_vertici);
						Cell1DsExtrema = Cell1DsConverter(L_s_g, mappa_vertici, mappa_lati);	
					} 
				}
				else if (q == 5) {
					int F = 20;
					F_s_g = 20 * (3 * b * b + 3 * b);
					V_s_g = 12 + 30 * (2 * b - 1) + 20 * (3 * b * b / 2.0 - 3 * b / 2.0 + 1);
					L_s_g = 30 * (2 * b) + 20 * (9 * b * b / 2.0 + 3 * b / 2.0);
					int id_vertice = 0;
					int id_lato = 0;
					int id_faccia = 0;
					map<array<int,3> , int> mappa_vertici;
					map<pair<array<int,3>, array<int,3>>, int> mappa_lati;
					map<array<array<int, 3>, 3>, int> mappa_facce;
					for(int i = 0; i < F; i++) {
						int id_A = iCell2DsVertices[i][0];
						int id_B = iCell2DsVertices[i][1]; 
						int id_C = iCell2DsVertices[i][2];
						Vector3d A = iCell0DsCoordinates[id_A];
						Vector3d B = iCell0DsCoordinates[id_B];
						Vector3d C = iCell0DsCoordinates[id_C];
						vector<Vector3d> points = punti_triangolazione_II(A, B, C, b);
						vector<Vector3d> punti_n_n = punti_triangolazione_II_n_n(A, B, C, b);
						if(!file_vertici_II(points, mappa_vertici, id_vertice, s_g_Cell0Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						};

						if (!file_lati_II(punti_n_n, mappa_lati, mappa_vertici, id_lato, b, s_g_Cell1Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						}; 
						if (!file_facce_II(punti_n_n, mappa_facce,mappa_lati, mappa_vertici, id_faccia, b, s_g_Cell2Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						}; 
						if(!file_poliedro(F_s_g, V_s_g, L_s_g, s_g_Cell3Ds))
						{
							cerr << "errore nella compilazione del file" << endl;
							return 1;
						};

						Cell0DsCoordinates = Cell0DsConverter(V_s_g, mappa_vertici);
						Cell1DsExtrema = Cell1DsConverter(L_s_g, mappa_vertici, mappa_lati);	
					} 
				}
			}				
		}
		else{
			cerr << "Input non valido" << endl;
		}
	    vector<double> ShortPath_v(V_s_g, 0.0);
		vector<double> ShortPath_l(L_s_g, 0.0);
		if (argc == 7){
			int start = stoi(argv[5]);
			int end = stoi(argv[6]);
			// Controllo validità input
			if (start < 0 || start >= V_s_g){
				cerr << "ERRORE: id del vertice di partenza fuori dal range: [" << 0 << ", " << V_s_g -1 << "]" << endl;
				return 1;
			}
			if (end < 0 || end >= V_s_g){
				cerr << "ERRORE: id del vertice di arrivo fuori dal range: [" << 0 << ", " << V_s_g -1 << "]" << endl;
				return 1;
			}
			vector<vector<int>> lista_ad(V_s_g);
			vector<vector<double>> pesi(V_s_g);

			// Costruzione lista di adiacenza
			for (int i = 0; i < L_s_g; ++i) {  // Entrambi i versi perchè non orientato
				int from = Cell1DsExtrema(0, i);
				int to   = Cell1DsExtrema(1, i);
				Vector3d p_from = Cell0DsCoordinates.col(from);
				Vector3d p_to = Cell0DsCoordinates.col(to);
				
				double peso = (p_from - p_to).norm();
				
				lista_ad[from].push_back(to);
				lista_ad[to].push_back(from);
				pesi[from].push_back(peso);
				pesi[to].push_back(peso);
			}
			
			vector<int> cammino_min = dijkstra(V_s_g, lista_ad, pesi, start, end) ;
			cout << "Il cammino minimo è: ";
			for (size_t k = 0; k < cammino_min.size(); k++){
				cout << cammino_min[k] << " ";
				ShortPath_v[cammino_min[k]] = 1.0;
			}
			cout << endl;

			for (size_t i = 0; i < cammino_min.size() -1 ; ++i) {
				int u = cammino_min[i];
				int v = cammino_min[i + 1];

				for (int j = 0; j < Cell1DsExtrema.cols(); ++j) {
					int from = Cell1DsExtrema(0, j);
					int to   = Cell1DsExtrema(1, j);

					if ((from == u && to == v) || (from == v && to == u)) {
						ShortPath_l[j] = 1.0;
						break;
					}
				}
			}
					
		}
		
		Gedim::UCDUtilities utilities;
		{
			
			vector<Gedim::UCDProperty<double>> cell0Ds_properties(1);

			cell0Ds_properties[0].Label = "Marker";
			cell0Ds_properties[0].UnitLabel = "-";
			cell0Ds_properties[0].NumComponents = 1;
			
			cell0Ds_properties[0].Data = ShortPath_v.data();

				
			utilities.ExportPoints("./Cell0Ds.inp",
								   Cell0DsCoordinates,
								   cell0Ds_properties);

		}

		{

			vector<Gedim::UCDProperty<double>> cell1Ds_properties(1);
			
			cell1Ds_properties[0].Label = "Marker";
			cell1Ds_properties[0].UnitLabel = "-";
			cell1Ds_properties[0].NumComponents = 1;

			cell1Ds_properties[0].Data = ShortPath_l.data();
				
			utilities.ExportSegments("./Cell1Ds.inp",
									 Cell0DsCoordinates,
									 Cell1DsExtrema,
									 {},
									 cell1Ds_properties);
		}
	}
	else {
		cerr << "Input non valido" << endl;
	}
}
