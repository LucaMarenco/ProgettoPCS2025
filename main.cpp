/*#include <iostream>
#include "TriangularMesh.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace TriangularLibrary;

int main()
{
    TriangularMesh mesh;

    if(!ImportMesh(mesh))
    {
        cerr << "file not found" << endl;
        return 1;
    }

    /// Per visualizzare online le mesh:
    /// 1. Convertire i file .inp in file .vtu con https://meshconverter.it/it
    /// 2. Caricare il file .vtu su https://kitware.github.io/glance/app/

    Gedim::UCDUtilities utilities;
    {
        vector<Gedim::UCDProperty<double>> cell0Ds_properties(1);


        cell0Ds_properties[0].NumComponents = 1;

        vector<double> cell0Ds_marker(mesh.NumCell0Ds, 0.0);


        utilities.ExportPoints("./Cell0Ds.inp",
                               mesh.Cell0DsCoordinates);

    }

    {

        vector<Gedim::UCDProperty<double>> cell1Ds_properties(1);

        cell1Ds_properties[0].NumComponents = 1;

        vector<double> cell1Ds_marker(mesh.NumCell1Ds, 0.0);

        utilities.ExportSegments("./Cell1Ds.inp",
                                 mesh.Cell0DsCoordinates,
                                 mesh.Cell1DsExtrema,
                                 {});
        utilities.ExportSegments("./Cell1Ds.inp",
                                 mesh.Cell0DsCoordinates,
                                 mesh.Cell1DsExtrema,
                                 {});
    }
    
	{

        vector<Gedim::UCDProperty<double>> cell1Ds_properties(1);

        cell1Ds_properties[0].NumComponents = 1;

        vector<double> cell1Ds_marker(mesh.NumCell1Ds, 0.0);

        utilities.ExportSegments("./Cell1Ds.inp",
                                 mesh.Cell0DsCoordinates,
                                 mesh.Cell1DsExtrema,
                                 {});
        utilities.ExportSegments("./Cell1Ds.inp",
                                 mesh.Cell0DsCoordinates,
                                 mesh.Cell1DsExtrema,
                                 {});
    }*/
	
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
#include "Utils.hpp"

using namespace std;
using namespace Eigen;

int main()
{
	int p=3;
	int b=2;
	int c=0;   
	int q=5;
	
	/*vector<vector<int>> cells = {
        {0, 3, 3, 0, 1, 2, 0, 3, 1},
        {1, 3, 3, 0, 1, 3, 0, 4, 2},   // Questo è il Cell2 del tetraedro di partenza
        {2, 3, 3, 0, 3, 2, 2, 5, 1},
        {3, 3, 3, 1, 2, 3, 3, 5, 4}
    };*/
	
	/*vector<vector<int>> cells = {
        {0, 3, 3, 0, 1, 2, 0, 4, 1},
        {1, 3, 3, 0, 2, 3, 1, 5, 2},   // Questo è il Cell2 dell'ottaedro di partenza
        {2, 3, 3, 0, 3, 4, 2, 6, 3},
        {3, 3, 3, 0, 4, 1, 3, 7, 0},
		{4, 3, 3, 5, 2, 1, 9, 4, 8},
		{5, 3, 3, 5, 3, 2, 10, 5, 9},
		{6, 3, 3, 5, 4, 3, 11, 6, 10},
		{7, 3, 3, 5, 1, 4, 8, 7, 11}
    };*/
	
	vector<vector<int>> cells = {
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

    for (const auto& cell : cells) {
        int face_id = cell[0]; // ID della faccia
        int num_vertices = cell[1]; // Numero di vertici

        // Estrai i vertici dalla riga
        vector<int> vertices(cell.begin() + 3, cell.begin() + 3 + num_vertices);

        // Aggiungi alla mappa
        tCell2DsVertices[face_id] = vertices;
    }
		
	/*vector<pair<int, Vector3d>> vertices = {
        {0, { 0.57735027,  0.57735027,  0.57735027}},
        {1, {-0.57735027, -0.57735027,  0.57735027}},    // Questo è il Cell0 del tetraedro di partenza
        {2, {-0.57735027,  0.57735027, -0.57735027}},
        {3, { 0.57735027, -0.57735027, -0.57735027}}
    };*/
	
	/*vector<pair<int, Vector3d>> vertices = {
        {0, {0,  0, 1}},
        {1, {0, -1, 0}},    // Questo è il Cell0 dell'ottaedro di partenza
        {2, { 1, 0, 0}},
        {3, { 0, 1, 0}},
		{4, {-1, 0, 0}},
		{5, {0, 0, -1}}
    };*/
	
	vector<pair<int, Vector3d>> vertices = {
        {0, {0.52573111, 0.85065081, 0.00000000}},
		{1, {0.00000000, 0.52573111, 0.85065081}},
        {2, {-0.52573111, 0.85065081, 0.00000000}},
		{3, {0.85065081, 0.00000000, 0.52573111}},		// Questo è il Cell0 dell'icosaedro di partenza
        {4, {0.85065081, 0.00000000, -0.52573111}},
		{5, {0.00000000, 0.52573111, -0.85065081}},
		{6, {-0.85065081, 0.00000000, -0.52573111}},
		{7, {-0.85065081, 0.00000000, 0.52573111}},
		{8, {0.52573111, -0.85065081, 0.00000000}},
        {9, {-0.52573111, -0.85065081, 0.00000000}},
		{10, {0.00000000, -0.52573111, 0.85065081}},
		{11, {0.00000000, -0.52573111, -0.85065081}}
    };

    map<int, Vector3d> tCell0DsCoordinates; // Mappa che associa l'ID di un vertice alle sue coordinate

    for (const auto& [id, coords] : vertices) {
        tCell0DsCoordinates[id] = coords;
    }


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
				int F_s_g = 4 * T;
				int V_s_g = 2 * T + 2;
				int L_s_g = 6 * T;
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
				
				
				
				MatrixXd Cell0DsCoordinates(3, V_s_g);
				for (const auto& [coord,id] : mappa_vertici){
					for (int d = 0; d < 3; ++d) {
						Cell0DsCoordinates(d,id) = coord[d];
					}
			    }

			}
			
			if(q==4){
				int F = 8;
				int F_s_g = 8 * T;
				int V_s_g = 4 * T + 2;
				int L_s_g = 12 * T;
				int id_vertice = 0;
				int id_lato = 0;
				int id_faccia = 0;
				map<array<int,3> , int> mappa_vertici;    
				map<pair<array<int,3>, array<int,3>>, int> mappa_lati;
				map<int, pair<Vector3i, Vector3i>> mappa_facce;
				//Siamo nel caso ottaedro:
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
				
				
				
				MatrixXd Cell0DsCoordinates(3, V_s_g);
				for (const auto& [coord,id] : mappa_vertici){
					for (int d = 0; d < 3; ++d) {
						Cell0DsCoordinates(d,id) = coord[d];
					}
			    }
				
			}
			
			if(q==5){
				int F = 20;
				int F_s_g = 20 * T;
				int V_s_g = 10 * T + 2;
				int L_s_g = 30 * T;
				int id_vertice = 0;
				int id_lato = 0;
				int id_faccia = 0;
				map<array<int,3> , int> mappa_vertici;    
				map<pair<array<int,3>, array<int,3>>, int> mappa_lati;
				map<int, pair<Vector3i, Vector3i>> mappa_facce;
				//Siamo nel caso icosaedro:
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
				
				
				
				MatrixXd Cell0DsCoordinates(3, V_s_g);
				for (const auto& [coord,id] : mappa_vertici){
					for (int d = 0; d < 3; ++d) {
						Cell0DsCoordinates(d,id) = coord[d];
					}
			    }
				
			}
		}
	}
}




/*funzione che ci permette di trasformare le mappe in vettori/matrici, come richiesto per ucd

void ConvertMapToExportData(const std::map<std::array<int, 3>, int>& input_map,
                            Eigen::MatrixXd& points){

		std::map<std::array<int, 3>, int> mappa_vertici;

		const int num_points = mappa_vertici.size();
		Eigen::MatrixXd points(3, num_points);


		int i = 0;
		for (const auto& [coord, id] : mappa_vertici) {
			points(i, 0) = static_cast<double>(coord[0]);
			points(i, 1) = static_cast<double>(coord[1]);
			points(i, 2) = static_cast<double>(coord[2]);
			

			// Puoi aggiungere qui eventuali proprietà se ne hai
			++i;
		}
	}

// Chiamata finale
ucdUtilities.ExportPoints("output.ucd", points, properties, materials);

*/

