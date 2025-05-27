#include <iostream>
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
    }
	
/*Consideriamo di partire che abbiamo già dei file con le celle che rappresentano i solidi platonici con p = 3, 
quindi tetraedro, ottaedro e icosaedro; e consideriamo anche di aver già aperto i file in lettura 
nel programma (presi in input).
Dati in input p,q,b,c bisogna innanzitutto, perché sia un solido geodetico, avere p = 3.
Faccio una funzione void (perché non ritorna niente) che, dati in input i 4 interi p,q,b e c,
crea 4 file txt che rappresentano il poliedro geodetico che si ottiene a partire da un certo solido platonico:
*/









// Correzioni Luca


// Da mettere prima del main
std::array<double, 3> to_array(const Eigen::Vector3d& v) {
    return {v[0], v[1], v[2]};
}

int main()
{
	int p=3;
	int b=2;
	int c=0;
	int q=3;
	
	std::vector<std::vector<int>> cells = {
        {0, 3, 3, 0, 1, 2, 0, 3, 1},
        {1, 3, 3, 0, 1, 3, 0, 4, 2},
        {2, 3, 3, 0, 3, 2, 2, 5, 1},
        {3, 3, 3, 1, 2, 3, 3, 5, 4}
    };

    // Mappa che associa l'ID della faccia al vettore dei suoi vertici
    std::map<int, std::vector<int>> tCell2DsVertices;

    for (const auto& cell : cells) {
        int face_id = cell[0]; // ID della faccia
        int num_vertices = cell[1]; // Numero di vertici

        // Estrai i vertici dalla riga
        std::vector<int> vertices(cell.begin() + 3, cell.begin() + 3 + num_vertices);

        // Aggiungi alla mappa
        tCell2DsVertices[face_id] = vertices;
    }
		
	std::vector<std::pair<int, Vector3d>> vertices = {
        {0, { 0.57735027,  0.57735027,  0.57735027}},
        {1, {-0.57735027, -0.57735027,  0.57735027}},
        {2, {-0.57735027,  0.57735027, -0.57735027}},
        {3, { 0.57735027, -0.57735027, -0.57735027}}
    };

    std::map<int, Vector3d> tCell0DsCoordinates;

    for (const auto& [id, coords] : vertices) {
        tCell0DsCoordinates[id] = coords;
    }


	if (p == 3) {
		if (b >= 1 && c == 0 ){
			int T = b * b + b * c + c * c; 
			ofstream s_g_Cell0Ds("s_g_Cell0Ds.txt");
			ofstream s_g_Cell1Ds("s_g_Cell1Ds.txt");
			ofstream s_g_Cell2Ds("s_g_Cell2Ds.txt");
			ofstream //s_g_Cell3Ds("//s_g_Cell3Ds.txt");
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
				//Siamo nel caso tetraedro:
				//in questo caso il programma deve anche restituirmi un poliedro di Goldberg di classe I:
				for(int i = 0; i < F; i++) {
					int id_A = tCell2DsVertices[i][0];
					int id_B = tCell2DsVertices[i][1]; 
					int id_C = tCell2DsVertices[i][2];
					Vector3d A = tCell0DsCoordinates[id_A];
					Vector3d B = tCell0DsCoordinates[id_B];
					Vector3d C = tCell0DsCoordinates[id_C];	
					vector<Vector3d> points;
					for (int k = 0; k <= b; k++) {
						for (int j = 0; j <= b - k; j++) {
							double c_A = 1.0 - (double)k / b - (double)j / b;
							double c_B = (double)k / b;
							double c_C = (double)j / b;
							Vector3d P = c_A * A + c_B * B + c_C * C;
							points.push_back(P / P.norm());
						}
					}
					
					for(int z = 0; z < points.size(); z++) {
						 auto key = to_array(points[z]);    
						if (mappa_vertici.find(key) == mappa_vertici.end()) {
							mappa_vertici[key] = id_vertice;
							
							double eps = 1e-5;
							s_g_Cell0Ds << id_vertice << " " << key[0]*eps << " " << key[1]*eps << " " << key[2]*eps << "\n";
							id_vertice++;	
							
						}
					}
					
					int d = 0;
					for (int f = 0; f <= b; f++) {   
						for (int j = d; j < d + b - f; j++) {
							auto key_j = to_array(points[j]);
							auto key_j1 = to_array(points[j+1]);
							auto key_jb = to_array(points[j+b+1-f]);
							
							if (mappa_lati.find({min({key_j,key_jb}), max({key_j,key_jb})}) != mappa_lati.end()){}
							else{ mappa_lati.insert({{min({key_j,key_jb}), max({key_j,key_jb})},id_lato});
							s_g_Cell1Ds << id_lato << " " <<mappa_vertici[min({key_j,key_jb})] << " "<<  mappa_vertici[max({key_j,key_jb})]  <<"\n";
							id_lato++;
							}
							if (mappa_lati.find({min({key_j,key_j1}), max({key_j,key_j1})}) != mappa_lati.end()){}
							else{ mappa_lati.insert({{min({key_j,key_j1}), max({key_j,key_j1})} ,id_lato});
							s_g_Cell1Ds << id_lato << " " << mappa_vertici[min({key_j,key_jb})]  << " " << mappa_vertici[max({key_j,key_jb})] << "\n";
							
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
					
					map<int, pair<Vector3i, Vector3i>> mappa_facce;

					Vector3i id_vertici_faccia;
					Vector3i id_lati_faccia;
					for (int i = 0; i < b; i++) {   
						for (int j = d; j < d + b - i; j++){
							auto key_j = to_array(points[j]);
							auto key_j1 = to_array(points[j+1]);
							auto key_jb = to_array(points[j+b+1-i]);
											
							if(i==0){
								id_vertici_faccia = Vector3i{mappa_vertici[key_j], mappa_vertici[key_j1], mappa_vertici[key_jb]}; // insieme id vertici per ogni faccia
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
					
					for (const auto& [key, value] : mappa_facce)
					{
						s_g_Cell2Ds << key << " " << "3 3 " << value.first.transpose()<< " " << value.second.transpose() << endl;
					}
					
				}
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
				
				
				MatrixXd Cell0DsCoordinates(3, V_s_g);
				for (const auto& [coord,id] : mappa_vertici){
					for (int d = 0; d < 3; ++d) {
						Cell0DsCoordinates(d,id) = coord[d];
					}
			    }
				for (int d = 0; d < 3; ++d) {
					std::cout << "Riga " << d << ": ";
					for (int id = 0; id < V_s_g; ++id) {
						std::cout << Cell0DsCoordinates(d,id) << " ";
					}
					std::cout << "\n";
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