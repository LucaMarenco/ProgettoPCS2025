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

const double EPS = 1e-12;

// Crea una chiave "intera" da coordinate (x, y, z) con arrotondamento
tuple<int, int, int> key_from_point(double x, double y, double z) {
    return make_tuple(
        static_cast<int>(std::round(x / EPS)),
        static_cast<int>(std::round(y / EPS)),
        static_cast<int>(std::round(z / EPS))
    );
}

class GestoreVertici {
private:
    map<tuple<int, int, int>, int> mappa_vertici;  // chiave → ID
    vector<array<double, 3>> lista_vertici;         // ID → punto
    int next_id = 0;

public:
    // Ritorna l'ID del punto (x, y, z), lo aggiunge se è nuovo
    int ottieni_o_aggiungi(double x, double y, double z) {
        auto chiave = key_from_point(x, y, z);

        auto it = mappa_vertici.find(chiave);
        if (it != mappa_vertici.end()) {
            return it->second; // Punto già esistente
        }

        int id = next_id++;
        mappa_vertici[chiave] = id;
        lista_vertici.push_back({x, y, z});
        return id;
    }

    // Accesso alla lista dei vertici ordinata per ID
    const vector<array<double, 3>>& get_vertici() const {
        return lista_vertici;
    }

    int numero_vertici() const {
        return lista_vertici.size();
    }
};


void crea_poliedro_geodetico(int p, q, b, c) {
	if (p == 3) {
		if (b >= 1 && c == 0) {
			int T = b * b + b * c + c * c; 
			ofstream s_g_Cell0Ds("s_g_Cell0Ds.txt");
			ofstream s_g_Cell1Ds("s_g_Cell1Ds.txt");
			ofstream s_g_Cell2Ds("s_g_Cell2Ds.txt");
			ofstream s_g_Cell3Ds("s_g_Cell3Ds.txt");
			s_g_Cell0Ds << "Id" << "x" << "y" << "z" << "\n";
			s_g_Cell1Ds << "Id" << "start_vertex" << "end_vertex" << "\n"; 
			s_g_Cell2Ds << "Id" << "num_vertices" << "num_edges" << "vertices" << "edges" << "\n";
			s_g_Cell3Ds << "Id" << "num_vertices" << "num_edges" << "num_faces" << "vertices" << "edges" << "faces" << "\n";
			if (q == 3) {
				int F = 4;
				
				// Possiamo creare una funzione da qua in giù fino a q == 4
				int id_vertice = 0;
				int id_lato = 0;
				map<array<int,3> , int> mappa_vertici;    
				map<int, Vector3d> local_id_to_position;                    /// A ID associo coordinate o a coordinate associo ID ?????
				
				map<pair<Vector3d,Vector3d>, int> lati_map;
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
							points.push_back(P / sqrt(P[0] * P[0] + P[1] * P[1] + P[2] * P[2]));  //già normalizzate
						}
					}
					for(int l_id = 0; l_id < points.size(); l_id++) {
						local_id_to_position.insert({l_id, points[l_id]});
						auto [iterator, inserted] = mappa_vertici.insert({points[id_vertice], id_vertice}); 						/// stesso problema di sopra
						if(inserted)
							id_vertice++;
					}
					
					int d = 0;
					for (int i = 0; i < b; i++) {   
						for (int j = d; j < d + b - i; j++) {
							lati_map.insert({pair<local_id_to_position[j], local_id_to_position[j+1]>, id_lato});
							id_lato++;
							lati_map.insert({pair<local_id_to_position[j], local_id_to_position[j+b+1-i]>, id_lato});
							id_lato++;
							lati_map.insert({pair<local_id_to_position[j+1], local_id_to_position[j+b+1-i]>, id_lato});
							id_lato++;							
						}
						d = d + b + 1 - i;	
					}
					local_id_to_position.clear();                   

                }

			}
			if (q == 4) {
				int F = 8;
				int id_vertice = 0;
				int id_lato = 0;
				//Siamo nel caso ottaedro:
				//in questo caso il programma deve anche restituirmi un poliedro di Goldberg di classe I:
				
				//uso la stessa funzione del caso q==3
				
			}
			
			if (q == 5) {
				int F = 20;
				int id_vertice = 0;
				int id_lato = 0;
				//Siamo nel caso icosaedro:
				//in questo caso il programma deve anche restituirmi un poliedro di Goldberg di classe I:
				
				//uso la stessa funzione del caso q==3
			} 
		}
		
		if(b==c && b >=1){
			
		}
	}


    return 0;
}
*/
