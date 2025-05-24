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
				int F_s_g = 4 * T;
				int V_s_g = 2 * T + 2;
				int L_s_g = 6 * T;
				int V_p_o_f = (b + 1) * (b + 2) / 2;
				int L_p_o_f = 1.5 * b * (b + 1);
				int id_vertice = 0;
				int id_lato = 0;
				int id_faccia = 0;
				map<array<double,3> , int> mappa_vertici;    
				map<pair<array<double, 3> , array<double, 3>> , int> mappa_lati;
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
							points.push_back(P / sqrt(P[0] * P[0] + P[1] * P[1] + P[2] * P[2]));
						}
					}
					s_g_Cell2Ds << Id_faccia << V_p_o_f << L_p_o_f;
					for(int i = 0; i < points.size(); i++) {
						if (mappa_vertici.find(points[i]) != mappa_vertici.end()) {
							
						} else {
							mappa_vertici.insert({points[i], id_vertice});
							s_g_Cell0Ds << id_vertice << points[i][0] << points[i][1] << points[i][2] << "\n";
							id_vertice++;
						}
					}
					
					int d = 0;
					for (int i = 0; i < b; i++) {   
						for (int j = d; j < d + b - i; j++) {
							if (mappa_lati.find({points[j],points[j+b+1-i]} != mappa_lati.end()){}
							elseif(mappa_lati.find({points[j+b+1-i],points[j]} != mappa_lati.end()){}
							else{ mappa_lati.insert({points[j],points[j+b+1-i]},id_lato);
							s_g_Cell1Ds << Id_lato << mappa_vertici[points[j]] << mappa_vertici[points[j+b+1-i]] << "\n";
							id_lato++;}
							if (mappa_lati.find({points[j],points[j+1]} != mappa_lati.end()){}
							elseif(mappa_lati.find({points[j+1],points[j]} != mappa_lati.end()){}
							else{ mappa_lati.insert({points[j],points[j+1]},id_lato);
							s_g_Cell1Ds << Id_lato <<  mappa_vertici[points[j]] << mappa_vertici[points[j+1]] << "\n";
							id_lato++;}
							if (mappa_lati.find({points[j+1],points[j+b+1-i]} != mappa_lati.end()){}
							elseif(mappa_lati.find({points[j+1],points[j+b+1-i]} != mappa_lati.end()){}
							else{ mappa_lati.insert({points[j+1],points[j+b+1-i]},id_lato);
							s_g_Cell1Ds << id_lato <<  mappa_vertici[points[j+1]] << mappa_vertici[points[j+b+1-i]]<< "\n";
							id_lato++;}	
							s_g_Cell2Ds << mappa_vertici[points[j]] << mappa_vertici[points[j + 1]] << mappa_vertici[points[j+b+1-i]] << mappa_lati[
						}
						d = d + b + 1 - i;	
					}	
					
                    int d = 0;
					
					////
					
					for (int i = 0; i < b; i++) {   
						for (int j = d; j < d + b - i; j++){
							if( i==0){
								Vector3d id_vertici_faccia = [ mappa_vertici[points[j]], mappa_vertici[points[j+1]], mappa_vertici[points[j+b+1-i]]]; // insieme id vertici per ogni faccia  // controllare struttura
								Vector3d id_lati_faccia = [ mappa_lati[{points[j],points[j+1]}], mappa_lati[{points[j+1],points[j+b+1-i]}], mappa_lati[{points[j],points[j+b+1-i]}]];
							    mappa_facce.insert({id_faccia, {id_vertici_faccia, id_lati_faccia}};
								id_faccia++;
							}
							else {
								Vector3d id_vertici_faccia = [ mappa_vertici[points[j]], mappa_vertici[points[j+1]], mappa_vertici[points[j+b+1-i]]]; // insieme id vertici per ogni faccia  // controllare struttura
								Vector3d id_lati_faccia = [ mappa_lati[{points[j],points[j+1]}], mappa_lati[{points[j+1],points[j+b+1-i]}], mappa_lati[{points[j],points[j+b+1-i]}]];
								mappa_facce.insert({id_faccia, {id_vertici_faccia, id_lati_faccia}};
								id_faccia++;
								Vector3d id_vertici_faccia = [ mappa_vertici[points[j]], mappa_vertici[points[j+1]], mappa_vertici[points[j-b-1+i]]]; // insieme id vertici per ogni faccia  // controllare struttura
								Vector3d id_lati_faccia = [ mappa_lati[{points[j],points[j+1]}], mappa_lati[{points[j+1],points[j-b-1+i]}], mappa_lati[{points[j],points[j-b-1+i]}]];
								mappa_facce.insert({id_faccia, {id_vertici_faccia, id_lati_faccia}};
								id_faccia++;
							}
						}
					}
							
                }

			}
	if (q == 4) {
		Siamo nel caso ottaedro:
				
	}	 
	if (q == 5) {
		Siamo nel caso icosaedro:
				
	} 
}
*/

    return 0;
}