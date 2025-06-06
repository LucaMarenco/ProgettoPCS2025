#include "test_poliedro_geodetico.hpp"
#include <gtest/gtest.h>

TEST(Triangolazione_II_Test, NumeroPunti) {
	Vector3d A( 0.57735027,  0.57735027,  0.57735027);
	Vector3d B(-0.57735027, -0.57735027,  0.57735027);
	Vector3d C(-0.57735027,  0.57735027, -0.57735027);
	int b = 2;
	EXPECT_EQ(punti_triangolazione_II_n_n(A,B,C,b).size(), 16);
}  

TEST(Triangolazione_I_II_Test, Normalizzazione) {
	Vector3d A( 0.57735027,  0.57735027,  0.57735027);
	Vector3d B(-0.57735027, -0.57735027,  0.57735027);
	Vector3d C(-0.57735027,  0.57735027, -0.57735027);
	int b = 2;
	auto points = punti_triangolazione(A,B,C,b);
	auto points_II = punti_triangolazione_II(A,B,C,b);
	for(int i = 0; i < points.size(); i++) {
	EXPECT_NEAR(points[i].norm(), 1.0, 1e-6);
	}
	for(int i = 0; i < points_II.size(); i++) {
	EXPECT_NEAR(points_II[i].norm(), 1.0, 1e-6);
	}
}  

TEST(Id_vertici_test, Valore_id) {
	Vector3d A( 0.57735027,  0.57735027,  0.57735027);
	Vector3d B(-0.57735027, -0.57735027,  0.57735027);
	Vector3d C(-0.57735027,  0.57735027, -0.57735027);
	int b = 2;
	int v = 6;
	vector<Vector3d> points = punti_triangolazione(A,B,C,b);
	map<array<int,3> , int> mappa_vertici;
	int id_vertice = 0;
	ofstream s_g_Cell0Ds("s_g_Cell0Ds.txt");
	if
	
}

TEST(Calcola_intersezione_test, Intersezione) {
	Vector3d A(0, 0, 0);
	Vector3d B(2, 2, 0);
	Vector3d C(0, 2, 0);
	Vector3d D(2, 0, 0);
	Vector3d expected(1, 1, 0);
	EXPECT_NEAR((Calcola_intersezione(A, B, C, D) - expected).norm(), 0.0, 1e-6);
}

TEST(Punti_lungo_i_lati_e_lungo_i_lati_n_n_test, NumeroPunti_e_Normalizzazione) {
	int b = 2;
	Vector3d A(0.57735027,  0.57735027,  0.57735027);
	Vector3d B(-0.57735027, -0.57735027,  0.57735027);    
	Vector3d C(-0.57735027,  0.57735027, -0.57735027);
	auto points = punti_lungo_i_lati(A,B,C,b);
	auto points_II = punti_lungo_i_lati_n_n(A,B,C,b);
	EXPECT_EQ(points.size(), 12);
	EXPECT_NEAR(points_II.size(), 12);
	for(int i = 0; i < points_II.size(); i++) {
		EXPECT_NEAR(points_II[i].norm(), 1.0, 1e-6);
	}
}

TEST(Trova_punti_vicini, Numero_punti_e_vicinanza) {
	Vector3d punto(-0.85065081, 0.00000000, 0.52573111);
	vector<Vector3d> points = {{ 0.57735027,  0.57735027,  0.57735027}, {-0.57735027, -0.57735027,  0.57735027},
	{-0.57735027,  0.57735027, -0.57735027}, { 0.57735027, -0.57735027, -0.57735027}, {0.52573111, 0.85065081, 0.00000000},
	{-0.52573111, 0.85065081, 0.00000000}, {0.52573111, -0.85065081, 0.00000000}, {-0.52573111, -0.85065081, 0.00000000},
	{0.85065081, 0.00000000, 0.52573111}, {-0.85065081, 0.00000000, 0.52573111}};
	EXPECT_EQ(trova_punti_vicini(points).size(), 6);
	
}