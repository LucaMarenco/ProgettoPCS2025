#include "test_poliedro_geodetico.hpp"
#include <gtest/gtest.h>

TEST(Triangolazione_I_Test, NumeroPunti) {
	Vector3d A( 0.57735027,  0.57735027,  0.57735027);
	Vector3d B(-0.57735027, -0.57735027,  0.57735027);
	Vector3d C(-0.57735027,  0.57735027, -0.57735027);
	int b = 2;
	EXPECT_EQ(punti_triangolazione(A,B,C,b).size(), 6);
	EXPECT_EQ(punti_triangolazione_II(A,B,C,b).size(), 16);
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