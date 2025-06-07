#pragma once
#include <gtest/gtest.h>

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
#include "Utils.hpp"

using namespace std;

TEST(Triangolazione_Test, NumeroPunti) {
	Vector3d A( 0.57735027,  0.57735027,  0.57735027);
	Vector3d B(-0.57735027, -0.57735027,  0.57735027);
	Vector3d C(-0.57735027,  0.57735027, -0.57735027);
	int b = 2;
	EXPECT_EQ(punti_triangolazione(A,B,C,b).size(), 6);
	EXPECT_EQ(punti_triangolazione_II_n_n(A,B,C,b).size(), 16);
}  

TEST(Triangolazione_Test, Normalizzazione) {
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

TEST(Calcola_intersezione_test, SiIntersezione) {
	Vector3d A(0, 0, 0);
	Vector3d B(2, 2, 0);
	Vector3d C(0, 2, 0);
	Vector3d D(2, 0, 0);
	Vector3d expected(1, 1, 0);
	auto intersect = calcola_intersezione(A, B, C, D);
	ASSERT_TRUE(intersect.has_value());
	EXPECT_TRUE(intersect->isApprox(expected, 1e-6));
}

TEST(Calcola_intersezione_test, NoIntersezione) {
	Vector3d A(0, 0, 0);
	Vector3d B(1, 0, 0);
	Vector3d C(0, 1, 0);
	Vector3d D(1, 2, 0);
	auto intersect = calcola_intersezione(A, B, C, D);
	EXPECT_FALSE(intersect.has_value());
}

TEST(Punti_lungo_i_lati_e_lungo_i_lati_n_n_test, NumeroPunti) {
	int b = 2;
	Vector3d A(0.57735027,  0.57735027,  0.57735027);
	Vector3d B(-0.57735027, -0.57735027,  0.57735027);    
	Vector3d C(-0.57735027,  0.57735027, -0.57735027);
	
	EXPECT_EQ(punti_lungo_i_lati(b,A,B,C).size(), 12);
	EXPECT_EQ(punti_lungo_i_lati_normalizzati(b,A,B,C).size(), 12);
}

TEST(Punti_lungo_i_lati_n_n_test, Normalizzazione) {
	int b = 2;
	Vector3d A(0.57735027,  0.57735027,  0.57735027);
	Vector3d B(-0.57735027, -0.57735027,  0.57735027);    
	Vector3d C(-0.57735027,  0.57735027, -0.57735027);
	auto points_II = punti_lungo_i_lati_normalizzati(b,A,B,C);
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
	EXPECT_EQ(trova_punti_vicini(punto, points).size(), 6);
	
}

TEST(Id_vertici_test, Valore_id) {
	Vector3d A( 0.57735027,  0.57735027,  0.57735027);
	Vector3d B(-0.57735027, -0.57735027,  0.57735027);
	Vector3d C(-0.57735027,  0.57735027, -0.57735027);
	int b = 2;
	int v = 6;
	vector<Vector3d> points = {A, B, C};
	map<array<int,3> , int> mappa_vertici;
	int id_vertice = 0;
	ostringstream s_g_Cell0Ds;
	
	bool result = file_vertici(points, mappa_vertici, id_vertice, s_g_Cell0Ds);
	
	EXPECT_TRUE(result);
	
	EXPECT_EQ(mappa_vertici.size(), 3);
    EXPECT_EQ(id_vertice, 3);
	
	string output = s_g_Cell0Ds.str();
	
	EXPECT_NE(output.find("0 0.577 0.577 0.577"), string::npos);
	
}