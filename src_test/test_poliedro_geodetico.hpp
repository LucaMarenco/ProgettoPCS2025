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
	EXPECT_EQ(punti_triangolazione_II(A,B,C,b).size(), 16);
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

TEST(Punti_lungo_i_lati_test, NumeroPunti) {
	int b = 2;
	Vector3d A(0.57735027,  0.57735027,  0.57735027);
	Vector3d B(-0.57735027, -0.57735027,  0.57735027);    
	Vector3d C(-0.57735027,  0.57735027, -0.57735027);
	
	EXPECT_EQ(punti_lungo_i_lati(b,A,B,C).size(), 12);
}

TEST(TrovapuntiviciniTest, Numero_punti_e_vicinanza) {
	Vector3d punto(-0.85065081, 0.00000000, 0.52573111);
	vector<Vector3d> points = {{ 0.57735027,  0.57735027,  0.57735027}, {-0.57735027, -0.57735027,  0.57735027},
	{-0.57735027,  0.57735027, -0.57735027}, { 0.57735027, -0.57735027, -0.57735027}, {0.52573111, 0.85065081, 0.00000000},
	{-0.52573111, 0.85065081, 0.00000000}, {0.52573111, -0.85065081, 0.00000000}, {-0.52573111, -0.85065081, 0.00000000},
	{0.85065081, 0.00000000, 0.52573111}, {-0.85065081, 0.00000000, 0.52573111}};
	size_t k = 6;
	
	vector<Vector3d> vicini = trova_k_punti_vicini(punto, points, k);
	EXPECT_EQ(vicini.size(), 6);
	
	for (size_t i = 1; i < vicini.size(); ++i) {
        double d_prec = (vicini[i - 1] - vicini[i]).norm();
        for (size_t j = i + 1; j < vicini.size(); ++j) {
            double d_altro = (vicini[i - 1] - vicini[j]).norm();
            EXPECT_LE(d_prec, d_altro + 1e-9);
        }
	}
}

TEST(FileVerticiTest, PuntiCorretti) {
	vector<Vector3d> points = {Vector3d(0.57735027,  0.57735027,  0.57735027), Vector3d(-0.57735027, -0.57735027,  0.57735027), Vector3d(-0.57735027,  0.57735027, -0.57735027)};
	for(size_t i = 0; i < points.size(); i++){
		EXPECT_NEAR(points[i].norm(), 1.0, 1e-6);
	}
	
	map<array<int,3> , int> mappa_vertici;
	int id_vertice = 0;
	ostringstream s_g_Cell0Ds;
	bool duale = false;
	
	bool result = file_vertici(points, mappa_vertici, id_vertice, s_g_Cell0Ds, duale);
	
	EXPECT_TRUE(result);
	
	EXPECT_NE(mappa_vertici.size(), 0);
    EXPECT_NE(id_vertice, 0);
	
	string output = s_g_Cell0Ds.str();
	for (int i = 0; i < id_vertice; ++i) {
        string expected_line = to_string(i);
        EXPECT_NE(output.find(expected_line), string::npos);
	}
	
}

TEST(FileLatiTest, LatiCorretti) {
	vector<Vector3d> points = {Vector3d(0.57735027,  0.57735027,  0.57735027), Vector3d(-0.57735027, -0.57735027,  0.57735027), Vector3d(-0.57735027,  0.57735027, -0.57735027)};
	
	map<array<int,3> , int> mappa_vertici;
	map<pair<array<int,3>, array<int,3>>, int> mappa_lati;
	int id_lato = 0;
	ostringstream s_g_Cell1Ds;
	int b = 2;
	bool duale = false;
	
	bool result = file_lati(points, mappa_lati, mappa_vertici, id_lato, b, s_g_Cell1Ds, duale);
	
	EXPECT_TRUE(result);
	
	EXPECT_NE(mappa_lati.size(), 0);
    EXPECT_NE(id_lato, 0);
	
	string output = s_g_Cell1Ds.str();
	for (int i = 0; i < id_lato; ++i) {
        string expected_line = to_string(i);
        EXPECT_NE(output.find(expected_line), string::npos);
	}
	
}

TEST(FileLati_II_Test, LatiCorretti) {
	vector<Vector3d> punti_unici = {
		Vector3d(0.577,0.577,0.577), Vector3d(0,0,1), Vector3d(-0.577,-0.577,0.577),
		Vector3d(-1,0,0), Vector3d(-0.577,0.577,-0.577), Vector3d(0,1,0), Vector3d(-0.577,0.577,0.577)
	};
	
	map<array<int,3> , int> mappa_vertici;
	map<pair<array<int,3>, array<int,3>>, int> mappa_lati;
	int id_lato = 0;
	ostringstream s_g_Cell1Ds;
	int b = 1;
	bool duale = false;
	
	bool result = file_lati_II(punti_unici, mappa_lati, mappa_vertici, id_lato, b, s_g_Cell1Ds, duale);
	
	EXPECT_TRUE(result);
	
	EXPECT_NE(mappa_lati.size(), 0);
    EXPECT_NE(id_lato, 0);
	
	string output = s_g_Cell1Ds.str();
	for (int i = 0; i < id_lato; ++i) {
        string expected_line = to_string(i);
        EXPECT_NE(output.find(expected_line), string::npos);
	}
	
}

TEST(FileFacceTest, FacceCorrette) {
	vector<Vector3d> points = {Vector3d(0.57735027,  0.57735027,  0.57735027), Vector3d(-0.57735027, -0.57735027,  0.57735027), Vector3d(-0.57735027,  0.57735027, -0.57735027)};
	
	map<array<int,3> , int> mappa_vertici;
	map<pair<array<int,3>, array<int,3>>, int> mappa_lati;
	map<int, pair<Vector3i, Vector3i>> mappa_facce;
	int id_faccia = 0;
	int b = 2;
	
	bool result = file_facce(points, mappa_facce, mappa_lati, mappa_vertici, id_faccia, b);
	
	EXPECT_TRUE(result);
	
	EXPECT_NE(mappa_facce.size(), 0);
    EXPECT_NE(id_faccia, 0);
	
}

TEST(FileFacce_II_Test, FacceCorrette) {
    vector<Vector3d> punti_unici = {
        Vector3d(0.577,0.577,0.577), Vector3d(0,0,1), Vector3d(-0.577,-0.577,0.577),
        Vector3d(-1,0,0), Vector3d(-0.577,0.577,-0.577), Vector3d(0,1,0), Vector3d(-0.577,0.577,0.577)
    };

    map<array<int,3> , int> mappa_vertici;
    map<pair<array<int,3>, array<int,3>>, int> mappa_lati;
    map<array<int, 3>, int> mappa_facce_2;
    map<int, pair<Vector3i, Vector3i>> mappa_facce;
    int id_faccia = 0;
    ostringstream s_g_Cell2Ds;
    int b = 0;
    bool duale = false;

    for (size_t i = 0; i < punti_unici.size(); ++i) {
        mappa_vertici[to_array(punti_unici[i].normalized())] = i;
    }

    mappa_lati[{to_array(punti_unici[0].normalized()), to_array(punti_unici[1].normalized())}] = 0;
    mappa_lati[{to_array(punti_unici[1].normalized()), to_array(punti_unici[0].normalized())}] = 0;

    mappa_lati[{to_array(punti_unici[1].normalized()), to_array(punti_unici[2].normalized())}] = 1;
    mappa_lati[{to_array(punti_unici[2].normalized()), to_array(punti_unici[1].normalized())}] = 1;

    mappa_lati[{to_array(punti_unici[2].normalized()), to_array(punti_unici[0].normalized())}] = 2;
    mappa_lati[{to_array(punti_unici[0].normalized()), to_array(punti_unici[2].normalized())}] = 2;

    bool result = file_facce_II(punti_unici, mappa_facce_2, mappa_lati, mappa_vertici, id_faccia, b, s_g_Cell2Ds, mappa_facce, duale);

    EXPECT_TRUE(result);
    EXPECT_NE(mappa_facce.size(), 0);
    EXPECT_NE(mappa_facce_2.size(), 0);
    EXPECT_NE(id_faccia, 0);

    string output = s_g_Cell2Ds.str();
    for (int i = 0; i < id_faccia; ++i) {
        string expected_line = to_string(i);
        EXPECT_NE(output.find(expected_line), string::npos);
    }
}

TEST(FilePoliedroTest, PoliedroCorretto) {
	int V = 4, L = 6, F = 4;
	ostringstream s_g_Cell3Ds;
	
	bool result = file_poliedro(F, V, L, s_g_Cell3Ds);
	
	EXPECT_TRUE(result);
	string output = s_g_Cell3Ds.str();
	EXPECT_EQ(output, "0 4 6 4 0 1 2 3 0 1 2 3 4 5 0 1 2 3 ");
	
}

TEST(CellConverterTest, ConversioneCorretta) {
	int V = 4, L = 6;
	map<array<int,3> , int> mappa_vertici;
	map<pair<array<int,3>, array<int,3>>, int> mappa_lati;
	
	EXPECT_NE(Cell0DsConverter(V, mappa_vertici).size(), 0);
	EXPECT_NE(Cell1DsConverter(L, mappa_vertici, mappa_lati).size(), 0);
}

TEST(GetLatoIDTest, LatiPresenti) {
    array<int,3> a = {1, 2, 3};
    array<int,3> b = {4, 5, 6};
    array<int,3> c = {7, 8, 9};

    map<pair<array<int,3>, array<int,3>>, int> mappa_lati;
    mappa_lati[{a, b}] = 0;

    EXPECT_EQ(getLatoID(a, b, mappa_lati), 0);
    EXPECT_EQ(getLatoID(b, a, mappa_lati), 0);

    EXPECT_EQ(getLatoID(a, c, mappa_lati), -1);
    EXPECT_EQ(getLatoID(c, b, mappa_lati), -1);
}

TEST(FileVerticiDualeTest, PuntiCorretti) {
	vector<Vector3d> points = {
        Vector3d(1, 0, 0),
        Vector3d(0, 1, 0),
        Vector3d(0, 0, 1)
    };
	
	map<array<int,3> , int> mappa_vertici;
	for (size_t i = 0; i < points.size(); ++i) {
        mappa_vertici[to_array(points[i].normalized())] = i;
    }

	map<int, pair<Vector3i, Vector3i>> mappa_facce = {
        {0, {Vector3i(0, 1, 2), Vector3i(0, 1, 2)}},
        {1, {Vector3i(1, 2, 0), Vector3i(1, 2, 0)}},
        {2, {Vector3i(2, 0, 1), Vector3i(2, 0, 1)}}
    };
		
	int F_s_g = 1;
    map<array<int,3> , int> mappa_vertici_duale;
	
	ostringstream s_g_Cell0Ds;
	
	auto result = file_vertici_duale(F_s_g, mappa_facce, mappa_vertici, mappa_vertici_duale, s_g_Cell0Ds);
	
	EXPECT_NE(mappa_vertici_duale.size(), 0);
	EXPECT_NE(result.first.size(), 0);
	EXPECT_NE(result.second.size(), 0);
	
	string output = s_g_Cell0Ds.str();
	istringstream iss(output);
    string line;
    int line_count = 0;

    while (getline(iss, line)) {
        string expected_prefix = to_string(line_count) + " ";
        EXPECT_TRUE(line.find(expected_prefix) == 0) << "Riga non valida: " << line;
        ++line_count;
    }
	
}

TEST(FileLatiDualeTest, LatiCorretti) {
    map<int, Vector3d> mappa_baricentri = {
        {0, Vector3d(0, 0, 0)},
        {1, Vector3d(1, 0, 0)},
        {2, Vector3d(0, 1, 0)}
    };

    map<int, pair<Vector3i, Vector3i>> mappa_facce = {
        {0, {Vector3i(0, 1, 2), Vector3i(0, 1, 2)}},
        {1, {Vector3i(0, 1, 2), Vector3i(0, 1, 2)}},
        {2, {Vector3i(0, 1, 2), Vector3i(0, 1, 2)}}
    };

    int L_s_g = 3;

    map<array<int,3>, int> mappa_vertici_duale;
    for (const auto& [id, bar] : mappa_baricentri) {
        mappa_vertici_duale[to_array(bar)] = id;
    }

    int id_lato = 0;
    ostringstream s_g_Cell1Ds;

    auto result = file_lati_duale(L_s_g, mappa_baricentri, mappa_facce, mappa_vertici_duale, id_lato, s_g_Cell1Ds);

    EXPECT_NE(result.size(), 0);
    EXPECT_NE(id_lato, 0);

    string output = s_g_Cell1Ds.str();
    for (int i = 0; i < id_lato; ++i) {
        string expected_line = to_string(i);
        EXPECT_NE(output.find(expected_line), string::npos);
	}
}

TEST(FileFacceDualeTest, FacceCorrette) {
    map<int, Vector3d> mappa_baricentri = {
        {0, Vector3d(0, 0, 0)},
        {1, Vector3d(1, 0, 0)},
        {2, Vector3d(0, 1, 0)}
    };

    map<int, pair<Vector3i, Vector3i>> mappa_facce = {
        {0, {Vector3i(0, 1, 2), Vector3i(0, 1, 2)}},
        {1, {Vector3i(0, 1, 2), Vector3i(0, 1, 2)}},
        {2, {Vector3i(0, 1, 2), Vector3i(0, 1, 2)}}
    };

    map<pair<array<int,3>, array<int,3>>, int> mappa_lati_duale = {
        {{to_array(mappa_baricentri[0]), to_array(mappa_baricentri[1])}, 0},
        {{to_array(mappa_baricentri[1]), to_array(mappa_baricentri[2])}, 1},
        {{to_array(mappa_baricentri[2]), to_array(mappa_baricentri[0])}, 2}
    };

    map<array<int,3>, int> mappa_vertici = {
        {to_array(mappa_baricentri[0]), 0},
        {to_array(mappa_baricentri[1]), 1},
        {to_array(mappa_baricentri[2]), 2}
    };

    int id_faccia_duale = 0;
    ostringstream s_g_Cell2Ds;

    auto result = file_facce_duale(mappa_facce, mappa_lati_duale, mappa_vertici, mappa_baricentri, id_faccia_duale, s_g_Cell2Ds);

    EXPECT_NE(result.size(), 0);
    EXPECT_NE(id_faccia_duale, 0);

    string output = s_g_Cell2Ds.str();
    for (int i = 0; i < id_faccia_duale; ++i) {
        string expected_line = to_string(i);
        EXPECT_NE(output.find(expected_line), string::npos);
    }
}

TEST(DijkstraTest, CamminoMinimoSemplice) {
    int n = 5;
    vector<vector<int>> adiac_nodi = {
        {1, 2},    
        {2, 3},    
        {3},       
        {4},       
        {}         
    };

    vector<vector<double>> adiac_pesi = {
        {2.0, 4.0},    
        {1.0, 7.0},    
        {3.0},         
        {1.0},         
        {}             
    };

    vector<int> expected_path = {0, 1, 2, 3, 4};
    vector<int> result = dijkstra(n, adiac_nodi, adiac_pesi, 0, 4);

    EXPECT_EQ(result, expected_path);
}