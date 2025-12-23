#ifndef MOD_SINGULAR_TREATMENT_H
#define MOD_SINGULAR_TREATMENT_H

#pragma once

#include <array>
#include <vector>
#include <iostream>
#include"Mod_Global_Init.h"
#include "Mod_Type.h"
#include "Mod_Function.h"

struct FKNResult {
	double P_BUF;
	Vector DE_P_BUF;
	std::array<Vector, 3> L_BUF;
	std::array<Vector, 3> CS_BUF;
};

FKNResult F_K_N(std::array<Point, 3>& Point_v_IN, Point& Point_t_IN);

std::vector<double> F_K_1_N(std::array<Point, 3>& Point_v_IN, Point& Point_t_IN, int N_in);

std::vector<std::array<double, 3>> F_K_2_N(std::array<Point, 3>& Point_v_IN, Point& Point_t_IN,
	int Point_U_IN, int N_in);

#endif



