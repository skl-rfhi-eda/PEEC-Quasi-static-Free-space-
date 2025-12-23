#ifndef MOD_SAVE_CIRCUIT
#define MOD_SAVE_CIRCUIT

#pragma once

#include <complex>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <map>
#include <set>

#include"Mod_Global_Init.h"
#include "Mod_Type.h"
#include "Mod_Function.h"
#include "Mod_MKL_Interface.h"

void Save_PEEC_Circuit();
int Save_Netlist();
#endif