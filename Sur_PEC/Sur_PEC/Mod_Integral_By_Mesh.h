#ifndef MOD_INTEGRAL_BY_MESH_H
#define MOD_INTEGRAL_BY_MESH_H

#pragma once

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>

#include"Mod_Global_Init.h"
#include "Mod_Type.h"
#include "Mod_Function.h"   
#include "Gaussion.h"
#include "Mod_Singular_Treatment.h"

Ele_Mesh_Int Ele_Cal(int IND_A, int IND_B);
Ele_Mesh_Int Trans_Ele(int IND_A, int IND_B);

Ele_Mesh_Int Ele_Integral_T(int IND_A, int IND_B);
Ele_Mesh_Int Ele_Integral_F(int IND_A, int IND_B);



#endif   // MOD_INTEGRAL_BY_MESH_H
