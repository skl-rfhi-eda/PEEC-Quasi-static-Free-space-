#ifndef MOD_ELEMENT_BUILT_H
#define MOD_ELEMENT_BUILT_H

#pragma once

#include <iostream>
#include <vector>


#include"Mod_Global_Init.h"
#include "Mod_Type.h"
#include "Mod_Function.h"
#include "Gaussion.h"
#include "Mod_Integral_By_Mesh.h"

// 主调度函数
void Element_Built();

// 内部实现函数
void Ini_Element_Sz();
void Build_MS_Ele();
void Re_Construct_Ele();
//void Built_Branch_L();
//void Built_Node_P();

#endif // MOD_ELE_BUILT_H