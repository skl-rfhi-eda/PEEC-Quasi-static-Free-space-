#ifndef MOD_CONNECT_FUNCTION_H
#define MOD_CONNECT_FUNCTION_H

#pragma once

#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>  
#include <cstdlib>  // For std::exit
#include "Mod_Function.h"



// 将面网格转换为线网格
Mesh Trans_L(const Mesh& ms_in, int adj);

// 将临时网格转换为PEC表面三角形
Mesh Trans_MS_PEC(const Temp_Mesh& ms_in);

// 从四面体网格中提取一个三角形面
Mesh Trans_MS_TETRA(const Temp_Mesh& ms_in, int un_p);

// 判断两个表面是否相同(共享顶点数量)
int Same_Sur(const Mesh& ms_a, const Mesh& ms_b);

// 计算两个表面的邻接关系
void Adj_Sur(const Mesh& ms_a, const Mesh& ms_b, int& adj_a, int& adj_b, int& cm_p);

// 计算三角形表面与线表面的公共顶点数量及邻接编码
void Adj_Sur_L(const Mesh& ms_a, const Mesh& ms_b, int& adj_a, int& adj_b, int& cm_p);

// 设置网格的基本属性
void Set_Up_Mesh(Mesh& ms_in);

// 判断两个四面体网格是否相邻
void Judge_Tetra_C(Temp_Mesh& ms_a, Temp_Mesh& ms_b, int& err);

// 计算一个与平面垂直且指向外部的法向量
Vector Oth_V(const Point* cp, const Point& up);

// 在所有PEC表面中查找包含指定点的端口节点
int Port_Node(Point& p_in);


#endif // MOD_CONNECT_FUNCTION_H