#pragma once
#ifndef MOD_GLOBAL_INIT_H
#define MOD_GLOBAL_INIT_H

#pragma once

#include <complex>
#include <array>
#include <vector>
#include <cstdint>
#include <string>
#include <cmath>



#include "Mod_Type.h"

// 常用物理常数（const变量在头文件中直接定义）
const double PI = 4.0 * std::atan(1.0);          
const double E0 = 8.85418782e-12;                // 真空介电常数
const double U0 = 4.0 * PI * 1e-7;               // 真空磁导率
const double IMPD = std::sqrt(U0 / E0);           // 自由波阻抗S
const double IMPD_P2 = IMPD * IMPD;               // 自由波阻抗的平方
const double SPEED = 1.0 / std::sqrt(E0 * U0);    // 光速

//相关路径
extern std::string PATH;      
extern std::string INPUT_FILE;
extern std::string SET_FILE;	
extern std::string DIELECTRIC_FILE;  						                         
extern std::string MAP_PATH;  

// 类型相关的常量数组
constexpr std::array<int, 15> TYP_LENGTH = { 2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 18, 14, 1 };
constexpr std::array<int, 15> TYP_FACE = { 0, 5, 6, 4, 6, 6, 5, 0, 1, 1, 4, 6, 6, 5, 0 }; 

/*
elm - type
defines the geometrical type of the n - th element :

1
2 - node line.

2
3 - node triangle.

3
4 - node quadrangle.

4
4 - node tetrahedron.

5
8 - node hexahedron.

6
6 - node prism.

7
5 - node pyramid.

8
3 - node second order line(2 nodes associated with the vertices and 1 with the edge).

9
6 - node second order triangle(3 nodes associated with the vertices and 3 with the edges).

10
9 - node second order quadrangle(4 nodes associated with the vertices, 4 with the edges and 1 with the face).

11
10 - node second order tetrahedron(4 nodes associated with the vertices and 6 with the edges).

12
27 - node second order hexahedron(8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume).

13
18 - node second order prism(6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces).

14
14 - node second order pyramid(5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face).

15
1 - node point.
*/



// S参数结果数据结构
constexpr std::array<int, 2> S_LENGTH = { 3, 9 };

// 物理属性数组维度常量
const int P_Y = 10;                               // PHYZ数组的列数
const int PM_X = 100;                             // PHYZ数组的行数


const int UN_S_CONS = 1e9;                      
const int MAX_PORT = 4;                         // 最大端口数量

extern std::array<int, 15> NUMBER_OF_TYPE;        // 每种类型的计数

// 点数据
extern std::vector<Point> PT;                     // 点集合
extern int N_P;                                   // 点的数量

// 临时网格数据
extern std::vector<Temp_Mesh> T_MS;               // 临时网格集合
extern int T_MS_N;                                // 临时网格的数量

// 端口数据
extern std::vector<Port> PORT_DATA;               // 端口数据集合
extern int N_PORT;                                // 端口数量

// 网格和积分线计数
extern int NODE_N;   

//!DIE & PEC & LINE & PORT MESH
extern int N_S;                                  
extern int DS_N;                                  
extern int PEC_N;                                 
extern int L_N;                                  

extern std::vector<Mesh> Sur_M;                   // 表面网格集合

// 连接关系计数
extern int CN_N;                                  // 总连接数量
extern int C_DS_N;                                // DS连接数量
extern int C_PEC_N;                               // PEC连接数量

// 连接关系数据
extern std::vector<int> DS_PEC_Branch;            // DS-PEC分支
extern std::vector<std::vector<int>> CN;          
extern std::vector<std::vector<int>> CN_P;        

// 重构后对应数据值
extern std::vector<std::vector<Ele_BL>> B_IND;    // 电感矩阵
extern std::vector<std::vector<Ele_BC>> B_CAP;    
extern std::vector<std::vector<Ele_P>> N_CAP;     //电容矩阵
extern std::vector<std::vector<Ele_CS>> B_CS;     
extern std::vector<std::vector<Ele_CS>> N_CS;     

// 积分相关数据
extern std::vector<std::vector<double>> PDD_N;    
extern std::vector<std::vector<double>> PD0_N;    
extern std::vector<std::vector<double>> LD0_N;    
extern std::vector<std::vector<Vector>> PV_DD;   
extern std::vector<std::vector<Vector>> PV_D0;    
extern std::vector<std::vector<Vector>> LV_D0;    

// 电阻矩阵
extern std::vector<std::vector<double>> B_RR;     // 边界电阻矩阵

// 单元网格积分数据
extern std::vector<std::vector<Ele_Mesh_Int>> MS_ELE;  // 网格单元积分


extern std::vector<std::vector<Ele_BC>> TC_CAP;   
extern std::vector<std::vector<Ele_CS>> TC_CS;    
extern std::vector<std::vector<Ele_P>> N_TC_CAP;  
extern int C_TC_N;                               
extern int TC_N;                                 


extern std::vector<std::vector<int>> TC_B2M;      
extern std::vector<std::vector<int>> TC_B2M_P;    
extern std::vector<int> TC_N2M;                   
extern std::vector<int> TC_M2N;                   

// 介电常数索引数组
extern std::array<int, 30> BLK_D;                 // 块介电常数索引

// 子问题和块计数
extern int N_SP;                                  // 子问题数量
extern int N_BLOCK;                               // 块数量

// 物理属性数组
extern std::array<std::array<int, P_Y>, PM_X> PHYZ;  // 物理属性数组
extern std::array<int, PM_X> DIE;                 // 介电常数数组

// S参数数据
extern std::vector<std::vector<std::vector<double>>> S_PAR_ABS;  // S参数幅度
extern std::vector<std::vector<std::vector<double>>> S_PAR_PHA;  // S参数相位

// 阻抗和导纳参数
extern std::vector<std::vector<std::vector<std::complex<double>>>> Z_PT;  // 阻抗参数
extern std::vector<std::vector<std::vector<std::complex<double>>>> Y_PT;  // 导纳参数


extern int INFI_G;
extern double Z_G;                              // 全局阻抗
extern int SLV_W;                               // 求解器选择
extern int Swap_Log;                            // 交换日志
extern long long FS, FE;                              // 起始\结束频率
extern int N_FP;                                // 频率点数
extern double PRECISION;                        // 计算精度
extern int DIM;                                 // 维度
extern std::vector<std::complex<double>> DIELECTRIC; // 介电常数数组

extern std::vector<double> CF;
extern double W;
extern double LAMDA_E;
extern double MAX_L;

extern ElementInfo ELM_TYPE;
extern int Solver_SET;
extern std::vector<std::vector<double>> LL, PP;
extern std::vector<std::vector<double>> A_E, A_ET;
extern int M_SZ;
#endif    // Global
