#ifndef GAUSSION_H
#define GAUSSION_H

#pragma once
#include <vector>
#include <array>
#include <stdexcept>
#include <cmath>
#include "Mod_Type.h"


// 1D高斯积分权重（AK[阶数][点数]，阶数0-9对应1-10阶）
extern std::array<std::array<double, 10>, 10> AK;
// 1D高斯积分点坐标（NK[阶数][点数]，阶数0-9对应1-10阶）
extern std::array<std::array<double, 10>, 10> NK;

// 归一化后的高斯点坐标（U_NK = (NK + 1)/2，范围映射到[0,1]）
extern std::array<std::array<double, 10>, 10> U_NK;
// 高斯点坐标的平方根（S_NK = sqrt(U_NK)）
extern std::array<std::array<double, 10>, 10> S_NK;
// 高斯点坐标的立方根（Q_NK = U_NK^(1/3)）
extern std::array<std::array<double, 10>, 10> Q_NK;

// 二维权重数组（T_AK[阶数].W[行][列]，阶数0-9对应1-10阶）
extern std::array<W_NK, 10> T_AK;

// 三角形单元高斯积分参数
extern std::array<std::array<std::array<double, 3>, 25>, 10> SG;  // 坐标（阶数0-9，点数0-24，x/y/z分量0-2）
extern std::array<std::array<double, 25>, 10> SW;                 // 权重（阶数0-9，点数0-24）
extern std::array<int, 10> SN;                                     // 各阶点数（阶数0-9对应1-10阶的点数）

// 四面体单元高斯积分参数
extern std::array<std::array<std::array<double, 4>, 15>, 5> VG;   // 坐标（阶数0-4，点数0-14，4个分量）
extern std::array<std::array<double, 15>, 5> VW;                  // 权重（阶数0-4，点数0-14）
extern std::array<int, 7> VN;                                      // 各阶点数（阶数0-4对应1-5阶的点数，后2位补0）

// 矩形单元高斯积分参数
extern std::array<int, 10> RN;                                     // 各阶点数（阶数0-9对应1-10阶，为1D点数的平方）
extern std::vector<std::vector<std::array<double, 2>>> RG;
extern std::vector<std::vector<double>> RW;


// 动态数组
extern std::vector<std::vector<std::vector<double>>> SG_C;  // 三角形动态坐标数组
extern std::vector<std::vector<double>> SW_C;               // 三角形动态权重数组
extern std::vector<int> SN_C;                               // 三角形动态点数数组



int ARR_MULTI(int n);

void initialize_arrays();

void Ini_Gaussion();


#endif  // GAUSSION_H