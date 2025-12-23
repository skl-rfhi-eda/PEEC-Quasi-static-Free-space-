#ifndef MOD_FUNCTION_H
#define MOD_FUNCTION_H

#pragma once

#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>  
#include <chrono>
#include <thread>

#include"Mod_Global_Init.h"
#include "Mod_Type.h"

#define MAX(a, b) ((a) >= (b) ? (a) : (b))
#define MIN(a, b) ((a) <= (b) ? (a) : (b))

double ABSS(std::complex<double> s);  //计算复数的幅度并以db为单位返回
double Phase(std::complex<double> z); //计算复数的相位角

int TRANS_IN_LU(int IND_A, int IND_B, int SZ_A);
int TRANS_IN_SQ(int IND_A, int IND_B, int SZ_A);

void PRT(int i, int T, int A);

double V_Dot(const Vector& V1, const Vector& V2);  //两个三维向量点乘

Vector V_Add(const Vector& V1, const Vector& V2); //两个三维向量相加

Vector V_Sub(const Vector& V1, const Vector& V2); //两个三维向量相减

Vector V_Mul(const double s, const Vector& V1); //三维向量乘法

Vector V_Div(const double s, const Vector& V1); //三维向量除法

Vector V_Cross(const Vector& V1, const Vector& V2);  //叉乘

double  Trian_Area(const Vector& P_A, const Vector& P_B, const Vector& P_C);  //计算三角形面积

double Distance(const Vector& V1, const Vector& V2);  //计算两个三维点的欧几里得距离

Vector Oth_S(std::array<Point, 2>& CP, Point& UP);  //计算两个点连线的正交向量

Vector Othognal_Vector(std::array<Point, 3>& CP, Point& UP); //面法向量

double Tetrahedral_Volum(std::array<Point, 4>& V);  //计算四面体体积,要改

int P_In_S(Point& P, std::array<Point, 3> MS);  //判断点是否在三角形内

//读取机器时间
double Get_Time();

//计算并打印时间差
void Time_Diff(double T_S, double T_E);

#endif // 

