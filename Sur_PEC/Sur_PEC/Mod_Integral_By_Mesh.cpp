#include "Mod_Integral_By_Mesh.h"


// @brief Trans_Ele 函数：根据条件选择计算或转换元素间的相互作用矩阵。
 
Ele_Mesh_Int Trans_Ele(int IND_A, int IND_B) {
	Ele_Mesh_Int result;
	if ((Sur_M[IND_A].M_TYP[0] != Sur_M[IND_A].M_TYP[1]) && Sur_M[IND_A].IS_PEC == 1) {
		std::cout << "HAHAHA" << std::endl;
		result = Ele_Cal(IND_A, IND_B);
	}
	else {

		result = MS_ELE[IND_B][IND_A];

		// 转置L矩阵
		for (int i = 0; i < 3; ++i) {
			for (int j = i; j < 3; ++j) {
				std::swap(result.L[i][j], result.L[j][i]);
			}
		}
		// 转置C矩阵
		for (int i = 0; i < 3; ++i) {
			for (int j = i; j < 3; ++j) {
				std::swap(result.C[i][j], result.C[j][i]);
			}
		}
		// 转置CS矩阵
		for (int i = 0; i < 3; ++i) {
			for (int j = i; j < 3; ++j) {
				std::swap(result.CS[i][j], result.CS[j][i]);
			}
		}
	}
	return result;
}


/**
 * @brief Ele_Cal 函数：根据两个元素间的距离，选择不同的积分方法来计算相互作用。
 * @param IND_A 整数，第一个元素的索引 。
 * @param IND_B 整数，第二个元素的索引 。
 * @return ELE_MESH_INT 计算得到的元素相互作用数据。
 */
Ele_Mesh_Int Ele_Cal(int IND_A, int IND_B) {
	Ele_Mesh_Int result;
	double Dist = 0.0;

	Dist = Distance(Sur_M[IND_A].MID_CO, Sur_M[IND_B].MID_CO);

	if (Dist < 2 * MAX_L) {
		result = Ele_Integral_T(IND_A, IND_B);
	}
	else {
		// 远场计算
		result = Ele_Integral_F(IND_A, IND_B);
	}

	if (IND_A == IND_B) {
		// 自作用时，CS矩阵清零
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				result.CS[i][j] = 0.0;
			}
		}
	}

	return result;
}


/*
 * @brief Ele_Integral_T 子程序：执行近场积分计算。
 * @param IND_A 整数，第一个元素的索引; param IND_B 整数，第二个元素的索引。
 * @param MS_Ele 引用类型，用于存储输出的计算结果。
 */
Ele_Mesh_Int Ele_Integral_T(int IND_A, int IND_B) {
	std::array< Point, 3> V_PA;
	std::array< Point, 3> V_PB;
	Point  U_PA;
	std::array < Vector, 3> V_A;

	std::array<double, 3> XA;

	Ele_Mesh_Int MS_Ele;   //函数返回输出值
	FKNResult result_F_K_N;   //F_K_N函数返回值

	std::array<Vector, 3> LC_SINGU, CS_SINGU;
	Vector DE_P_SINGU;
	double P_SINGU;
	

	std::array<std::array<double, 3>, 3> L_TEMP{ 0 }, C_TEMP{ 0 }, CS_TEMP{ 0 };
	double DE_P_TEMP = 0.0;
	double P_TEMP = 0.0;
	std::array<std::array<double, 3>, 3> TEMP_NT{ 0 }; 

	int NA = 9;  
	int NB = 5;

	// 初始化网格输出结构
	const Mesh& MS_A = Sur_M[IND_A];   
	const Mesh& MS_B = Sur_M[IND_B];   

	for (int i = 0; i < 3; ++i) {
		V_PA[i] = PT[MS_A.P[i]];   
		V_PB[i] = PT[MS_B.P[i]];   
	}

	int num_gauss_points = SN[NA];
	for (int GA = 0; GA < num_gauss_points; ++GA) {
		for (int i = 0; i < 3; ++i) {
			XA[i] = V_PA[0].X[i] * SG[NA][GA][0] +
				V_PA[1].X[i] * SG[NA][GA][1] +
				V_PA[2].X[i] * SG[NA][GA][2];
		}
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				V_A[i][j] = XA[j] - V_PA[i].X[j];
				//if (V_A[i][j] < XA[j] / 1e5)
				//	V_A[i][j] = 0;
			}
		}
		for (int i = 0; i < 3; ++i) {
			U_PA.X[i] = XA[i]; 
		}
		// 调用核心计算函数 F_K_N
		result_F_K_N = F_K_N(V_PB, U_PA);

		//_F_K_N函数返回结果的赋值操作
		CS_SINGU = result_F_K_N.CS_BUF;
		LC_SINGU = result_F_K_N.L_BUF;
		DE_P_SINGU = result_F_K_N.DE_P_BUF;
		P_SINGU = result_F_K_N.P_BUF;
		
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				L_TEMP[i][j] = L_TEMP[i][j] + SW[NA][GA] * V_Dot(V_A[i], LC_SINGU[j]);
			}
		}
		P_TEMP = P_TEMP +  SW[NA][GA] * P_SINGU;
		
	}

	// 根据面积进行归一化,计算P_TEMP
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			L_TEMP[i][j] /= MS_B.AREA;
		}
	}
	P_TEMP /= MS_B.AREA;

	// 计算 C_TEMP
	double dot_product = V_Dot(MS_A.NOR_S, MS_B.NOR_S);
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			C_TEMP[i][j] = L_TEMP[i][j] * dot_product;
		}
	}

	// 最终赋值给输出参数
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			MS_Ele.L[i][j] = L_TEMP[i][j] / (16.0 * PI);    //计算 MS_Ele.L
		}
	}
	MS_Ele.P = P_TEMP / (4.0 * PI);   //计算 MS_Ele.P

	return MS_Ele;
	
}


/**
 * @brief Ele_Integral_F 函数：执行远场积分计算。
 * @param IND_A 整数，第一个（源）元素的索引。
 * @param IND_B 整数，第二个（场）元素的索引。
 * @param MS_Ele 引用类型，用于存储输出的计算结果。
 */
Ele_Mesh_Int Ele_Integral_F(int IND_A, int IND_B) {

	// 获取单元顶点坐标
	std::array<Point, 3> V_PA;
	std::array<Point, 3> V_PB;
	Point U_PA;
	std::array<double, 3> XA;    
	std::array<double, 3> XB;   

	std::array<std::array<double, 3>, 3> L_TEMP{ 0 };
	double P_TEMP = 0.0;
	double WT;
	int NA, NB, GA, GB;
	double DIST;
	double R;
	Vector V_R;

	Ele_Mesh_Int MS_Ele;

	// 获取网格单元数据
	Mesh& MS_A = Sur_M[IND_A];
	Mesh& MS_B = Sur_M[IND_B];

	for (int i = 0; i < 3; ++i) {
		V_PA[i] = PT[MS_A.P[i]];
		V_PB[i] = PT[MS_B.P[i]];
	}

	// 根据距离动态确定高斯积分的阶数
	DIST = Distance(MS_A.MID_CO, MS_B.MID_CO);   // 计算两个单元的中心点距离
	NA = std::max(0, static_cast<int>(5.0 - 0.5 * DIST / LAMDA_E)-1);  // 网格A高斯积分阶数
	NB = std::max(NA - 1, 0);

	// 双重高斯积分循环
	int num_gauss_points_A = SN[NA];  // 计算源单元上的积分点数
	int num_gauss_points_B = SN[NB];  // 计算场单元上的积分点数

	for (GA = 0; GA < num_gauss_points_A; ++GA) {
		XA[0] = V_PA[0].X[0] * SG[NA][GA][0] + V_PA[1].X[0] * SG[NA][GA][1] + V_PA[2].X[0] * SG[NA][GA][2];
		XA[1] = V_PA[0].X[1] * SG[NA][GA][0] + V_PA[1].X[1] * SG[NA][GA][1] + V_PA[2].X[1] * SG[NA][GA][2];
		XA[2] = V_PA[0].X[2] * SG[NA][GA][0] + V_PA[1].X[2] * SG[NA][GA][1] + V_PA[2].X[2] * SG[NA][GA][2];

		// 计算源单元上的基函数向量 V_A（计算网格A顶点到积分点XA的向量V_A）
		std::array<Vector, 3> V_A;
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				V_A[i][j] = XA[j] - V_PA[i].X[j];
			}
		}

		for (GB = 0; GB < num_gauss_points_B; ++GB) {
			XB[0] = V_PB[0].X[0] * SG[NB][GB][0] + V_PB[1].X[0] * SG[NB][GB][1] + V_PB[2].X[0] * SG[NB][GB][2];
			XB[1] = V_PB[0].X[1] * SG[NB][GB][0] + V_PB[1].X[1] * SG[NB][GB][1] + V_PB[2].X[1] * SG[NB][GB][2];
			XB[2] = V_PB[0].X[2] * SG[NB][GB][0] + V_PB[1].X[2] * SG[NB][GB][1] + V_PB[2].X[2] * SG[NB][GB][2];

			// 计算场单元上的基函数向量 V_B
			std::array<Vector, 3> V_B;
			for (int i = 0; i < 3; ++i) {
				for (int j = 0; j < 3; ++j) {
					V_B[i][j] = XB[j] - V_PB[i].X[j];
				}
			}

			// 计算两积分点之间的距离 R
			R = Distance(XA, XB);
			if (R != 0.0) {
				V_R[0] = (XB[0] - XA[0]) / R;
				V_R[1] = (XB[1] - XA[1]) / R;
				V_R[2] = (XB[2] - XA[2]) / R;
			}

			// 计算当前积分点对的权重
			WT = (SW[NA][GA] * SW[NB][GB] / R);

			// 累加标量势 P
			P_TEMP += WT;

			// 累加矢量势 L
			for (int i = 0; i < 3; ++i) {
				for (int j = 0; j < 3; ++j) {
					L_TEMP[i][j] += WT * V_Dot(V_A[i], V_B[j]);
				}
			}
		}
	}

	// 最终缩放并赋值给输出结构体
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			MS_Ele.L[i][j] = L_TEMP[i][j] / (16.0 * PI);
		}
	}
	MS_Ele.P = P_TEMP / (4.0 * PI);

	return MS_Ele;
}   // Ele_Integral_F


