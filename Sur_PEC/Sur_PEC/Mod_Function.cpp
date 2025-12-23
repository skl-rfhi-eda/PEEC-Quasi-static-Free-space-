#include "Mod_Function.h"

//const double PI = 4.0 * std::atan(1.0);          // 圆周率


//向量运算部分函数

/**
* @brief 计算正交平面的向量,计算两个向量点积。
* @param V1 第一个向量。V2 第二个向量。
*/
double V_Dot(const Vector& V1, const Vector& V2) {
	double V_Dot = V1[0] * V2[0] + V1[1] * V2[1] + V1[2] * V2[2];
	return V_Dot;
}  //   V_Dot

Vector V_Add(const Vector& V1, const Vector& V2) {
	Vector result;
	for (int i = 0; i < 3; ++i) {
		result[i] = V1[i] + V2[i];
	}
	return result;
}

Vector V_Sub(const Vector& V1, const Vector& V2) {
	Vector result;
	for (int i = 0; i < 3; ++i) {
		result[i] = V1[i] - V2[i];
	}
	return result;
}

Vector V_Mul(const double s, const Vector& V1) {
	Vector result;
	for (int i = 0; i < 3; ++i) {
		result[i] = s * V1[i];
	}
	return result;
}

Vector V_Div(const double s, const Vector& V1) {
	Vector result;
	for (int i = 0; i < 3; ++i) {
		result[i] = V1[i] / s;
	}
	return result;
}

/**
* 计算两个三维向量的叉积（向量积）
* @param V1 第一个输入向量
* @param V2 第二个输入向量
* @return 叉积结果向量（V1 × V2）
*/
Vector V_Cross(const Vector& V1, const Vector& V2) {
	Vector result;
	for (int i = 0; i < 3; ++i) {
		int idx1 = (i + 1) % 3;
		int idx2 = (i + 2) % 3;
		result[i] = V1[idx1] * V2[idx2] - V1[idx2] * V2[idx1];
	}
	return result;
}


//计算两点间距离
double Distance(const Vector& V1, const Vector& V2) {
	double Dis = 0.0;
	int i;
	for (i = 0; i < 3; ++i) {
		Dis += (V1[i] - V2[i]) * (V1[i] - V2[i]);
	}
	Dis = std::sqrt(Dis);
	return Dis;
}// Distance

/**
* 计算三角形的面积
* @param P_A， P_B，P_C 三角形的三个点
* @return Triangle_Area三角形面积
*/
double  Trian_Area(const Vector& P_A, const Vector& P_B, const Vector& P_C) {
	std::array<double, 3> length;

	// 使用 Distance 函数计算三条边的长度
	length[0] = Distance(P_A, P_B); // 边 a (在 P_A 和 P_B 之间)
	length[1] = Distance(P_A, P_C); // 边 b (在 P_A 和 P_C 之间)
	length[2] = Distance(P_B, P_C); // 边 c (在 P_B 和 P_C 之间)

	// 计算半周长(semi - perimeter)
	double ind = (length[0] + length[1] + length[2]) / 2.0;

	// 应用海伦公式: Area = sqrt(s * (s-a) * (s-b) * (s-c))
	double Area_Squared = ind * (ind - length[0]) * (ind - length[1]) * (ind - length[2]);

	double Triangle_Area = std::sqrt(Area_Squared);
	return Triangle_Area;
}

double ABSS(std::complex<double> s) {
	double ABSS = 20.0 * std::log10(abs(s));
	return ABSS;
}// abss

//计算一个复数的相位（相位角），并以度为单位返回结果
double Phase(std::complex<double> z) {
	double a = z.real();
	double b = z.imag();
	const double eps = 1.0e-8;
	double phase_val = 0.0;
	// 情况1：复数接近原点（实部和虚部都接近0），相位定义为0°
	if (std::abs(a) < eps && std::abs(b) < eps) {
		phase_val = 0.0;
		return phase_val;
	}
	// 情况2：虚部接近0（点在实轴上）
	if (std::abs(b) < eps && a > 0.0) {
		phase_val = 0.0;
		return phase_val;
	}
	// 情况3：虚部接近0且实部为负，相位为180°
	if (std::abs(b) < eps && a < 0.0) {
		phase_val = 180.0;
		return phase_val;
	}
	// 情况4：实部接近0且虚部为正（点在正虚轴上），相位为90°
	if (std::abs(a) < eps && b > 0.0) {
		phase_val = 90.0;
		return phase_val;
	}
	// 情况5：在负虚轴上 (a=0, b<0)，相位为-90°
	if (std::abs(a) < eps && b < 0.0) {
		phase_val = -90.0;
		return phase_val;
	}// 
	//情况6：第一象限(a > 0, b > 0)
	if (a > 0.0 && b > 0.0) {
		phase_val = std::atan(b / a) * 180.0 / PI;
		return phase_val;
	}
	// 情况7：第二象限 (a<0, b>0)
	if (a < 0.0 && b > 0.0) {
		phase_val = 180.0 - std::atan(std::abs(b / a)) * 180.0 / PI;
		return phase_val;
		//phase_val = -90.0;
	}
	// 情况8：第三象限 (a<0, b<0)
	if (a < 0.0 && b < 0.0) {
		phase_val = -180.0 + std::atan(std::abs(b / a)) * 180.0 / PI;
		return phase_val;
	}
	// 情况9：第四象限 (a>0, b<0)
	if (a > 0.0 && b < 0.0) {
		phase_val = -std::atan(std::abs(b / a)) * 180.0 / PI;
		return phase_val;
	}
	return 0.0;
}// Phase

//矩阵转置函数
/**
* 计算下三角矩阵中的索引转换,将二维矩阵的下三角部分（包括对角线按行优先顺序映射到一维数组的索引转换函数。(从对角线开始存储)
* @param IND_A 行索引（原矩阵）
* @param IND_B 列索引（原矩阵）
* @param SZ_A 矩阵的大小（维度）
* @return 转换后的线性索引
*/
//上三角每行存储个数按等差数列存储计算（ 第0行：SZ_A ， 第一行：SZ_A -1，……， 第i行：SZ_A -i）
int TRANS_IN_LU(int IND_A, int IND_B, int SZ_A) {
	return (SZ_A + SZ_A - IND_A + 1) * (IND_A) / 2 + IND_B - IND_A;
}

/**
* 计算方阵中的索引转换
* @param IND_A 行索引（原矩阵）
* @param IND_B 列索引（原矩阵）
* @param SZ_A 矩阵的大小（维度）
* @return 转换后的线性索引
*/
int TRANS_IN_SQ(int IND_A, int IND_B, int SZ_A) {
	// 计算方阵按行存储时的线性索引
	return SZ_A * (IND_A)+IND_B;
}


/**
* 功能：根据条件打印空格或点，用于输出格式化
* @param i 循环索引
* @param t 总迭代次数
* @param A 控制空格数量的参数
*/
void PRT(int i, int T, int A) {
	int j;
	if (i == 1) {
		for (j = 0; j <= 5 * A; ++j) {
			std::cout << " " << std::flush;
		}
	}
	// 当i是(1 + t/20)的倍数时，打印一个点
	if (i % (1 + static_cast<int>(T / 20.0)) == 0) {
		std::cout << "." << std::flush;
	}
	if (i == T) {
		std::cout << " " << std::endl; //输出换行符
	}
}// PRT

//面法向量
Vector Othognal_Vector(std::array<Point, 3>& CP, Point& UP) {
	Vector V_1 = V_Sub(CP[1].X, CP[0].X);
	Vector V_2 = V_Sub(CP[2].X, CP[0].X);
	Vector normal_vec = V_Cross(V_1, V_2);
	normal_vec = V_Div(std::sqrt(V_Dot(normal_vec, normal_vec)), normal_vec);
	return normal_vec;
}

//边法向量(给定一条边和一个点)
Vector Oth_S(std::array<Point, 2>& CP, Point& UP) {
	Vector CP_vector = V_Sub(CP[1].X, CP[0].X);
	Vector UP_vector = V_Sub(UP.X, CP[0].X);
	Vector triangle_normal = V_Cross(UP_vector, CP_vector);
	Vector in_plane_normal = V_Cross(triangle_normal, CP_vector);
	return V_Div(std::sqrt(V_Dot(in_plane_normal, in_plane_normal)), in_plane_normal);
}

/**
	* @brief  计算由四个点定义的四面体的体积。
	* @param V 包含四个顶点的数组。
	* @return 返回四面体的体积。
	*/
double Tetrahedral_Volum(std::array<Point, 4>& V) {
	int i, j;
	double temp = 0.0;
	std::array<std::array<double, 3>, 3> ind;
	double volume = 0.0;
	// 或者直接使用三重积 = v1 · (v2 x v3) double triple_product = V_Dot(v1, V_Cross(v2, v3)); V=abs(triple_product) / 6.0;
	for (i = 0; i < 3; ++i) {
		for (j = 0; j < 3; ++j) {
			ind[i][j] = V[i + 1].X[j] - V[0].X[j]; //计算基于V[0].X[j]三条矢量边
		}
	}
	for (i = 0; i < 3; ++i) {
		double positive = ind[i % 3][0] * ind[(i + 1) % 3][1] * ind[(i + 2) % 3][2];  // 行列式正项
		double negative = ind[i % 3][0] * ind[(i + 2) % 3][1] * ind[(i + 1) % 3][2];  // 行列式负项
		temp = (positive - negative);
		volume += temp;
	}
	volume = std::fabs(volume) / 6.0;
	return volume;

}// Tetrahedral_Volum

/**
* 判断三维点P是否在由3个顶点组成的三角形MS所在平面内且位于三角形内部
* @param P 待判断的点（三维）
* @param MS 三角形的3个顶点数组（MS[0]~MS[2]，对应原Fortran的MS(1)~MS(3)）
* @return: 1：点在三角形内部；0：点不在三角形内部（或不在平面内）
*/
int P_In_S(Point& P, std::array<Point, 3> MS) {
	Vector V0, V1, V2;
	double INV, U, V, dot00, dot01, dot02, dot11, dot12, vol;
	std::array<Point, 4> temp_p;
	//std::array<int, 4> V_IND;
	int p_in_side = 0;

	temp_p[0] = P;

	for (int j = 0; j < 3; ++j) {
		temp_p[j + 1] = MS[j];
	}
	vol = Tetrahedral_Volum(temp_p);

	for (int j = 0; j < 3; ++j) {
		V0[j] = temp_p[3].X[j] - temp_p[1].X[j];
		V1[j] = temp_p[2].X[j] - temp_p[1].X[j];
		V2[j] = temp_p[0].X[j] - temp_p[1].X[j];
	}

	dot00 = V_Dot(V0, V0);
	dot01 = V_Dot(V0, V1);
	dot02 = V_Dot(V0, V2);
	dot11 = V_Dot(V1, V1);
	dot12 = V_Dot(V1, V2);

	// 计算并求解 u, v(重心坐标)
	INV = 1.0 / (dot00 * dot11 - dot01 * dot01);
	U = (dot11 * dot02 - dot01 * dot12) * INV;
	//
	if (U > -1e-10 && U < 1 + 1e-10) {
		V = (dot00 * dot12 - dot01 * dot02) * INV;
		if (V > -1e-10 && V < 1 + 1e-10) {
			if (U + V < 1+ 1e-10 && vol < 1e-20) {
				p_in_side = 1;
			}
		}
	}
	return p_in_side;
}// P_In_S

double Get_Time() {
	auto now = std::chrono::system_clock::now();  // time_point
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
		now.time_since_epoch());  // 转换为秒
	return static_cast<double>(duration.count());
}

void Time_Diff(double t1, double t2) {
	std::cout << " Finished. Using:" <<std::fabs(t2 - t1)<< "ms" <<std::endl;
	std::cout << std::endl;
}




