#include "Mod_Singular_Treatment.h"

FKNResult F_K_N(std::array<Point, 3>& Point_v_IN, Point& Point_t_IN)
{
	// 返回值
	FKNResult result;

	std::array<double, 3> R, R_0, L, S_P, S_M, T_0, BELTA;
	Vector Normal_S_vec, U_vec, V_vec;
	std::array<double, 2> K_1_N = { 0.0 };
	std::array<std::array<double, 2>, 3> I_0;

	std::array<Vector, 3> M_vec, K_2_N, K_4_N, r_q;
	std::array<Vector, 2> I_M = { 0.0 };
	Vector K_3_N;
	std::array<Point, 3> Point_v;
	std::array<Point, 2> PN;
	Point Point_t, P_T, PP;


	Point_v[0].X = { 0.0,0.0,0.0 };
	Point_v[1].X = V_Sub(Point_v_IN[1].X, Point_v_IN[0].X);
	Point_v[2].X = V_Sub(Point_v_IN[2].X, Point_v_IN[0].X);

	Point_t.X = V_Sub(Point_t_IN.X, Point_v_IN[0].X);
	Normal_S_vec = Othognal_Vector(Point_v, Point_t);

	for (int i = 0; i < 3; ++i) {
		R[i] = Distance(Point_v[i].X, Point_t.X);
		L[i] = Distance(Point_v[(i + 2) % 3].X, Point_v[(i + 1) % 3].X);
	}

	U_vec = V_Div(L[2], Point_v[1].X);
	V_vec = V_Cross(Normal_S_vec, U_vec);

	double W0 = V_Dot(V_Sub(Point_t.X, Point_v[0].X), Normal_S_vec);

	double U0 = V_Dot(V_Sub(Point_t.X, Point_v[0].X), U_vec);

	double V0 = V_Dot(V_Sub(Point_t.X, Point_v[0].X), V_vec);

	double U3 = V_Dot(V_Sub(Point_v_IN[2].X, Point_v_IN[0].X), V_Sub(Point_v_IN[1].X, Point_v_IN[0].X)) / L[2];

	double V3 = 2.0 * Trian_Area(Point_v[0].X, Point_v[1].X, Point_v[2].X) / L[2];

	S_M[0] = -((L[2] - U0) * (L[2] - U3) + V0 * V3) / L[0];
	S_M[1] = -(U3 * (U3 - U0) + V3 * (V3 - V0)) / L[1];
	S_M[2] = -U0;

	S_P = V_Add(S_M, L);

	T_0[0] = (V0 * (U3 - L[2]) + V3 * (L[2] - U0)) / L[0];
	T_0[1] = (U0 * V3 - V0 * U3) / L[1];
	T_0[2] = V0;

	for (int i = 0; i < 3; ++i) R_0[i] = std::sqrt(T_0[i] * T_0[i] + W0 * W0);

	for (int i = 0; i < 3; ++i)
	{
		I_0[i][0] = std::log(MAX(R[(i + 2) % 3] + S_P[i], LAMDA_E * 1.0E-6) /
			MAX(R[(i + 1) % 3] + S_M[i], LAMDA_E * 1.0E-6));
	}

	for (int i = 0; i < 3; ++i)
	{
		I_0[i][1] = (S_P[i] * R[(i + 2) % 3] - S_M[i] * R[(i + 1) % 3] + R_0[i] * R_0[i] * I_0[i][0]) / 2;
	}

	if (std::fabs(W0) < L[0] * 1.0e-3)
	{
		for (int i = 0; i < 3; ++i) K_1_N[1] += T_0[i] * I_0[i][0];
	}
	else {
		for (int i = 0; i < 3; ++i)
		{
			double a1 = std::atan2(T_0[i] * S_P[i], R_0[i] * R_0[i] + std::fabs(W0) * R[(i + 2) % 3]);
			double a2 = std::atan2(T_0[i] * S_M[i], R_0[i] * R_0[i] + std::fabs(W0) * R[(i + 1) % 3]);
			BELTA[i] = a1 - a2;
			if (std::isnan(BELTA[i])) BELTA[i] = 0.0;
		}
		K_1_N[0] = (BELTA[0] + BELTA[1] + BELTA[2]) / std::fabs(W0);
		for (int i = 0; i < 3; ++i) K_1_N[1] += T_0[i] * I_0[i][0];
		K_1_N[1] -=  W0 * W0 * K_1_N[0];
	}

	// 计算边法向量 M_vec
	for (int i = 0; i < 3; ++i) {

		PN[0] = Point_v[(i + 2) % 3];
		PN[1] = Point_v[(i + 1) % 3];

		PP = Point_v[i];
		M_vec[i] = Oth_S(PN, PP);
	}

	// I_M
	for (int j = 0; j < 2; ++j) {
		I_M[j] = { 0,0,0 };
		for (int i = 0; i < 3; ++i) {
			I_M[j] = V_Add(I_M[j], V_Mul(I_0[i][j], M_vec[i]));

		}
	}
	Vector ROLL = V_Sub(Point_t.X, V_Mul(W0, Normal_S_vec));

	// K_2_N
	for (int i = 0; i < 3; ++i) {
		K_2_N[i] = V_Add(V_Mul(K_1_N[1], (V_Sub(ROLL, Point_v[i].X))), I_M[1]);
	}

	for (int i = 0; i < 3; ++i)
	{
		r_q[i] = V_Sub(Point_t.X, Point_v[i].X);
	}

	K_3_N = V_Sub(V_Mul(-W0 * K_1_N[0], Normal_S_vec), I_M[0]);
	
	for (int i = 0; i < 3; ++i) 
	{
		K_4_N[i] = V_Cross(r_q[i], K_3_N);
	}

	result.P_BUF = K_1_N[1];
	result.DE_P_BUF = K_3_N;

	for (int i = 0; i < 3; ++i) 
	{
		result.L_BUF[i] = K_2_N[i];
		result.CS_BUF[i] = K_4_N[i];
	}

	return result;
}

std::vector<double> F_K_1_N(std::array<Point, 3>& Point_v_IN, Point& Point_t_IN, int N_in) 
{
	int N = (N_in + 1) / 2;

	std::array<Point, 3>Point_v{};
	Point Point_t, P_t;
	std::array<double, 3> R, L, T_0, S_M, S_P, R_0;
	std::vector<std::vector<double>> I_0(3, std::vector<double>(N + 2, 0.0));

	std::array<double, 3> BELTA{};
	std::vector<double> K_1_N(N + 2, 0.0);
	std::vector<double> BUFF(N + 1, 0.0);

	Point_v[0].X = { 0.0,0.0,0.0 };
	Point_v[1].X = V_Sub(Point_v_IN[1].X, Point_v_IN[0].X);
	Point_v[2].X = V_Sub(Point_v_IN[2].X, Point_v_IN[0].X);
	Point_t.X = V_Sub(Point_t_IN.X, Point_v_IN[0].X);

	Vector Normal_S_vec = Othognal_Vector(Point_v,Point_t);

	if (V_Dot(Normal_S_vec, Point_t.X) < -1e-30)
	{
		Normal_S_vec = V_Mul(-1.0, Normal_S_vec);
		Point a = Point_v[2];
		Point_v[2] = Point_v[1];
		Point_v[1] = a;
	}

	for (int i = 0; i < 3; ++i) 
	{
		R[i] = Distance(Point_v[i].X, Point_t.X);
		L[i] = Distance(Point_v[(i + 2) % 3].X, Point_v[(i + 1) % 3].X);
	}

	Vector U_vec = V_Div(L[2], Point_v[1].X);
	Vector V_vec = V_Cross(Normal_S_vec, U_vec);

	double W0 = V_Dot(V_Sub(Point_t.X, Point_v[0].X), Normal_S_vec);

	double U0 = V_Dot(V_Sub(Point_t.X, Point_v[0].X), U_vec);

	double V0 = V_Dot(V_Sub(Point_t.X, Point_v[0].X), V_vec);

	double U3 = V_Dot(V_Sub(Point_v_IN[2].X, Point_v_IN[0].X), V_Sub(Point_v_IN[1].X, Point_v_IN[0].X)) / L[2];

	double V3 = 2.0 * Trian_Area(Point_v[0].X, Point_v[1].X, Point_v[2].X) / L[2];

	S_M[0] = -((L[2] - U0) * (L[2] - U3) + V0 * V3) / L[0];
	S_M[1] = -(U3 * (U3 - U0) + V3 * (V3 - V0)) / L[1];
	S_M[2] = -U0;

	S_P[0] = S_M[0] + L[0];
	S_P[1] = S_M[1] + L[1];
	S_P[2] = S_M[2] + L[2];

	T_0[0] = (V0 * (U3 - L[2]) + V3 * (L[2] - U0)) / L[0];
	T_0[1] = (U0 * V3 - V0 * U3) / L[1];
	T_0[2] = V0;

	for (int i = 0; i < 3; ++i) R_0[i] = std::sqrt(T_0[i] * T_0[i] + W0 * W0);

	for (int i = 0; i < 3; ++i) 
	{
		I_0[i][0] = std::log(MAX(R[(i + 2) % 3] + S_P[i], LAMDA_E * 1.0e-6) /
			MAX(R[(i + 1) % 3] + S_M[i], LAMDA_E * 1.0e-6));
	}

	for (int j = 1; j <= N; ++j) 
	{
		for (int i = 0; i < 3; ++i) 
		{
			double term1 = S_P[i] * std::pow(R[(i + 2) % 3], j * 2 - 1);
			double term2 = S_M[i] * std::pow(R[(i + 1) % 3], j * 2 - 1);
			I_0[i][j] = (term1 - term2 + (j * 2 - 1) * R_0[i] * R_0[i] * I_0[i][j - 1]) / (j * 2);
		}
	}

	if (std::fabs(W0) < L[0] * 1.0e-3) 
	{
		for (int j = 0; j <= N; ++j) 
		{
			double sum = 0.0;
			for (int i = 0; i < 3; ++i) sum += T_0[i] * I_0[i][j];
			K_1_N[j + 1] = sum / (j * 2 + 1);
		}
	}
	else {
		for (int i = 0; i < 3; ++i) 
		{
			double a1 = std::atan2(T_0[i] * S_P[i], R_0[i] * R_0[i] + std::fabs(W0) * R[(i + 2) % 3]);
			/*std::cout << a1 << std::endl;*/
			double a2 = std::atan2(T_0[i] * S_M[i], R_0[i] * R_0[i] + std::fabs(W0) * R[(i + 1) % 3]);
			BELTA[i] = a1 - a2;
			if (std::isnan(BELTA[i])) BELTA[i] = 0.0;
		}
		K_1_N[0] = (BELTA[0] + BELTA[1] + BELTA[2]) / std::fabs(W0);
		for (int j = 0; j <= N; ++j) 
		{
			double sum = 0.0;
			for (int i = 0; i < 3; ++i) sum += T_0[i] * I_0[i][j];
			K_1_N[j + 1] = (sum + (j * 2 - 1) * (W0 * W0) * K_1_N[j]) / (j * 2 + 1);
		}
	}

	for (int j = 0; j <= N; ++j) BUFF[j] = K_1_N[j + 1];
	return BUFF;
}

std::vector<std::array<double, 3>> F_K_2_N(std::array<Point, 3>& Point_v_IN, Point& Point_t_IN,	int Point_U_IN,	int N_in) 
{
	int N = (N_in + 1) / 2;

	std::array<Point, 3> Point_v;
	std::array<Point, 2> PN;
	Point Point_t, Point_U, P_t, PP;
	std::array<double, 3> R, L, S_P, R_0, S_M, T_0, BELTA;

	std::array<Vector, 3> M_vec{};

	std::vector<std::vector<double>> I_0(3, std::vector<double>(N + 2, 0.0));
	std::vector<double> K_1_N(N + 2, 0.0);
	std::vector<std::array<double, 3>> K_2_N(N + 1, { 0.0, 0.0, 0.0 });
	std::vector<std::array<double, 3>> I_M(N + 2, { 0.0, 0.0, 0.0 });
	std::vector<std::array<double, 3>> BUFF(N + 1, { 0.0, 0.0, 0.0 });

	// 相对坐标
	Point_v[0].X = { {0.0, 0.0, 0.0} };
	Point_v[1].X = V_Sub(Point_v_IN[1].X, Point_v_IN[0].X);
	/*std::cout << Point_v[1].V[0]<<std::endl;*/
	Point_v[2].X = V_Sub(Point_v_IN[2].X, Point_v_IN[0].X);
	Point_t.X = V_Sub(Point_t_IN.X, Point_v_IN[0].X);
	Point_U.X = V_Sub(Point_v_IN[Point_U_IN].X, Point_v_IN[0].X);

	Vector Normal_S_vec = Othognal_Vector(Point_v, Point_t);
	if (V_Dot(Normal_S_vec, Point_t.X) > 0.0)
	{
		Normal_S_vec = V_Mul(-1.0, Normal_S_vec);
		P_t = Point_v[2];
		Point_v[2] = Point_v[1];
		Point_v[1] = P_t;
	}

	// 计算 R 和 L
	for (int i = 0; i < 3; ++i) {
		R[i] = Distance(Point_v[i].X, Point_t.X);
		L[i] = Distance(Point_v[(i + 2) % 3].X, Point_v[(i + 1) % 3].X);
	}
	/*std::cout << R[1] << std::endl;
	std::cout << L[1] << std::endl;*/

	// 基向量 U, V
	Vector U_vec = V_Div(L[2], Point_v[1].X);
	Vector V_vec = V_Cross(Normal_S_vec, U_vec);

	double W0 = V_Dot(V_Sub(Point_t.X, Point_v[0].X), Normal_S_vec);
	double U0 = V_Dot(V_Sub(Point_t.X, Point_v[0].X), U_vec);
	double V0 = V_Dot(V_Sub(Point_t.X, Point_v[0].X), V_vec);

	double U3 = V_Dot(V_Sub(Point_v[2].X, Point_v[0].X), V_Sub(Point_v[1].X, Point_v[0].X)) / L[2];
	double V3 = 2.0 * Trian_Area(Point_v[0].X, Point_v[1].X, Point_v[2].X) / L[2];

	/*std::cout << U3 << std::endl;*/

	S_M[0] = -((L[2] - U0) * (L[2] - U3) + V0 * V3) / L[0];
	S_M[1] = -(U3 * (U3 - U0) + V3 * (V3 - V0)) / L[1];
	S_M[2] = -U0;

	/*std::cout << S_M[2] << std::endl;*/

	S_P = V_Add(S_M, L);

	/*std::cout << S_P[1] << std::endl;*/

	T_0[0] = (V0 * (U3 - L[2]) + V3 * (L[2] - U0)) / L[0];
	T_0[1] = (U0 * V3 - V0 * U3) / L[1];
	T_0[2] = V0;

	/*std::cout << T_0[0] << std::endl;*/

	for (int i = 0; i < 3; ++i) R_0[i] = std::sqrt(T_0[i] * T_0[i] + W0 * W0);

	// I_0 初始化
	for (int i = 0; i < 3; ++i) {
		I_0[i][0] = std::log(MAX(R[(i + 2) % 3] + S_P[i], LAMDA_E * 1.0e-6) / MAX(R[(i + 1) % 3] + S_M[i], LAMDA_E * 1.0e-6));
	}

	for (int j = 1; j <= N + 1; ++j) {
		for (int i = 0; i < 3; ++i) {
			double term1 = S_P[i] * std::pow(R[(i + 2) % 3], j * 2 - 1);
			double term2 = S_M[i] * std::pow(R[(i + 1) % 3], j * 2 - 1);
			I_0[i][j] = (term1 - term2 + (j * 2 - 1) * R_0[i] * R_0[i] * I_0[i][j - 1]) / (j * 2);
		}
	}

	//std::cout << I_0[2][2];

	// K_1_N
	if (std::fabs(W0) < L[0] * 1.0e-3) {
		for (int j = 1; j <= N + 1; ++j) {
			double sum = 0.0;
			for (int i = 0; i < 3; ++i) sum += T_0[i] * I_0[i][j - 1];
			K_1_N[j] = sum / (j * 2 - 1);
		}
	}
	else {
		for (int i = 0; i < 3; ++i) {
			double num = T_0[i] * S_P[i];
			double den = R_0[i] * R_0[i] + std::fabs(W0) * R[(i + 2) % 3];
			double a1 = std::atan2(num, den);
			num = T_0[i] * S_M[i];
			den = R_0[i] * R_0[i] + std::fabs(W0) * R[(i + 1) % 3];
			double a2 = std::atan2(num, den);
			BELTA[i] = a1 - a2;
			if (std::isnan(BELTA[i])) BELTA[i] = 0.0;
		}
		K_1_N[0] = (BELTA[0] + BELTA[1] + BELTA[2]) / std::fabs(W0);
		for (int j = 0; j <= N; ++j) {
			double sum = 0.0;
			for (int i = 0; i < 3; ++i) sum += T_0[i] * I_0[i][j];
			K_1_N[j + 1] = (sum + (j * 2 - 1) * (W0 * W0) * K_1_N[j]) / (j * 2 + 1);
		}
	}

	// 计算边法向量 M_vec
	for (int i = 0; i < 3; ++i) {

		PN[0] = Point_v[(i + 2) % 3];
		PN[1] = Point_v[(i + 1) % 3];

		PP = Point_v[i];
		M_vec[i] = Oth_S(PN, PP);
	}

	// I_M
	for (int j = 0; j <= N + 1; ++j) {
		for (int i = 0; i < 3; ++i) {
			for (int k = 0; k < 3; ++k)
				I_M[j][k] += M_vec[i][k] * I_0[i][j];
		}
	}

	Vector ROLL = V_Sub(Point_t.X, V_Mul(W0, Normal_S_vec));

	// K_2_N
	for (int j = 0; j <= N; ++j) {
		for (int k = 0; k < 3; ++k) {
			K_2_N[j][k] = K_1_N[j + 1] * (ROLL[k] - Point_U.X[k]) + I_M[j + 1][k] / (j * 2 + 1);
		}
		BUFF[j] = K_2_N[j];
	}

	return BUFF;
}






