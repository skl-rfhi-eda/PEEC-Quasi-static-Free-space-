#include"Mod_Connection.h"

void Creat_CNN();
void Creat_PEC_C();
void Creat_Mesh_Data();
void Trans_Sur_M();
void Sort_In_Mesh();
void Creat_PEC_Mesh();
void Creat_Port();


std::vector<Temp_Mesh> T_PEC, T_LMS, T_TETRA;
int N_PEC_T = 0, N_TETRA_T = 0, N_LMS_T = 0;

std::vector<Mesh> M_PEC, M_LMS, M_DS;
int PEC_T = 0, DS_T = 0, LMS_T = 0;

std::vector<std::vector<int>> CN_TEMP, CN_P_TEMP;

void Creat_Connection()
{
	double TIME_S, TIME_E;
	std::cout << "P2: Setting up connection.." << std::endl;
	TIME_S = Get_Time();

	T_PEC.resize(T_MS_N);
	T_LMS.resize(T_MS_N);
	T_TETRA.resize(T_MS_N);

	Sort_In_Mesh();

	M_PEC.resize(2 * N_PEC_T);
	M_LMS.resize(2 * N_LMS_T + N_PEC_T);
	M_DS.resize(4 * N_TETRA_T);

	Creat_Mesh_Data();
	Creat_CNN();
	Creat_Port();

	TIME_E = Get_Time();
	Time_Diff(TIME_S, TIME_E);
	std::cout << "TOTAL MESH:" << N_S << "  DS MESH:" << DS_N << "  PEC MESH:" << PEC_N << "  LINE MESH:" << L_N << std::endl;
	std::cout << "TOTAL C:" << CN_N << "  PEC C:" << C_PEC_N << "  DS C:" << C_DS_N << std::endl;
	std::cout << std::endl;
}

void Creat_CNN()
{
	CN_N = 0;
	CN_TEMP.assign(3 * N_S, std::vector<int>(2, 0));
	CN_P_TEMP.assign(3 * N_S, std::vector<int>(2, 0));

	Creat_PEC_C();

	CN.resize(CN_N, std::vector<int>(2, 0));
	CN_P.resize(CN_N, std::vector<int>(2, 0));

	for (int i = 0; i < CN_N; ++i)
	{
		CN[i] = CN_TEMP[i];
		CN_P[i] = CN_P_TEMP[i];
	}
}

void Creat_PEC_C()
{
	int ADJ_A, ADJ_B, COM_P;
	C_PEC_N = 0;
	for(int i = 0; i<PEC_N;++i)
		for (int j = PEC_N + 1; j < PEC_N + L_N; ++j)
		{
			Adj_Sur_L(Sur_M[i], Sur_M[j], ADJ_A, ADJ_B, COM_P);
			if (COM_P == 2)
			{
				C_PEC_N ++;
				CN_N ++;
				CN_TEMP[CN_N-1][0] = i;
				CN_TEMP[CN_N-1][1] = j;

				CN_P_TEMP[CN_N-1][0] = ADJ_A;
				CN_P_TEMP[CN_N-1][1] = 0;

				Sur_M[i].ADJ_NODE[ADJ_A] = j+1;
			}
		}
	for(int i=0;i<PEC_N;++i)
		for (int j = i + 1; j < PEC_N; j++)
		{
			ADJ_A = 6;
			ADJ_B = 6;
			COM_P = 0;

			Adj_Sur(Sur_M[i], Sur_M[j], ADJ_A, ADJ_B, COM_P);
			if(COM_P == 2)
				if (Sur_M[i].ADJ_NODE[ADJ_A] + Sur_M[j].ADJ_NODE[ADJ_B] == 0)
				{
					C_PEC_N++;
					CN_N++;
					CN_TEMP[CN_N-1][0] = i;
					CN_TEMP[CN_N-1][1] = j;

					CN_P_TEMP[CN_N-1][0] = ADJ_A;
					CN_P_TEMP[CN_N-1][1] = ADJ_B;

					Sur_M[i].ADJ_NODE[ADJ_A] = j;
					Sur_M[j].ADJ_NODE[ADJ_B] = i;
				}
		}
}

void Creat_Mesh_Data()
{

	Creat_PEC_Mesh();

	PEC_N = PEC_T;
	L_N = LMS_T;
	N_S = PEC_N + L_N;

	Sur_M.resize(N_S);
	Trans_Sur_M();
}

void Trans_Sur_M()
{
	int INDEX;
	std::array<Point, 2> P_TEMP;
	NODE_N = 0;

	for (int i = 0; i < PEC_T; ++i)
	{
		NODE_N++;
		M_PEC[i].NODE_NUMBER = NODE_N;
	}

	for (int i = 0; i < LMS_T; ++i)
	{
		NODE_N++;
		M_LMS[i].NODE_NUMBER = NODE_N;
	}

	INDEX = 0;

	for (int i = 0; i < PEC_T; ++i)
	{
		INDEX++;
		Sur_M[INDEX-1] = M_PEC[i];
		Sur_M[INDEX - 1].ADJ_NODE = {0};
	}

	for (int i = 0; i < LMS_T; ++i)
	{
		INDEX++;
		Sur_M[INDEX - 1] = M_LMS[i];
		Sur_M[INDEX - 1].ADJ_NODE = { 0 };
	}

	for (int i = 0; i < N_S; ++i)
		for (int j = 0; j < 3; ++j) {
			MAX_L = MAX(Sur_M[i].L[j], MAX_L);
		}
}

void Sort_In_Mesh()
{
	for (int i = 0; i < T_MS_N; ++i)
	{
		if (T_MS[i].TYP == 4)
		{
			N_TETRA_T++;
			T_TETRA[N_TETRA_T-1] = T_MS[i];
		}
		else if (T_MS[i].TYP == 2)
		{
			N_PEC_T++;
			T_PEC[N_PEC_T-1] = T_MS[i];
		}
		else if (T_MS[i].TYP == 1)
		{
			N_LMS_T++;
			T_LMS[N_LMS_T-1] = T_MS[i];
		}
	}
}

void Creat_PEC_Mesh()
{
	PEC_T = 0;
	for (int i = 0; i < N_PEC_T; i++)
	{
		PEC_T++;
		M_PEC[PEC_T-1] = Trans_MS_PEC(T_PEC[i]);

	}
}

bool L_Mesh_EXIST(int& PA, int& PB, int& IS_PEC)
{
	int P_L, P_S;
	P_S = MIN(PA, PB);
	P_L = MAX(PA, PB);
	for (int i = 0; i < LMS_T; ++i)
		if (P_S == M_LMS[i].P[0] && P_L == M_LMS[i].P[1])
			if (IS_PEC == M_LMS[i].IS_PEC)
				return 1;
	return 0;
}

void Creat_Port()
{
	std::array<Temp_Mesh, 10> Port_Mesh;
	N_PORT = 0;
	for(int i=0;i<T_MS_N;++i)
		if (T_MS[i].TYP == 1)
		{
			N_PORT++;
			Port_Mesh[N_PORT - 1] = T_MS[i];
		}
	PORT_DATA.resize(N_PORT);
	for (int i = 0; i < N_PORT; ++i)
	{
		PORT_DATA[i].N[0] = Port_Node(PT[Port_Mesh[i].P[0]]);
		PORT_DATA[i].N[1] = Port_Node(PT[Port_Mesh[i].P[1]]);
		PORT_DATA[i].Zin = T_MS[i].M_TYP;

	}
}