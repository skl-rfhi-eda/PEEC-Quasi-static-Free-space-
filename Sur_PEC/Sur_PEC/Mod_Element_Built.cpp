#include "Mod_Element_Built.h"


// 构建矩阵元素的主函数
void Element_Built() {
    std::cout << "P3: Building matrix element ...\n";

    double TIME_S = Get_Time();

    Ini_Element_Sz();

    double ms_ele_mem = static_cast<double>(N_S) * N_S * sizeof(Ele_Mesh_Int) / 1e6;
    double b_ind_mem = static_cast<double>(CN_N) * C_PEC_N * sizeof(Ele_BL) / 1e6;
    double n_cap_mem = static_cast<double>(N_S) * N_S * sizeof(Ele_P) / 1e6;
    double b_cap_mem = static_cast<double>(CN_N) * CN_N * sizeof(Ele_BC) / 1e6;
    double b_cs_mem = static_cast<double>(CN_N) * sizeof(Ele_CS) / 1e6;

    std::cout << "     Memory usage: MS: " << static_cast<int>(ms_ele_mem) << " MB.\n";
    std::cout << "                 : LL: " << std::setw(5) << static_cast<int>(b_ind_mem) << " MB,  CC:  "
        << std::setw(5) << static_cast<int>(b_cap_mem) << " MB,  PP: " << std::setw(5) << static_cast<int>(n_cap_mem)
        << " MB,    CS " << std::setw(5) << static_cast<int>(b_cs_mem) << " MB.\n";

    // 核心计算流程
    Build_MS_Ele();
    Re_Construct_Ele();

    // MS_ELE.clear(); // 如果不需要，可以清空

    double TIME_E = Get_Time();

    // 计算并打印耗时
    Time_Diff(TIME_S, TIME_E);
}


// 分配内存
void Ini_Element_Sz() {

    B_IND.resize(CN_N, std::vector<Ele_BL>(C_PEC_N));
    N_CAP.resize(N_S, std::vector<Ele_P>(N_S));
    MS_ELE.resize(N_S, std::vector<Ele_Mesh_Int>(N_S));

}


void Build_MS_Ele() {
    std::cout << "\n    Build Mesh Integral...\n";

    int COUNT = 0;
    int TOTAL =  PEC_N * PEC_N;

    for (int i = 0; i < PEC_N; ++i) {
        for (int j = 0; j < PEC_N; ++j) {
            COUNT++;
            MS_ELE[i][j] = Ele_Cal(i, j);
            PRT(COUNT, TOTAL, 1); 
        }
    }
}


void Re_Construct_Ele() {
    std::cout << "    Construct Potential Coefficient...\n";
    for (int i = 0; i < PEC_N; ++i) {
        for (int j = 0; j < PEC_N; ++j) {
            N_CAP[i][j].VAL = MS_ELE[i][j].P;
            N_CAP[i][j].E[0] = DIELECTRIC[Sur_M[j].M_TYP[0]].real();
            N_CAP[i][j].E[1] = DIELECTRIC[Sur_M[j].M_TYP[1]].real();

            if (Sur_M[j].IS_PEC == 1) {
                N_CAP[i][j].E[0] *= 2.0;
                N_CAP[i][j].E[1] *= 2.0;
                N_CAP[i][j].E[0] = -N_CAP[i][j].E[0];
            }
        }
    }

    std::cout << "              Inductive Coupling...\n";
    for (int i = 0; i < CN_N; ++i) {
        for (int j = 0; j < C_PEC_N; ++j) {
            double TEMP = 0.0;

            for (int ii = 0; ii < 2; ++ii) {
                for (int jj = 0; jj < 2; ++jj) {

                    int cn_p_i = CN_P[i][ii];
                    int cn_p_j = CN_P[j][jj];

                    if (cn_p_i != -1 || cn_p_j != -1) {
                        int PM = ((ii + jj ) % 2 == 0) ? 1 : -1;
                        
                        TEMP += PM * MS_ELE[CN[i][ii]][CN[j][jj]].L[cn_p_i][cn_p_j];
                    }
                }
            }
            B_IND[i][j].VAL = TEMP;

        }
    }
}



//void Built_Branch_L() {
//    std::cout << "    Branch L calculating...\n";
//    for (int i = 0; i < CN_N; ++i) {
//        for (int j = 0; j < C_PEC_N; ++j) {
//            B_IND[i][j] = B_L_Cal(i, j);
//
//        }
//        PRT(i + 1, CN_N, 1);
//    }
//
//    // --- 检查矩阵性质 ---
//    for (int i = 0; i < C_PEC_N; ++i) {
//        for (int j = i + 1; j < C_PEC_N; ++j) {
//            if (B_IND[i][j].VAL * B_IND[j][i].VAL > B_IND[i][i].VAL * B_IND[j][j].VAL) {
//                std::cout << "Warning: Matrix property violated at (" << i << "," << j << ")\n";
//                std::cout << B_IND[i][j].VAL << "," << B_IND[j][i].VAL << std::endl;
//                std::cout << B_IND[i][i].VAL << "," << B_IND[j][j].VAL << std::endl;
//               
//            }
//        }
//    }
//}


//void Built_Node_P() {
//    std::cout << "    Node PP calculating...\n";
//    for (int i = 0; i < DS_N + PEC_N; ++i) {
//        for (int j = 0; j < DS_N + PEC_N; ++j) {
//            N_CAP[i][j] = N_PP_Cal(i, j);
//        }
//        PRT(i + 1, DS_N + PEC_N, 1);
//    }
//}

