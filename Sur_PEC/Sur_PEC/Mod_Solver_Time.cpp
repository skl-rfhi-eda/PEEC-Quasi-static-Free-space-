#include"Mod_Solver_Time.h"

int S_BL[3];
int L_BL[3];
int N_TPT;
int Source_Sz = 1;
double T_STEP;
std::vector<std::vector<double>> PP_AET;
std::vector<std::vector<double>> VI, D_VI, VI_TEMP, D_VI_TEMP;
std::vector<std::vector<double>> A_M, B_M;
std::vector<std::vector<std::vector<double>>> x_t;
std::vector<std::vector<std::complex<double>>> y_t;

using namespace std;
using Matrix = std::vector<std::vector<double>>;

bool invertMatrix(const Matrix& A, Matrix& inv)
{
    int n = static_cast<int>(A.size());
    if (n == 0) return false;
    if (A[0].size() != static_cast<size_t>(n)) return false;

    // 构造 [A | I] 增广矩阵
    Matrix aug(n, std::vector<double>(2 * n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j)
            aug[i][j] = A[i][j];
        aug[i][i + n] = 1.0;               // 单位矩阵
    }

    // 高斯-约当消元（带部分主元）
    for (int p = 0; p < n; ++p) {
        // 找第 p 列绝对值最大的主元
        int pivot = p;
        for (int i = p + 1; i < n; ++i)
            if (std::abs(aug[i][p]) > std::abs(aug[pivot][p]))
                pivot = i;

        if (std::abs(aug[pivot][p]) < 1e-12)   // 奇异
            return false;

        // 交换行
        std::swap(aug[p], aug[pivot]);

        // 归一化主元行
        double div = aug[p][p];
        for (int j = 0; j < 2 * n; ++j)
            aug[p][j] /= div;

        // 消去其它行
        for (int i = 0; i < n; ++i) {
            if (i == p) continue;
            double factor = aug[i][p];
            for (int j = 0; j < 2 * n; ++j)
                aug[i][j] -= factor * aug[p][j];
        }
    }

    // 提取右半部分为逆矩阵
    inv.assign(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            inv[i][j] = aug[i][j + n];

    return true;
}

// 读取矩阵文件为二维 vector
vector<vector<double>> readMatrix(const string& filename) {
    ifstream fin(filename);
    vector<vector<double>> matrix;
    string line;

    if (!fin.is_open()) {
        cerr << "无法打开文件：" << filename << endl;
        exit(1);
    }

    while (getline(fin, line)) {
        istringstream iss(line);
        vector<double> row;
        double value;
        while (iss >> value)
            row.push_back(value);
        if (!row.empty())
            matrix.push_back(row);
    }

    return matrix;
}

// 比较两个矩阵是否相同（允许一定浮点误差）
bool compareMatrix(const vector<vector<double>>& A,
    const vector<vector<double>>& B,
    double rel_tol = 1e-7) {
    if (A.size() != B.size())
        return false;
    for (size_t i = 0; i < A.size(); ++i) {
        if (A[i].size() != B[i].size())
            return false;
        for (size_t j = 0; j < A[i].size(); ++j) {
            double a = A[i][j];
            double b = B[i][j];

            // 处理两个数都为零的情况
            if (a == 0.0 && b == 0.0)
                continue;

            // 计算相对误差
            double scale = fmax(fabs(a), fabs(b));
            double rel_error = fabs(a - b) / scale;

            if (rel_error > rel_tol)
                return false;
        }
    }
    return true;
}

void saveMatrixToTxt(const std::string& filename, const std::vector<std::vector<double>>& matrix)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return;
    }

    file << std::fixed << std::setprecision(15); // 设置小数点精度
    for (const auto& row : matrix) {
        for (size_t j = 0; j < row.size(); ++j) {
            file << row[j];
            if (j < row.size() - 1) file << " ";
        }
        file << "\n";
    }

    file.close();
    std::cout << "矩阵已保存到文件: " << filename << std::endl;
}


void Ini_Element() {
    std::cout << "Ini_Element..." << std::endl;

    NODE_N = N_S - DS_N; //有效节点数
    L_BL[0] = NODE_N;
    L_BL[1] = C_PEC_N; //PEC导体块
    L_BL[2] = N_PORT; //端口块

    S_BL[0] = 0;
    S_BL[1] = S_BL[0] + L_BL[0];
    S_BL[2] = S_BL[1] + L_BL[1];
    /*cout << S_BL[0] << "," << S_BL[1] << "," << S_BL[2] << endl;*/

    M_SZ = S_BL[2] + L_BL[2]; // 总未知维数
    

    //分配A_E、A_ET
    A_E.assign(C_PEC_N + N_PORT, std::vector<double>(NODE_N, 0.0));
    A_ET.assign(NODE_N, std::vector<double>(C_PEC_N + N_PORT, 0.0));

    //用CN构造A_E、A_ET
    for (int i = 0; i < C_PEC_N; ++i) {
        for (int j = 0; j < 2; ++j) {
            int idx = CN[i][j] - DS_N;
            if (idx < NODE_N && idx >= 0) {
                A_E[i][idx] = std::pow(-1.0, j + 2);
                A_ET[idx][i] = std::pow(-1.0, j + 1);
            }
        }
    }

    for (int i = 0; i < N_PORT; ++i) {
        int idx1 = PORT_DATA[i].N[0] - DS_N;
        int idx2 = PORT_DATA[i].N[1] - DS_N;
        if (idx1 < NODE_N && idx1 >= 0) {
            A_E[i + C_PEC_N][idx1] = 1.0;
            A_ET[idx1][i + C_PEC_N] = -1.0;
        }
        if (idx2 < NODE_N && idx2 >= 0) {
            A_E[i + C_PEC_N][idx2] = -1.0;
            A_ET[idx2][i + C_PEC_N] = 1.0;
        }
    }

    /*saveMatrixToTxt("C:\\Users\\Administrator\\Desktop\\Old_ver\\A_E.txt", A_E);
    saveMatrixToTxt("C:\\Users\\Administrator\\Desktop\\Old_ver\\A_ET.txt", A_ET);*/

    //分配VI，D_VI等
    VI.assign(M_SZ, std::vector<double>(Source_Sz, 0.0));
    D_VI.assign(M_SZ, std::vector<double>(Source_Sz, 0.0));
    VI_TEMP.assign(M_SZ, std::vector<double>(Source_Sz, 0.0));
    D_VI_TEMP.assign(M_SZ, std::vector<double>(Source_Sz, 0.0));

};

void Ini_Source() {
    double RF_T, RAISE_S, RAISE_E;
    int R_S, R_E;
    double Amplitude = 1.0;

    N_TPT = 60; //时间步数
    double T_E = 5e-11; //总时间
    T_STEP = T_E / N_TPT; //时间步长

    x_t.assign(N_TPT + 1,
        std::vector<std::vector<double>>(Source_Sz,
            std::vector<double>(Source_Sz, 0.0)));

    y_t.assign(N_TPT + 1, std::vector<std::complex<double>>(N_PORT, { 0.0,0.0 }));
    //上升时间窗口+保持
    RF_T = 1.2e-11;
    RAISE_S = 0;
    RAISE_E = RF_T + RAISE_S;

    R_S = static_cast<int>(RAISE_S / T_STEP);
    R_E = R_S + static_cast<int>(RF_T / T_STEP) + 1;

    for (int i = R_S + 1; i <= R_E; ++i) {
        double frac = static_cast<double>(i - R_S) / (R_E - R_S);
        double window = (1.0 - std::cos(PI * frac)) * 0.5 * Amplitude;
        for (int a = 0; a < Source_Sz; ++a) {
            for (int b = 0; b < Source_Sz; ++b) {
                x_t[i][a][b] = window;
            }
        }
    }

    for (int i = R_E + 1; i <= N_TPT; ++i) {
        for (int a = 0; a < Source_Sz; ++a) {
            for (int b = 0; b < Source_Sz; ++b) {
                x_t[i][a][b] = Amplitude;
            }
        }
    }


};

void Build_PP() {

    //初始化PP矩阵
    PP.assign(NODE_N, std::vector<double>(NODE_N, 0.0));

    //对角线元素设为1
    for (int i = PEC_N; i < NODE_N; ++i) {
        PP[i][i] = 1.0;
    }

    if (DS_N != 0) {

    }
    else {
        //计算电容矩阵
        for (int i = 0; i < PEC_N; ++i) {
            for (int j = 0; j < PEC_N; ++j) {
                PP[i][j] = ((1.0 / N_CAP[i][j].E[1] - 1.0 / N_CAP[i][j].E[0]) * N_CAP[i][j].VAL) / E0;
            }
        }
    }
    /*saveMatrixToTxt("C:\\Users\\Administrator\\Desktop\\Old_ver\\PP.txt", PP);*/

    // 计算PP_AET = PP * A_ET
    PP_AET.assign(NODE_N, std::vector<double>(C_PEC_N + N_PORT, 0.0));
    PP_AET = Product_M(PP, A_ET);


};

void Build_LL() {
    LL.assign(C_PEC_N, std::vector<double>(C_PEC_N, 0.0));

    if (DS_N != 0) {

    }
    else {
        for (int i = 0; i < C_PEC_N; ++i) {
            for (int j = 0; j < C_PEC_N; ++j) {
                LL[i][j] = U0 * B_IND[i][j].VAL;
            }
        }
    }
    /*saveMatrixToTxt("C:\\Users\\Administrator\\Desktop\\Old_ver\\LL.txt", LL);*/
};

void Build_Circuit() {

    Build_PP();

    Build_LL();

};



void Build_CP_M() {
    // 分配临时矩阵
    std::vector<std::vector<double>> E_TEMP(M_SZ, std::vector<double>(M_SZ, 0.0));
    std::vector<std::vector<double>> A_TEMP(M_SZ, std::vector<double>(M_SZ, 0.0));
    std::vector<std::vector<double>> B_TEMP(M_SZ, std::vector<double>(Source_Sz, 0.0));

    /*cout << S_BL[0] << "," << S_BL[1] << "," << S_BL[2] << endl;*/
    // 构建E_TEMP矩阵
    for (int i = 0; i < PEC_N; ++i) {
        E_TEMP[i][i] = 1.0;
    }

    // 添加电感矩阵
    for (int i = 0; i < C_PEC_N; ++i) {
        for (int j = 0; j < C_PEC_N; ++j) {
            E_TEMP[i + S_BL[1]][j + S_BL[1]] = LL[i][j];
        }
    }
    /*saveMatrixToTxt("C:\\Users\\Administrator\\Desktop\\Old_ver\\E_TEMP.txt", E_TEMP);*/

    // 构建A_TEMP矩阵
    for (int i = 0; i < NODE_N; ++i) {
        for (int j = 0; j < C_PEC_N + N_PORT; ++j) {
            A_TEMP[i][j + S_BL[1]] = PP_AET[i][j];
            A_TEMP[j + S_BL[1]][i] = A_E[j][i];
        }
    }

    // 添加对角项
    for (int i = 0; i < C_PEC_N; ++i) {
        A_TEMP[i + S_BL[1]][i + S_BL[1]] = -1e-10;
    }

    // 添加端口阻抗
    for (int i = 0; i < N_PORT; ++i) {
        A_TEMP[i + S_BL[2]][i + S_BL[2]] = -PORT_DATA[i].Zin;
    }

    /*saveMatrixToTxt("C:\\Users\\Administrator\\Desktop\\Old_ver\\A_TEMP.txt", A_TEMP);*/
    // 构建源项
    B_TEMP[M_SZ - 1][Source_Sz - 1] = -1.0;

    // 计算E_HA = E_TEMP - 0.5 * T_STEP * A_TEMP

    std::vector<std::vector<double>> E_HA(M_SZ, std::vector<double>(M_SZ, 0.0));
    for (int i = 0; i < M_SZ; ++i) {
        for (int j = 0; j < M_SZ; ++j) {
            E_HA[i][j] = E_TEMP[i][j] - 0.5 * T_STEP * A_TEMP[i][j];
        }
    }
    /*saveMatrixToTxt("C:\\Users\\Administrator\\Desktop\\Old_ver\\E_HA.txt", E_HA);*/

    // 计算逆矩阵
    bool a = Inverse_M(E_HA, 0);
    /*bool a = invertMatrix(E_HA, E_HA);*/

    /*saveMatrixToTxt("C:\\Users\\Administrator\\Desktop\\Old_ver\\E_HAI.txt", E_HA);*/

    // 计算A_M = E_HA^(-1) * A_TEMP
    A_M.assign(M_SZ, std::vector<double>(M_SZ, 0.0));
    A_M = Product_M(E_HA, A_TEMP);

    // 计算B_M = E_HA^(-1) * B_TEMP
    B_M.assign(M_SZ, std::vector<double>(Source_Sz, 0.0));
    B_M = Product_M(E_HA, B_TEMP);
    /*saveMatrixToTxt("C:\\Users\\Administrator\\Desktop\\Old_ver\\B_M.txt", B_M);

    saveMatrixToTxt("C:\\Users\\Administrator\\Desktop\\Old_ver\\A_M.txt", A_M);*/

    
};

void Update_VI(const std::vector<std::vector<double>>& S_in) {

    std::vector<std::vector<double>> A_M_BUF(M_SZ, std::vector<double>(Source_Sz, 0.0));
    std::vector<std::vector<double>> B_M_BUF(M_SZ, std::vector<double>(Source_Sz, 0.0));
    // 保存当前值
    D_VI_TEMP = D_VI;
    VI_TEMP = VI;

    // 预测步
    for (int i = 0; i < M_SZ; ++i) {
        for (int j = 0; j < Source_Sz; ++j) {
            VI[i][j] += 0.5 * T_STEP * D_VI[i][j];
        }
    }

    // 计算A_M * VI
    for (int i = 0; i < M_SZ; ++i) {
        for (int j = 0; j < Source_Sz; ++j) {
            for (int k = 0; k < M_SZ; ++k) {
                A_M_BUF[i][j] += A_M[i][k] * VI[k][j];
            }
        }
    }

    // 计算B_M * S_in
    for (int i = 0; i < M_SZ; ++i) {
        for (int j = 0; j < Source_Sz; ++j) {
            for (int k = 0; k < Source_Sz; ++k) {
                B_M_BUF[i][j] += B_M[i][k] * S_in[k][j];
            }
        }
    }

    // 更新D_VI
    for (int i = 0; i < M_SZ; ++i) {
        for (int j = 0; j < Source_Sz; ++j) {
            D_VI[i][j] = A_M_BUF[i][j] + B_M_BUF[i][j];
        }
    }

    // 校正步
    for (int i = 0; i < M_SZ; ++i) {
        for (int j = 0; j < Source_Sz; ++j) {
            VI[i][j] = VI_TEMP[i][j] + 0.5 * T_STEP * (D_VI_TEMP[i][j] + D_VI[i][j]);
        }
    }
};

void Time_Solver() {
    std::cout << "Start solving..." << std::endl;

    Ini_Element();
    Ini_Source();
    Build_Circuit();
    Build_CP_M();


    //时间步循环
    std::ofstream outFile(PATH + "output\\time.txt");
    for (int i = 1; i <= N_TPT; ++i) {
        Update_VI(x_t[i]);

        // 输出结果
        outFile << i * T_STEP << ", "
            << x_t[i][0][0] << ", "
            << VI[PORT_DATA[0].N[0] - DS_N ][0] - VI[PORT_DATA[0].N[1] - DS_N][0] << ", "
            << VI[PORT_DATA[1].N[0] - DS_N ][0] - VI[PORT_DATA[1].N[1] - DS_N ][0] << std::endl;
    }
    outFile.close();

    std::cout << "Finish solving." << std::endl;
}