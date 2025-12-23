#include "Mod_MKL_Interface.h" 


// ===================================================================
//特征值分解 (Eigenvalue Decomposition)
// ===================================================================
int Eigenvalue_Decomposition(
    const std::vector<std::vector<double>>& Matrix_Input,
    std::vector<std::complex<double>>& Eigen_Value,
    std::vector<std::vector<double>>& V_L,
    std::vector<std::vector<double>>& V_R)
{
    int N;
	N = static_cast<int>(Matrix_Input.size());
    // 检查输入参数合法性
    if (N <= 0) {
        throw std::invalid_argument("Eigen_Decomposition: 矩阵维度N必须为正整数");
    }
    if (Matrix_Input.size() != static_cast<size_t>(N) || Matrix_Input[0].size() != static_cast<size_t>(N)) {
        throw std::invalid_argument("Eigen_Decomposition: 输入矩阵维度与N不匹配");
    }

    // 设置MKL参数
    char jobl = 'V';  //控制是否计算左特征向量
    char jobr = 'V';  //控制是否计算右特征向量
    int lda = N;  //输入矩阵的阶数
    int ldvl = N;
    int ldvr = N;
    int info = 0;  //函数的执行状态

    // 分配内存
    std::vector<double> A_col_major(static_cast<size_t>(N) * N);
    std::vector<double> WR(N), WI(N);
    std::vector<double> VL_col_major(static_cast<size_t>(N) * N);
    std::vector<double> VR_col_major(static_cast<size_t>(N) * N);

    // 将C++的行主序矩阵转换为Fortran的列主序
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A_col_major[j * N + i] = Matrix_Input[i][j];
        }
    }

    // 查询最佳工作空间大小
    int lwork = -1;
    double work_query;
    dgeev_(&jobl, &jobr, &N, A_col_major.data(), &lda,
        WR.data(), WI.data(), VL_col_major.data(), &ldvl,
        VR_col_major.data(), &ldvr, &work_query, &lwork, &info);

    // 分配实际工作空间并执行特征值分解
    lwork = static_cast<int>(work_query);
    std::vector<double> work(lwork);
    dgeev_(&jobl, &jobr, &N, A_col_major.data(), &lda,
        WR.data(), WI.data(), VL_col_major.data(), &ldvl,
        VR_col_major.data(), &ldvr, work.data(), &lwork, &info);

    if (info != 0) {

        std::cout << "警告: Eigen_Decomposition 的 dgeev_ 执行失败，错误码=" << info << std::endl;
    }

    // 整理特征值结果
    Eigen_Value.resize(N);
    for (int i = 0; i < N; ++i) {
        Eigen_Value[i] = std::complex<double>(WR[i], WI[i]);
    }

    // 将列主序的特征向量转换回行主序的 vector<vector>
    V_L.assign(N, std::vector<double>(N));
    V_R.assign(N, std::vector<double>(N));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            V_L[i][j] = VL_col_major[j * N + i];
            V_R[i][j] = VR_col_major[j * N + i];
        }
    }
    return info;
}

