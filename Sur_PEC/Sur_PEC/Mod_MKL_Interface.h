#ifndef MOD_MKL_INTERFACE_H
#define MOD_MKL_INTERFACE_H

#pragma once

#include <vector>
#include <complex>
#include <string>
#include <stdexcept>
#include <iostream>
#include <type_traits> 
#include "mkl.h"     


using cplx64 = std::complex<double>;


int Eigenvalue_Decomposition(
    const std::vector<std::vector<double>>& Matrix_Input,
    std::vector<std::complex<double>>& Eigen_Value,
    std::vector<std::vector<double>>& V_L,
    std::vector<std::vector<double>>& V_R
);


// ===================================================================
//                 通用线性方程求解 (模板接口)
// ===================================================================
template<typename T>
int Solve_M(
    std::vector<std::vector<T>>& mat_A,
    std::vector<std::vector<T>>& mat_B,
    int mat_type)
{
    // --- 维度检查与推断 ---
    if (mat_A.empty() || mat_A.size() != mat_A[0].size()) {
        throw std::invalid_argument("Solve_M: 系数矩阵 A 必须为非空方阵。");
    }
    if (mat_B.empty() || mat_A.size() != mat_B.size()) {
        throw std::invalid_argument("Solve_M: 矩阵 A 和 B 的行数必须相等。");
    }
    lapack_int N = mat_A.size();
    lapack_int NRHS = mat_B[0].size();
    if (NRHS == 0) return 0; // 没有方程组需要求解

    // --- 将二维输入“扁平化”为一维 ---
    std::vector<T> a_1d(N * N);
    std::vector<T> b_1d(N * NRHS);
    for (lapack_int i = 0; i < N; ++i) for (lapack_int j = 0; j < N; ++j) a_1d[i * N + j] = mat_A[i][j];
    for (lapack_int i = 0; i < N; ++i) for (lapack_int j = 0; j < NRHS; ++j) b_1d[i * NRHS + j] = mat_B[i][j];

    // --- 调用 LAPACKE ---
    lapack_int info = 0;
    std::vector<lapack_int> ipiv(N);

    // 根据类型选择正确的 LAPACKE 例程
    if constexpr (std::is_same_v<T, double>) {
        double* a_ptr = a_1d.data();
        double* b_ptr = b_1d.data();
        if (mat_type == 0) { // 通用矩阵
            info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, N, N, a_ptr, N, ipiv.data());
            if (info == 0) info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', N, NRHS, a_ptr, N, ipiv.data(), b_ptr, NRHS);
        }
        else if (mat_type == 1) { // 对称矩阵
            info = LAPACKE_dsytrf(LAPACK_ROW_MAJOR, 'U', N, a_ptr, N, ipiv.data());
            if (info == 0) info = LAPACKE_dsytrs(LAPACK_ROW_MAJOR, 'U', N, NRHS, a_ptr, N, ipiv.data(), b_ptr, NRHS);
        }
        else {
            throw std::invalid_argument("Solve_M: 无效的 mat_type。");
        }
    }
    else if constexpr (std::is_same_v<T, cplx64>) {
        auto* a_ptr = reinterpret_cast<MKL_Complex16*>(a_1d.data());
        auto* b_ptr = reinterpret_cast<MKL_Complex16*>(b_1d.data());
        if (mat_type == 0) { // 通用矩阵
            info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, N, N, a_ptr, N, ipiv.data());
            if (info == 0) info = LAPACKE_zgetrs(LAPACK_ROW_MAJOR, 'N', N, NRHS, a_ptr, N, ipiv.data(), b_ptr, NRHS);
        }
        else if (mat_type == 1) { // 对称矩阵
            info = LAPACKE_zsytrf(LAPACK_ROW_MAJOR, 'U', N, a_ptr, N, ipiv.data());
            if (info == 0) info = LAPACKE_zsytrs(LAPACK_ROW_MAJOR, 'U', N, NRHS, a_ptr, N, ipiv.data(), b_ptr, NRHS);
        }
        else {
            throw std::invalid_argument("Solve_M: 无效的 mat_type。");
        }
    }
    else {
        static_assert(std::is_same_v<T, double> || std::is_same_v<T, cplx64>, "Solve_M only supports double and std::complex<double>.");
    }

    // --- 如果成功，将一维结果写回二维矩阵 B ---
    if (info == 0) {
        // 将解矩阵 X 写回 mat_B
        for (lapack_int i = 0; i < N; ++i) {
            for (lapack_int j = 0; j < NRHS; ++j) {
                mat_B[i][j] = b_1d[i * NRHS + j];
            }
        }
        // (可选) 如果需要，也可以将分解后的 A 写回 mat_A
        for (lapack_int i = 0; i < N; ++i) {
            for (lapack_int j = 0; j < N; ++j) {
                mat_A[i][j] = a_1d[i * N + j];
            }
        }
    }

    return info;
}


// ===================================================================
//                 通用矩阵乘法 (模板接口)
// ===================================================================
template<typename T>
std::vector<std::vector<T>> Product_M(
    const std::vector<std::vector<T>>& a,
    const std::vector<std::vector<T>>& b
) {
    // --- 维度检查 ---
    if (a.empty() || a[0].empty() || b.empty() || b[0].empty() || a[0].size() != b.size()) {
        throw std::invalid_argument("Product_M: 矩阵维度不合法或不匹配。");
    }
    int m = a.size();
    int k = a[0].size();
    int n = b[0].size();

    // --- 扁平化 ---
    std::vector<T> a_1d(m * k);
    std::vector<T> b_1d(k * n);
    for (int i = 0; i < m; ++i) for (int j = 0; j < k; ++j) a_1d[i * k + j] = a[i][j];
    for (int i = 0; i < k; ++i) for (int j = 0; j < n; ++j) b_1d[i * n + j] = b[i][j];

    std::vector<T> c_1d(m * n);

    // --- 根据类型调用对应的 MKL 函数 (编译时选择) ---
    if constexpr (std::is_same<T, double>::value) {
        const double alpha = 1.0;
        const double beta = 0.0;
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            m, n, k,
            alpha, a_1d.data(), k,
            b_1d.data(), n,
            beta, c_1d.data(), n);
    }
    else if constexpr (std::is_same<T, std::complex<double>>::value) {
        const std::complex<double> alpha = { 1.0, 0.0 };
        const std::complex<double> beta = { 0.0, 0.0 };
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            m, n, k,
            &alpha, reinterpret_cast<const MKL_Complex16*>(a_1d.data()), k,
            reinterpret_cast<const MKL_Complex16*>(b_1d.data()), n,
            &beta, reinterpret_cast<MKL_Complex16*>(c_1d.data()), n);
    }
    else {
        throw std::invalid_argument("Product_M: 不支持的矩阵元素类型。");
    }

    mkl_free_buffers();

    // --- 转换回二维矩阵 ---
    std::vector<std::vector<T>> c(m, std::vector<T>(n));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            c[i][j] = c_1d[i * n + j];
        }
    }
    return c;
} //


// ===================================================================
//                 通用矩阵求逆 (模板接口)
// ===================================================================
template<typename T>
int Inverse_M(
    std::vector<std::vector<T>>& A,
    int method_key)
{
    if (A.empty() || A.size() != A[0].size()) {
        throw std::invalid_argument("Inverse_M: 矩阵必须为方阵。");
    }
    lapack_int N = A.size();

    // 1. 将二维输入扁平化为一维
    std::vector<T> A_1d(N * N);
    for (lapack_int i = 0; i < N; ++i) {
        for (lapack_int j = 0; j < N; ++j) {
            A_1d[i * N + j] = A[i][j];
        }
    }
    lapack_int lda = N;
    lapack_int info = 0;
    std::vector<lapack_int> ipiv(N);

    // 2. 根据类型调用对应的 LAPACKE 函数 (编译时选择)
    if constexpr (std::is_same_v<T, double>) {
        if (method_key == 0) { // 一般矩阵
            info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, N, N, A_1d.data(), lda, ipiv.data());
            if (info == 0) info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, N, A_1d.data(), lda, ipiv.data());
        }
        else if (method_key == 1) { // 对称矩阵
            char uplo = 'U';
            info = LAPACKE_dsytrf(LAPACK_ROW_MAJOR, uplo, N, A_1d.data(), lda, ipiv.data());
            if (info == 0) info = LAPACKE_dsytri(LAPACK_ROW_MAJOR, uplo, N, A_1d.data(), lda, ipiv.data());
        }
        else {
            throw std::invalid_argument("Inverse_M: 无效的 method_key。");
        }
    }
    else if constexpr (std::is_same_v<T, cplx64>) {
        auto* a_ptr = reinterpret_cast<MKL_Complex16*>(A_1d.data());
        if (method_key == 0) { // 一般矩阵
            info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, N, N, a_ptr, lda, ipiv.data());
            if (info == 0) info = LAPACKE_zgetri(LAPACK_ROW_MAJOR, N, a_ptr, lda, ipiv.data());
        }
        else if (method_key == 1) { // 对称矩阵
            char uplo = 'U';
            info = LAPACKE_zsytrf(LAPACK_ROW_MAJOR, uplo, N, a_ptr, lda, ipiv.data());
            if (info == 0) info = LAPACKE_zsytri(LAPACK_ROW_MAJOR, uplo, N, a_ptr, lda, ipiv.data());
        }
        else {
            throw std::invalid_argument("Inverse_M: 无效的 method_key。");
        }
    }
    else {
        static_assert(std::is_same_v<T, double> || std::is_same_v<T, cplx64>, "Inverse_M only supports double and std::complex<double>.");
    }

    // 3. 如果成功，将一维结果写回二维矩阵
    if (info == 0) {
        for (lapack_int i = 0; i < N; ++i) {
            for (lapack_int j = 0; j < N; ++j) {
                A[i][j] = A_1d[i * N + j];
            }
        }
        // 对于对称情况，补全下三角
        if (method_key == 1) {
            for (int i = 0; i < N; ++i) {
                for (int j = i + 1; j < N; ++j) {
                    A[j][i] = A[i][j];
                }
            }
        }
    }

    return info;
}//

#endif // 