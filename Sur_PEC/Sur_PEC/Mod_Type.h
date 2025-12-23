#ifndef MOD_TYPE_H
#define MOD_TYPE_H

#pragma once

#include <array>
#include <vector>
#include <complex>

// 精度映射
using f64 = double;
using cplx64 = std::complex<double>;

using Vector = std::array<double, 3>;

struct Point {
    int number{ 0 };
    Vector X{ 0.0, 0.0, 0.0 };
};

struct Port {
    std::array<int, 2> N{ };
    double Zin{};
};

struct Temp_Mesh {
    // !P_TYP store the physical block number
    // !M_typ store the dielectric constant
    std::array<int, 8> P{ 0 };
    int TYP;
    int M_TYP;
    int P_TYP;
    std::array<int, 5> ADJ_NODE{ 0 };
};

struct Mesh {
    std::array<int, 2> P_TYP{ };
    std::array<int, 2> M_TYP{ };
    int IS_PEC{ 0 };
    std::array<int, 3> P{ };
    int NODE_NUMBER{ };

    std::array<double, 3> MID_CO{ };
    std::array<std::array<double, 3>, 3> B_MID_CO{ }; 
    double AREA{ };

    std::array<int, 3> ADJ_NODE{ };

    Vector NOR_S{ };
    std::array<Vector, 3> L_N{};
    std::array<double, 3> L{ };

    int TC_PEC{ };
    int TC_DIE{ };

    std::array<double, 2> ER{ 1.0, 1.0 };
};

struct Ele_Mesh_Int {
    std::array<std::array<double, 3>, 3> L{ };  
    std::array<std::array<double, 3>, 3> C{ };
    std::array<std::array<double, 3>, 3> CS{ };
    double P{ };
    double DE_P{ };
};

struct Ele_BL { //存储电感积分对应值
    double VAL{ };
    std::array<double, 2> U{ };
};

struct Ele_BC {
    double VAL{};
    std::array<double, 2> E{ };
};

struct Ele_P { //存储电容积分值
    double VAL{ };
    std::array<double, 2> E{ };
};

struct Ele_CS {
    double VAL{ };
};

struct W_NK {
    std::vector<std::vector<double>> W; 
};

struct ElementInfo {
    int type{2};  // 1:一维, 2:三角形, 3:矩形, 4:四面体
    int order{5}; // 积分阶数 ( 1 到 10等)
};




#endif // MOD_TYPE_H