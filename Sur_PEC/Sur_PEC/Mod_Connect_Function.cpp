#include "Mod_Connect_Function.h"

	
void Test_Node() {
    std::cout << "Mesh Data Testing " << N_S << std::endl;
    
    auto printMesh = [](int index, const Mesh& m) {
        printf("NUM %4d   PT %4d%4d%4d  PEC, %4d    M_TYP %4d%4d    P_TYP %4d%4d   Node Num %4d\n",
            index, m.P[0], m.P[1], m.P[2], m.IS_PEC,
            m.M_TYP[0], m.M_TYP[1], m.P_TYP[0], m.P_TYP[1], m.NODE_NUMBER);
        };

    std::cout << "DS Surface " << DS_N << std::endl;
    for (int i = 0; i < DS_N; i += 10) {
        printMesh(i, Sur_M[i]);
    }

    std::cout << "PEC Surface " << PEC_N << std::endl;
    for (int i = 0; i < PEC_N; i += 10) {
        printMesh(i, Sur_M[i]);
    }

    std::cout << "Line mesh Surface " << L_N << std::endl;
    int start_index = PEC_N;  // 从PEC_N开始（0-based）
    for (int i = 0; i < L_N; ++i) {  // 逐个处理，步长为1
        printMesh(start_index + i, Sur_M[PEC_N + i]);
    }

    system("pause");
}  //  Test_Point 

// 打印连接信息
void Test_Branch() {
    std::cout << "Connection Testing " << CN_N << std::endl;
    auto printConnection = [](int index, int node1, int node2, int vert1, int vert2, int nodeNum1, int nodeNum2) {
        printf("NUM %4d   NODE %4d%4d  VERT, %4d%4d    NODE NUMBER %4d%4d\n",
            index, node1, node2, vert1, vert2, nodeNum1, nodeNum2);
		};

    std::cout << "PEC Connection " << C_PEC_N << std::endl;
    for (int i = 0; i < C_PEC_N; i += 40) {
        int idx = i * 2;  // CN是2列数组
        printConnection(i + 1, CN[i][0], CN[i][1], CN_P[i][0], CN_P[i][1],
            Sur_M[CN[i][0]].NODE_NUMBER, Sur_M[CN[i][1]].NODE_NUMBER);
    }

    std::cout << "DS Connection " << C_DS_N << std::endl;
    int index = 0;
    for (int i = C_PEC_N; i < CN_N; i += 40) {
        index++;
        printConnection(i + 1, CN[i][0], CN[i][1], CN_P[i][0], CN_P[i][1],
            Sur_M[CN[i][0]].NODE_NUMBER, Sur_M[CN[i][1]].NODE_NUMBER);
    }

    system("pause");
}

// 将面网格转换为线网格
/**
* @brief 从一个三角形网格中提取一条边，并将其作为线段网格返回。
* @param ms_in 输入的三角形网格对象。
* @param adj   指定要提取哪条边。在C++中，adj=0表示提取顶点P[0]对面的边(P[1]-P[2])，
*              adj=1表示提取P[1]对面的边(P[2]-P[0])，adj=2表示提取P[2]对面的边(P[0]-P[1])。
* @return Mesh 返回一个代表线段的Mesh对象。其中P[0]和P[1]存储线段的两个顶点索引，
*/
Mesh Trans_L(const Mesh& ms_in, int adj) {
    //创建一个新的 Mesh 对象用于存储结果
    Mesh line_mesh;

    // 计算顶点索引 
    int idx1 = ms_in.P[(adj + 1) % 3];
    int idx2 = ms_in.P[(adj + 2) % 3];

    // 排序顶点
    line_mesh.P[0] = std::min(idx1, idx2);
    line_mesh.P[1] = std::max(idx1, idx2);
    line_mesh.P[2] = -1;
        

    // 复制类型信息
    line_mesh.M_TYP = ms_in.M_TYP;
    line_mesh.P_TYP = ms_in.P_TYP;
    line_mesh.IS_PEC = ms_in.IS_PEC;

    // 计算长度(存储在AREA字段)
    line_mesh.AREA = Distance(PT[line_mesh.P[0]].X, PT[line_mesh.P[1]].X);

    return line_mesh;
}

/**
* @brief 将一个临时网格转换为PEC（理想电导体）表面三角形。
* @param ms_in 输入的临时网格，函数将使用其前三个顶点。
* @return result_mesh 返回一个配置好的、代表PEC表面的Mesh对象。
*/
Mesh Trans_MS_PEC(const Temp_Mesh& ms_in) {
    void Set_Up_Mesh(Mesh & ms_in);
	Vector Oth_V(const Point * cp, const Point & up); // 函数声明
    Mesh result_mesh;  //创建返回结果对象
    Point infi_p = { -1e3, -1e3, -1e3 };   //义一个“无穷远点”作为法向量计算的参考点

    //复制顶点索引
    for (int i = 0; i < 3; ++i) {
        result_mesh.P[i] = ms_in.P[i];
    }
    

    // 设置类型信息
    result_mesh.M_TYP[0] = ms_in.M_TYP;
    result_mesh.M_TYP[1] = ms_in.M_TYP;
    result_mesh.P_TYP[0] = ms_in.P_TYP;
    result_mesh.P_TYP[1] = 0;

    result_mesh.IS_PEC = 1;  // 标记为PEC边界

    // 计算法向量
    Point cp[3] = { PT[result_mesh.P[0]], PT[result_mesh.P[1]], PT[result_mesh.P[2]] };
    result_mesh.NOR_S = Oth_V(cp, infi_p);

    // 设置其他网格属性
    Set_Up_Mesh(result_mesh);

    return result_mesh;
}

/**
* @brief 从一个四面体网格中提取一个三角形面。
* @param ms_in 输入的四面体网格 (Temp_Mesh)。
* @param un_p  要排除的顶点的索引（0-based），剩下的三个点组成一个面。
* @return Mesh 返回一个代表四面体一个面的Mesh对象。
*/
Mesh Trans_MS_TETRA(const Temp_Mesh& ms_in, int un_p) {
	void Set_Up_Mesh(Mesh & ms_in);
    Vector Oth_V(const Point * cp, const Point & up);  // 函数声明
    // 创建返回结果对象
    Mesh result_mesh;
    int loc = 0;

    // 复制除un_p外的顶点
    for (int i = 0; i < 4; ++i) {
        if (i != un_p ) {  
            result_mesh.P[loc] = ms_in.P[i];
            loc++;
        }
    }

    // 设置类型信息
    result_mesh.M_TYP[0] = ms_in.M_TYP;
    result_mesh.M_TYP[1] = 1;
    result_mesh.P_TYP[0] = ms_in.P_TYP;
    result_mesh.P_TYP[1] = 0;

    result_mesh.IS_PEC = 0;  // 非PEC边界

    // 计算法向量
    Point cp[3] = { PT[result_mesh.P[0]], PT[result_mesh.P[1]], PT[result_mesh.P[2]] };
    result_mesh.NOR_S = Oth_V(cp, PT[ms_in.P[un_p ]]);  // 

    // 设置其他网格属性
    Set_Up_Mesh(result_mesh);

    return result_mesh;
}

// 判断两个表面是否相同(共享顶点数量)
int Same_Sur(const Mesh& ms_a, const Mesh& ms_b) {
    int count = 0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (ms_a.P[i] == ms_b.P[j]) {
                count++;
            }
        }
    }
    return count;
}

// 计算两个表面的邻接关系
void Adj_Sur(const Mesh& ms_a, const Mesh& ms_b, int& adj_a, int& adj_b, int& cm_p) {
    cm_p = 0;
    adj_a = 6;
    adj_b = 6;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (ms_a.P[i] == ms_b.P[j]) {
                cm_p++;
                adj_a /= (i + 1);   //
                adj_b /= (j + 1);   //
            }
        }
    }
    adj_a--;
	adj_b--;
}

/**
* 功能：计算三角形表面与线表面的公共顶点数量及邻接编码,三角形-线表面邻接判断
*  @ms_a      输入：三角形表面（3个顶点）
*  @ms_b      输入：线表面（2个顶点，使用P[0]和P[1]）
*  @adj_a     输出：三角形表面的邻接编码（基于0基索引计算）
*  @adj_b     输出：线表面的邻接编码
*  @cm_p      输出：公共顶点的数量
*/
void Adj_Sur_L(const Mesh& ms_a, const Mesh& ms_b, int& adj_a, int& adj_b, int& cm_p) {
    cm_p = 0;
    adj_a = 6;
    adj_b = 6;   //6?2

    // 如果PEC属性不同，则不相邻
    if (ms_a.IS_PEC != ms_b.IS_PEC) {
        return;
    }

    for (int i = 0; i < 3; ++i) {   //遍历三角形ms_a的所有顶点 (i = 0, 1, 2)
        for (int j = 0; j < 2; ++j) {   //遍历线段ms_b的两个顶点 (j = 0, 1)
            if (ms_a.P[i] == ms_b.P[j]) {
                cm_p++;
                adj_a /= (i + 1);  // 避免除以0
                adj_b /= (j + 1);
            }
        }
    }
}

// 设置网格的基本属性
void Set_Up_Mesh(Mesh& ms_in) {
    // 获取三个顶点的引用
    Point& p1 = PT[ms_in.P[0]];
    Point& p2 = PT[ms_in.P[1]];
    Point& p3 = PT[ms_in.P[2]];
        
    // 计算重心坐标
    for (int i = 0; i < 3; ++i) {
        ms_in.MID_CO[i] = (p1.X[i] +p2.X[i] +p3.X[i]) / 3.0;
    }

    // 计算面积
    ms_in.AREA = Trian_Area(p1.X,p2.X,p3.X);

    // 计算边信息
    for (int j = 0; j < 3; ++j) {
            
        std::array<Point, 2> cp = {
        PT[ms_in.P[(j + 1) % 3]],
        PT[ms_in.P[(j + 2) % 3]]
        };
        //Point& cp1 = PT[ms_in.P[(j + 1) % 3]]; // 边的第一个顶点
        //Point& cp2 = PT[ms_in.P[(j + 2) % 3]]; // 边的第二个顶点
        //std::array<Point, 2> cp = { cp1, cp2 };
        Point up = PT[ms_in.P[j]];   //// 当前顶点索引

        // 计算边的法向分量和长度
        ms_in.L_N[j] = Oth_S(cp, up);   //Oth_S传入的参数为非常量引用，Point cp[2] 和 std::array<Point, 2>不等价
        ms_in.L[j] = Distance(cp[0].X, cp[1].X);

        // 计算边界中点坐标
        for (int i = 0; i < 3; ++i) {
            ms_in.B_MID_CO[j][i] = (PT[ms_in.P[0]].X[i] +
                PT[ms_in.P[1]].X[i] +
                PT[ms_in.P[2]].X[i] -
                PT[ms_in.P[j]].X[i] / 3.0) * 3.0 / 8.0;
        }
    }
}

/**
* @brief 判断两个四面体网格是否相邻
* @param ms_a      输入/输出：第一个四面体网格
* @param ms_b      输入/输出：第二个四面体网格
* @param err       输出：返回状态码
*                  - 0: 不相邻或材料类型不同
*                  - 1: 相邻且成功标记
*                  - -1: 相邻但已被标记过
*/
void Judge_Tetra_C(Temp_Mesh& ms_a, Temp_Mesh& ms_b, int& err) {
    err = 0;

    // 如果类型不同，则不相邻
    if (ms_a.M_TYP != ms_b.M_TYP) {
        err = 0;
        return;
    }

    int temp = 0;
    int va = 24;
    int vb = 24;

    // 比较两个四面体的所有顶点
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (ms_a.P[i] == ms_b.P[j]) {
                temp++;
                va /= (i + 1);  // 
                vb /= (j + 1);
            }
        }
    }

    // 如果有3个共同顶点，则标记为相邻
    if (temp == 3) {
        // 检查并标记第一个四面体的邻接信息
        if (ms_a.ADJ_NODE[va-1] != 0) {  
            err = -1;
        }
        else {
            ms_a.ADJ_NODE[va-1] = 1;
            err = 1;
        }
        // 检查并标记第二个四面体的邻接信息
        if (ms_b.ADJ_NODE[vb-1] != 0) {  
            err = -1;
        }
        else {
            ms_b.ADJ_NODE[vb-1] = 1;
            err = 1;
        }
    }
}

/**
    * @brief 计算一个与平面垂直且指向外部的法向量
    * @param cp    输入：平面上的3个点（cp[0], cp[1], cp[2]）
    * @param up    输入：平面外的一个参考点
    * @return      返回归一化的法向量
    */
Vector Oth_V(const Point* cp, const Point& up) {
    Vector result{};
    std::array<std::array<double, 2>,3> T;
        

    // 计算边向量
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
            T[j][i] = cp[i + 1].X[j] - cp[0].X[j];
        }
    }

    // 计算法向量
    for (int i = 0; i < 3; ++i) {
        int idx1 = (i % 3);
        int idx2 = ((i + 1) % 3);
        result[i] = T[idx1][0] * T[idx2][1] - T[idx2][0] * T[idx1][1];
    }
    /*  result.V[0] = T[1][0] * T[2][1] - T[2][0] * T[1][1];
    result.V[1] = T[2][0] * T[0][1] - T[0][0] * T[2][1];
    result.V[2] = T[0][1] * T[1][0] - T[1][1] * T[0][0];*/

    // 归一化
    double length = std::sqrt(result[0] * result[0] +
        result[1] * result[1] +
        result[2] * result[2]);

    //  计算临时变量temp,调整法向量方向
    double temp = 0.0;
    for (int i = 0; i < 3; ++i) {
        temp += result[i] * (up.X[i] - cp[0].X[i]);
    }

    double sign = -temp / (length * std::abs(temp));
    for (int i = 0; i < 3; ++i) {
        result[i] *= sign;
    }

    return result;
} //  Oth_V

/**
* @brief 在所有 PEC 表面中查找包含指定点的端口节点。
* @param P_in 输入的三维空间点，用于在 PEC 表面中进行定位。
*
* @return int 返回匹配的表面网格在全局数组 SUR_M 中的索引。

*/
int Port_Node( Point& p_in) {
    int port_node = -1;  // 初始化找到的端口节点索引为 0，表示尚未找到。

    for (int i = 0; i < PEC_N; ++i) {
        // 构建一个包含三个 Point 对象的数组，代表这个三角形表面。
        std::array<Point, 3> pts = {
            PT[Sur_M[i].P[0]],
            PT[Sur_M[i].P[1]],
            PT[Sur_M[i].P[2]]
        };
        // P_In_S 返回 1 表示“是”，返回其他值表示“否”
        if (P_In_S(p_in, pts) == 1) {
            if (port_node != -1) {
                std::cout << "Error: one port find two matched node!" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            port_node = i ;  // 找到匹配的表面，记录下它的索引 
        }
    }

    if (port_node == -1) {
        std::cout << "Error: one port do not find matched node!" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    return port_node;
}  // Port_Node




