#include"Mod_Global_Init.h"

//相关路径
std::string PATH;                 //路径
std::string INPUT_FILE;           //输入网格文件
std::string SET_FILE;	          //设置文件
std::string DIELECTRIC_FILE;      //介电常数文件							                         
std::string MAP_PATH;      //输出S参数的路径

std::array<int, 15> NUMBER_OF_TYPE;        // 每种类型的计数

// 点数据
std::vector<Point> PT;                     // 点集合
int N_P;                                   // 点的数量

// 临时网格数据
std::vector<Temp_Mesh> T_MS;               // 临时网格集合
int T_MS_N;                                // 临时网格的数量

// 端口数据
std::vector<Port> PORT_DATA;               // 端口数据集合
int N_PORT;                                // 端口数量

// 网格和积分线计数
int NODE_N;                                // 节点数量

//!DIE & PEC & LINE & PORT MESH
int N_S;                                   // 表面总数
int DS_N;                                  // DS表面数量
int PEC_N;                                 // PEC表面数量
int L_N;                                   // 线网格数量

std::vector<Mesh> Sur_M;                   // 表面网格集合

//!DIE & PEC& PORT CONNECTION
int CN_N;                                  
int C_DS_N;                                
int C_PEC_N;                               

// 连接关系数据
std::vector<int> DS_PEC_Branch;            // DS-PEC分支
std::vector<std::vector<int>> CN;         // CN索引：存储的是构成基函数的两个三角单元的全局索引
std::vector<std::vector<int>> CN_P;       // CN_P索引：存储的是与基函数相关的自由顶点在该三角形内的局部索引


// 边界和电容数据结构
std::vector<std::vector<Ele_BL>> B_IND;   //电感矩阵
std::vector<std::vector<Ele_BC>> B_CAP;    
std::vector<std::vector<Ele_P>> N_CAP;    //电容矩阵
std::vector<std::vector<Ele_CS>> B_CS;    
std::vector<std::vector<Ele_CS>> N_CS;     

// 积分相关数据
std::vector<std::vector<double>> PDD_N;    // 点-点双积分
std::vector<std::vector<double>> PD0_N;    // 点-0积分
std::vector<std::vector<double>> LD0_N;    // 线-0积分
std::vector<std::vector<Vector>> PV_DD;    // 点向量-点双积分
std::vector<std::vector<Vector>> PV_D0;    // 点向量-0积分
std::vector<std::vector<Vector>> LV_D0;    // 线向量-0积分

// 电阻矩阵
std::vector<std::vector<double>> B_RR;     // 边界电阻矩阵

// 单元网格积分数据
std::vector<std::vector<Ele_Mesh_Int>> MS_ELE;  // 网格单元积分

std::vector<std::vector<Ele_BC>> TC_CAP;   
std::vector<std::vector<Ele_CS>> TC_CS;    
std::vector<std::vector<Ele_P>> N_TC_CAP;  
int C_TC_N;                                
int TC_N;                                 

// 传输线映射关系
std::vector<std::vector<int>> TC_B2M;      
std::vector<std::vector<int>> TC_B2M_P;    
std::vector<int> TC_N2M;                   
std::vector<int> TC_M2N;                  

// 介电常数索引数组
std::array<int, 30> BLK_D;                 // 块介电常数索引

// 子问题和块计数
int N_SP;                                  // 子问题数量
int N_BLOCK;                               // 块数量

// 物理属性数组
std::array<std::array<int, P_Y>, PM_X> PHYZ;  // 物理属性数组
std::array<int, PM_X> DIE;                 // 介电常数数组

// S参数数据
std::vector<std::vector<std::vector<double>>> S_PAR_ABS;  // S参数幅度
std::vector<std::vector<std::vector<double>>> S_PAR_PHA;  // S参数相位

// 阻抗和导纳参数
std::vector<std::vector<std::vector<std::complex<double>>>> Z_PT;  // 阻抗参数
std::vector<std::vector<std::vector<std::complex<double>>>> Y_PT;  // 导纳参数

int INFI_G = 0;
double Z_G = 0.0;                        // 全局阻抗
int SLV_W = 0;                           // 求解器选择
int Swap_Log = 0;                        // 交换日志
long long FS, FE;                        // 起始\结束频率
int N_FP;                                // 频率点数
double PRECISION;                        // 计算精度
int DIM;                                 // 维度

std::vector<double> CF;
double W{ 0.0f };


double MAX_L{ 0.01 };
double LAMDA_E{ 1e5 };
ElementInfo ELM_TYPE;

std::vector<std::complex<double>> DIELECTRIC; // 介电常数数组
int Solver_SET{ 0 };

std::vector<std::vector<double>> LL, PP;
std::vector<std::vector<double>> A_E, A_ET;
int M_SZ;