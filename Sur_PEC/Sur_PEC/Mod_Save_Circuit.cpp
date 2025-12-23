#include "Mod_Save_Circuit.h"

void Save_PEEC_Circuit()
{
    int N_C = PEC_N;
    int N_L = C_PEC_N;

    // 构造路径
    std::string Cir_Save_Path = PATH + "output/CirNetList.txt";

    // --- 处理 CC 矩阵（real(PP_00)） ---
    std::vector<std::vector<double>> CC(N_C, std::vector<double>(N_C, 0.0));
    for (int i = 0; i < N_C; ++i)
        for (int j = 0; j < N_C; ++j)
            CC[i][j] = PP[i][j];


    // --- 写出到文件 ---
    std::ofstream ofs(Cir_Save_Path);
    if (!ofs.is_open()) {
        // 无法打开文件 — 在项目中应有统一错误处理，这里简单返回
        return;
    }

    // 写 N_PORT
    ofs << "N_PORT " << N_PORT << "\n";

    // 写每个 port
    for (int idx = 0; idx < N_PORT; ++idx) {
        std::ostringstream ss_nodes;
        ss_nodes << std::setw(8) << PORT_DATA[idx].N[0]
            << std::setw(8) << PORT_DATA[idx].N[1];
        ofs << "PO_" << (idx + 1) << "  " << ss_nodes.str() << "\n";
    }

    int N_ELEMENT = static_cast<int>((N_C * (N_C + 1) * 0.5) + N_L * N_L);
    ofs << "N_ELEMENT " << N_ELEMENT << "\n";

    // 写 LL 对角与 LM（I!=J）
    for (int i = 0; i < N_L; ++i) {
        // LL_i_i
        std::ostringstream ss;
        ss << std::setw(10) << CN[i][0]
            << std::setw(10) << CN[i][1]
            << std::scientific << std::uppercase << std::setprecision(10) << std::setw(18)
            << std::real(LL[i][i]);
        ofs << "LL_" << (i + 1) << "_" << (i + 1) << "  " << ss.str() << "\n";

        for (int j = 0; j < N_L; ++j) {
            if (i == j) continue;
            std::ostringstream ss2;
            ss2 << std::setw(10) << CN[i][0]
                << std::setw(10) << CN[i][1]
                << std::setw(10) << CN[j][0]
                << std::setw(10) << CN[j][1]
                << std::scientific << std::uppercase << std::setprecision(10) << std::setw(18)
                << std::real(LL[i][j]);
            ofs << "LM_" << (i + 1) << "_" << (j + 1) << "  " << ss2.str() << "\n";
        }
    }

    // CC 矩阵
    for (int i = 0; i < N_C; ++i) {
        for (int j = i; j < N_C; ++j) {
            std::ostringstream ss;
            if (i != j) {
                ss << std::setw(10) << (i + 1) << std::setw(10) << (j + 1)
                    << std::scientific << std::uppercase << std::setprecision(10) << std::setw(18)
                    << -CC[i][j];
            }
            else {
                double row_sum = 0.0;
                for (int k = 0; k < N_C; ++k) {
                    row_sum += CC[i][k];
                }
                ss << std::setw(10) << (i + 1) << std::setw(10) << 0
                    << std::scientific << std::uppercase << std::setprecision(10) << std::setw(18)
                    << row_sum;
            }
            ofs << "CC_" << (i + 1) << "_" << (j + 1) << "  " << ss.str() << "\n";
        }
    }

    ofs.close();
}

static void parse_indices(const std::string& name, int& a, int& b) {
    // name = LL_12_12 or LM_3_8
    size_t p1 = name.find('_');
    size_t p2 = name.find('_', p1 + 1);
    a = std::stoi(name.substr(p1 + 1, p2 - p1 - 1));
    b = std::stoi(name.substr(p2 + 1));
}

int Save_Netlist() {

    struct Inductor {
        int n1;
        int n2;
        double L;
    };

    struct Capacitor {
        int n1;
        int n2;
        double C;
    };

    struct Mutual {
        int i;
        int j;
        double M;
    };

    std::string inputFile = PATH + "output/CirNetList.txt";
    std::string outputFile = PATH + "output/CirNetList.sp";

    //const double K_CUTOFF = 5e-4;
    const double R_REF = 1e15;
    const double R_LLOSS = 1e-9;

    std::ifstream fin(inputFile);
    if (!fin) {
        std::cerr << "ERROR: cannot open " << inputFile << std::endl;
        return 1;
    }

    std::map<int, Inductor> inductors;
    std::vector<Mutual> mutuals;
    std::vector<Capacitor> capacitors;
    std::set<std::pair<int, int>> used_pairs;
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '*') continue;

        std::istringstream iss(line);
        std::string name;
        iss >> name;

        if (name.rfind("LL_", 0) == 0) {
            int idx, dummy;
            parse_indices(name, idx, dummy);
            Inductor ind;
            iss >> ind.n1 >> ind.n2 >> ind.L;
            inductors[idx] = ind;
        }
        else if (name.rfind("LM_", 0) == 0) {
            Mutual m;
            parse_indices(name, m.i, m.j);
            int tmp;
            iss >> tmp >> tmp >> tmp >> tmp >> m.M;
            mutuals.push_back(m);
        }
        else if (name.rfind("CC_", 0) == 0) {
            Capacitor c;
            iss >> c.n1 >> c.n2 >> c.C;
            capacitors.push_back(c);
        }
    }
    fin.close();

    std::ofstream fout(outputFile);
    if (!fout) {
        std::cerr << "ERROR: cannot create " << outputFile << std::endl;
        return 1;
    }

    fout << "*=========================================================\n";
    fout << "* Stable PEEC netlist for PrimeSim HSPICE\n";
    fout << "* Includes DC reference + inductor loss\n";
    fout << "*=========================================================\n\n";

    fout << ".option post=2 accurate measform=3\n";
    fout << ".option s\n\n";

    fout << "*---------------------------\n";
    fout << "* Port definitions\n";
    fout << "*---------------------------\n";
    for (int i = 0; i < N_PORT; ++i) {
        fout << "P"<<i+1<<" " << PORT_DATA[i].N[0] << " " << PORT_DATA[i].N[1] << " PORT="<<i+1<<" Z0=" << PORT_DATA[i].Zin << "\n";
    }

    fout << "*---------------------------\n";
    fout << "* DC reference resistors\n";
    fout << "*---------------------------\n";
    for (int i = 0; i < N_PORT; ++i) {
        fout << "RREF" << i + 1 << " " << PORT_DATA[i].N[0] << " 0 " << R_REF << "\n";
    }

    fout << std::scientific << std::setprecision(6);

    fout << "*---------------------------\n";
    fout << "* Self inductances + series loss\n";
    fout << "*---------------------------\n";
    for (const auto& kv : inductors) {
        int id = kv.first;
        const auto& ind = kv.second;

        std::string mid = "nL_" + std::to_string(id);

        fout << "RL" << id << " "
            << ind.n1 + 1 << " "
            << mid << " "
            << R_LLOSS << "\n";

        fout << "L" << id << " "
            << mid << " "
            << ind.n2 + 1 << " "
            << ind.L << "\n";
    }

    fout << "\n*---------------------------\n";
    fout << "* Mutual couplings\n";
    fout << "*---------------------------\n";

    for (const auto& m : mutuals) {

        int i = m.i;
        int j = m.j;

        // 统一顺序，确保 i < j
        if (i == j) continue;
        if (i > j) std::swap(i, j);

        // 去重：已经生成过就跳过
        if (used_pairs.count({ i, j })) continue;
        used_pairs.insert({ i, j });

        double Li = inductors[i].L;
        double Lj = inductors[j].L;
        double kij = m.M / std::sqrt(Li * Lj);

        fout << "K" << i << "_" << j
            << " L" << i
            << " L" << j
            << " " << kij << "\n";
    }

    fout << "\n*---------------------------\n";
    fout << "* Capacitances (CC)\n";
    fout << "*---------------------------\n";
    int cid = 1;
    for (const auto& c : capacitors) {
        fout << "C" << cid++
            << " " << c.n1
            << " " << c.n2
            << " " << c.C << "\n";
    }

    fout << "\n*---------------------------\n";
    fout << "* AC / S-parameter analysis\n";
    fout << "*---------------------------\n";
    fout << ".AC LIN "<< N_FP <<" "<<FS<<" "<<FE<<"\n";
    fout << ".LIN sparcalc=1 filename=CirNetLis format=touchstone dataformat=db\n";
    fout << ".PRINT LIN S11 S21 S12 S22\n";
    fout << ".end\n";

    fout.close();
    std::cout << "Stable HSPICE netlist generated: "
        << outputFile << std::endl;
    return 0;
}