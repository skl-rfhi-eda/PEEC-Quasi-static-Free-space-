# PEEC-Quasi-static-Free-space-
This is the first open source version of PEEC of our group. Its a quasi-static model for free-space electromagnetic problem. It gives equivalent circuit model conductors in a free-space. Gmsh is required as input (old version). .sp file and S-parameter simulation can be the output. The GUI are developping and will come soon.
## Project Overview

The Quasi-Static Free Space PEEC Electromagnetic Solver Engine is an electromagnetic simulation tool developed based on the Partial Element Equivalent Circuit (PEEC) method. It is suitable for simulating and analyzing the electromagnetic characteristics of electronic components such as integrated circuits, RF devices, metamaterials, and array antennas.

> Release Date: 2025
> Version: V1.0

----------

## 1. Functional Overview

This software completes the electromagnetic simulation process through the following six core modules:

1.  **Initialization**
    – Generates necessary parameter libraries, such as Gaussian integration databases, to support subsequent calculations.

2.  **Data Import**
    – Reads input information required for simulation, including the physical model of the electronic device, material parameters, and simulation settings.

3.  **Data Preprocessing**
    – Performs topological analysis on the physical model, extracts grid connectivity relationships, and generates information required for the equivalent circuit.

4.  **Equivalent Circuit Extraction**
    – Generates and outputs the partial element equivalent circuit model by calculating discretized integrals over the elements.

5.  **Circuit Solution and Analysis**
    – Solves for key electromagnetic characteristics in the circuit, such as unknown voltages, currents, and charge distributions, using the Modified Nodal Analysis method in the frequency domain.

6.  **Performance Evaluation**
    – Calculates target parameters including port scattering parameters (S-parameters), impedance parameters (Z-parameters), and port voltage/current distributions.

***
## 2. Model Preparation (Optional)

This software requires a mesh file with the `".msh"` extension, which can be generated using **GMSH version 2.2**.

GMSH is an open-source finite element mesh generation tool that can discretize complex geometric models into elements such as triangles and tetrahedra. It is widely used in various simulation computations. This software currently primarily supports mesh data based on triangular elements.

For use, rename the generated mesh file to **model.msh**. The following example uses a C-shaped resonator to illustrate how to generate the required model file (based on GMSH-2.2).

First, define the physical model using the following format:

***Example: C-shaped Resonator***

Define Points
```
Point(1) = {20, 0-20, 1.8, lc};
...
Point(26) = {10.5, 43-20, 1.8 , lc};
```
Define Lines
```
Line(2) = {1, 6};
...
Line(16) = {7, 6};
Line(23) = {17, 16};
Line(24) = {15, 18};

Line(51) = {13, 23};
...
Line(62) = {20, 22};
```
Define Line Loops
```
Line Loop(63) = {14, 8, 9, 10, 11, 12, -53, 54, 55, 56, -51};
Line Loop(65) = {15, -3, -16, -11};
Line Loop(67) = {2, 3, 4, 5, 6, 7, 57, 60, 61, 62, 59};
```
Define Surfaces
```
Plane Surface(64) = {63};
Plane Surface(66) = {65};
Plane Surface(68) = {67};
```
***
**It must be emphasized: To ensure correct simulation by this software, excitation lines must be added and their impedance specified via line properties. The endpoints of excitation lines should avoid coinciding with mesh vertices; they should be located on mesh boundaries or inside elements. For example, to define lines numbered 23 and 24 as excitation lines with an impedance of 120Ω, the corresponding GMSH command is as follows:**
```
Physical Line(120) = {24, 23};
```
***
The initial property for all surfaces is set to 1, representing a perfect electric conductor (PEC). Subsequent algorithm versions will be extended to support heterogeneous material analysis. Definition example:
```
Physical Surface(1) = {64, 81, 93, 66, 68, 69};
```

**After completing the above setup, use the GMSH software to generate the corresponding `.msh` file.**
***
The physical model is converted into a mesh model containing physical information. The following lists the content of the mesh file for the C-shaped resonator to illustrate the main format of the input mesh data (based on GMSH-2.2).

***1. Version Number***
```
$MeshFormat
2.2 0 8
$EndMeshFormat
```
***2. Node Information***
```
$Nodes
1050
1 20 -20 1.8
2 0 -20 1.8
...
1049 18.49999999999983 -16.35000000000062 -1.8
1050 18.50000000000017 16.35000000000063 -1.8
$EndNodes
```
***3. Element Information***
```
$Elements
1634
1 1 2 120 23 17 16
2 1 2 120 24 15 18
...
1632 2 2 1 93 377 1012 1050
1633 2 2 1 93 262 1049 1013
1634 2 2 1 93 378 1050 1014
$EndElements
```
***Note: A line must not be segmented into multiple parts during meshing; it must remain a single, complete element.***

## 3. Usage Instructions

This section introduces the software's features and usage method through an application case involving the simulation and analysis of a C-shaped resonator (Figure 1).

<div align="center">
 <img width="1951" height="1002" alt="5" src="https://github.com/user-attachments/assets/3135430b-5170-4661-8fa2-70ae4df0dfa6" />
</div>

<p align="center">Figure 1: Schematic Diagram of the C-shaped Resonator</p>

### 3.1 Data Import

The software reads three types of input data: mesh information of the physical structure, material information, and simulation settings. The following describes the directory structure of the code.

```
Sur_PEC/

├── Data/        # Data storage directory

    ├── model.msh         # Model file for simulation

    ├── Output/         # Output results directory

        ├── CirNetList.txt   # Circuit file

        ├── CirNetList.sp   # SPICE netlist

        ├── map.txt     # Port scattering parameters

        ├── Z_in.txt     # Port impedance

        └── Y_in.txt     # Port admittance

    ├── DIELECTRIC.txt  # Dielectric parameters file

    └── set.txt         # Configuration file

├── libiomp5md.dll      # Dynamic link library

└── Sur_PEC.exe         # Simulation executable program
```

#### 3.1.1 **model.msh** - Model File

-   **Purpose**: Specifies the mesh file required for simulation.
-   **Note**: Rename the model file to be simulated as `model.msh` and place it in the `Data` folder to ensure the simulation runs correctly and smoothly.

***
#### 3.1.2 **DIELECTRIC.txt** - Material Information Document

-   **Purpose**: Defines the dielectric material parameters used in the simulation.
-   **Explanation**: This file is a preset for future support of heterogeneous material simulation. **In the current version, no modification is needed.**

***
#### 3.1.3 **set.txt** - Main Configuration File

-   **Purpose**: Sets the main parameters for the simulation.
-   **Important Note**: **This file must be modified according to specific simulation requirements.**
-   **Six simulation parameters will be entered. The main parameters that need modification include**:
    -   **Solver_SET** (Frequency Domain/Time Domain Solver, where 0 indicates frequency domain solving, 1 indicates time domain solving)
    -   **FS** (Simulation Start Frequency)
    -   **FE** (Simulation End Frequency)
    -   **N_FP** (Number of Frequency Points to Solve)
    -   **DIM** (Scaling Dimension)
    -   **ER** (Dielectric Constant)

```
Solver_SET    0
FS 		5e8
FE 		5e9
N_FP 	200
DIM		1e3
ER 		1.0
```

As shown in the example configuration above, it represents frequency domain solving, with a frequency range from 0.5 GHz to 5 GHz for the electromagnetic device. Within this range, 200 frequency points are evenly selected as calculation frequencies. The model is defined as having millimeter-scale dimensions (1m = 1000mm), and the dielectric constant is 1.

***Note: The time domain solver in this version is a demo version. Please use the frequency domain solver and set Solver_SET to 0 by default.***

***

### 3.2 Equivalent Circuit Extraction

As shown in Figure 2, based on the Partial Element Equivalent Circuit theory, the software converts the model's mesh and connection network into a partial element equivalent circuit composed of distributed capacitances and inductances.

<div align="center">
 <img width="1020" height="333" alt="2" src="https://github.com/user-attachments/assets/56da100f-5a2b-4352-8649-7c579ab2114d" />
</div>

<p align="center">Figure 2: Conceptual Diagram of Model Discretization Translated into a Partial Element Equivalent Circuit</p>

This section generates two matrices through integration: the potential coefficient matrix (P) and the inductance matrix (L).

### 3.3 Circuit Solution and Analysis

This section utilizes self-developed code to analyze the extracted equivalent circuit in the frequency domain using the Modified Nodal Analysis (MNA) method, solving for information such as currents, potentials, and charges within the circuit. The core equation is as follows:

<div align="center">
 <img width="442" height="133" alt="3" src="https://github.com/user-attachments/assets/b4328c20-57eb-4fd8-8d61-177e7fbad914" />
</div>

<p align="center">Core Equation</p>

In this example, a total of 200 frequency points between 0.5GHz and 5GHz are selected for calculation. During the solving process, an equivalent circuit document is also generated containing information such as: number of ports, port nodes, number of elements, equivalent element names, connecting nodes, and coupling strength. The content is as follows:

```
N_PORT 2            Number of Ports                                              Port Information
PO_1      1369     264
PO_2       815    1092        Port Nodes

N_ELEMENT 6260928
LL_1_1           0         2  8.7224591179E-11                            Inductance Information
LM_1_2           0         2         0         3 -2.1744683326E-11
LM_1_3           0         2         1        29  8.6677865715E-13
...Equivalent Element Name               Connecting Nodes                  Element Value
LM_2220_2218        1621      1622      1619      1631 -9.2443493820E-13
LM_2220_2219        1621      1622      1620      1623 -3.8458394158E-11

CC_1_1           1         0  4.2690843441E-14                            Capacitance Information
CC_1_2           1         2 -2.1754451536E-18
...
CC_1631_1632        1631      1632 -3.2329462645E-20
CC_1632_1632        1632         0  3.0259206314E-14
```

***Additionally, this software provides a SPICE netlist. For higher accuracy, you may import the `.sp` file into HSPICE (or compatible SPICE simulators) for independent solving. The figure below, for reference, illustrates results generated by HuaDa SPICE.***

<div align="center">
 <img width="9648" height="6840" alt="6" src="https://github.com/user-attachments/assets/480179d3-d487-4a66-940e-e32773fbfb68" />
</div>

<p align="center">Figure 3: relevant results generated by HuaDa SPICE</p>

### 3.4 Performance Evaluation

This section performs post-processing on the obtained results, generating the following three documents stored in the `Output` folder:

-   **Port Scattering Parameters (S-parameters)**: Stored in `Sur_PEC\Data\Output\map.txt`
    ```
       Frequency (Hz)  Real(S11)(dB)     Imag(S11)    Real(S22)(dB)     Imag(S22)    Real(S12)(dB)    Imag(S12)     Real(S21)(dB)    Imag(S21)
       5.225000e+08,  -1.482463e+01,  3.962179e+01,  -1.453806e-01,  1.296218e+02,  -1.454241e-01,  1.296218e+02,  -1.482463e+01,  3.962182e+01
       5.450000e+08,  -1.446840e+01,  3.743193e+01,  -1.580369e-01,  1.274319e+02,  -1.580843e-01,  1.274319e+02,  -1.446840e+01,  3.743195e+01
       5.675000e+08,  -1.412683e+01,  3.523903e+01,  -1.712264e-01,  1.252390e+02,  -1.712779e-01,  1.252390e+02,  -1.412683e+01,  3.523901e+01
       ...
    ```
-   **Port Impedance**: Stored in `Sur_PEC\Data\Output\Z_in.txt`
    ```
       Frequency (Hz)   Real(Z11)       Imag(Z11)       Real(Z22)      Imag(Z22)      Real(Z12)       Imag(Z12)      Real(Z21)       Imag(Z21)
       5.150000e+08,  0.000000e+00,  -1.355504e+02,  -0.000000e+00,  2.078720e+02,  -0.000000e+00,  2.078706e+02,  0.000000e+00,  -1.355489e+02
       5.300000e+08,  0.000000e+00,  -1.288976e+02,  -0.000000e+00,  2.038118e+02,  -0.000000e+00,  2.038103e+02,  0.000000e+00,  -1.288961e+02
       5.450000e+08,  0.000000e+00,  -1.225001e+02,  -0.000000e+00,  2.000577e+02,  -0.000000e+00,  2.000562e+02,  0.000000e+00,  -1.224985e+02
       5.600000e+08,  0.000000e+00,  -1.163339e+02,  -0.000000e+00,  1.965888e+02,  -0.000000e+00,  1.965872e+02,  0.000000e+00,  -1.163323e+02
       ...
    ```
-   **Port Admittance**: Stored in `Sur_PEC\Data\Output\Y_in.txt`
    ```
       Frequency (Hz)   Real(Y11)       Imag(Y11)       Real(Y22)       Imag(Y22)     Real(Y12)       Imag(Y12)      Real(Y21)      Imag(Y21)
       5.150000e+08,  0.000000e+00,  -5.457594e-03,  0.000000e+00,  -8.369530e-03,  0.000000e+00,  -8.369473e-03,  0.000000e+00,  -5.457653e-03
       5.300000e+08,  0.000000e+00,  -5.171451e-03,  0.000000e+00,  -8.177151e-03,  0.000000e+00,  -8.177091e-03,  0.000000e+00,  -5.171511e-03
       5.450000e+08,  0.000000e+00,  -4.896669e-03,  0.000000e+00,  -7.996965e-03,  0.000000e+00,  -7.996903e-03,  0.000000e+00,  -4.896731e-03
       5.600000e+08,  0.000000e+00,  -4.632269e-03,  0.000000e+00,  -7.828027e-03,  0.000000e+00,  -7.827964e-03,  0.000000e+00,  -4.632333e-03
       ...
    ```

***
## 4. Important Notes

1.  Ensure all necessary DLL files are located in the runtime directory.
2.  The input format only supports GMSH mesh files (`.msh`). It is recommended to use GMSH for generation and ensure the mesh consists of triangular elements.
3.  Material files must be written strictly according to the specified format to avoid parsing errors.
4.  Be mindful of memory usage and computation time when simulating large models.

## 5. Recommended Visualization Tools

-   Mesh files (`.msh`) are recommended to be opened with **GMSH.exe**. GMSH.exe will directly generate interactive 3D graphics. (GMSH is included in the folder).
-   Parameters such as current, voltage, port scattering parameters, and impedance parameters (`.txt`) are recommended to be visualized using data analysis tools like **Excel** or **Origin**.

## 6. References

[1] Jiang Y, Wu K L. Quasi-static surface-PEEC modeling of electromagnetic problem with finite dielectrics[J]. IEEE Transactions on Microwave Theory and Techniques, 2018, 67(2): 565-576.

## License and Acknowledgments

This software is developed for research and educational purposes. Copyright belongs to the authors.

Feedback and collaboration are welcome to jointly advance the development of domestic EDA tools.
