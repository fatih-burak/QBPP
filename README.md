# Mixed-integer programming formulations for the quadratic bin packing problem
This github page contains the C++ implementations of the formulations and hybrid approach studied in the paper **"Mixed-integer programming formulations for the quadratic bin packing problem: linearization, hybridization, and enhancements"** by Fatih Burak Akçay, Maxence Delorme and Renata Sotirov.

The **quadratic bin packing problem (QBPP)** extends the classical bin packing problem by considering both a fixed cost for each used bin and pairwise interaction costs between items assigned to the same bin. 
Here we provide the source code of the models implemented in the paper. 

All implementations use the commercial **Gurobi Optimizer** to solve the optimization models.

## Structure

The code is organized into seven main folders, each corresponding to a method discussed in the paper.

| Folder | Main model file | Description |
| --- | --- | --- |
| `1. BQP_S/` | `BQP_S_R.cpp` | Corresponding to $\mathrm{F^S_{BQP}}$ |
| `2. FP_S/` | `FP_S_p_R.cpp` | Corresponding to $\mathrm{F^{S+}_{FP}}$ |
| `3. GW_S/` | `GW_S_p_R.cpp` | Corresponding to $\mathrm{F^{S+}_{GW}}$ |
| `4. RLT_S/` | `F_S_RLT1_R.cpp` | Corresponding to $\mathrm{F^S_{RLT1}}$ |
| `5. G_S/` | `G_S_p_R_3_CONT.cpp` | Corresponding to $\mathrm{F^{S+}_{G}}$| 
| `6. P/` | `PatternILP.cpp` | Corresponding to $\mathrm{F_P}$|
| `7. Hybrid/` | `H_P_GW_I2_bothL.cpp` | Corresponding to $\mathrm{F^{S+}_{HI}}$ $`I_{10}^{(2)}\text{-}D^\downarrow`$ with global+cumulative valid inequalities|

### Files within each folder

Each folder follows a common structure:

- `main.cpp` : Main file for running the corresponding method.
- `main.h` : Header associated with `main.cpp`.
- `helper_functions.cpp` and `helper_functions.h` : Instance loading, result reporting and other shared helper functions.
- `model.h` : Declaration of the model implemented in the method-specific `.cpp` file listed above.

## Instance format

We did not generate the instances ourselves but we still provide them under the folder INSTANCES.
They come from Chagas et al. (2026), https://onlinelibrary.wiley.com/doi/full/10.1111/itor.70237 and can also be downloaded from https://github.com/regor-unimore/Quadratic-Bin-Packing-Problem. 
