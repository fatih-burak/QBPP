#include "helper_functions.h"
#include "model.h"

// FP formulation.
void FP_S_p_R(Instance& inst, Solution& sol) {
    try {
        const double totalStart = getCPUTime();

        const int n = inst.n;
        const int W = inst.W;
        const int c = inst.c;
        const vector<int>& w = inst.w;
        const vector<vector<int>>& d = inst.d;

        GRBEnv env(true);
        // env.set("OutputFlag", "0");
        env.start();

        GRBModel model(env);
        model.set(GRB_StringAttr_ModelName, "F_S_plus");

        vector<vector<GRBVar>> x(n);
        for (int i = 0; i < n; i++) {
            x[i].reserve(i+1);
            for (int k = 0; k <= i; k++) {
                x[i].push_back(model.addVar(0.0,1.0,0.0,GRB_BINARY,"x_" + std::to_string(i) + "_" + std::to_string(k)));
            }
        }

        vector<vector<GRBVar>> z(n, vector<GRBVar>(n));
        int nbPositivePairs = 0;
        int nbNegativePairs = 0;
        int nbZeroPairs = 0;
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                if (d[i][j] > 0) {
                    nbPositivePairs++;
                    z[i][j] = model.addVar(0.0,GRB_INFINITY,0.0,GRB_CONTINUOUS,"z_" + std::to_string(i) + "_" + std::to_string(j));
                }
                else if (d[i][j] < 0) {
                    nbNegativePairs++;
                    z[i][j] = model.addVar(0.0,GRB_INFINITY,0.0,GRB_CONTINUOUS,"z_" + std::to_string(i) + "_" + std::to_string(j));
                }
                else {
                    nbZeroPairs++;
                }
            }
        }
        for (int i = 0; i < n; i++) {
            GRBLinExpr assignment = 0.0;
            for (int k = 0; k <= i; k++) {
                assignment += x[i][k];
            }
            model.addConstr(assignment == 1.0);
        }
        for (int k = 0; k < n; k++) {
            GRBLinExpr binWeight = 0.0;
            for (int i = k; i < n; i++) {
                binWeight += w[i] * x[i][k];
            }
            model.addConstr(binWeight <= W * x[k][k]);
        }
        for (int i = 1; i < n; i++) {
            for (int k = 0; k < i; k++) {
                model.addConstr(x[i][k] <= x[k][k]);
            }
        }
         
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                if (d[i][j] <= 0) continue;
                for (int k = 0; k <= i; k++) {
                    model.addConstr(x[i][k] + x[j][k] - x[k][k] <= z[i][j]);
                }
            }
        }
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                if (d[i][j] >= 0) continue;
                for (int k = 0; k <= i; k++) {
                    model.addConstr(z[i][j] <= 1.0 - x[i][k] + x[j][k]);
                }
            }
        }
        for (int i = 1; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                if (d[i][j] >= 0) continue;
                for (int k = 0; k < i; k++) {
                    model.addConstr(z[i][j] <= 1.0 + x[i][k] - x[j][k]);
                }
            }
        }
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                if (d[i][j] >= 0) continue;
                for (int k = i + 1; k <= j; k++) {
                    model.addConstr(z[i][j] <= 1.0 - x[j][k]);
                }
            }
        }
        GRBLinExpr objective = 0.0;
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                if (d[i][j] != 0) {
                    objective += d[i][j] * z[i][j];
                }
            }
        }
        for (int k = 0; k < n; k++) {
            objective += c * x[k][k];
        }
        model.setObjective(objective, GRB_MINIMIZE);
        model.update();

        sol.nbVars = model.get(GRB_IntAttr_NumVars);
        sol.nbConstrs = model.get(GRB_IntAttr_NumConstrs);
        sol.nbNZs = model.get(GRB_IntAttr_NumNZs);

        cout << "Model size:\n";
        cout << "Variables = " << sol.nbVars << '\n';
        cout << "Constraints = " << sol.nbConstrs << '\n';
        cout << "Nonzeros = " << sol.nbNZs << '\n';

        model.getEnv().set(GRB_DoubleParam_MIPGap, 0.0);
        model.getEnv().set(GRB_IntParam_Threads, 1);
        model.getEnv().set(GRB_DoubleParam_TimeLimit, inst.timeLimit);
        model.getEnv().set(GRB_IntParam_Seed, 1);
        model.optimize();

        sol.timeT = getCPUTime() - totalStart;
        cout << "Solve time (seconds): " << sol.timeT << '\n';

        const int status = model.get(GRB_IntAttr_Status);
        const int solutionCount = model.get(GRB_IntAttr_SolCount);

        if (status == GRB_OPTIMAL) {
            sol.opt = 1;
            sol.valOpt = model.get(GRB_DoubleAttr_ObjVal);
            sol.valBest = model.get(GRB_DoubleAttr_ObjVal);
            sol.valBound = model.get(GRB_DoubleAttr_ObjBound);

            cout << "Optimal solution found. Objective = " << sol.valOpt << '\n';
        }
        else if (status == GRB_TIME_LIMIT) {
            sol.opt = 0;
            sol.valBound = model.get(GRB_DoubleAttr_ObjBound);

            cout << "Time limit reached.\n";
            cout << "Best bound = " << sol.valBound << '\n';

            if (solutionCount > 0) {
                sol.valBest = model.get(GRB_DoubleAttr_ObjVal);
                cout << "Best solution = " << sol.valBest << '\n';
            }
            else {
                cout << "No feasible solution was found.\n";
            }
        }
        else if (status == GRB_INFEASIBLE) {
            sol.opt = 0;
            cout << "The model is infeasible.\n";
        }
        else {
            sol.opt = 0;
            cout << "Optimization ended with status "<< status << ".\n";
        }

        if (solutionCount > 0) {
            cout << "\nUsed bins and assigned items:\n";
            for (int k = 0; k < n; k++) {
                if (x[k][k].get(GRB_DoubleAttr_X) > 0.5) {
                    cout << "Bin " << k << " used: ";
                    for (int i = k; i < n; i++) {
                        if (x[i][k].get(GRB_DoubleAttr_X) > 0.5) {
                            cout << "Item " << i << ' ';
                        }
                    }
                    cout << '\n';
                }
            }
            cout << "\nX =\n";

            for (int i = 0; i < n; i++) {
                for (int k = 0; k < n; k++) {
                    double value = 0.0;
                    if (k <= i) {
                        value = x[i][k].get(GRB_DoubleAttr_X);
                    }
                    cout << std::fixed << std::setprecision(2)<< value << '\t';
                }
                cout << '\n';
            }
            cout << "\nZ =\n";
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    double value = 0.0;
                    if (i < j && d[i][j] != 0) {
                        value = z[i][j].get(GRB_DoubleAttr_X);
                    }
                    else if (j < i && d[j][i] != 0) {
                        value = z[j][i].get(GRB_DoubleAttr_X);
                    }
                    cout << std::fixed << std::setprecision(2)<< value << '\t';
                }
                cout << '\n';
            }
        }
        const double relaxStart = getCPUTime();
        GRBModel relaxModel = model.relax();
        relaxModel.set(GRB_StringAttr_ModelName,"F_S_plus_FP_sign_reduced_relaxation");
        relaxModel.getEnv().set(GRB_DoubleParam_TimeLimit,inst.timeLimit);
        relaxModel.getEnv().set(GRB_IntParam_Threads, 1);

        relaxModel.optimize();
        sol.timeR = getCPUTime() - relaxStart;
        const int statusR = relaxModel.get(GRB_IntAttr_Status);

        if (statusR == GRB_OPTIMAL) {
            sol.optR = 1;
            sol.valRel = relaxModel.get(GRB_DoubleAttr_ObjVal);

            cout << "\nContinuous linear relaxation value = "<< sol.valRel << '\n';
            cout << "Relaxation solve time (seconds) = "<< sol.timeR << '\n';
        }
        else {
            sol.optR = 0;
            cout << "\nContinuous linear relaxation ended with status "<< statusR << ".\n";
        }
    }
    catch (const GRBException& e) {
        std::cerr << "Gurobi error code = "<< e.getErrorCode() << '\n'<< e.getMessage() << '\n';
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
    }
    catch (...) {
        std::cerr << "Unknown exception during optimization.\n";
    }
}
