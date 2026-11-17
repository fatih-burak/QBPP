#include "helper_functions.h"
#include "model.h"

// --------------------------------------------------
// RLT1 formulation F^S_RLT1
// --------------------------------------------------
void F_S_RLT1_R(Instance& inst, Solution& sol) {
    try {
        double totalStart = getCPUTime();
        int n = inst.n;
        int B = inst.B;
        int W = inst.W;
        int c = inst.c;
        vector<int> w = inst.w;
        vector<vector<int>> d = inst.d;

        GRBEnv env = GRBEnv(true);
        env.start();

        GRBModel model = GRBModel(env);
        double startTimeMain = model.get(GRB_DoubleAttr_Runtime);

        vector<vector<GRBVar>> x(n, vector<GRBVar>(B));
        for (int i = 0; i < n; i++) {
            for (int k = 0; k <= i; k++) {
                x[i][k] = model.addVar(0.0, 1.0,0.0,GRB_BINARY,"x_" + to_string(i) + "_" + to_string(k));
            }
        }

        vector<vector<vector<GRBVar>>> ghat(n,vector<vector<GRBVar>>(n, vector<GRBVar>(B)));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int maxK = (i < j ? i : j);
                for (int k = 0; k <= maxK; k++) {
                    ghat[i][j][k] = model.addVar(0.0,1.0,0.0,GRB_CONTINUOUS,"ghat_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k));
                }
            }
        }

        for (int i = 0; i < n; i++) {
            GRBLinExpr assign = 0;
            for (int k = 0; k <= i; k++) {
                assign += x[i][k];
            }
            model.addConstr(assign == 1);
        }

        for (int k = 0; k < B; k++) {
            GRBLinExpr cap = 0;
            for (int i = k; i < n; i++) {
                cap += w[i] * x[i][k];
            }
            model.addConstr(cap <= W * x[k][k]);
        }

        for (int i = 0; i < n; i++) {
            for (int k = 0; k < i; k++) {
                model.addConstr(x[i][k] <= x[k][k]);
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int maxK = (i < j ? i : j);
                for (int k = 0; k <= maxK; k++) {
                    model.addConstr(ghat[i][j][k] <= x[i][k]);
                }
            }
        }

        for (int i = 0; i < n; i++) {
            for (int j = i + 1;j < n; j++) {
                for (int k = 0; k <= i; k++) {
                    model.addConstr(ghat[i][j][k] == ghat[j][i][k]);
                }
            }
        }

        for (int k = 0; k < B; k++) {
            for (int j = k + 1; j < n; j++) {
                GRBLinExpr lhs = 0;
                for (int i = k; i < n; i++) {
                    lhs += w[i] * ghat[i][j][k];
                }
                model.addConstr(lhs <= W * ghat[k][j][k]);
            }
        }

        for (int k = 0; k < B; k++) {
            for (int j = k + 1; j < n; j++) {
                GRBLinExpr lhs = 0;
                for (int i = k; i < n; i++) {
                    lhs += w[i] * (x[i][k] - ghat[i][j][k]);
                }
                model.addConstr(lhs <= W *(x[k][k] - ghat[k][j][k]));
            }
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int maxKExclusive = (i < j ? i : j);
                for (int k = 0; k < maxKExclusive; k++) {
                    model.addConstr(ghat[i][j][k] <= ghat[j][k][k]);
                }
            }
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int maxKExclusive = (i < j ? i : j);
                for (int k = 0; k < maxKExclusive; k++) {
                    model.addConstr(x[i][k] - x[k][k]<= ghat[i][j][k] - ghat[j][k][k]);
                }
            }
        }

        for (int i = 0; i < n; i++) {
            for (int k = 0; k <= i; k++) {
                model.addConstr(x[i][k] == ghat[i][i][k]);
            }
        }

        // ghat >= 0 is imposed through each variable's lower bound.
        GRBLinExpr obj = 0;
        for (int k = 0; k < B; k++) {
            for (int i = k; i < n; i++) {
                for (int j = k; j < n; j++) {
                    obj += 0.5 * d[i][j] * ghat[i][j][k];
                }
            }
        }
        for (int k = 0; k < B; k++) {
            obj += c * x[k][k];
        }
        model.setObjective(obj, GRB_MINIMIZE);
        model.update();

        sol.nbVars = model.get(GRB_IntAttr_NumVars);
        sol.nbConstrs = model.get(GRB_IntAttr_NumConstrs);
        sol.nbNZs = model.get(GRB_IntAttr_NumNZs);

        cout << "Model size:\n";
        cout << "Variables = " << sol.nbVars << '\n';
        cout << "Constraints = " << sol.nbConstrs << '\n';
        cout << "Nonzeros = " << sol.nbNZs << '\n';

        model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
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

            for (int k = 0; k < B; k++) {
                if (x[k][k].get(GRB_DoubleAttr_X) > 0.5) {
                    cout << "Bin " << k << " used: ";
                    for (int i = k; i < n; i++) {
                        if (x[i][k].get(GRB_DoubleAttr_X) > 0.5) {
                            cout << "Item " << i << " ";
                        }
                    }
                    cout << endl;
                }
            }

            // Print X 
            cout << "\nX matrix:\n";
            for (int i = 0; i < n; i++) {
                for (int k = 0; k < B; k++) {
                    if (k <= i) {
                        double val = x[i][k].get(GRB_DoubleAttr_X);
                        cout << (val > 0.5 ? 1 : 0) << "\t";
                    }
                    else {
                        cout << "0\t";
                    }
                }
                cout << "\n";
            }

            // Print ghat matrices 
            cout << "\nGhat matrices:\n";
            for (int k = 0; k < B; k++) {
                cout << "\nGhat matrix for bin k = " << k << ":\n";
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        if (k <= i && k <= j) {
                            double val = ghat[i][j][k].get(GRB_DoubleAttr_X);
                            cout << (val > 0.5 ? 1 : 0) << "\t";
                        }
                        else {
                            cout << "0\t";
                        }
                    }
                    cout << "\n";
                }
            }
        }
        else if (status == GRB_TIME_LIMIT) {
            sol.opt = 0;
            sol.valBound = model.get(GRB_DoubleAttr_ObjBound);

            if (solutionCount > 0) {
                sol.valBest = model.get(GRB_DoubleAttr_ObjVal);
            }

            cout << "Time limit reached." << endl;
            if (solutionCount > 0) {
                cout << "Best solution = " << sol.valBest<< ", Best bound = " << sol.valBound << endl;
            }
            else {
                cout << "No incumbent solution found. Best bound = " << sol.valBound << endl;
            }
        }
        else {
            sol.opt = 0;
            cout << "Optimization stopped with status = " << status << endl;
        }

        // LP RELAXATION
        GRBModel relaxModel = model.relax();
        relaxModel.getEnv().set(GRB_DoubleParam_TimeLimit, inst.timeLimit);

        double startTimeR = relaxModel.get(GRB_DoubleAttr_Runtime);
        relaxModel.optimize();

        vector<vector<GRBVar>> xR(n, vector<GRBVar>(B));
        vector<vector<vector<GRBVar>>> ghatR(n,vector<vector<GRBVar>>(n, vector<GRBVar>(B)));

        // Recover x var
        for (int i = 0; i < n; i++) {
            for (int k = 0; k <= i; k++) {
                xR[i][k] = relaxModel.getVarByName("x_" + to_string(i) + "_" + to_string(k));
            }
        }

        //Recover ghat var
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int maxK = (i < j ? i : j);
                for (int k = 0; k <= maxK; k++) {
                    ghatR[i][j][k] = relaxModel.getVarByName("ghat_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k));
                }
            }
        }

        double endTimeR = relaxModel.get(GRB_DoubleAttr_Runtime);
        sol.timeR = endTimeR - startTimeR;

        int statusR = relaxModel.get(GRB_IntAttr_Status);
        if (statusR == GRB_OPTIMAL) {
            sol.optR = 1;
            sol.valRel = relaxModel.get(GRB_DoubleAttr_ObjVal);

            // Print X  (LP)
            cout << "\nX matrix (LP relaxation):\n";
            for (int i = 0; i < n; i++) {
                for (int k = 0; k < B; k++) {
                    if (k <= i) {
                        double val = xR[i][k].get(GRB_DoubleAttr_X);
                        cout << fixed << setprecision(2) << val << "\t";
                    }
                    else {
                        cout << fixed << setprecision(2) << 0.0 << "\t";
                    }
                }
                cout << "\n";
            }

            // Print ghat matrices (LP)
            cout << "\nGhat matrices (LP relaxation):\n";
            for (int k = 0; k < B; k++) {
                cout << "\nGhat matrix for bin k = " << k << ":\n";
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        if (k <= i && k <= j) {
                            double val = ghatR[i][j][k].get(GRB_DoubleAttr_X);
                            cout << fixed << setprecision(2) << val << "\t";
                        }
                        else {
                            cout << fixed << setprecision(2) << 0.0 << "\t";
                        }
                    }
                    cout << "\n";
                }
            }
        }
        else {
            sol.optR = 0;
        }
    }
    catch (GRBException e) {
        cerr << "Error code = " << e.getErrorCode() << endl;
        cerr << e.getMessage() << endl;
    }
    catch (const exception& e) {
        cerr << "Exception during optimization: " << e.what() << endl;
    }
    catch (...) {
        cerr << "Exception during optimization" << endl;
    }
}
