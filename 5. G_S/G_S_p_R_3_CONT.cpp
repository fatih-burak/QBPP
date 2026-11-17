#include "helper_functions.h"
#include "model.h"

// Glover formulation F^S_{G,+} 
void G_S_p_R(Instance& inst, Solution& sol) {
    try {
        const double totalStart = getCPUTime();

        const int n = inst.n;
        const int W = inst.W;
        const int c = inst.c;
        const vector<int>& w = inst.w;
        const vector<vector<int>>& d = inst.d;

        vector<vector<bool>> compatible(n);
        for (int i = 0; i<n; i++) {
            compatible[i].resize(i+1, true);
            for (int k = 0; k<i; k++) {
                compatible[i][k] = (w[i]+w[k] <= W);
            }
        }

        const double boundsStart = getCPUTime();
        GRBEnv boundsEnv(true);
        boundsEnv.set("OutputFlag", "0");
        boundsEnv.start();
        // Solve one knapsacktype bound problem.

        auto solveBoundKnapsack = [&](int item,int capacity,bool useNegativeTerms,const std::string& modelName) -> double {
            vector<int> candidates;
            candidates.reserve(n - item - 1);
            for (int j = item+1; j<n; j++) {
                if ((useNegativeTerms && d[item][j]<0) ||(!useNegativeTerms && d[item][j] > 0)) {
                    candidates.push_back(j);
                }
            }
            if (candidates.empty()) return 0.0;        

            GRBModel boundModel(boundsEnv);
            boundModel.set(GRB_StringAttr_ModelName, modelName);
            boundModel.getEnv().set(GRB_IntParam_Threads, 1);
            boundModel.getEnv().set(GRB_IntParam_Seed, 1);
            boundModel.getEnv().set(GRB_DoubleParam_MIPGap, 0.0);

            vector<GRBVar> y;
            y.reserve(candidates.size());

            for (int j : candidates) {
                y.push_back(boundModel.addVar(0.0,1.0,0.0,GRB_BINARY,"y_"+std::to_string(j)));
            }

            GRBLinExpr knapsackWeight = 0.0;
            GRBLinExpr knapsackObjective = 0.0;

            for (std::size_t p = 0; p<candidates.size(); p++) {
                const int j = candidates[p];
                knapsackWeight += w[j] * y[p];
                knapsackObjective += d[item][j] * y[p];
            }

            boundModel.addConstr(knapsackWeight <= capacity);
            boundModel.setObjective(knapsackObjective,useNegativeTerms ? GRB_MINIMIZE : GRB_MAXIMIZE);
            boundModel.optimize();

            const int boundStatus = boundModel.get(GRB_IntAttr_Status);
            if (boundStatus != GRB_OPTIMAL) {
                throw std::runtime_error("Gurobi did not solve bound subproblem "+modelName +" to proven optimality. Status = " +std::to_string(boundStatus));
            }
            double value = boundModel.get(GRB_DoubleAttr_ObjVal);
            if (std::abs(value) <= 1e-9) value = 0.0;
            return value;
        };

        vector<vector<double>> lTilde(n - 1);
        vector<vector<double>> uInTilde(n - 1);
        vector<vector<double>> uOutTilde(n - 1);

        for (int i = 0; i<n - 1; i++) {
            lTilde[i].assign(i+1, 0.0);
             uInTilde[i].assign(i+1, 0.0);
            uOutTilde[i].assign(i+1, 0.0);

            for (int k = 0; k <= i; k++) {
                // When x[i][k] = 1, item i is in bin k. If i > k, item k is also in the bin.
                const int inCapacity =W - w[i] - ((i > k) ? w[k] : 0);
                //here compatible is actually same as checking whether inCapcity<0 or not.
                if (compatible[i][k]) {
                    lTilde[i][k] = solveBoundKnapsack(i,inCapacity,true,"lower_bound_"+std::to_string(i)+"_" +std::to_string(k));
                    uInTilde[i][k] = solveBoundKnapsack(i,inCapacity,false,"upper_in_bound_"+std::to_string(i)+"_" +std::to_string(k));
                }
                else {
                    // x[i][k] and t[i][k] will be fixed to zero.
                    lTilde[i][k] = 0.0;
                    uInTilde[i][k] = 0.0;
                }
                if (k == i) {
                    // For the diagonal, W - w[k] = W - w[i], so the out and in upper-bound subproblems equaj
                    uOutTilde[i][k] = uInTilde[i][k];
                }
                else {
                    //When x[i][k] = 0 but bin k is open, item k still occupies capacity. Item i does not occupy cap
                    uOutTilde[i][k] = solveBoundKnapsack(i,W - w[k],false,"upper_out_bound_"+std::to_string(i)+"_" +std::to_string(k));
                }

                if (lTilde[i][k] > 1e-7 ||uInTilde[i][k]<-1e-7 ||uOutTilde[i][k]<-1e-7) {
                    throw std::runtime_error("BIIIP: Inconsistent signs in item-bin bounds for (i,k) = (" +std::to_string(i)+","+std::to_string(k)+").");
                }

                if (compatible[i][k] &&uInTilde[i][k] > uOutTilde[i][k]+1e-7) {
                    throw std::runtime_error("BIIIP: uIn exceeds uOut for compatible pair (i,k) = (" +std::to_string(i)+","+std::to_string(k)+").");
                }
            }
        }
        const double boundsTime = getCPUTime() - boundsStart;
        cout << "Item-bin-dependent knapsack bounds:\n";
        cout << "item\tbin\tcompatible\tlTilde\tuInTilde\tuOutTilde\n";
        for (int i = 0; i<n - 1; i++) {
            for (int k = 0; k <= i; k++) {
                cout << i << '\t'<< k << '\t'<< (compatible[i][k] ? 1 : 0) << '\t'<< lTilde[i][k] << '\t'<< uInTilde[i][k] << '\t'<< uOutTilde[i][k] << '\n';
            }
        }
        cout << "Bound computation time (seconds) = "<< boundsTime << '\n';
        sol.timePatternGen = boundsTime;

        GRBEnv env(true);
        // env.set("OutputFlag", "0");
        env.start();

        GRBModel model(env);
        model.set(GRB_StringAttr_ModelName,"F_S_G_plus_R_item_bin_ILP_bounds");

        vector<vector<GRBVar>> x(n);
        for (int i = 0; i<n; i++) {
            x[i].reserve(i+1);
            for (int k = 0; k <= i; k++) {
                const double upperBound = compatible[i][k] ? 1.0 : 0.0;
                x[i].push_back(model.addVar(0.0,upperBound,0.0,GRB_BINARY,"x_"+std::to_string(i)+"_" +std::to_string(k)));
            }
        }
        vector<vector<GRBVar>> t(n - 1);
        for (int i = 0; i<n - 1; i++) {
            t[i].reserve(i+1);
            for (int k = 0; k <= i; k++) {
                t[i].push_back(model.addVar(lTilde[i][k],uInTilde[i][k],0.0,GRB_CONTINUOUS,"t_"+std::to_string(i)+"_" +std::to_string(k)));
            }
        }
        for (int i = 0; i<n; i++) {
            GRBLinExpr assignment = 0.0;
            for (int k = 0; k <= i; k++) {
                assignment += x[i][k];
            }
            model.addConstr(assignment == 1.0);
        }

        for (int k = 0; k<n; k++) {
            GRBLinExpr binWeight = 0.0;
            for (int i = k; i<n; i++) {
                binWeight += w[i] * x[i][k];
            }
            model.addConstr(binWeight <= W * x[k][k]);
        }
        for (int i = 1; i<n; i++) {
            for (int k = 0; k<i; k++) {
                model.addConstr(x[i][k] <= x[k][k]);
            }
        }
        for (int i = 0; i<n - 1; i++) {
            for (int k = 0; k <= i; k++) {
                model.addConstr(lTilde[i][k] * x[i][k] <= t[i][k]);
            }
        }
        for (int i = 1; i<n - 1; i++) {
            for (int k = 0; k<i; k++) {
                GRBLinExpr dissimilarity = 0.0;
                for (int j = i+1; j<n; j++) {
                    dissimilarity += d[i][j] * x[j][k];
                }
                model.addConstr(dissimilarity -uOutTilde[i][k] * (1.0 - x[i][k]) <= t[i][k]);
            }
        }
        for (int i = 0; i<n - 1; i++) {
            GRBLinExpr diagonalDissimilarity = 0.0;
            for (int j = i+1; j<n; j++) {
                diagonalDissimilarity += d[i][j] * x[j][i];
            }
            model.addConstr(diagonalDissimilarity <= t[i][i]);
        }

        GRBLinExpr objective = 0.0;
        for (int i = 0; i<n -1; i++) {
            for (int k = 0; k <= i; k++) {
                objective += t[i][k];
            }
        }
        for (int k = 0; k<n; k++) {
            objective += c*x[k][k];
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
        cout << "Total preprocessing and solve time (seconds): "<< sol.timeT << '\n';

        const int status = model.get(GRB_IntAttr_Status);
        const int solutionCount = model.get(GRB_IntAttr_SolCount);

        if (status == GRB_OPTIMAL) {
            sol.opt = 1;
            sol.valOpt = model.get(GRB_DoubleAttr_ObjVal);
            sol.valBest = model.get(GRB_DoubleAttr_ObjVal);
            sol.valBound = model.get(GRB_DoubleAttr_ObjBound);

            cout << "Optimal solution found. Objective = "<< sol.valOpt << '\n';
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
            for (int k = 0; k<n; k++) {
                if (x[k][k].get(GRB_DoubleAttr_X) > 0.5) {
                    cout << "Bin " << k << " used: ";
                    for (int i = k; i<n; i++) {
                        if (x[i][k].get(GRB_DoubleAttr_X) > 0.5) {
                            cout << "Item " << i << ' ';
                        }
                    }
                    cout << '\n';
                }
            }

            cout << "\nX =\n";
            for (int i = 0; i<n; i++) {
                for (int k = 0; k<n; k++) {
                    double value = 0.0;
                    if (k <= i) {
                        value = x[i][k].get(GRB_DoubleAttr_X);
                    }
                    cout << std::fixed << std::setprecision(2)<< value << '\t';
                }
                cout << '\n';
            }
            cout << "\nT =\n";
            for (int i = 0; i<n; i++) {
                for (int k = 0; k<n; k++) {
                    double value = 0.0;
                    if (i<n - 1 && k <= i) {
                        value = t[i][k].get(GRB_DoubleAttr_X);
                    }
                    cout << std::fixed << std::setprecision(2)<< value << '\t';
                }
                cout << '\n';
            }
        }
        const double relaxStart = getCPUTime();

        GRBModel relaxModel = model.relax();

        relaxModel.set(GRB_StringAttr_ModelName,"F_S_G_plus_R_item_bin_ILP_bounds_relaxation");
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
        std::cerr << "Gurobi error code = " << e.getErrorCode()<< '\n' << e.getMessage() << '\n';
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
    }
    catch (...) {
        std::cerr << "Unknown exception during optimization.\n";
    }
}
