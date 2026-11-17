#include "helper_functions.h"
#include "model.h"

void BQP_S_R(Instance& inst, Solution& sol) {
    try {
        double totalStart = getCPUTime();

        const int n = inst.n;
        const int W = inst.W;
        const int c = inst.c;

        const std::vector<int>& w = inst.w;
        const std::vector<std::vector<int>>& d = inst.d;

        GRBEnv env(true);
        // env.set("LogFile", "qbpp.log");
        // env.set("OutputFlag", "0");
        env.start();

        GRBModel model(env);

        std::vector<std::vector<GRBVar>> x(n); //x[i][k] exists only for 0 <= k <= i < n

        for (int i = 0; i < n; i++) {
            x[i].reserve(i + 1);

            for (int k = 0; k <= i; k++) {
                x[i].push_back(model.addVar(0, 1, 0, GRB_BINARY, "x_" + std::to_string(i) + "_" + std::to_string(k)));
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

        for (int k = 0; k < n; ++k) {
            for (int i = k + 1; i < n; ++i) {
                model.addConstr(x[i][k] <= x[k][k]);
            }
        }

        // WE did not create variables x[i][k] with i < k soo constraints x[i][k] = 0 are unnecessary.

        GRBQuadExpr objective = 0;

        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                for (int k = 0; k <= i; k++) {
                    objective += d[i][j] *x[i][k] * x[j][k];
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

        model.getEnv().set(GRB_DoubleParam_MIPGap, 0.0);
        model.getEnv().set(GRB_IntParam_Threads, 1);
        model.getEnv().set(GRB_DoubleParam_TimeLimit, inst.timeLimit);
        model.getEnv().set(GRB_IntParam_Seed, 1);
        // model.getEnv().set(GRB_IntParam_NonConvex, 2); //not necassary in the new versions of gurobi. gurobi automatically adjusts

        model.optimize();

        sol.timeT = getCPUTime() - totalStart;
        std::cout << "Solve time (seconds): " << sol.timeT << '\n';

        const int status = model.get(GRB_IntAttr_Status);
        const int solutionCount = model.get(GRB_IntAttr_SolCount);

        if (status == GRB_OPTIMAL) {
            sol.opt = 1;
            sol.valOpt = model.get(GRB_DoubleAttr_ObjVal);
            sol.valBest = model.get(GRB_DoubleAttr_ObjVal);
            sol.valBound = model.get(GRB_DoubleAttr_ObjBound);

            std::cout << "Optimal solution found. Objective = " << sol.valOpt << '\n';

            // Print used bins and assigned items
            for (int k = 0; k < n; k++) {
                if (x[k][k].get(GRB_DoubleAttr_X) > 0.5) {
                    std::cout << "Bin " << k << " used: ";
                    for (int i = k; i < n; i++) {
                        if (x[i][k].get(GRB_DoubleAttr_X) > 0.5) {
                            std::cout << "Item " << i << ' ';
                        }
                    }
                    std::cout << '\n';
                }
            }

            std::cout << "X =\n";

            for (int i = 0; i < n; i++) {
                for (int k = 0; k < n; k++) {
                    const double value = (k <= i) ? x[i][k].get(GRB_DoubleAttr_X) : 0.0;

                    std::cout << std::fixed << std::setprecision(2) << value << '\t';
                }
                std::cout << '\n';
            }
        }
        else if (status == GRB_TIME_LIMIT) {
            std::cout << "Time limit reached.\n";
            sol.opt = 0;
            sol.valBound = model.get(GRB_DoubleAttr_ObjBound);
            std::cout << "Best bound = " << sol.valBound << '\n';

            if (solutionCount > 0) {
                sol.valBest = model.get(GRB_DoubleAttr_ObjVal);

                std::cout << "Best solution = " << sol.valBest << '\n';
            }
            else {
                std::cout << "No feasible solution was found.\n";
            }
        }
        else {
            sol.opt = 0;
            sol.valBound = model.get(GRB_DoubleAttr_ObjBound);

            std::cout << "Optimization ended with status " << status << '\n';
        }

        // CONTINUOUS QUADRATIC RELAXATION
        /*
        double relaxStart = getCPUTime();
        GRBModel relaxModel = model.relax();

        relaxModel.getEnv().set(GRB_DoubleParam_TimeLimit, inst.timeLimit);
        relaxModel.getEnv().set(GRB_IntParam_Threads, 1);
        relaxModel.getEnv().set(GRB_DoubleParam_MIPGap, 0.0);

        relaxModel.optimize();

        sol.timeR = getCPUTime() - relaxStart;
        
        const int statusR = relaxModel.get(GRB_IntAttr_Status);
        if (statusR == GRB_OPTIMAL) {
            sol.optR = 1;
            sol.valRel = relaxModel.get(GRB_DoubleAttr_ObjVal);
            std::cout << "Continuous quadratic relaxation value = " << sol.valRel << '\n';
        }
        else {
            sol.optR = 0;
            std::cout << "Continuous quadratic relaxation ended with status "  << statusR << '\n';
        }
        */
    }
    catch (const GRBException& e) {
        std::cerr << "Gurobi error code = " << e.getErrorCode() << '\n' << e.getMessage() << '\n';
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
    }
    catch (...) {
        std::cerr << "Unknown exception during optimization!!\n";
    }
}