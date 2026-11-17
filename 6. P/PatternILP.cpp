#include "helper_functions.h"
#include "model.h"

struct Pattern {  //what a pattern represents
    vector<int> items;  //indices of the items in the pattern
    double weight = 0.0; // total weight of those items
    double cost = 0.0; //fixed cost of using this in a bin plus the pairwise costs between items in the group
};


// generate all feasible patterns 
static void generatePatternsDFS(const Instance& inst, int pos, vector<int>& currentItems, double currentWeight, double currentPairCost, vector<Pattern>& patterns, vector<vector<int>>& itemToPatterns) {
    const double eps = 1e-9;

    if (pos == inst.n) { //when this is true, all items are considered
        if (!currentItems.empty()) { //if the selected subset is not empty, we create a pattern
            Pattern p;
            p.items = currentItems;
            p.weight = currentWeight;
            p.cost = inst.c + currentPairCost;
            int patternIndex = static_cast<int>(patterns.size()); // This pattern will be inserted at this index.
            patterns.push_back(p);

            // Record that every selected item belongs to this pattern.
            for (int item : currentItems) {
                itemToPatterns[item].push_back(patternIndex);
            }
        }
        return;
    }
    // For each item, make exactly two choices: skip it or include it.
    // Case 1: do not include item-pos (branch-1). recursion moves to the next itm without selecting the current item
    generatePatternsDFS(inst, pos + 1, currentItems, currentWeight, currentPairCost, patterns, itemToPatterns);

    // Case 2: include item-pos if capacity allows (branch-2)
    double newWeight = currentWeight + inst.w[pos];

    if (newWeight <= inst.W + eps) {
        double addedPairCost = 0.0;

        for (int i : currentItems) {
            addedPairCost += inst.d[min(i, pos)][max(i, pos)];
        }
        //backtracking
        currentItems.push_back(pos); //add the item

        generatePatternsDFS(inst, pos + 1, currentItems, newWeight, currentPairCost + addedPairCost, patterns, itemToPatterns); // explore ALL subsets containing it recursively
        //when all subsets containing it finishes, return to the state where this item was not included..
        currentItems.pop_back(); //remove it before exploring other branches
    }
}

void generateAllFeasiblePatterns(const Instance& inst, vector<Pattern>& patterns, vector<vector<int>>& itemToPatterns) {
    //vector<Pattern> patterns;
    vector<int> currentItems;
    itemToPatterns.resize(inst.n);

    //check any item is heavier than W. This should not happen as it creates infeasibility..
    for (int i = 0; i < inst.n; ++i) {
        if (inst.w[i] > inst.W + 1e-9) {
            throw runtime_error("Item " + to_string(i) + "has weight larger than bin capacity.");
        }
    }

    generatePatternsDFS(inst, 0, currentItems, 0.0, 0.0, patterns, itemToPatterns);
}

void PatternILP(Instance& inst, Solution& sol) {
    try {
        double totalStart = getCPUTime();

        int n = inst.n;
        cout << "\nGenerating all feasible patterns for PatternILP..." << endl;
        double patternStart = getCPUTime();
        vector<Pattern> patterns;
        vector<vector<int>> itemToPatterns;
        generateAllFeasiblePatterns(inst, patterns, itemToPatterns);
        double patternGenerationTime = getCPUTime() - patternStart;
        cout << "Number of feasible patterns generated: " << patterns.size() << endl;
        cout << "Pattern generation time (seconds): "<< patternGenerationTime << endl;
        sol.timePatternGen = patternGenerationTime; //save
        sol.nbPatterns = patterns.size(); //save
        if (patterns.empty()) throw runtime_error("No feasible patterns generated.");

        for (int i = 0; i < n; i++) {
            if (itemToPatterns[i].empty()) throw runtime_error("Item " + to_string(i) + " is not contained in any feasible pattern.");
        }

        GRBEnv env = GRBEnv(true);
        // env.set("OutputFlag", "0");
        env.start();

        GRBModel model = GRBModel(env);

        double startTimeMain = model.get(GRB_DoubleAttr_Runtime);

        // lambda[p] = 1 if pattern p is selected
        int P = static_cast<int>(patterns.size());

        vector<GRBVar> lambda(P);
        for (int p = 0; p < P; p++) {
            lambda[p] = model.addVar(0.0, 1.0, patterns[p].cost, GRB_BINARY, "lambda_" + to_string(p));
        }

        for (int i = 0; i < n; i++) {
            GRBLinExpr cover = 0.0;
            for (int p : itemToPatterns[i]) {
                cover += lambda[p];
            }
            model.addConstr(cover == 1.0, "cover_" + to_string(i));
        }

        // Objective coefficients were already aded in addVar.
        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        model.update();
        // ----------------------------------------------------------
        // INTEGER MODEL SIZE
        // ----------------------------------------------------------
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
        model.set(GRB_IntParam_Presolve, 0);
	    model.getEnv().set(GRB_IntParam_Sifting, 0);

        model.optimize();

        double endTimeMain = model.get(GRB_DoubleAttr_Runtime);

        // This includes model optimize time.
        // Pattern generation time is included below using getCPUTime().
        sol.timeT = getCPUTime() - totalStart;

        cout << "Gurobi solve time ILP (seconds): "  << endTimeMain - startTimeMain << endl;
        cout << "Total PatternILP time incl. pattern generation (seconds): " << sol.timeT << endl;

        int status = model.get(GRB_IntAttr_Status);
        int solCount = model.get(GRB_IntAttr_SolCount);

        if (status == GRB_OPTIMAL) {
            sol.opt = 1;
            sol.valOpt = model.get(GRB_DoubleAttr_ObjVal);
            sol.valBest = model.get(GRB_DoubleAttr_ObjVal);
            sol.valBound = model.get(GRB_DoubleAttr_ObjBound);

            cout << "Optimal PatternILP solution found. Objective = " << sol.valOpt<< endl;
        }
        else if (status == GRB_TIME_LIMIT) {
            sol.opt = 0;

            if (solCount > 0) {
                sol.valBest = model.get(GRB_DoubleAttr_ObjVal);
                sol.valOpt = sol.valBest;
                cout << "Best solution = " << sol.valBest << endl;
            }

            sol.valBound = model.get(GRB_DoubleAttr_ObjBound);
            cout << "Best bound = " << sol.valBound << endl;

            cout << "PatternILP time limit reached." << endl;
                             
        }
        else if (status == GRB_INFEASIBLE) {
            sol.opt = 0;
            cout << "PatternILP is infeasible." << endl;
        }
        else {
            sol.opt = 0;

            if (solCount > 0) {
                sol.valBest = model.get(GRB_DoubleAttr_ObjVal);
                sol.valOpt = sol.valBest;
            }

            sol.valBound = model.get(GRB_DoubleAttr_ObjBound);

            cout << "PatternILP ended with status " <<status << endl;
        }
        if (solCount > 0) {
            cout << "\nSelected patterns:" << endl;

            for (int p = 0; p < P; ++p) {
                if (lambda[p].get(GRB_DoubleAttr_X) > 0.5) {
                    cout << "Pattern " << p
                        << " | cost = " << patterns[p].cost
                        << " | weight = " << patterns[p].weight
                        << " | items = { ";
                    for (int item : patterns[p].items) cout << item << " ";
                    cout << "}" << endl;
                }
            }
        }

        // Create the LP relaxation.
        GRBModel relaxModel = model.relax();

        relaxModel.getEnv().set(GRB_DoubleParam_TimeLimit, inst.timeLimit);
        relaxModel.getEnv().set(GRB_IntParam_Threads, 1);
        //relaxModel.getEnv().set(GRB_IntParam_OutputFlag, 0);

        double startTimeRel = relaxModel.get(GRB_DoubleAttr_Runtime);

        relaxModel.optimize();

        double endTimeRel = relaxModel.get(GRB_DoubleAttr_Runtime);

        sol.timeR = endTimeRel - startTimeRel;

        int statusR = relaxModel.get(GRB_IntAttr_Status);
        if (statusR == GRB_OPTIMAL) {
            sol.optR = 1;
            sol.valRel = relaxModel.get(GRB_DoubleAttr_ObjVal);

            cout << "Pattern LP relaxation value = " << sol.valRel << endl;
        }
        else {
            sol.optR = 0;
            cout << "Pattern LP relaxation ended with status " << statusR << endl;
        }
           
    }
    catch (GRBException e) {
        cerr << "Error code = " << e.getErrorCode() << endl;
        cerr << e.getMessage() << endl;
    }
    catch (exception& e) {
        cerr << "Exception during PatternILP: "
            << e.what() << endl;
    }
    catch (...) {
        cerr << "Exception during PatternILP optimization" << endl;
    }
}