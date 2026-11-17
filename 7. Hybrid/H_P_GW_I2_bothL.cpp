#include "helper_functions.h"
#include "model.h"

// hybrid pattern + Glover-Woolsey formulation with an item subset S.

//AUXILIARY BPP RESULT
struct PrefixBPPResult {
    // L[m-1] is the exact BPP for reordered items 0,...,m-1 whenever exact[m-1] == true.
    vector<int> L;
    vector<bool> exact;
    vector<double> timeByPrefix;
    double totalTime = 0.0;
};

// COMPUTE EXACT L_m VALUES FOR Sbar
PrefixBPPResult computePrefixBPPBounds(const Instance& inst,const vector<int>& rw,const int nSbar) {
    PrefixBPPResult result;

    result.L.assign(nSbar, 0);
    result.exact.assign(nSbar, false);
    result.timeByPrefix.assign(nSbar, 0.0);

    const int W = inst.W;
    const double totalAuxStart = getCPUTime();
    cout << "AUXILIARY BPPs FOR CUMULATIVE BIN-COUNT INEQUALITIES\n";

    int previousExactLB = 0;
    for (int m = 1; m <= nSbar; m++) {
        const double prefixStart = getCPUTime();
        // L_m >= L_{m-1}. If the previous prefix was solved exactly. we use it as an additional valid lower bound for the auxiliary BPP.
        const int knownLB = previousExactLB;

        GRBEnv env(true);
        env.set("OutputFlag", "0");
        env.start();
        GRBModel auxModel(env);
        auxModel.set(GRB_StringAttr_ModelName,"Auxiliary_BPP_" + std::to_string(m));

        vector<GRBVar> u(m);
        for (int b = 0; b < m; b++) {
            u[b] = auxModel.addVar(0.0,1.0,0.0,GRB_BINARY,"u_" + std::to_string(b));
        }
        vector<vector<GRBVar>> a(m,vector<GRBVar>(m));

        for (int i = 0; i < m; i++) {
            for (int b = 0; b < m; b++) {
                a[i][b] = auxModel.addVar(0.0,1.0,0.0,GRB_BINARY,"a_" + std::to_string(i) +"_" + std::to_string(b));
            }
        }
        for (int i = 0; i < m; i++) {
            GRBLinExpr assignment = 0.0;
            for (int b = 0; b < m; b++) {
                assignment += a[i][b];
            }
            auxModel.addConstr(assignment == 1.0);
        }
        for (int b = 0; b < m; b++) {
            GRBLinExpr binWeight = 0.0;
            for (int i = 0; i < m; i++) {
                binWeight += rw[i] * a[i][b];
            }
            auxModel.addConstr(binWeight <= W * u[b]);
        }
        for (int i = 0; i < m; i++) {
            for (int b = 0; b < m; b++) {
                auxModel.addConstr(a[i][b] <= u[b]);
            }
        }
        for (int b = 0; b < m-1; b++) {
            auxModel.addConstr(u[b] >= u[b + 1]);
        }

        GRBLinExpr usedBins = 0.0;
        for (int b = 0; b < m; b++) {
            usedBins += u[b];
        }
        auxModel.addConstr(usedBins >= knownLB);
        auxModel.setObjective(usedBins,GRB_MINIMIZE);
        auxModel.update();

        auxModel.getEnv().set(GRB_DoubleParam_MIPGap,0.0);
        auxModel.getEnv().set(GRB_IntParam_Threads,1);
        auxModel.getEnv().set(GRB_IntParam_Seed,1);
        auxModel.getEnv().set(GRB_DoubleParam_TimeLimit,inst.timeLimit);

        auxModel.optimize();
        result.timeByPrefix[m-1] =getCPUTime() - prefixStart;
        const int status =auxModel.get(GRB_IntAttr_Status);

        if (status == GRB_OPTIMAL) {
            result.exact[m-1] = true;
            result.L[m - 1] =static_cast<int>(std::llround(auxModel.get(GRB_DoubleAttr_ObjVal)));
            previousExactLB =result.L[m-1];
            cout<< "m = " << m<< " | L_m = " << result.L[m - 1]<< " | time = " << result.timeByPrefix[m - 1]<< " s\n";
        }
        else {
            result.exact[m-1] = false;
        }
    }
    result.totalTime =getCPUTime() - totalAuxStart;
    cout<< "Total prefix-BPP preprocessing time = "<< " seconds\n";

    return result;
}

struct GlobalBPPResult {
    int L = 0;
    bool exact = false;
    double time = 0.0;
    int status = 0;
};

// COMPUTE EXACT GLOBAL L_n
GlobalBPPResult computeGlobalBPPBound(const Instance& inst,const vector<int>& weights,const int knownLBFromPrefixes) {
    GlobalBPPResult result;
    const double start = getCPUTime();
    const int n = inst.n;
    const int W = inst.W;

    const int knownLB =knownLBFromPrefixes;
    cout << "AUXILIARY BPP FOR GLOBAL BIN-COUNT INEQUALITY\n";
    cout << "Known lower bound givrn to global BPP = "<< knownLB << '\n';
    GRBEnv env(true);
    env.set("OutputFlag", "0");
    env.start();

    GRBModel auxModel(env);

    auxModel.set(GRB_StringAttr_ModelName, "BPP_global_Ln");
    vector<GRBVar> u(n);
    for (int b = 0; b < n; b++) {
        u[b] = auxModel.addVar(0.0,1.0,0.0,GRB_BINARY,"u_" + std::to_string(b));
    }
    vector<vector<GRBVar>> a(n,vector<GRBVar>(n));
    for (int i = 0; i < n; i++) {
        for (int b = 0; b < n; b++) {
            a[i][b] = auxModel.addVar(0.0,1.0,0.0,GRB_BINARY,"a_" + std::to_string(i) +"_" + std::to_string(b));
        }
    }
    for (int i = 0; i < n; i++) {
        GRBLinExpr assignment = 0.0;
        for (int b = 0; b < n; b++) {
            assignment += a[i][b];
        }
        auxModel.addConstr(assignment == 1.0);
    }

    for (int b = 0; b < n; b++) {
        GRBLinExpr binWeight = 0.0;
        for (int i = 0; i < n; i++) {
            binWeight += weights[i] * a[i][b];
        }
        auxModel.addConstr(binWeight <= W * u[b]);
    }

    for (int i = 0; i < n; i++) {
        for (int b = 0; b < n; b++) {
            auxModel.addConstr(a[i][b] <= u[b]);
        }
    }
    for (int b = 0; b < n - 1; b++) {
        auxModel.addConstr(u[b] >= u[b+1]);
    }

    GRBLinExpr usedBins = 0.0;
    for (int b = 0; b < n; b++) {
        usedBins += u[b];
    }
    auxModel.addConstr(usedBins >= knownLB);
    auxModel.setObjective(usedBins,GRB_MINIMIZE);
    auxModel.update();

    auxModel.getEnv().set(GRB_DoubleParam_MIPGap,0.0);
    auxModel.getEnv().set(GRB_IntParam_Threads,1);
    auxModel.getEnv().set(GRB_IntParam_Seed,1);
    auxModel.getEnv().set(GRB_DoubleParam_TimeLimit,inst.timeLimit);
    auxModel.optimize();

    result.time =getCPUTime() - start;

    result.status =auxModel.get(GRB_IntAttr_Status);
    if (result.status == GRB_OPTIMAL) {
        result.exact = true;
        result.L =static_cast<int>(std::llround(auxModel.get(GRB_DoubleAttr_ObjVal)));
        cout << "Exact global L_n = "<< result.L<< '\n';
    }
    else {
        result.exact = false;
        cout << "Global auxiliary BPP was not solved to optimality.\n";
        cout << "Exact global L_n cut will be skipped.\n";
    }
    cout << "Global-BPP preprocessing time = "<< result.time<< " seconds\n";
    cout << "====================================================\n\n";
    return result;
}
// PATTERN STRUCTURE
struct SubsetPattern {
    vector<int> items;
    int weight = 0;
    double cost = 0.0;
};

// ENUMERATE EVERY FEASIBLE NONEMPTY SUBSET OF S
void enumerateSPatternsRec(const Instance& inst,const vector<int>& reorderedWeights,const vector<vector<int>>& reorderedD,const vector<int>& sItems,const int pos,const int currentWeight,const double currentInteraction,vector<int>& currentItems,vector<SubsetPattern>& patterns) {
    if (!currentItems.empty()) {
        SubsetPattern p;
        p.items =currentItems;
        p.weight =currentWeight;
        p.cost =static_cast<double>(inst.c)+ currentInteraction;
        patterns.push_back(p);
    }
    for (int r = pos;r < static_cast<int>(sItems.size());r++) {
        const int item =sItems[r];
        const int newWeight =currentWeight+ reorderedWeights[item];
        if (newWeight > inst.W) {
            continue;
        }
        double addedInteraction = 0.0;
        for (const int other : currentItems) {
            const int i =std::min(item, other);
            const int j =std::max(item, other);
            addedInteraction +=reorderedD[i][j];
        }
        currentItems.push_back(item);
        enumerateSPatternsRec(inst,reorderedWeights,reorderedD,sItems,r+1,newWeight,currentInteraction+ addedInteraction,currentItems,patterns);
        currentItems.pop_back();
    }
}
// GENERATE S-ONLY PATTERNS
vector<SubsetPattern> generateSPatterns(const Instance& inst,const vector<int>& reorderedWeights,const vector<vector<int>>& reorderedD,const vector<int>& sItems) {
    vector<SubsetPattern> patterns;
    vector<int> currentItems;
    currentItems.reserve(static_cast<int>(sItems.size()));
    enumerateSPatternsRec(inst,reorderedWeights,reorderedD,sItems,0,0,0.0,currentItems,patterns);
    return patterns;
}
// HYBRID S-PATTERN/GW MODEL

void H_GW_S_p_R_sign(Instance& inst, Solution& sol, const int ANCHORS) {
    try {
        const double totalStart = getCPUTime();
        const int n =inst.n;
        const int W =inst.W;
        const int c =inst.c;
        const vector<int>& w =inst.w;
        const vector<vector<int>>& d =inst.d;

        vector<long long> interactionStrength(n, 0);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    continue;
                }
                const int a =std::min(i, j);
                const int b =std::max(i, j);
                interactionStrength[i] +=std::abs(d[a][b]);
            }
        }
        // Order from weakest interaction to strongest.
        vector<int> interactionOrder(n);
        for (int i = 0; i < n; i++) {
            interactionOrder[i] = i;
        }
        std::sort(interactionOrder.begin(),interactionOrder.end(),[&](const int a, const int b) {
                if (interactionStrength[a]!= interactionStrength[b]) {
                    return interactionStrength[a]< interactionStrength[b];
                }
                return a < b;
            }
        );
        cout << "\nItems ordered by interaction strength:\n";
        for (int r = 0; r < n; r++) {
            const int item =interactionOrder[r];
            cout<< "Item "<< item<< " -> interaction strength = "<< interactionStrength[item];
            if (r < ANCHORS) {
                cout << "  [SELECTED]";
            }
            cout << '\n';
        }
        cout << '\n';
        // Put ANCHORS weakest-interaction items into S.
        vector<bool> originalInS(n, false);
        for (int r = 0; r < ANCHORS; ++r) {
            originalInS[interactionOrder[r]] = true;
        }
        // INTERNAL REORDERING:I2-DEGREE-DESC
        // Sbar first, ordered by descending number of nonzero interactions. S afterwards
        vector<int>nonzeroDegree(n, 0);
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (d[i][j] != 0) {
                    nonzeroDegree[i]++;
                    nonzeroDegree[j]++;
                }
            }
        }
        // collect Sbar.
        vector<int> sbarItems;
        sbarItems.reserve(n - ANCHORS);
        for (int old = 0; old < n; ++old) {
            if (!originalInS[old]) {
                sbarItems.push_back(old);
            }
        }
        // put most-connected Sbar items first.
        std::sort(sbarItems.begin(),sbarItems.end(),[&](const int a, const int b) {
                if (nonzeroDegree[a]!= nonzeroDegree[b]) {
                    return nonzeroDegree[a] > nonzeroDegree[b];
                }
                return a < b;
            }
        );
        // Build reordered index
        vector<int> newToOld;
        newToOld.reserve(n);
        // Sbar first.
        for (const int old : sbarItems) {
            newToOld.push_back(old);
        }
        // S afterwards.
        for (int old = 0; old < n; old++) {
            if (originalInS[old]) {

                newToOld.push_back(old);
            }
        }
        vector<int>oldToNew(n, -1);
        for (int i = 0; i < n; i++) {
            oldToNew[newToOld[i]] = i;
        }
        const int nSbar =n - ANCHORS;
        // Reordered weights and interactions.
        vector<int>rw(n, 0);
        vector<vector<int>>rd(n,vector<int>(n, 0));
        for (int i = 0; i < n; i++) {
            rw[i] =w[newToOld[i]];
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    rd[i][j] = 0;
                }
                else {
                    const int oldI =newToOld[i];
                    const int oldJ =newToOld[j];
                    const int a =std::min(oldI,oldJ);
                    const int b =std::max(oldI,oldJ);
                    rd[i][j] =d[a][b];
                }
            }
        }
        // COMPUTE EXACT BPP VALUES L_m FOR Sbar 
        const PrefixBPPResult prefixBPP =computePrefixBPPBounds(inst,rw,nSbar);
        int largestExactPrefixLB = 0;
        for (int m = 1; m <= nSbar; m++) {
            if (prefixBPP.exact[m-1]) {
                largestExactPrefixLB =std::max(largestExactPrefixLB,prefixBPP.L[m-1]);
            }
        }
        const GlobalBPPResult globalBPP = computeGlobalBPPBound(inst,rw,largestExactPrefixLB);
        sol.timePrepr = prefixBPP.totalTime + globalBPP.time;

        // S items occupy last ANCHORS positions.
        vector<int> sItems;
        sItems.reserve(ANCHORS);
        for (int i = nSbar;i < n;i++) {
            sItems.push_back(i);
        }
        cout<< "\nHybrid S-pattern / GW preprocessing:\n";
        cout<< "ANCHORS = "<< ANCHORS<< '\n';
        cout<< "S (original item indices): ";
        if (ANCHORS == 0) {
            cout << "none";
        }
        else {
            for (int old = 0;old < n;old++) {
                if (originalInS[old]) {
                    cout<< old<< ' ';          
                }
            }
        }
        cout << '\n';
        cout<< "Sbar (original item indices): ";
        if (nSbar == 0) {
            cout << "none";
        }
        else {
            for (int old = 0;old <n;old++) {
                if (!originalInS[old]) {
                    cout<< old<< ' ';
                }
            }
        }
        cout << '\n';
        cout<< "Internal order new->old: ";
        for (int i = 0; i < n; i++) {
            cout<< i<< "->"<< newToOld[i]<< ' ';
        }
        cout << "\n\n";

        // PATTERN ENUMERATION OVER S ONLY
        const double patternStart =getCPUTime();

        vector<SubsetPattern> patterns =generateSPatterns(inst,rw,rd,sItems);
        const double patternTime =getCPUTime()- patternStart;
        const int numPatterns =static_cast<int>(patterns.size());

        sol.timePatternGen =patternTime;
        sol.nbPatterns =numPatterns;

        vector<vector<int>>patternIdsByItem(n);
        for (int p = 0;p < numPatterns;p++) {
            for (const int item:patterns[p].items) {
                patternIdsByItem[item].push_back(p);
            }
        }
        cout<< "Enumerated S-only patterns = "<< numPatterns<< '\n';
        cout<< "Pattern enumeration time (seconds) = "<< patternTime<< "\n\n";
        // GUROBI MODEL
        GRBEnv env(true);
        env.start();
        GRBModel model(env);
        model.set(GRB_StringAttr_ModelName, "F_S_pattern_GW_S_plus_sign_reduced");

        // Pattern variables.
        vector<GRBVar>lambda(numPatterns);
        for (int p = 0;p < numPatterns;p++) {
            lambda[p] =model.addVar(0.0,1.0,0.0,GRB_BINARY,"lambda_"+ std::to_string(p));
        }
        // x[i][k] exists iff k is in Sbar and k <= i
        vector<vector<GRBVar>> x(n,vector<GRBVar>(nSbar));
        vector<vector<bool>> hasX(n,vector<bool>(nSbar,false));

        int nbXVars = 0;
        for (int i = 0; i < n; i++) {
            const int kMax =std::min(i,nSbar - 1);
            for (int k = 0;k <= kMax;k++) {
                x[i][k] =model.addVar(0.0,1.0,0.0,GRB_BINARY,"x_"+ std::to_string(i)+ "_"+ std::to_string(k));
                hasX[i][k] =true;
                nbXVars++;
            }
        }
        // g[i][j][k] exists only for k < i < j, k in Sbar,rd[i][j] != 0
        vector<vector<vector<GRBVar>>> g(n,vector<vector<GRBVar>>(n,vector<GRBVar>(nSbar)));
        vector<vector<vector<bool>>>hasG(n,vector<vector<bool>>(n,vector<bool>(nSbar,false)));

        int nbPositivePairs = 0;
        int nbNegativePairs = 0;
        int nbZeroPairs = 0;
        int nbGVars = 0;

        for (int i = 0;i < n-1;i++) {
            for (int j = i+1;j < n;j++) {
                if (rd[i][j] > 0) {
                    nbPositivePairs++;
                }
                else if (rd[i][j] < 0) {
                    nbNegativePairs++;
                }
                else {
                    nbZeroPairs++;
                    continue;
                }
                const int kMax =std::min(i - 1,nSbar - 1);
                for (int k = 0;k <= kMax;k++) {
                    g[i][j][k] =model.addVar(0.0,1.0,0.0,GRB_CONTINUOUS,"g_"+ std::to_string(i)+ "_"+ std::to_string(j)+ "_"+ std::to_string(k));
                    hasG[i][j][k] =true;
                    nbGVars++;
                }
            }
        }

        for (int i = 0; i < n; i++) {
            GRBLinExpr assignment =0.0;
            if (i >= nSbar) {
                for (const int p :patternIdsByItem[i]) {
                    assignment +=lambda[p];
                }
            }
            const int kMax =std::min(i,nSbar-1);
            for (int k = 0;k <= kMax;k++) {
                assignment +=x[i][k];
            }
            model.addConstr(assignment == 1.0);
        }
        for (int k = 0;k < nSbar;k++) {
            GRBLinExpr binWeight =0.0;
            for (int i = k;i < n;i++) {
                binWeight +=rw[i]* x[i][k];
            }
            model.addConstr(binWeight<=W * x[k][k]);
        }
        for (int i = 0; i < n; i++) {
            const int kMax =std::min(i - 1,nSbar - 1);
            for (int k = 0;k <= kMax;k++) {
                model.addConstr(x[i][k]<=x[k][k]);
            }
        }

        int nbCumulativeCuts = 0;
        int lastAddedM = 0;
        int lastAddedL = 0;

        for (int m = 1; m <= nSbar; m++) {
            if (!prefixBPP.exact[m - 1]) {
                cout<< "Skipping cumulative cut for m = "<< m<< " because exact L_m was not obtained.\n";
                continue;
            }
            const int Lm =prefixBPP.L[m - 1];
            // m=1 is already enforced by the assignment x[0][0] = 1.
            if (m == 1) {
                lastAddedM = 1;
                lastAddedL = Lm;
                continue;
            }
            // If L_m has not increased to the strongest exact cut already available the new inequality is redundant
            if (Lm <= lastAddedL) {
                continue;
            }
            GRBLinExpr cumulativeBins = 0.0;
            for (int k = 0; k < m; k++) {
                cumulativeBins += x[k][k];
            }
            model.addConstr(cumulativeBins >= Lm,"cumulative_bin_count_" +std::to_string(m));

            nbCumulativeCuts++;
            lastAddedM = m;
            lastAddedL = Lm;
            cout<< "Added cumulative cut for m = "<< m<< ": sum_{k=0}^{"<< m - 1<< "} x[k][k] >= "<< Lm<< '\n';
        }

        cout<< "Number of cumulative bin-count cuts added = "<< nbCumulativeCuts<< "\n\n";
        // GLOBAL BIN-COUNT VALID INEQUALITY
        bool globalCutAdded = false;

        if (globalBPP.exact) {
            GRBLinExpr globalBinCount = 0.0;
            for (int p = 0; p < numPatterns; p++) {
                globalBinCount += lambda[p];
            }
            for (int k = 0; k < nSbar; k++) {
                globalBinCount += x[k][k];
            }
            model.addConstr(globalBinCount >= globalBPP.L);
            globalCutAdded = true;
            cout<< "Added global bin-count cut: "<< "sum(lambda_p) + sum(x_kk) >= "<< globalBPP.L<< '\n';
        }
        else {
            cout<< "Global bin-count cut skipped because exact L_n was not obtained.\n";
        }
        cout<< "Global bin-count cut added = "<< (globalCutAdded ? 1 : 0)<< "\n\n";

        for (int i = 0;i < n - 1;i++) {
            for (int j = i + 1;j < n;j++) {
                if (rd[i][j] <= 0) {
                    continue;
                }
                const int kMax =std::min(i - 1,nSbar - 1);
                for (int k = 0;k <= kMax;k++) {
                    if (!hasG[i][j][k]) {
                        continue;
                    }
                    model.addConstr(x[i][k]+ x[j][k]- x[k][k]<=g[i][j][k]);
                }
            }
        }
        for (int i = 0;i < n - 1;i++) {
            for (int j = i + 1;j < n;j++) {
                if (rd[i][j] >= 0) {
                    continue;
                }
                const int kMax =std::min(i - 1,nSbar - 1);
                for (int k = 0;k <= kMax;k++) {
                    if (!hasG[i][j][k]) {
                        continue;
                    }
                    model.addConstr(g[i][j][k]<=x[i][k]);
                    model.addConstr(g[i][j][k]<=x[j][k]);
                }
            }
        }
        GRBLinExpr objective =0.0;
        // S-only pattern costs
        for (int p = 0;p < numPatterns;p++) {
            objective +=patterns[p].cost* lambda[p];
        }

        // Eliminated  terms  g[i][j][i] = x[j][i] for i in Sbar.
        for (int i = 0;i < nSbar;i++) {
            for (int j = i + 1;j < n;j++) {
                if (rd[i][j] == 0) {
                    continue;
                }
                objective +=rd[i][j]* x[j][i];
            }
        }
        // Remaining g terms
        for (int i = 0;i < n - 1;i++) {
            for (int j = i + 1;j < n;j++) {
                if (rd[i][j] == 0) {
                    continue;
                }
                const int kMax =std::min(i - 1,nSbar - 1);
                for (int k = 0;k <= kMax;k++) {
                    if (!hasG[i][j][k]) {
                        continue;
                    }
                    objective +=rd[i][j]* g[i][j][k];
                }
            }
        }
        // GW fixed costs
        for (int k = 0;k < nSbar;k++) {
            objective +=c * x[k][k];
        }
        model.setObjective(objective,GRB_MINIMIZE);
        model.update();
        // INTEGER MODEL SIZE
        cout<< "GW statistics after internal reordering:\n";
        cout<< "GW representatives |Sbar| = "<< nSbar<< '\n';
        cout<< "x variables = "<< nbXVars<< '\n';
        cout<< "Positive item pairs = "<< nbPositivePairs<< '\n';
        cout<< "Negative item pairs = "<< nbNegativePairs<< '\n';
        cout<< "Zero item pairs     = "<< nbZeroPairs<< '\n';
        cout<< "Remaining g variables = "<< nbGVars<< '\n';

        sol.nbVars =model.get(GRB_IntAttr_NumVars);
        sol.nbConstrs =model.get(GRB_IntAttr_NumConstrs);
        sol.nbNZs =model.get(GRB_IntAttr_NumNZs);
        cout<< "Model size:\n";
        cout<< "Variables = "<< sol.nbVars<< '\n';
        cout<< "Constraints = "<< sol.nbConstrs<< '\n';
        cout<< "Nonzeros = "<< sol.nbNZs<< '\n';

        model.getEnv().set(GRB_DoubleParam_MIPGap,0.0);
        model.getEnv().set(GRB_IntParam_Threads,1);
        model.getEnv().set(GRB_DoubleParam_TimeLimit,inst.timeLimit);
        model.getEnv().set(GRB_IntParam_Seed,1);

        model.optimize();

        sol.timeT =getCPUTime()- totalStart;

        const int status =model.get(GRB_IntAttr_Status);
        const int solutionCount =model.get(GRB_IntAttr_SolCount);

        if (status == GRB_OPTIMAL) {
            sol.opt = 1;
            sol.valOpt =model.get(GRB_DoubleAttr_ObjVal);
            sol.valBest =model.get(GRB_DoubleAttr_ObjVal);
            sol.valBound =model.get(GRB_DoubleAttr_ObjBound);
            cout<< "Optimal solution found. Objective = "<< sol.valOpt<< '\n';
        }

        else if (status == GRB_TIME_LIMIT) {
            sol.opt = 0;
            sol.valBound =model.get(GRB_DoubleAttr_ObjBound);
            cout<< "Time limit reached.\n";
            cout<< "Best bound = "<< sol.valBound<< '\n';
            if (solutionCount > 0) {
                sol.valBest =model.get(GRB_DoubleAttr_ObjVal);
                cout<< "Best solution = "<< sol.valBest<< '\n';
            }
            else {
                cout<< "No feasible solution was found.\n";
            }
        }
        else if (status == GRB_INFEASIBLE) {
            sol.opt = 0;
            cout<< "The model is infeasible.\n";
        }
        else if (status == GRB_UNBOUNDED) {
            sol.opt = 0;
            cout<< "The model is unbounded.\n";
        }
        else if (status == GRB_INF_OR_UNBD) {
            sol.opt = 0;
            cout<< "The model is infeasible or unbounded.\n";
        }
        else {
            sol.opt = 0;
            cout<< "Optimization ended with status "<< status<< ".\n";
        }

        //CONTINUOUS LINEAR RELAXATION
        const double relaxStart =getCPUTime();
        GRBModel relaxModel =model.relax();
        relaxModel.set(GRB_StringAttr_ModelName,"F_S_pattern_GW_S_plus_sign_reduced_relaxation");
        relaxModel.getEnv().set(GRB_DoubleParam_TimeLimit,inst.timeLimit);
        relaxModel.getEnv().set(GRB_IntParam_Threads,1);
        relaxModel.optimize();
        sol.timeR =getCPUTime()- relaxStart;
        const int statusR =relaxModel.get(GRB_IntAttr_Status);

        if (statusR == GRB_OPTIMAL) {
            sol.optR =1;
            sol.valRel =relaxModel.get(GRB_DoubleAttr_ObjVal);
            cout<< "\nContinuous linear relaxation value = "<< sol.valRel<< '\n';

            cout<< "Relaxation solve time (seconds) = "<< sol.timeR<< '\n';
        }
        else {
            sol.optR =0;
        }
    }

    catch (const GRBException& e) {
        std::cerr<< "Gurobi error code = "<< e.getErrorCode()<< '\n'<< e.getMessage()<< '\n';
    }
    catch (const std::exception& e) {
        std::cerr<< "Error: "<< e.what()<< '\n';
    }
    catch (...) {
        std::cerr<< "Unknown exception during optimization.\n";
    }
}