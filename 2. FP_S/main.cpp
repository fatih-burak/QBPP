#include "main.h"
#include <iostream>
#include <string>
using namespace std;

int main(int argc, char* argv[]) {
    try {
        const string instancePath = argv[1];
        const string resultsPath = argv[2];

        double timeLimit = 3600.0;

        NamedInstance namedInstance = loadChagasInstance(instancePath);

        Instance inst = namedInstance.inst;
        inst.timeLimit = timeLimit;

        Solution sol(namedInstance.name);

        cout << "Instance   : " << namedInstance.name << '\n';
        cout << "Input path : " << namedInstance.path << '\n';
        cout << "Time limit : " << inst.timeLimit << " seconds\n";


        FP_S_p_R(inst, sol);

        printSolutionSummary(sol);
        printInfo(resultsPath, sol);
        return 0;
    }
    catch (const GRBException& e) {
        cerr << "Gurobi error code = "<< e.getErrorCode() << '\n';
        cerr << e.getMessage() << '\n';
        return 1;
    }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << '\n';
        return 1;
    }
    catch (...) {
        cerr << "Unknown exception during execution.\n";
        return 1;
    }
}