#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

#include "gurobi_c++.h"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <regex>
#include <functional>


using namespace std;

struct Instance
{
    int n = 0;
    int B = 0;                    // number of bins
    int W = 0;
    int c = 0;                    // fixed cost per used bin
    vector<vector<double>> d;     // interaction matrix
    vector<double> w;             // weights
    double timeLimit = 1200.0;    // time limit in seconds
};

struct Solution
{
    string name = "instance";
    int opt = 0;
    int optR = 0;

    double valOpt = 0.0;
    double valRel = 0.0;
    double valRootRel = 0.0;
    double valBest = 0.0;
    double valBound = 0.0;

    double timeT = 0.0;
    double timeR = 0.0;
    double timePatternGen = 0.0;
    int nbPatterns = 0;

    // Original integer model size
    int nbVars = 0;
    int nbConstrs = 0;
    int nbNZs = 0;

    Solution() = default;
    Solution(const string& name_val) : name(name_val) {}
};

class MessageCallback : public GRBCallback {
private:
    Solution* solData;

public:
    explicit MessageCallback(Solution* s);

protected:
    void callback() override;
};


Instance loadInstance(const string& path);
void sort_decreasing(vector<double>& w, vector<vector<double>>& d);
void printInfo(const string& outputPath, const Solution& sol);

struct ChagasInstanceInfo {
    string name;
    char signal = '?';   // '+', '-', '*'
    int items = 0;       // 25, 30, 35, 40, 45
    int density = 0;     // 25, 50, 75
    int mu = 0;          // 6, 10, 20
    int id = 0;          // 1, 2, 3, 4, 5
};

struct ChagasFilter {
    vector<char> signals;     // empty means all
    vector<int> items;        // empty means all
    vector<int> densities;    // empty means all
    vector<int> mus;          // empty means all
    vector<int> ids;          // empty means all
};

struct NamedInstance {
    string name;
    Instance inst;
    string path;
    ChagasInstanceInfo info;
};

bool parseChagasInstanceName(const string& name, ChagasInstanceInfo& info);

NamedInstance loadChagasInstance(const string& instanceFolderOrFilePath);

vector<NamedInstance> loadChagasInstancesFromFolder(const string& rootFolder, const ChagasFilter& filter = ChagasFilter());

// Other helpers
void printInstance(const Instance& inst);
void printSolutionSummary(const Solution& sol);

// CPU timer 
double getCPUTime();

#endif