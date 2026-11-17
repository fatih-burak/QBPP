#include "helper_functions.h"
#include <filesystem>
#include <stdexcept>


double getCPUTime() {
    return static_cast<double>(clock()) / CLOCKS_PER_SEC;
}

Instance loadInstance(const string& path) {
    Instance inst;

    ifstream in(path);
    if (!in) {
        cerr << "Error opening file: " << path << endl;
        exit(1);
    }

    in >> inst.n >> inst.W >> inst.c;
    inst.B = inst.n;

    inst.w.resize(inst.n);
    for (int i = 0; i < inst.n; i++) {
        in >> inst.w[i];
    }

    inst.d.assign(inst.n, vector<int>(inst.n, 0.0));
    for (int i = 0; i < inst.n; i++) {
        for (int j = 0; j < inst.n; j++) {
            in >> inst.d[i][j];
        }
    }

    in.close();
    return inst;
}

namespace fs = std::filesystem;

static bool containsChar(const vector<char>& values, char x) {
    return values.empty() || find(values.begin(), values.end(), x) != values.end();
}

static bool containsInt(const vector<int>& values, int x) {
    return values.empty() || find(values.begin(), values.end(), x) != values.end();
}

static bool matchesChagasFilter(const ChagasInstanceInfo& info, const ChagasFilter& filter) {
    return containsChar(filter.signals, info.signal)
        && containsInt(filter.items, info.items)
        && containsInt(filter.densities, info.density)
        && containsInt(filter.mus, info.mu)
        && containsInt(filter.ids, info.id);
}

bool parseChagasInstanceName(const string& name, ChagasInstanceInfo& info) {
    const string prefix = "QBPP_HJ";

    if (name.rfind(prefix, 0) != 0) {
        return false;
    }

    size_t pos = prefix.size();
    if (pos >= name.size()) {
        return false;
    }

    char signal = name[pos];
    if (signal != '+' && signal != '-' && signal != '*') {
        return false;
    }
    ++pos;

    if (pos >= name.size() || name[pos] != '_') {
        return false;
    }

    string rest = name.substr(pos + 1);
    vector<string> parts;
    string token;
    stringstream ss(rest);

    while (getline(ss, token, '_')) {
        parts.push_back(token);
    }

    if (parts.size() != 4) {
        return false;
    }

    try {
        info.name = name;
        info.signal = signal;
        info.items = stoi(parts[0]);
        info.density = stoi(parts[1]);
        info.mu = stoi(parts[2]);
        info.id = stoi(parts[3]);
    }
    catch (...) {
        return false;
    }

    return true;
}

static Instance readChagasInstanceFile(const string& filePath, string& instanceName) {
    ifstream in(filePath);
    if (!in) throw runtime_error("Could not open Chagas instance file: " + filePath);

    Instance inst;
    if (!(in >> instanceName)) throw runtime_error("Could not read instance name from: " + filePath);

    if (!(in >> inst.n >> inst.W >> inst.c)) {
        throw runtime_error("Could not read n, W, c from: " + filePath);
    }

    inst.B = inst.n;
    inst.w.resize(inst.n);
    for (int i = 0; i < inst.n; i++) {
        if (!(in >> inst.w[i])) {
            throw runtime_error("Could not read weight " + to_string(i) + " from: " + filePath);
        }
    }

    inst.d.assign(inst.n, vector<int>(inst.n, 0.0));
    for (int i = 0; i < inst.n; i++) {
        for (int j = 0; j < inst.n; j++) {
            if (!(in >> inst.d[i][j])) {
                throw runtime_error("Could not read d[" + to_string(i) + "][" + to_string(j) + "] from: " + filePath);
            }
        }
    }
    return inst;
}

NamedInstance loadChagasInstance(const string& instanceFolderOrFilePath) {
    fs::path inputPath(instanceFolderOrFilePath);
    fs::path filePath;

    if (fs::is_directory(inputPath)) {
        string folderName = inputPath.filename().string();
        fs::path expectedFile = inputPath / (folderName + ".in");

        if (fs::exists(expectedFile) && fs::is_regular_file(expectedFile)) {
            filePath = expectedFile;
        }
        else {
            for (const auto& entry : fs::directory_iterator(inputPath)) {
                if (entry.is_regular_file() && entry.path().extension() == ".in") {
                    filePath = entry.path();
                    break;
                }
            }
        }
        if (filePath.empty()) {
            throw runtime_error("No .in file found in folder: " + inputPath.string());
        }
    }
    else filePath = inputPath;

    string instanceName;
    Instance inst = readChagasInstanceFile(filePath.string(), instanceName);

    ChagasInstanceInfo info;
    if (!parseChagasInstanceName(instanceName, info)) {
        throw runtime_error("Invalid Chagas instance name: " + instanceName);
    }

    NamedInstance named;
    named.name = instanceName;
    named.inst = inst;
    named.path = filePath.string();
    named.info = info;

    return named;
}

vector<NamedInstance> loadChagasInstancesFromFolder(const string& rootFolder, const ChagasFilter& filter) {
    fs::path root(rootFolder);
    if (!fs::exists(root)) {
        throw runtime_error("Root folder does not exist: " + rootFolder);
    }
    if (!fs::is_directory(root)) {
        throw runtime_error("Path is not a folder: " + rootFolder);
    }
    vector<NamedInstance> instances;
    for (const auto& entry : fs::directory_iterator(root)) {
        if (!entry.is_directory()) {
            continue;
        }
        string folderName = entry.path().filename().string();

        ChagasInstanceInfo folderInfo;
        if (!parseChagasInstanceName(folderName, folderInfo)) {
            cerr << "Skipping folder with invalid Chagas name: " << folderName << "\n";
            continue;
        }

        if (!matchesChagasFilter(folderInfo, filter)) continue;

        NamedInstance named = loadChagasInstance(entry.path().string());
        instances.push_back(named);
    }

    sort(instances.begin(), instances.end(), [](const NamedInstance& a, const NamedInstance& b) {
        if (a.info.signal != b.info.signal) return a.info.signal < b.info.signal;
        if (a.info.items != b.info.items) return a.info.items < b.info.items;
        if (a.info.density != b.info.density) return a.info.density < b.info.density;
        if (a.info.mu != b.info.mu) return a.info.mu < b.info.mu;
        return a.info.id < b.info.id;
        });

    return instances;
}

void sort_decreasing(vector<double>& w, vector<vector<double>>& d) {
    int n = static_cast<int>(w.size());

    vector<int> idx(n);
    for (int i = 0; i < n; i++) {
        idx[i] = i;
    }

    sort(idx.begin(), idx.end(), [&](int a, int b) {return w[a] > w[b];});

    vector<double> w_sorted(n);
    for (int i = 0; i < n; i++) {
        w_sorted[i] = w[idx[i]];
    }

    vector<vector<double>> d_sorted(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            d_sorted[i][j] = d[idx[i]][idx[j]];
        }
    }
    w = move(w_sorted);
    d = move(d_sorted);
}

void printInstance(const Instance& inst) {
    cout << "n = " << inst.n << ", W = " << inst.W << ", c = " << inst.c << "\n";
    cout << "weights:\n";
    for (double wi : inst.w) {
        cout << wi << "\t";
    }
    cout << "\n";

    cout << "interaction matrix:\n";
    for (int i = 0; i < inst.n; i++) {
        for (int j = 0; j < inst.n; j++) {
            cout << inst.d[i][j] << "\t";
        }
        cout << "\n";
    }
    cout << "\n";
}

void printSolutionSummary(const Solution& sol) {
    cout << "\n===== SOLUTION SUMMARY =====\n";
    cout << "Model      : " << sol.name << "\n";
    cout << "Optimal    : " << sol.opt << "\n";
    cout << "Obj value  : " << sol.valOpt << "\n";
    cout << "Best value : " << sol.valBest << "\n";
    cout << "Best bound : " << sol.valBound << "\n";
    cout << "Solve time : " << sol.timeT << "\n";
    cout << "Relax opt  : " << sol.optR << "\n";
    cout << "Relax val  : " << sol.valRel << "\n";
    cout << "Relax time : " << sol.timeR << "\n";
    //cout << "Root rel   : " << sol.valRootRel << "\n";

void printInfo(const string& outputPath, const Solution& sol) {
    ofstream file(outputPath.c_str(), ios::out | ios::app);

    file << fixed << setprecision(2)
        << sol.name << "\t"  << sol.opt << "\t"<< sol.timeT << "\t"  << sol.valOpt << "\t" << sol.valBest << "\t" << sol.valBound << "\t"<< sol.optR << "\t" << sol.timeR << "\t" << sol.valRel << "\t"
        << sol.nbVars << "\t"<< sol.nbConstrs << "\t" << sol.nbNZs<< endl;

    file.close();
}