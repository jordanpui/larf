#pragma once

#include "global.h"

class TdmDB;
class Edge;
class Troncon;
class XdrVar;
class TimingGraph;

class Memorization;
class TdmLegalize;
class WireData;

class TdmLegalize {
public:
    TdmLegalize();
    void solve(bool writeDB = true, vector<double> *result = NULL);

private:
    TimingGraph *_timingGraph;
    vector<XdrVar *> &_optXdrVars;

    map<Troncon *, vector<int>> _tronconToXdrVar;
    map<Troncon *, int> _tronconToIdx;
    map<XdrVar *, int> _varToOptIdx;

    vector<WireData> _wireData;
    Setting::LgMethod _flow;

    double legalizeTroncon(const WireData &data, Memorization &memorization) const;

    void updateWireData();
};

class OneWireData {
public:
    OneWireData(XdrVar *var) : _var(var) {}

    bool isForward() const;
    void setVal(double val);
    double getVal() const;

    double _k;
    double _b;
    XdrVar *_var;
};

class WireData {
public:
    Troncon *_troncon;
    vector<OneWireData *> _vars;
    unsigned _numForwardVars;
    double _totDisp;
    double _maxDisp;

    int calcEndIdx(unsigned idx, int choice) const;
    int calcEndIdx(unsigned idx, int choice, const vector<pair<int, int>> &choiceRanges) const;
    int getMinWireRequire(unsigned idx) const;
    int getMinWireRequire(unsigned idx, const vector<pair<int, int>> &choiceRanges) const;
    void calcDisp(const Memorization &memorization);
    void sortVars(Setting::LgMethod flow);

    double legalizeTronconDisp(unsigned idx, int p, Memorization &memorization) const;
    double legalizeTronconDispPrune(unsigned idx, int p, Memorization &memorization) const;
    double legalizeTronconDisp(unsigned idx,
                               int p,
                               Memorization &memorization,
                               const vector<pair<int, int>> &choiceRanges) const;
    double legalizeTronconDispPrune(unsigned idx,
                                    int p,
                                    Memorization &memorization,
                                    const vector<pair<int, int>> &choiceRanges,
                                    const vector<int> &minWire) const;
    double legalizeTronconMaxDisp(Memorization &memorization) const;

    void recoverSol(const Memorization &memorization);
    void dumpSol(vector<double> &result, const map<XdrVar *, int> &varToOptIdx, const Memorization &memorization) const;

private:
    void sortByDisp();
};

class Memorization {
public:
    Memorization(int n, int p) {
        _bestEndIdx.resize(n, vector<int>(p, -1));
        _bestChoice.resize(n, vector<int>(p, -1));
        _bestCost.resize(n, vector<double>(p, DBL_MAX));
        _optimized.resize(n, vector<bool>(p, false));
    }

    void saveBest(int idx, int p, double cost, int choice, int endIdx) {
        _bestEndIdx[idx][p - 1] = endIdx;
        _bestCost[idx][p - 1] = cost;
        _bestChoice[idx][p - 1] = choice;
        _optimized[idx][p - 1] = true;
    }

    double getBestCost(int idx, int p) const { return _bestCost[idx][p - 1]; }
    int getBestEndIdx(int idx, int p) const { return _bestEndIdx[idx][p - 1]; }
    int getBestChoice(int idx, int p) const { return _bestChoice[idx][p - 1]; }
    bool isOptimized(int idx, int p) const { return _optimized[idx][p - 1]; }

private:
    vector<vector<int>> _bestEndIdx;
    vector<vector<int>> _bestChoice;
    vector<vector<double>> _bestCost;
    vector<vector<bool>> _optimized;
};
