#pragma once

#include "global.h"
#include "gurobi_c++.h"

class TdmDB;
class XdrVar;
class TimingGraph;
class Troncon;

class TdmLpSolver {
public:
    TdmLpSolver(bool useLP);
    void solve();

protected:
    int _nChoice = 1600 / 8 + 1;
    int _nVar;

    TimingGraph *_timingGraph;
    vector<XdrVar *> &_optXdrVars;

    vector<GRBVar> _xdrVar;
    vector<GRBVar> _gateVar;
    vector<GRBVar> _usageVar;
    vector<vector<GRBVar>> _extraXdrVar;

    GRBEnv _env;
    GRBModel _model;

    void getResult();

    void addOneChoiceConstraint();
    void addLimitConstraint();
    void addTimingEdgeConstraint();

private:
    bool _useLP;
    map<Troncon *, vector<int>> _tronconToXdrVar;
    const int _timeLimit = 10000;
    const int _moreChoiceIdx = 2;

    void init();
    void addExactLimitConstraint();
    void genILPInitSol();
    int getXdrChoice(int extraXdrVarIdx);
};
