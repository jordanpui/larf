#pragma once

#include "global.h"
#include "gurobi_c++.h"

class TdmDB;
class XdrVar;
class TimingGraph;
class Troncon;

class TdmRefineLP {
public:
    TdmRefineLP(bool genContSol);
    static void printUsageError();
    void solve();

private:
    vector<vector<GRBVar>> _xdrVar;
    vector<vector<GRBVar>> _extraXdrVar;
    vector<GRBVar> _gateVar;
    vector<GRBVar> _usageVar;

    GRBEnv _env;
    GRBModel _model;

    TimingGraph *_timingGraph;
    vector<XdrVar *> &_optXdrVars;
    map<Troncon *, vector<int>> _tronconToXdrVar;
    int _nVar;

    vector<pair<int, int>> _choiceRanges;

    int getChoiceNum(int idx) const { return _choiceRanges[idx].second - _choiceRanges[idx].first + 1; }
    GRBVar &getXdrVar(int idx, int choiceIdx) { return _xdrVar[idx][choiceIdx - _choiceRanges[idx].first]; }
    int getXdrChoice(int extraXdrVarIdx);
    bool withinXdrChoiceRange(int idx, int choice);

    // TODO: will have bug if _choiceBnd<1
    // TODO: now adding some useless variable
    const int _timeLimit = 10000;
    const int _choiceBnd = 8;
    const int _moreChoiceIdx = 2;

    bool _genContSol;

    void init();
    void genILPInitSol();
    void getResult();

    void addOneChoiceConstraint();
    void addLimitConstraint();
    void addExactLimitConstraint();
    void addTimingEdgeConstraint();
};
