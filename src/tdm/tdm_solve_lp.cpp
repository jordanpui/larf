#include "tdm_solve_lp.h"
#include "tdm_db.h"
#include "tdm_net.h"
#include "timing_graph.h"

void TdmLpSolver::addOneChoiceConstraint() {
    // only one choice
    for (int i = 0; i < _nVar; i++) {
        GRBLinExpr expr = GRBLinExpr();
        for (int j = 0; j < _nChoice; j++) expr += _xdrVar[i * _nChoice + j];
        if (_useLP && _moreChoiceIdx > 0) {
            for (int j = 0, sz = _extraXdrVar[i].size(); j < sz; j++) {
                if (getXdrChoice(j) <= TdmDB::getXdrChoice(_nChoice - 1)) expr += _extraXdrVar[i][j];
            }
        }

        _model.addConstr(expr == 1);
    }
}

void TdmLpSolver::addTimingEdgeConstraint() {
    vector<int> xdrToLPIdx(_timingGraph->getNumNodes(), -1);
    for (int i = 0; i < _nVar; i++) xdrToLPIdx[_optXdrVars[i]->_id] = i;

    _model.addConstr(_gateVar[_timingGraph->getSource()->_id] == 0);

    for (int i = 0, gateVarNum = _timingGraph->getNumNodes(); i < gateVarNum; i++) {
        for (auto edge : _timingGraph->getNode(i)->_drivers) {
            GRBLinExpr expr;
            expr += _gateVar[i];
            expr -= _gateVar[edge->_driver->_id];
            if (edge->_net && edge->_net->isInterNet()) {
                int LPIdx = xdrToLPIdx[edge->_net->_xdrVar->_id];
                if (LPIdx != -1) {
                    for (int j = 0; j < _nChoice; j++) {
                        expr -= _xdrVar[LPIdx * _nChoice + j] * edge->_tdmCoef * TdmDB::getXdrChoice(j);
                    }
                    if (_useLP && _moreChoiceIdx > 0) {
                        for (int j = 0, sz = _extraXdrVar[LPIdx].size(); j < sz; j++) {
                            if (getXdrChoice(j) <= TdmDB::getXdrChoice(_nChoice - 1))
                                expr -= edge->_tdmCoef * getXdrChoice(j) * _extraXdrVar[LPIdx][j];
                        }
                    }
                    _model.addConstr(expr >= edge->_constDelay);
                } else {
                    _model.addConstr(expr >= edge->_constDelay + edge->_tdmCoef);
                }
            } else {
                _model.addConstr(expr >= edge->_constDelay);
            }
        }
    }
}

void TdmLpSolver::addLimitConstraint() {
    for (auto &pair : _tronconToXdrVar) {
        GRBLinExpr expr;
        for (auto i : pair.second) {
            for (int j = 0; j < _nChoice; j++) {
                double usage = 1.0 / TdmDB::getXdrChoice(j);
                expr += _xdrVar[i * _nChoice + j] * usage;
            }

            if (_useLP && _moreChoiceIdx > 0) {
                for (int j = 0, sz = _extraXdrVar[i].size(); j < sz; j++) {
                    if (getXdrChoice(j) <= TdmDB::getXdrChoice(_nChoice - 1))
                        expr += _extraXdrVar[i][j] * 1.0 / getXdrChoice(j);
                }
            }
        }
        _model.addConstr(expr <= pair.first->_limit);
    }
}

void TdmLpSolver::addExactLimitConstraint() {
    int cnt = 0;
    for (auto &pair : _tronconToXdrVar) {
        for (int j = 0; j < _nChoice; j++) {
            GRBLinExpr expr1, expr2;
            for (auto i : pair.second) {
                if (_optXdrVars[i]->isForward()) {
                    expr1 += _xdrVar[i * _nChoice + j];
                } else {
                    expr2 += _xdrVar[i * _nChoice + j];
                }
            }
            _model.addConstr(expr1 <= _usageVar[cnt * _nChoice * 2 + 2 * j] * TdmDB::getXdrChoice(j));
            _model.addConstr(expr2 <= _usageVar[cnt * _nChoice * 2 + 2 * j + 1] * TdmDB::getXdrChoice(j));
        }

        GRBLinExpr expr;
        for (int i = 0; i < _nChoice; i++) {
            expr += _usageVar[cnt * _nChoice * 2 + 2 * i] + _usageVar[cnt * _nChoice * 2 + 2 * i + 1];
        }
        _model.addConstr(expr <= pair.first->_limit);
        cnt++;
    }
}

void TdmLpSolver::init() {
    unsigned sinkId = _timingGraph->getSink()->_id;

    if (_useLP) {
        for (unsigned i = 0, sz = _xdrVar.size(); i < sz; i++) _xdrVar[i] = _model.addVar(0.0, 1, 0, GRB_CONTINUOUS);
        if (_moreChoiceIdx > 0) {
            for (unsigned i = 0, sz1 = _extraXdrVar.size(); i < sz1; i++)
                for (int j = 0, sz2 = _extraXdrVar[i].size(); j < sz2; j++)
                    _extraXdrVar[i][j] = _model.addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS);
        }
        for (unsigned i = 0, sz = _gateVar.size(); i < sz; i++)
            _gateVar[i] = _model.addVar(0.0, GRB_INFINITY, i == sinkId, GRB_CONTINUOUS);
        addLimitConstraint();
    } else {
        for (unsigned i = 0, sz = _xdrVar.size(); i < sz; i++) _xdrVar[i] = _model.addVar(0.0, 1, 0, GRB_BINARY);
        for (unsigned i = 0, sz = _gateVar.size(); i < sz; i++)
            _gateVar[i] = _model.addVar(0.0, GRB_INFINITY, i == sinkId, GRB_CONTINUOUS);
        for (unsigned i = 0, sz = _usageVar.size(); i < sz; i++)
            _usageVar[i] = _model.addVar(0.0, GRB_INFINITY, 0, GRB_INTEGER);
        addExactLimitConstraint();
    }

    addOneChoiceConstraint();
    addTimingEdgeConstraint();

    if (!_useLP) genILPInitSol();
}

void TdmLpSolver::genILPInitSol() {
    printlog(LOG_INFO, "ILP initial solution");
    int cnt = 0;
    for (auto &pair : _tronconToXdrVar) {
        vector<int> forwVars, backVars;
        for (auto i : pair.second) {
            XdrVar *var = _optXdrVars[i];
            if (var->isForward())
                forwVars.push_back(i);
            else
                backVars.push_back(i);
        }

        int limit = pair.first->_limit;
        int forwChoiceIdx = -1, backChoiceIdx = -1;
        forwChoiceIdx = TdmDB::getCeilChoiceIdx(pair.second.size() * 1.0 / limit);
        int forwUsage = ceil(forwVars.size() * 1.0 / TdmDB::getXdrChoice(forwChoiceIdx));
        backChoiceIdx = TdmDB::getCeilChoiceIdx(backVars.size() * 1.0 / (limit - forwUsage));
        int backUsage = ceil(backVars.size() * 1.0 / TdmDB::getXdrChoice(backChoiceIdx));

        for (auto i : forwVars) _xdrVar[i * _nChoice + forwChoiceIdx].set(GRB_DoubleAttr_Start, 1);
        for (auto i : backVars) _xdrVar[i * _nChoice + backChoiceIdx].set(GRB_DoubleAttr_Start, 1);

        printlog(LOG_INFO,
                 "troncon#%d: #forwVars=%lu, choice=%d, #backVars=%lu, choice=%d, forwUsage/backUsage/Limit=%d/%d/%d",
                 pair.first->_id,
                 forwVars.size(),
                 TdmDB::getXdrChoice(forwChoiceIdx),
                 backVars.size(),
                 TdmDB::getXdrChoice(backChoiceIdx),
                 forwUsage,
                 backUsage,
                 limit);

        for (int j = 0; j < _nChoice; j++) {
            double forwardUsage = 0, backwardUsage = 0;
            if (j == forwChoiceIdx) forwardUsage = forwVars.size();
            if (j == backChoiceIdx) backwardUsage = backVars.size();
            _usageVar[cnt * _nChoice * 2 + 2 * j].set(GRB_DoubleAttr_Start,
                                                      ceil(forwardUsage / TdmDB::getXdrChoice(j)));
            _usageVar[cnt * _nChoice * 2 + 2 * j + 1].set(GRB_DoubleAttr_Start,
                                                          ceil(backwardUsage / TdmDB::getXdrChoice(j)));
        }

        cnt++;
    }
}

void TdmLpSolver::solve() {
    log() << "==================== begin TDM analytical solving ====================" << endl;

    tdmDatabase.reportSol();

    init();

    _model.getEnv().set(GRB_DoubleParam_TimeLimit, _timeLimit);
    _model.getEnv().set(GRB_IntParam_Threads, setting.nThreads);
    // _model.getEnv().set(GRB_IntParam_OutputFlag, 0);

    _model.optimize();
    getResult();

    log() << "---------------- finish TDM analytical solving ----------------" << endl << endl;
}

void TdmLpSolver::getResult() {
    int status = _model.get(GRB_IntAttr_Status);
    if (status == GRB_OPTIMAL) {
        log() << "optimal, objective=" << _model.get(GRB_DoubleAttr_ObjVal) << endl;
    } else if (status == GRB_SUBOPTIMAL) {
        log() << "suboptimal, objective=" << _model.get(GRB_DoubleAttr_ObjVal) << endl;
    } else if (status == GRB_INF_OR_UNBD) {
        log() << "infeasible" << endl;
    } else if (status == GRB_TIME_LIMIT) {
        log() << "time out" << endl;
    } else {
        log() << "unknown return: " << status << endl;
        getchar();
    }

    for (int i = 0; i < _nVar; i++) {
        double val = 0;
        for (int j = 0; j < _nChoice; j++)
            val += _xdrVar[i * _nChoice + j].get(GRB_DoubleAttr_X) * TdmDB::getXdrChoice(j);

        if (_useLP && _moreChoiceIdx > 0) {
            for (int j = 0, sz = _extraXdrVar[i].size(); j < sz; j++) {
                if (getXdrChoice(j) <= TdmDB::getXdrChoice(_nChoice - 1))
                    val += _extraXdrVar[i][j].get(GRB_DoubleAttr_X) * getXdrChoice(j);
            }
        }

        _optXdrVars[i]->setVal(val);
    }

    tdmDatabase.reportSol();
    // tdmDatabase.reportTdmAssignment();
}

TdmLpSolver::TdmLpSolver(bool useLP) : _optXdrVars(tdmDatabase.getOptXdrVars()), _model(_env), _useLP(useLP) {
    _nVar = _optXdrVars.size();
    tdmDatabase.checkFeasibility(_nChoice);

    _timingGraph = tdmDatabase.getTimingGraph();
    for (int i = 0; i < _nVar; i++) {
        Troncon *troncon = _optXdrVars[i]->getNet()->_troncon;
        _tronconToXdrVar[troncon].push_back(i);
    }

    int xdrVarNum = _nVar * _nChoice;
    int gateVarNum = _timingGraph->getNumNodes();
    int usageVarNum = _tronconToXdrVar.size() * _nChoice * 2;

    _xdrVar.resize(xdrVarNum);
    _gateVar.resize(gateVarNum);
    _usageVar.resize(usageVarNum);

    if (_useLP && _moreChoiceIdx > 0) {
        _extraXdrVar.resize(_optXdrVars.size());
        for (int i = 0, sz = _optXdrVars.size(); i < sz; i++) _extraXdrVar[i].resize(6 + (_moreChoiceIdx - 1) * 7);
    }
}

int TdmLpSolver::getXdrChoice(int extraXdrVarIdx) {
    if (extraXdrVarIdx >= 6) {
        int tmp = extraXdrVarIdx - 6;
        int baseChoiceVal = TdmDB::getXdrChoice(tmp / 7 + 1);
        return baseChoiceVal + tmp % 7 + 1;
    } else {
        return extraXdrVarIdx + 2;
    }
}
