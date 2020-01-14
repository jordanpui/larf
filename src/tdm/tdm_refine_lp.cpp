#include "tdm_refine_lp.h"
#include "tdm_db.h"
#include "tdm_net.h"
#include "timing_graph.h"

void TdmRefineLP::getResult() {
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
        for (int j = _choiceRanges[i].first; j <= _choiceRanges[i].second; j++)
            val += getXdrVar(i, j).get(GRB_DoubleAttr_X) * TdmDB::getXdrChoice(j);
        if (_genContSol && _moreChoiceIdx > 0) {
            for (int j = 0, sz = _extraXdrVar[i].size(); j < sz; j++) {
                if (withinXdrChoiceRange(i, getXdrChoice(j)))
                    val += _extraXdrVar[i][j].get(GRB_DoubleAttr_X) * getXdrChoice(j);
            }
        }
        _optXdrVars[i]->setVal(val);
    }

    // for (int i = 0; i < _nVar; i++) {
    //     double val = 0;
    //     for (int j = _choiceRanges[i].first; j <= _choiceRanges[i].second; j++)
    //         val += getXdrVar(i, j).get(GRB_DoubleAttr_X) * 1.0 / TdmDB::getXdrChoice(j);
    //     if (_genContSol && _moreChoiceIdx > 0) {
    //         for (int j = 0, sz = _extraXdrVar[i].size(); j < sz; j++) {
    //             if (withinXdrChoiceRange(i, getXdrChoice(j)))
    //                 val += _extraXdrVar[i][j].get(GRB_DoubleAttr_X) * 1.0 / getXdrChoice(j);
    //         }
    //     }
    //     _optXdrVars[i]->setVal(1.0 / val);
    // }
}

TdmRefineLP::TdmRefineLP(bool genContSol)
    : _model(_env), _optXdrVars(tdmDatabase.getOptXdrVars()), _genContSol(genContSol) {
    _nVar = _optXdrVars.size();

    _timingGraph = tdmDatabase.getTimingGraph();
    for (int i = 0; i < _nVar; i++) {
        Troncon *troncon = _optXdrVars[i]->getNet()->_troncon;
        _tronconToXdrVar[troncon].push_back(i);
    }

    _choiceRanges.resize(_optXdrVars.size());
    for (int i = 0, sz = _optXdrVars.size(); i < sz; i++) {
        int centerIdx = TdmDB::getClosestChoiceIdx(_optXdrVars[i]->getVal());
        _choiceRanges[i].first = max(0, centerIdx - _choiceBnd);
        _choiceRanges[i].second = min(TdmDB::getNumChoices() - 1, centerIdx + _choiceBnd);
    }

    int gateVarNum = _timingGraph->getNumNodes();
    int usageVarNum = _tronconToXdrVar.size() * TdmDB::getNumChoices() * 2;

    _gateVar.resize(gateVarNum);
    _usageVar.resize(usageVarNum);

    _xdrVar.resize(_optXdrVars.size());
    for (int i = 0, sz = _optXdrVars.size(); i < sz; i++) _xdrVar[i].resize(getChoiceNum(i));

    if (_genContSol && _moreChoiceIdx > 0) {
        _extraXdrVar.resize(_optXdrVars.size());
        for (int i = 0, sz = _optXdrVars.size(); i < sz; i++)
            if (_choiceRanges[i].first < _moreChoiceIdx) _extraXdrVar[i].resize(6 + (_moreChoiceIdx - 1) * 7);
    }
}

void TdmRefineLP::addOneChoiceConstraint() {
    // only one choice
    for (int i = 0; i < _nVar; i++) {
        GRBLinExpr expr = GRBLinExpr();
        for (int j = _choiceRanges[i].first; j <= _choiceRanges[i].second; j++) expr += getXdrVar(i, j);

        if (_genContSol && _moreChoiceIdx > 0) {
            for (int j = 0, sz = _extraXdrVar[i].size(); j < sz; j++) {
                if (withinXdrChoiceRange(i, getXdrChoice(j))) expr += _extraXdrVar[i][j];
            }
        }

        _model.addConstr(expr == 1);
    }
}

void TdmRefineLP::addTimingEdgeConstraint() {
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
                    for (int j = _choiceRanges[LPIdx].first; j <= _choiceRanges[LPIdx].second; j++)
                        expr -= edge->_tdmCoef * TdmDB::getXdrChoice(j) * getXdrVar(LPIdx, j);

                    if (_genContSol && _moreChoiceIdx > 0) {
                        for (int j = 0, sz = _extraXdrVar[LPIdx].size(); j < sz; j++) {
                            if (withinXdrChoiceRange(LPIdx, getXdrChoice(j)))
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

void TdmRefineLP::addLimitConstraint() {
    for (auto &pair : _tronconToXdrVar) {
        GRBLinExpr expr;
        for (auto i : pair.second) {
            for (int j = _choiceRanges[i].first; j <= _choiceRanges[i].second; j++) {
                double usage = 1.0 / TdmDB::getXdrChoice(j);
                expr += getXdrVar(i, j) * usage;
            }
        }
        if (_genContSol && _moreChoiceIdx > 0) {
            for (auto i : pair.second) {
                for (int j = 0, sz = _extraXdrVar[i].size(); j < sz; j++) {
                    if (withinXdrChoiceRange(i, getXdrChoice(j))) expr += _extraXdrVar[i][j] * 1.0 / getXdrChoice(j);
                }
            }
        }
        _model.addConstr(expr <= pair.first->_limit);
    }
}

void TdmRefineLP::addExactLimitConstraint() {
    int cnt = 0;
    int nChoice = TdmDB::getNumChoices();
    for (auto &pair : _tronconToXdrVar) {
        for (int j = 0; j < nChoice; j++) {
            GRBLinExpr expr1, expr2;
            for (auto i : pair.second) {
                if (j >= _choiceRanges[i].first && j <= _choiceRanges[i].second) {
                    if (_optXdrVars[i]->isForward()) {
                        expr1 += getXdrVar(i, j);
                    } else {
                        expr2 += getXdrVar(i, j);
                    }
                }
            }
            _model.addConstr(expr1 <= _usageVar[cnt * nChoice * 2 + 2 * j] * TdmDB::getXdrChoice(j));
            _model.addConstr(expr2 <= _usageVar[cnt * nChoice * 2 + 2 * j + 1] * TdmDB::getXdrChoice(j));
        }

        GRBLinExpr expr;
        for (int i = 0; i < nChoice; i++) {
            expr += _usageVar[cnt * nChoice * 2 + 2 * i] + _usageVar[cnt * nChoice * 2 + 2 * i + 1];
        }
        _model.addConstr(expr <= pair.first->_limit);
        cnt++;
    }
}

void TdmRefineLP::init() {
    unsigned sinkId = _timingGraph->getSink()->_id;

    if (_genContSol) {
        for (unsigned i = 0, sz1 = _xdrVar.size(); i < sz1; i++)
            for (int j = 0, sz2 = _xdrVar[i].size(); j < sz2; j++)
                _xdrVar[i][j] = _model.addVar(0.0, 1, 0, GRB_CONTINUOUS);
        if (_moreChoiceIdx > 0) {
            for (unsigned i = 0, sz1 = _extraXdrVar.size(); i < sz1; i++)
                for (int j = 0, sz2 = _extraXdrVar[i].size(); j < sz2; j++)
                    _extraXdrVar[i][j] = _model.addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS);
        }
        for (unsigned i = 0, sz = _gateVar.size(); i < sz; i++)
            _gateVar[i] = _model.addVar(0.0, GRB_INFINITY, i == sinkId, GRB_CONTINUOUS);
        addLimitConstraint();
    } else {
        for (unsigned i = 0, sz1 = _xdrVar.size(); i < sz1; i++)
            for (int j = 0, sz2 = _xdrVar[i].size(); j < sz2; j++) _xdrVar[i][j] = _model.addVar(0.0, 1, 0, GRB_BINARY);
        for (unsigned i = 0, sz = _gateVar.size(); i < sz; i++)
            _gateVar[i] = _model.addVar(0.0, GRB_INFINITY, i == sinkId, GRB_CONTINUOUS);
        for (unsigned i = 0, sz = _usageVar.size(); i < sz; i++)
            _usageVar[i] = _model.addVar(0.0, GRB_INFINITY, 0, GRB_INTEGER);
        addExactLimitConstraint();
    }

    addOneChoiceConstraint();
    addTimingEdgeConstraint();

    if (!_genContSol) genILPInitSol();
}

void TdmRefineLP::genILPInitSol() {
    printlog(LOG_INFO, "ILP initial solution");
    for (auto &pair : _tronconToXdrVar) {
        vector<int> forwVars, backVars;
        for (auto i : pair.second) {
            XdrVar *var = _optXdrVars[i];
            int choiceIdx = TdmDB::getClosestChoiceIdx(var->getVal());
            int choice = TdmDB::getXdrChoice(choiceIdx);
            var->setVal(choice);
            getXdrVar(i, choiceIdx).set(GRB_DoubleAttr_Start, choice);
        }
    }

    int cnt = 0;
    int nChoice = TdmDB::getNumChoices();
    for (auto &pair : _tronconToXdrVar) {
        for (int j = 0; j < nChoice; j++) {
            double forwardUsage = 0, backwardUsage = 0;
            for (auto i : pair.second) {
                XdrVar *var = _optXdrVars[i];
                if (TdmDB::getXdrChoice(j) == var->getVal()) {
                    if (var->isForward()) {
                        forwardUsage++;
                    } else {
                        backwardUsage++;
                    }
                }
            }
            _usageVar[cnt * nChoice * 2 + 2 * j].set(GRB_DoubleAttr_Start, ceil(forwardUsage / TdmDB::getXdrChoice(j)));
            _usageVar[cnt * nChoice * 2 + 2 * j + 1].set(GRB_DoubleAttr_Start,
                                                         ceil(backwardUsage / TdmDB::getXdrChoice(j)));
        }

        cnt++;
    }
}

void TdmRefineLP::solve() {
    if (_genContSol) {
        log() << "==================== begin TDM refinement(lp cont) ====================" << endl;
    } else {
        log() << "==================== begin TDM refinement(ilp disc) ====================" << endl;
    }
    tdmDatabase.reportSol();

    init();

    _model.getEnv().set(GRB_DoubleParam_TimeLimit, _timeLimit);
    _model.getEnv().set(GRB_IntParam_Threads, setting.nThreads);
    // _model.getEnv().set(GRB_IntParam_OutputFlag, 0);

    _model.optimize();
    getResult();

    tdmDatabase.reportSol();
    log() << "---------------- finish TDM refinement ----------------" << endl << endl;
}

int TdmRefineLP::getXdrChoice(int extraXdrVarIdx) {
    if (extraXdrVarIdx >= 6) {
        int tmp = extraXdrVarIdx - 6;
        int baseChoiceVal = TdmDB::getXdrChoice(tmp / 7 + 1);
        return baseChoiceVal + tmp % 7 + 1;
    } else {
        return extraXdrVarIdx + 2;
    }
}

bool TdmRefineLP::withinXdrChoiceRange(int idx, int choice) {
    return choice >= TdmDB::getXdrChoice(_choiceRanges[idx].first) &&
           choice <= TdmDB::getXdrChoice(_choiceRanges[idx].second);
}

void TdmRefineLP::printUsageError() {
    const int choiceBnd = 5;
    const int moreChoiceIdx = 2;

    ofstream fout("../usageError.txt");
    for (int val = 1; val <= 1600; val++) {
        int floorChoiceIdx = TdmDB::getFloorChoiceIdx(val);
        int ceilChoiceIdx = TdmDB::getCeilChoiceIdx(val);
        double maxActualUsage = 0;

        for (int choiceBegIdx = max(0, ceilChoiceIdx - choiceBnd * 2); choiceBegIdx <= floorChoiceIdx; choiceBegIdx++) {
            int choiceEndIdx = min(choiceBegIdx + choiceBnd * 2, 200);

            GRBEnv env;
            GRBModel model(env);

            vector<GRBVar> vars(choiceEndIdx - choiceBegIdx + 1);
            vector<GRBVar> extraVars(6 + (moreChoiceIdx - 1) * 7);

            for (int i = 0, sz = vars.size(); i < sz; i++) vars[i] = model.addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS);
            for (int i = 0, sz = extraVars.size(); i < sz; i++)
                extraVars[i] = model.addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS);
            GRBVar obj = model.addVar(0.0, GRB_INFINITY, 1, GRB_CONTINUOUS);

            GRBLinExpr expr1;
            for (int i = 0, sz = vars.size(); i < sz; i++) expr1 += TdmDB::getXdrChoice(i + choiceBegIdx) * vars[i];
            if (choiceBegIdx < moreChoiceIdx) {
                int cnt = 0;
                for (int i = 0; i < moreChoiceIdx; i++) {
                    int varSize = i < 1 ? 6 : 7;
                    if (i >= choiceBegIdx && i < choiceEndIdx)
                        for (int j = 0; j < varSize; j++)
                            expr1 += (TdmDB::getXdrChoice(i) + j + 1) * extraVars[cnt + j];
                    cnt += varSize;
                }
            }
            model.addConstr(expr1 == val);

            GRBLinExpr expr2;
            for (int i = 0, sz = vars.size(); i < sz; i++) expr2 += vars[i];
            if (choiceBegIdx < moreChoiceIdx) {
                int cnt = 0;
                for (int i = 0; i < moreChoiceIdx; i++) {
                    int varSize = i < 1 ? 6 : 7;
                    if (i >= choiceBegIdx && i < choiceEndIdx)
                        for (int j = 0; j < varSize; j++) expr2 += extraVars[cnt + j];
                    cnt += varSize;
                }
            }
            model.addConstr(expr2 == 1);

            GRBLinExpr expr3;
            for (int i = 0, sz = vars.size(); i < sz; i++)
                expr3 += 1.0 / TdmDB::getXdrChoice(i + choiceBegIdx) * vars[i];
            if (choiceBegIdx < moreChoiceIdx) {
                int cnt = 0;
                for (int i = 0; i < moreChoiceIdx; i++) {
                    int varSize = i < 1 ? 6 : 7;
                    if (i >= choiceBegIdx && i < choiceEndIdx)
                        for (int j = 0; j < varSize; j++)
                            expr3 += 1.0 / (TdmDB::getXdrChoice(i) + j + 1) * extraVars[cnt + j];
                    cnt += varSize;
                }
            }
            model.addConstr(expr3 == obj);

            model.getEnv().set(GRB_IntParam_OutputFlag, 0);
            model.optimize();

            double actualUsage = 0;
            for (int i = 0, sz = vars.size(); i < sz; i++)
                actualUsage += vars[i].get(GRB_DoubleAttr_X) * 1.0 / TdmDB::getXdrChoice(i + choiceBegIdx);
            if (choiceBegIdx < moreChoiceIdx) {
                int cnt = 0;
                for (int i = 0; i < moreChoiceIdx; i++) {
                    int varSize = i < 1 ? 6 : 7;
                    if (i >= choiceBegIdx && i < choiceEndIdx)
                        for (int j = 0; j < varSize; j++)
                            actualUsage +=
                                1.0 / (TdmDB::getXdrChoice(i) + j + 1) * extraVars[cnt + j].get(GRB_DoubleAttr_X);
                    cnt += varSize;
                }
            }

            maxActualUsage = max(maxActualUsage, actualUsage);
            // log() << val << " " << actualUsage << " " << 1.0 / val << endl;
        }
        fout << maxActualUsage << " " << 1.0 / val << endl;
        // log() << val << " " << maxActualUsage << " " << 1.0 / val << endl;
    }

    fout.close();

    log() << "finish lp resource error calc...." << endl;
}
