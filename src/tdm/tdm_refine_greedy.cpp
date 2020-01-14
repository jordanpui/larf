#include "tdm_refine_greedy.h"
#include "tdm_db.h"
#include "timing_graph.h"
#include "tdm_net.h"

void TdmRefine::solve() {
    log() << "==================== begin TDM refinement(greedy) ====================" << endl;
    tdmDatabase.reportSol();

    tdmDatabase.updateTiming();
    auto timingGraph = tdmDatabase.getTimingGraph();
    vector<Edge *> criticalPath = timingGraph->getCriticalPath();

    int iter = 0;
    double bestCost = DBL_MAX;

    while (optimizePath(criticalPath)) {
        if (iter % 10 == 0) {
            if (bestCost <= tdmDatabase.getArrivalTime())
                break;
            else
                bestCost = tdmDatabase.getArrivalTime();
        }

        criticalPath = timingGraph->getCriticalPath();
        iter++;
    }

    tdmDatabase.reportSol();
    log() << "---------------- finish TDM refinement ----------------" << endl;
}

bool TdmRefine::optimizePath(vector<Edge *> &criticalPath) {
    // tdmDatabase.reportSol();

    unordered_set<XdrVar *> vars;
    for (auto edge : criticalPath) {
        if (!edge->isConstEdge()) vars.insert(edge->_net->getXdrVar());
    }

    unordered_map<Troncon *, vector<XdrVar *>> tronconToForwVars, tronconToBackVars;
    for (auto var : vars) {
        if (var->isForward())
            tronconToForwVars[tdmDatabase.getTroncon(var)].push_back(var);
        else
            tronconToBackVars[tdmDatabase.getTroncon(var)].push_back(var);
    }

    auto getRelateVars = [&](Troncon *troncon, vector<XdrVar *> &vars, bool isForward) {
        int sz = vars.size();

        auto &nets = troncon->getNets();
        for (auto net : nets) {
            if (net->isInterNet() && net->getXdrVar()->isForward() == isForward) {
                bool sameVar = false;
                for (int i = 0; i < sz; i++) {
                    if (net->getXdrVar() == vars[i]) {
                        sameVar = true;
                        break;
                    }
                }
                if (!sameVar) vars.push_back(net->getXdrVar());
            }
        }
    };

    for (auto &pair : tronconToForwVars) {
        int sz = pair.second.size();
        getRelateVars(pair.first, pair.second, true);
        if (optimizeTronconSort(pair.second, sz)) return true;
    }
    for (auto &pair : tronconToBackVars) {
        int sz = pair.second.size();
        getRelateVars(pair.first, pair.second, false);
        if (optimizeTronconSort(pair.second, sz)) return true;
    }

    return false;
}

bool TdmRefine::optimizeTroncon(vector<XdrVar *> &vars, int numOptVar) {
    for (int i = 0; i < numOptVar; i++) {
        double val = vars[i]->getVal();

        unordered_map<double, vector<XdrVar *>> candMap;
        for (unsigned j = numOptVar; j < vars.size(); j++)
            if (vars[j]->getVal() < val) candMap[vars[j]->getVal()].push_back(vars[j]);

        for (auto &pair : candMap) {
            int candVal = pair.first;

            for (auto var : pair.second) {
                bool canChange = true;
                vector<Edge *> relateEdges = tdmDatabase.getTimingGraph()->getEdges(var);
                for (auto edge : relateEdges) {
                    double slack = edge->_fanout->getSlack();
                    double diffDelay = edge->getDelay(val) - edge->getDelay();
                    if (slack < diffDelay) {
                        canChange = false;
                        break;
                    }
                }

                if (canChange) {
                    SwapHist swap1(vars[i], var);

                    if (hasSwap(swap1)) continue;

                    _hist.push_back(swap1);

                    vars[i]->setVal(candVal);
                    var->setVal(val);

                    // cout << val << " " << candVal << " " << var << endl;

                    double orgAT = tdmDatabase.getArrivalTime();
                    tdmDatabase.updateTiming();
                    if (tdmDatabase.getArrivalTime() > orgAT) {
                        vars[i]->setVal(val);
                        var->setVal(candVal);
                        tdmDatabase.updateTiming();
                        _hist.resize(_hist.size() - 1);
                    } else {
                        // cout << "succ: " << val << " " << candVal << " " << var << " " << vars[i] << endl;
                        return true;
                    }
                }
            }
        }
    }

    return false;
}

bool TdmRefine::optimizeTronconSort(vector<XdrVar *> &vars, int numOptVar) {
    for (int i = 0; i < numOptVar; i++) {
        double val = vars[i]->getVal();

        unordered_map<double, vector<XdrVar *>> candMap;
        for (unsigned j = numOptVar; j < vars.size(); j++)
            if (vars[j]->getVal() < val) candMap[vars[j]->getVal()].push_back(vars[j]);

        vector<pair<double, XdrVar *>> sortedCands;

        for (auto &pair : candMap) {
            for (auto var : pair.second) {
                bool canChange = true;
                double minResSlack = DBL_MAX;

                vector<Edge *> relateEdges = tdmDatabase.getTimingGraph()->getEdges(var);
                for (auto edge : relateEdges) {
                    double slack = edge->_fanout->getSlack();
                    double diffDelay = edge->getDelay(val) - edge->getDelay();
                    if (slack < diffDelay) {
                        canChange = false;
                        break;
                    }
                    minResSlack = min(minResSlack, slack - diffDelay);
                }

                if (canChange) sortedCands.emplace_back(-1 * minResSlack, var);
            }
        }

        sort(sortedCands.begin(), sortedCands.end());

        for (auto &pair : sortedCands) {
            int candVal = pair.second->getVal();

            SwapHist swap1(vars[i], pair.second);

            if (hasSwap(swap1)) continue;

            _hist.push_back(swap1);

            vars[i]->setVal(candVal);
            pair.second->setVal(val);

            // cout << val << " " << candVal << " " << pair.second << endl;

            double orgAT = tdmDatabase.getArrivalTime();
            tdmDatabase.updateTiming();
            if (tdmDatabase.getArrivalTime() > orgAT) {
                vars[i]->setVal(val);
                pair.second->setVal(candVal);
                tdmDatabase.updateTiming();
                _hist.resize(_hist.size() - 1);
            } else {
                // cout << "succ: " << val << " " << candVal << " " << pair.second << " " << vars[i] << endl;
                return true;
            }
        }
    }

    return false;
}

bool TdmRefine::hasSwap(const SwapHist &swap1) const {
    for (const auto &swap2 : _hist)
        if (swap1.isSame(swap2)) return true;
    return false;
}

SwapHist::SwapHist(XdrVar *u, XdrVar *v) {
    _u = u;
    _v = v;
    _uVal = u->getVal();
    _vVal = v->getVal();
}

bool SwapHist::isSame(const SwapHist &swap1) const {
    if (_u == swap1._u && _v == swap1._v && _uVal == swap1._uVal && _vVal == swap1._vVal) return true;
    if (_v == swap1._u && _u == swap1._v && _vVal == swap1._uVal && _uVal == swap1._vVal) return true;

    return false;
}