#include "tdm_leg.h"
#include "tdm_db.h"
#include "tdm_net.h"
#include "timing_graph.h"

TdmLegalize::TdmLegalize() : _optXdrVars(tdmDatabase.getOptXdrVars()) {
    _flow = setting.lg;
    _timingGraph = tdmDatabase.getTimingGraph();
    for (unsigned i = 0; i < _optXdrVars.size(); i++) {
        _varToOptIdx[_optXdrVars[i]] = i;
        Troncon *troncon = _optXdrVars[i]->getNet()->_troncon;
        _tronconToXdrVar[troncon].push_back(i);
    }

    int cnt = 0;
    for (auto &pair : _tronconToXdrVar) _tronconToIdx[pair.first] = cnt++;
}

void TdmLegalize::solve(bool writeDB, vector<double> *result) {
    timer::timer time;
    time.start();
    if (_flow == Setting::Lg_Disp) {
        log() << "=================== begin TDM legalization(Disp) ===================" << endl;
    } else if (_flow == Setting::Lg_MaxDisp) {
        log() << "=================== begin TDM legalization(MaxDisp) ===================" << endl;
    }

    tdmDatabase.updateTiming();

    double begAT = tdmDatabase.getArrivalTime();

    updateWireData();

    sort(_wireData.begin(), _wireData.end(), [&](const WireData &data1, const WireData &data2) {
        return data1._vars.size() > data2._vars.size();
    });

    double totalDisp = 0;
    double maxMaxDisp = 0;
    double avgMaxDisp = 0;

    std::mutex idx_mutex;
    int vIdx = 0;
    auto threadFunc = [&]() {
        while (true) {
            int idx;
            idx_mutex.lock();
            idx = vIdx++;
            idx_mutex.unlock();

            if (idx >= (int)_wireData.size()) break;

            auto data = _wireData[idx];

            double choiceVio = tdmDatabase.getChoiceVio(data._troncon);
            double limitVio = tdmDatabase.getContLimitVio(data._troncon);

            if (choiceVio == 0) {
                limitVio = tdmDatabase.getLimitVio(data._troncon);
                if (limitVio == 0) continue;
            }

            Memorization memorization(data._vars.size(), data._troncon->_limit);
            legalizeTroncon(data, memorization);

            data.calcDisp(memorization);
            if (writeDB) {
                data.recoverSol(memorization);
            } else {
                data.dumpSol(*result, _varToOptIdx, memorization);
            }

            idx_mutex.lock();
            totalDisp += data._totDisp;
            maxMaxDisp = max(maxMaxDisp, data._maxDisp);
            avgMaxDisp += data._maxDisp * data._vars.size();
            // printlog(LOG_INFO,
            //          "troncon#%d: cost=%.3f, usage=%.3f, #vars(fwd)=%lu(%d), disp:tot/avg/max=%.3f/%.3f/%.3f",
            //          data._troncon->_id,
            //          memorization.getBestCost(0, data._troncon->_limit),
            //          data._troncon->getUsage() * 1.0 / data._troncon->_limit,
            //          data._vars.size(),
            //          data._numForwardVars,
            //          data._totDisp,
            //          data._totDisp / data._vars.size(),
            //          data._maxDisp);
            idx_mutex.unlock();
        }
    };

    int nThreads = max(1, setting.nThreads);
    std::thread threads[nThreads];
    for (int i = 0; i < nThreads; i++) threads[i] = std::thread(threadFunc);
    for (int i = 0; i < nThreads; i++) threads[i].join();

    tdmDatabase.updateTiming();
    double endAT = tdmDatabase.getArrivalTime();

    printlog(LOG_INFO,
             "input/output: avgDisp=%.3f, maxMaxDisp=%.3f, avgMaxDisp=%.3f, arrival_time=%.3f/%.3f, time=%f",
             totalDisp / _optXdrVars.size(),
             maxMaxDisp,
             avgMaxDisp / _optXdrVars.size(),
             begAT,
             endAT,
             time.elapsed());

    log() << "---------------- finish TDM legalization ----------------" << endl;
}

void WireData::sortVars(Setting::LgMethod flow) {
    sortByDisp();

    _numForwardVars = 0;
    for (auto var : _vars) {
        if (var->isForward()) {
            _numForwardVars++;
        } else {
            break;
        }
    }
}

void WireData::sortByDisp() {
    sort(_vars.begin(), _vars.end(), [&](OneWireData *var1, OneWireData *var2) {
        if (var1->isForward() != var2->isForward()) return var1->isForward() > var2->isForward();
        return var1->getVal() < var2->getVal();
    });
}

void TdmLegalize::updateWireData() {
    _wireData.resize(_tronconToXdrVar.size());
    for (auto &pair : _tronconToXdrVar) {
        Troncon *troncon = pair.first;
        int idx = _tronconToIdx[troncon];

        _wireData[idx]._troncon = troncon;
        for (auto var : pair.second) _wireData[idx]._vars.push_back(new OneWireData(_optXdrVars[var]));
    }

    for (auto &data : _wireData) data.sortVars(_flow);
}

int WireData::calcEndIdx(unsigned idx, int choice) const {
    return _vars[idx]->isForward() ? min(idx + choice, _numForwardVars) : min(idx + choice, (unsigned)_vars.size());
}

int WireData::calcEndIdx(unsigned idx, int choice, const vector<pair<int, int>> &choiceRanges) const {
    int maxEndIdx = _vars[idx]->isForward() ? _numForwardVars : _vars.size();
    maxEndIdx = min(maxEndIdx, (int)idx + choice);

    for (int i = idx + 1; i < maxEndIdx; i++) {
        if (choice < TdmDB::getXdrChoice(choiceRanges[i].first) ||
            choice > TdmDB::getXdrChoice(choiceRanges[i].second)) {
            return i;
        }
    }

    return maxEndIdx;
}

int WireData::getMinWireRequire(unsigned idx) const {
    int minWireRequire = 0;
    int largestXdr = TdmDB::getXdrChoices().back();

    if (idx < _numForwardVars)
        minWireRequire = ceil(double(_numForwardVars - idx) / largestXdr) +
                         ceil(double(_vars.size() - _numForwardVars) / largestXdr);
    else
        minWireRequire = ceil(double(_vars.size() - idx) / largestXdr);

    return minWireRequire;
}

int WireData::getMinWireRequire(unsigned idx, const vector<pair<int, int>> &choiceRanges) const {
    // can be stored in advanced to avoid calculation
    int minWireRequire = 0;

    for (int i = idx, sz = choiceRanges.size(); i < sz;) {
        i = calcEndIdx(i, TdmDB::getXdrChoice(choiceRanges[i].second), choiceRanges);
        minWireRequire++;
    }
    return minWireRequire;
}

double TdmLegalize::legalizeTroncon(const WireData &data, Memorization &memorization) const {
    if (_flow == Setting::Lg_Disp) {
        // return data.legalizeTronconDisp(0, data._troncon->_limit, memorization);
        return data.legalizeTronconDispPrune(0, data._troncon->_limit, memorization);
    } else if (_flow == Setting::Lg_MaxDisp) {
        return data.legalizeTronconMaxDisp(memorization);
    } else {
        log() << "unkonw flow of legalization" << endl;
        return -1;
    }
}

double WireData::legalizeTronconDispPrune(unsigned idx, int p, Memorization &memorization) const {
    if (idx >= _vars.size()) return 0;
    if (p == 0) return -1;

    if (memorization.isOptimized(idx, p)) return memorization.getBestCost(idx, p);

    double bestCost = DBL_MAX;
    int bestChoice = -1;
    int bestEndIdx = -1;

    if (getMinWireRequire(idx) > p) {
        memorization.saveBest(idx, p, bestCost, bestChoice, bestEndIdx);
        return bestCost;
    }

    int candBeg = TdmDB::getClosestChoiceIdx(_vars[idx]->getVal()), candEnd = TdmDB::getXdrChoices().size();
    int prevEndIdx = idx;

    for (int c = candBeg; c < candEnd; c++) {
        int choice = TdmDB::getXdrChoice(c);

        unsigned minEndIdx = prevEndIdx + 1;
        unsigned endIdx = calcEndIdx(idx, choice);

        if (minEndIdx > endIdx) break;

        while (minEndIdx < endIdx) {
            if (TdmDB::getClosestChoice(_vars[minEndIdx]->getVal()) == choice) {
                minEndIdx++;
            } else {
                break;
            }
        }

        prevEndIdx = minEndIdx - 1;
        for (unsigned curEndIdx = minEndIdx; curEndIdx < endIdx; curEndIdx++) {
            if (_vars[curEndIdx]->getVal() <= choice) prevEndIdx = curEndIdx;
        }

        double curBest = DBL_MAX;
        for (unsigned curEndIdx = minEndIdx; curEndIdx <= endIdx; curEndIdx++) {
            double currentCost = 0;
            for (unsigned i = idx; i < curEndIdx; i++) currentCost += abs(_vars[i]->getVal() - choice);

            double cost = legalizeTronconDispPrune(curEndIdx, p - 1, memorization);
            double ret = cost == -1 ? DBL_MAX : currentCost + cost;

            if (currentCost > curBest) break;

            if (ret < curBest) {
                curBest = ret;

                if (ret < bestCost) {
                    bestCost = ret;
                    bestChoice = choice;
                    bestEndIdx = curEndIdx;
                }
            }
        }
    }

    memorization.saveBest(idx, p, bestCost, bestChoice, bestEndIdx);

    return bestCost;
}

double WireData::legalizeTronconDisp(unsigned idx, int p, Memorization &memorization) const {
    if (idx >= _vars.size()) return 0;
    if (p == 0) return -1;

    if (memorization.isOptimized(idx, p)) return memorization.getBestCost(idx, p);

    double bestCost = DBL_MAX;
    int bestChoice = -1;
    int bestEndIdx = -1;

    if (getMinWireRequire(idx) > p) {
        memorization.saveBest(idx, p, bestCost, bestChoice, bestEndIdx);
        return bestCost;
    }

    int candBeg = TdmDB::getClosestChoiceIdx(_vars[idx]->getVal()), candEnd = TdmDB::getXdrChoices().size();
    int prevEndIdx = idx;

    for (int c = candBeg; c < candEnd; c++) {
        int choice = TdmDB::getXdrChoice(c);

        unsigned minEndIdx = prevEndIdx + 1;
        unsigned endIdx = calcEndIdx(idx, choice);

        if (minEndIdx > endIdx) break;

        double curBest = DBL_MAX;
        for (unsigned curEndIdx = minEndIdx; curEndIdx <= endIdx; curEndIdx++) {
            double currentCost = 0;
            for (unsigned i = idx; i < curEndIdx; i++) currentCost += abs(_vars[i]->getVal() - choice);

            double cost = legalizeTronconDisp(curEndIdx, p - 1, memorization);
            double ret = cost == -1 ? DBL_MAX : currentCost + cost;

            if (currentCost > curBest) break;

            if (ret < curBest) {
                curBest = ret;

                if (ret < bestCost) {
                    bestCost = ret;
                    bestChoice = choice;
                    bestEndIdx = curEndIdx;
                }
            }
        }
    }

    memorization.saveBest(idx, p, bestCost, bestChoice, bestEndIdx);

    return bestCost;
}

double WireData::legalizeTronconDisp(unsigned idx,
                                     int p,
                                     Memorization &memorization,
                                     const vector<pair<int, int>> &choiceRanges) const {
    if (idx >= _vars.size()) return 0;
    if (p == 0) return -1;

    if (memorization.isOptimized(idx, p)) return memorization.getBestCost(idx, p);

    double bestCost = DBL_MAX;
    int bestChoice = -1;
    int bestEndIdx = -1;

    if (getMinWireRequire(idx, choiceRanges) > p) {
        memorization.saveBest(idx, p, bestCost, bestChoice, bestEndIdx);
        return bestCost;
    }

    int candBeg = choiceRanges[idx].first, candEnd = choiceRanges[idx].second;
    // int prevEndIdx = idx;

    for (int c = candBeg; c <= candEnd; c++) {
        int choice = TdmDB::getXdrChoice(c);

        // unsigned minEndIdx = prevEndIdx + 1;
        unsigned endIdx = calcEndIdx(idx, choice, choiceRanges);

        // if (minEndIdx > endIdx) break;

        // bool hasShift = false;
        // while (minEndIdx < endIdx) {
        //     if (TdmDB::getClosestChoice(minEndIdx) == choice) {
        //         minEndIdx++;
        //         hasShift = true;
        //     } else {
        //         break;
        //     }
        // }

        // prevEndIdx = hasShift ? minEndIdx : minEndIdx - 1;
        // for (unsigned curEndIdx = minEndIdx; curEndIdx <= endIdx; curEndIdx++) {
        //     if (_vars[curEndIdx - 1]->getVal() <= choice) prevEndIdx = curEndIdx;
        // }

        double curBest = DBL_MAX;
        // for (unsigned curEndIdx = minEndIdx; curEndIdx <= endIdx; curEndIdx++) {
        for (unsigned curEndIdx = idx + 1; curEndIdx <= endIdx; curEndIdx++) {
            double currentCost = 0;
            for (unsigned i = idx; i < curEndIdx; i++) currentCost += abs(_vars[i]->getVal() - choice);

            double cost = legalizeTronconDisp(curEndIdx, p - 1, memorization, choiceRanges);
            double ret = cost == -1 ? DBL_MAX : currentCost + cost;

            if (currentCost > curBest) break;

            if (ret < curBest) {
                curBest = ret;

                if (ret < bestCost) {
                    bestCost = ret;
                    bestChoice = choice;
                    bestEndIdx = curEndIdx;
                }
            }
        }
    }

    memorization.saveBest(idx, p, bestCost, bestChoice, bestEndIdx);

    return bestCost;
}

double WireData::legalizeTronconDispPrune(unsigned idx,
                                          int p,
                                          Memorization &memorization,
                                          const vector<pair<int, int>> &choiceRanges,
                                          const vector<int> &minWire) const {
    if (idx >= _vars.size()) return 0;
    // if (p == 0 || p < minWire[idx]) return -1;
    if (p == 0) return -1;

    if (memorization.isOptimized(idx, p)) return memorization.getBestCost(idx, p);

    double bestCost = DBL_MAX;
    int bestChoice = -1;
    int bestEndIdx = -1;

    if (getMinWireRequire(idx, choiceRanges) > p) {
        memorization.saveBest(idx, p, bestCost, bestChoice, bestEndIdx);
        return bestCost;
    }

    int candBeg = choiceRanges[idx].first, candEnd = choiceRanges[idx].second;
    int prevEndIdx = idx;

    for (int c = candBeg; c <= candEnd; c++) {
        if (c < TdmDB::getClosestChoiceIdx(_vars[idx]->getVal())) continue;
        int choice = TdmDB::getXdrChoice(c);

        unsigned minEndIdx = prevEndIdx + 1;
        unsigned endIdx = calcEndIdx(idx, choice, choiceRanges);

        if (minEndIdx > endIdx) break;

        while (minEndIdx < endIdx) {
            if (TdmDB::getClosestChoice(_vars[minEndIdx]->getVal()) == choice) {
                minEndIdx++;
            } else {
                break;
            }
        }

        prevEndIdx = minEndIdx - 1;
        for (unsigned curEndIdx = minEndIdx; curEndIdx < endIdx; curEndIdx++) {
            if (_vars[curEndIdx]->getVal() <= choice)
                prevEndIdx = curEndIdx;
            else
                break;
        }

        double curBest = DBL_MAX;
        for (unsigned curEndIdx = minEndIdx; curEndIdx <= endIdx; curEndIdx++) {
            // for (unsigned curEndIdx = idx + 1; curEndIdx <= endIdx; curEndIdx++) {
            double currentCost = 0;
            for (unsigned i = idx; i < curEndIdx; i++) currentCost += abs(_vars[i]->getVal() - choice);

            double cost = legalizeTronconDispPrune(curEndIdx, p - 1, memorization, choiceRanges, minWire);
            double ret = cost == -1 ? DBL_MAX : currentCost + cost;

            if (currentCost > curBest) break;

            if (ret < curBest) {
                curBest = ret;

                if (ret < bestCost) {
                    bestCost = ret;
                    bestChoice = choice;
                    bestEndIdx = curEndIdx;
                }
            }
        }
    }

    memorization.saveBest(idx, p, bestCost, bestChoice, bestEndIdx);

    return bestCost;
}

double WireData::legalizeTronconMaxDisp(Memorization &memorization) const {
    double maxMinDisp = 0;
    for (auto var : _vars) {
        maxMinDisp = max(maxMinDisp, abs(TdmDB::getClosestChoice(var->getVal()) - var->getVal()));
    }

    double bnd = maxMinDisp;
    double loBnd = maxMinDisp;
    double hiBnd = -1;

    vector<int> minWire(_vars.size(), 0);

    while (true) {
        vector<pair<int, int>> choiceRanges(_vars.size());
        for (int i = 0, sz = _vars.size(); i < sz; i++) {
            choiceRanges[i].first = TdmDB::getCeilChoiceIdx(_vars[i]->getVal() - bnd);
            choiceRanges[i].second = TdmDB::getFloorChoiceIdx(_vars[i]->getVal() + bnd);
        }

        if (getMinWireRequire(0, choiceRanges) > _troncon->_limit) {
            loBnd = bnd;
            if (hiBnd == -1) {
                bnd = bnd * 2;
            } else {
                bnd = (bnd + hiBnd) / 2;
            }
        } else {
            if (abs(loBnd - hiBnd) < 0.01) break;

            hiBnd = bnd;
            bnd = (loBnd + hiBnd) / 2;
        }
    }

    vector<pair<int, int>> choiceRanges(_vars.size());
    for (int i = 0, sz = _vars.size(); i < sz; i++) {
        choiceRanges[i].first = TdmDB::getCeilChoiceIdx(_vars[i]->getVal() - hiBnd);
        choiceRanges[i].second = TdmDB::getFloorChoiceIdx(_vars[i]->getVal() + hiBnd);
    }

    for (unsigned i = 0; i < _vars.size(); i++) {
        minWire[i] = getMinWireRequire(i, choiceRanges);
    }

    for (int i = 0, sz = _vars.size(); i < sz; i++) {
        double val = _vars[i]->getVal();
        assert(abs(val - TdmDB::getXdrChoice(choiceRanges[i].first)) <= hiBnd);
        assert(abs(val - TdmDB::getXdrChoice(choiceRanges[i].second)) <= hiBnd);
    }

    return legalizeTronconDispPrune(0, _troncon->_limit, memorization, choiceRanges, minWire);
    // return legalizeTronconDisp(0, _troncon->_limit, memorization, choiceRanges);
}

void WireData::recoverSol(const Memorization &memorization) {
    for (int i = 0, sz = _vars.size(), limit = _troncon->_limit; i < sz;) {
        int endIdx = memorization.getBestEndIdx(i, limit);
        int choice = memorization.getBestChoice(i, limit);
        for (int j = i; j < endIdx; j++) _vars[j]->setVal(choice);
        i = endIdx;
        limit--;
    }
}

void WireData::dumpSol(vector<double> &result,
                       const map<XdrVar *, int> &varToOptIdx,
                       const Memorization &memorization) const {
    for (int i = 0, sz = _vars.size(), limit = _troncon->_limit; i < sz;) {
        int endIdx = memorization.getBestEndIdx(i, limit);
        int choice = memorization.getBestChoice(i, limit);
        for (int j = i; j < endIdx; j++) {
            const auto iter = varToOptIdx.find(_vars[j]->_var);
            result[iter->second] = choice;
        }
        i = endIdx;
        limit--;
    }
}

void WireData::calcDisp(const Memorization &memorization) {
    _totDisp = 0;
    _maxDisp = 0;

    for (int i = 0, sz = _vars.size(), limit = _troncon->_limit; i < sz;) {
        int endIdx = memorization.getBestEndIdx(i, limit);
        int choice = memorization.getBestChoice(i, limit);

        for (int j = i; j < endIdx; j++) {
            double disp = abs(_vars[j]->getVal() - choice);
            _totDisp += disp;
            _maxDisp = max(_maxDisp, disp);
        }
        i = endIdx;
        limit--;
    }
}

bool OneWireData::isForward() const { return _var->isForward(); }

void OneWireData::setVal(double val) { _var->setVal(val); }

double OneWireData::getVal() const { return _var->getVal(); }
