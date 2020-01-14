#include "tdm_solve_lag.h"
#include "tdm_db.h"
#include "tdm_net.h"
#include "timing_graph.h"
#include "db/db.h"

void TdmLagSolver::solve() {
    log() << "==================== begin TDM analytical solving ====================" << endl;
    TimingGraph *timingGraph = _tdmLagData._timingGraph;

    ofstream convergeFile(db::database.bmName + ".curve");
    convergeFile << _nIter << endl;

    vector<double> bestSol;
    int bestIter = -1;
    double bestCost = DBL_MAX;

    vector<int> bestVals(_nIter, 0);
    for (int i = 0; i < _nIter; i++) {
        // log() << "======== iter " << i << " ========" << endl;

        if (i == 0)
            initLagMultiplier();
        else
            updateMultiplier(i);

        // reportLagMultiplier();

        solveLRS();

        double dualVal = computeDual(setting.computeDual);
        double primVal = timingGraph->getSinkAT();

        // log() << "primal:" << primVal << ", dual:" << dualVal << endl;
        convergeFile << primVal << " " << dualVal << endl;

        // tdmDatabase.reportSol();

        const int interval = 70;
        if (i - interval >= 0 && (bestVals[i - interval] - bestCost) < 1) {
            log() << "early break " << interval << " " << bestVals[i - interval] << " " << bestCost << endl;
            break;
        }

        if (timingGraph->getSinkAT() < bestCost) {
            bestIter = i;
            tdmDatabase.saveSol(bestSol);
            bestCost = timingGraph->getSinkAT();
        }
        bestVals[i] = bestCost;
    }

    convergeFile.close();

    double dual = computeDual(true);

    tdmDatabase.recoverSol(bestSol);
    timingGraph->updateArrivalTime();
    double primal = timingGraph->getSinkAT();
    dual = max(dual, computeDual(true));
    printlog(LOG_INFO,
             "recover solution from iter#%d, primal=%f, dual=%f, gap=%f",
             bestIter,
             primal,
             dual,
             (primal - dual) / dual);

    log() << "---------------- finish TDM analytical solving ----------------" << endl;
}

double TdmLagSolver::computeDual(bool computeDual) {
    if (!computeDual) return 0;
    TimingGraph *timingGraph = _tdmLagData._timingGraph;

    double result = 0;

    result += timingGraph->getSinkAT();

    for (auto &pair : _tdmLagData._tronconToXdrVar) {
        Troncon *troncon = pair.first;
        result += _tdmLagData.getLambda(troncon) * (troncon->getContUsage() - troncon->_limit);
    }

    for (int e = 0, sz = timingGraph->getNumEdges(); e < sz; e++) {
        Edge *edge = timingGraph->getEdge(e);
        double &mu = _tdmLagData.getMu(edge);

        result += mu * (edge->getArrivalTimeAlongEdge() - edge->_fanout->getArrivalTime());
    }

    return result;
}

void TdmLagSolver::initLagMultiplier() {
    TdmLagMultiplierInitializer initializer(_tdmLagData);
    initializer.run();
}

void TdmLagSolver::updateMultiplier(int iter) {
    TdmLagMultiplierUpdater updater(_tdmLagData);
    updater.run(iter);
}

void TdmLagSolver::solveLRS() {
    auto &optXdrVars = _tdmLagData._optXdrVars;
    auto *timingGraph = _tdmLagData._timingGraph;

    int cnt1 = 0, cnt2 = 0, cnt3 = 0;

    for (auto var : optXdrVars) {
        Troncon *troncon = tdmDatabase.getTroncon(var);
        double &lambda = _tdmLagData.getLambda(troncon);

        double sum = 0;
        for (auto edge : timingGraph->getEdges(var)) sum += edge->_tdmCoef * _tdmLagData.getMu(edge);

        if (sum != 0) {
            double newVal = sqrt(lambda / sum);
            if (newVal < 1 || newVal > _tdmLagData._maxChoice) {
                int nCrit = 0;
                for (auto edge : timingGraph->getEdges(var)) nCrit += edge->isCritical();
                // printlog(LOG_WARN,
                //          "OoBnd var: id=%d, mu=%.3f->%.3f, #edge(crit)=%lu(%d)",
                //          var->_id,
                //          var->getVal(),
                //          newVal,
                //          timingGraph->getEdges(var).size(),
                //          nCrit);
            }

            if (newVal < 1) {
                var->setVal(1);
                cnt1++;
            } else if (newVal > _tdmLagData._maxChoice) {
                var->setVal(_tdmLagData._maxChoice);
                cnt2++;
            } else {
                var->setVal(newVal);
            }
        } else {
            var->setVal(_tdmLagData._maxChoice);
            cnt3++;
        }
    }

    // printlog(LOG_INFO, "LRS: <lb=%d, >ub=%d, zero_mu=%d, total=%lu", cnt1, cnt2, cnt3, optXdrVars.size());

    timingGraph->updateArrivalTime();
}
