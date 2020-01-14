#include "tdm_solve_lag.h"
#include "tdm_db.h"
#include "tdm_net.h"
#include "timing_graph.h"

void TdmLagMultiplierInitializer::initMu() {
    // averaging the flow related to xdr edge
    TimingGraph *timingGraph = _tdmLagData._timingGraph;
    vector<vector<Node *>> &levels = timingGraph->getLevels();
    vector<double> nodePrecXdrEdge(timingGraph->getNumNodes(), 0), edgePrecXdrEdge(timingGraph->getNumEdges(), 0);
    for (int l = 0, sz = levels.size(); l < sz; l++) {
        vector<Node *> &level = levels[l];
        for (auto node : level) {
            for (auto driver : node->_drivers) {
                edgePrecXdrEdge[driver->_id] += timingGraph->isOptEdge(driver);
                nodePrecXdrEdge[node->_id] += edgePrecXdrEdge[driver->_id];
            }
            for (auto fanout : node->_fanouts) {
                edgePrecXdrEdge[fanout->_id] = nodePrecXdrEdge[node->_id] / node->_fanouts.size();
            }
        }
    }

    vector<vector<Node *>> &revLevels = timingGraph->getRevLevels();
    for (int l = 0, sz = revLevels.size(); l < sz; l++) {
        vector<Node *> &level = revLevels[l];

        for (auto node : level) {
            double fanoutSum = 0;
            if (node == timingGraph->getSink()) {
                fanoutSum = 1;
            } else {
                for (auto fanout : node->_fanouts) fanoutSum += _tdmLagData.getMu(fanout);
            }

            if (nodePrecXdrEdge[node->_id] == 0) {
                for (auto driver : node->_drivers) _tdmLagData.getMu(driver) = fanoutSum / node->_drivers.size();
            } else {
                for (auto driver : node->_drivers)
                    _tdmLagData.getMu(driver) = fanoutSum * edgePrecXdrEdge[driver->_id] / nodePrecXdrEdge[node->_id];
            }
        }
    }
}

void TdmLagMultiplierInitializer::initLambda() {
    // init to troncon limit
    auto &tronconToXdrVar = _tdmLagData._tronconToXdrVar;
    TimingGraph *timingGraph = _tdmLagData._timingGraph;
    auto &optXdrVars = _tdmLagData._optXdrVars;
    for (auto &pair : tronconToXdrVar) {
        Troncon *troncon = pair.first;
        double sum = 0;
        double maxTmpSum = 0;
        for (auto var : pair.second) {
            double tmpSum = 0;
            for (auto edge : timingGraph->getEdges(optXdrVars[var])) tmpSum += edge->_tdmCoef * _tdmLagData.getMu(edge);
            sum += sqrt(tmpSum);
            maxTmpSum = max(maxTmpSum, tmpSum);
        }
        double &lambda = _tdmLagData.getLambda(troncon);
        lambda = pow(sum / troncon->_limit, 2);
        lambda = max(lambda, maxTmpSum);
    }
}

void TdmLagMultiplierInitializer::run() {
    initMu();
    initLambda();

    assert(_tdmLagData.isLagMultiplierLegal());

    const auto &mu = _tdmLagData._mu;
    const auto &lambda = _tdmLagData._lambda;

    int cnt1 = 0, cnt2 = 0;
    for (unsigned i = 0; i < mu.size(); i++) cnt1 += mu[i] == 0;
    for (unsigned i = 0; i < lambda.size(); i++) cnt2 += lambda[i] == 0;

    printlog(LOG_INFO,
             "finish initialize LAG multiplier, (=0/tot)=mu:%d/%lu, lambda:%d/%lu",
             cnt1,
             mu.size(),
             cnt2,
             lambda.size());
}
