#include "tdm_solve_lag.h"
#include "tdm_db.h"
#include "tdm_net.h"
#include "timing_graph.h"

double &TdmLagData::getMu(Edge *edge) { return _mu[edge->_id]; }

double &TdmLagData::getLambda(Troncon *troncon) {
    auto iter = _tronconToIdx.find(troncon);
    assert(iter != _tronconToIdx.end());
    return _lambda[iter->second];
}

double TdmLagData::getMuVal(Edge *edge) const { return _mu[edge->_id]; }

double TdmLagData::getLambdaVal(Troncon *troncon) const {
    auto iter = _tronconToIdx.find(troncon);
    assert(iter != _tronconToIdx.end());
    return _lambda[iter->second];
}

void TdmLagData::reportLagMultiplier() {
    cout << "mu:";
    for (unsigned i = 0; i < _mu.size(); i++)
        if (_mu[i] != 0) cout << _mu[i] << " ";
    cout << endl;

    cout << "lambda:";
    // for (auto pair : _tronconToIdx)
    //     cout << pair.second << "," << _lambda[pair.second] << "," << _tronconToXdrVar[pair.first].size() << " ";
    for (unsigned i = 0; i < _lambda.size(); i++)
        if (_lambda[i] != 0) cout << _lambda[i] << " ";
    cout << endl;

    cout << endl;
}

double TdmLagData::getMuGrad(Edge *edge) { return edge->getArrivalTimeAlongEdge() - edge->_fanout->getArrivalTime(); }

double TdmLagData::getLambdaGrad(Troncon *troncon) { return troncon->getContUsage() - troncon->_limit; }

bool TdmLagData::isLagMultiplierLegal() {
    for (int i = 0, sz = _timingGraph->getNumNodes(); i < sz; i++) {
        Node *node = _timingGraph->getNode(i);

        double sum1 = 0, sum2 = 0;

        for (auto edge : node->_drivers) sum1 += _mu[edge->_id];
        for (auto edge : node->_fanouts) sum2 += _mu[edge->_id];
        if (node == _timingGraph->getSink()) {
            if (abs(sum1 - 1) > _epsilon) return false;
        } else if (node != _timingGraph->getSource()) {
            if (abs(sum1 - sum2) > _epsilon) return false;
        }
    }

    return true;
}

TdmLagData::TdmLagData() : _optXdrVars(tdmDatabase.getOptXdrVars()) {
    _timingGraph = tdmDatabase.getTimingGraph();
    for (unsigned i = 0; i < _optXdrVars.size(); i++) {
        Troncon *troncon = _optXdrVars[i]->getNet()->_troncon;
        _tronconToXdrVar[troncon].push_back(i);
    }

    int cnt = 0;
    for (auto &pair : _tronconToXdrVar) _tronconToIdx[pair.first] = cnt++;

    _mu.assign(_timingGraph->getNumEdges(), 0);
    _lambda.assign(_tronconToXdrVar.size(), 0);

    _maxChoice = tdmDatabase.getXdrChoices().back();
}