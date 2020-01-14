#pragma once

#include "global.h"

class TdmDB;
class Edge;
class Troncon;
class XdrVar;
class TimingGraph;
class TdmLagMultiplierUpdater;
class TdmLagMultiplierInitializer;
class Node;

class TdmLagData {
public:
    TdmLagData();

    map<Troncon *, vector<int>> _tronconToXdrVar;
    map<Troncon *, int> _tronconToIdx;
    TimingGraph *_timingGraph;
    vector<XdrVar *> &_optXdrVars;

    vector<double> _lambda;
    vector<double> _mu;
    vector<double> _lastMu;
    const double _epsilon = 0.00001;  // cannot be too small due to precision of cplex
    double _maxChoice;

    double &getMu(Edge *edge);
    double &getLambda(Troncon *troncon);
    double getMuVal(Edge *edge) const;
    double getLambdaVal(Troncon *troncon) const;

    void reportLagMultiplier();
    bool isLagMultiplierLegal();

    double getMuGrad(Edge *edge);
    double getLambdaGrad(Troncon *troncon);
};

class TdmLagSolver {
public:
    void solve();

private:
    const int _nIter = setting.lagIter;
    TdmLagData _tdmLagData;

    void initLagMultiplier();
    void updateMultiplier(int iter);

    void solveLRS();

    double computeDual(bool computeDual);
};

class TdmLagMultiplierUpdater {
public:
    TdmLagMultiplierUpdater(TdmLagData &tdmLagData);
    void run(int iter);

private:
    TdmLagData &_tdmLagData;

    vector<bool> _preserve;

    double getMuStepSize();
    double getLambdaStepSize();
    void getRatio(int iter);

    void updateMu();
    void updateLambda();

    void removeAccIssue(Node *node);
    void critFlow(Node *node, double driverSum, double fanoutSum);
    void decreaseFlow(double driverSum,
                       double fanoutSum,
                       vector<pair<Edge *, double>> &gradients,
                       int startIdx = 0,
                       Node *node = NULL);
    void increaseFlow(Node *node, double driverSum, double fanoutSum);
    void sinkFlow(double driverSum, double fanoutSum, vector<pair<Edge *, double>> &gradients);

    double _ratio = 0;

    const double _changeRate = 0.01;
    const double _baseRate = 0.2;
};

class TdmLagMultiplierInitializer {
public:
    TdmLagMultiplierInitializer(TdmLagData &tdmLagData) : _tdmLagData(tdmLagData) {}
    void run();

private:
    TdmLagData &_tdmLagData;

    void initMu();
    void initLambda();
};
