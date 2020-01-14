#pragma once

#include "global.h"

class Edge;
class XdrVar;

class SwapHist {
public:
    SwapHist() {}
    SwapHist(XdrVar *u, XdrVar *v);
    bool isSame(const SwapHist &swap1) const;

private:
    XdrVar *_u;
    XdrVar *_v;
    double _uVal;
    double _vVal;
};

class TdmRefine {
public:
    void solve();

private:
    bool optimizePath(vector<Edge *> &criticalPath);
    bool optimizeTroncon(vector<XdrVar *> &vars, int numOptVar);
    bool optimizeTronconSort(vector<XdrVar *> &vars, int numOptVar);
    bool hasSwap(const SwapHist &swap1) const;

    vector<SwapHist> _hist;
};