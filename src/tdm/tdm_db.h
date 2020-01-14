#pragma once

#include "global.h"

class TdmNet;
class TimingGraph;
class XdrVar;
class Troncon;

namespace db {
class Group;
}

class TdmDB {
public:
    void init(int nDevice, vector<db::Group> *groups);

    void updateTiming();
    double getArrivalTime() const;
    TimingGraph *getTimingGraph() const { return _timingGraph; }

    int getNumTroncon() const { return _nTroncon; }
    Troncon *getTroncon(int i) const { return getTroncon(_idxToTroncon[i].first, _idxToTroncon[i].second); }
    Troncon *getTroncon(int device1, int device2) const {
        return _troncons[min(device1, device2)][max(device1, device2)];
    }
    Troncon *getTroncon(XdrVar *xdrVar) const;
    Troncon *getTroncon(TdmNet *tdmNet) const;

    int getLimitVio() const;
    double getContLimitVio() const;
    double getChoiceVio() const;
    int getLimitVio(Troncon *troncon) const;
    double getContLimitVio(Troncon *troncon) const;
    double getChoiceVio(Troncon *troncon) const;

    vector<XdrVar *> &getXdrVars() { return _xdrVars; }
    vector<XdrVar *> &getOptXdrVars() { return _optXdrVars; }
    bool isOptVar(XdrVar *xdrVar) const;

    static vector<int> &getXdrChoices() { return _xdrChoices; }
    static int getXdrChoice(int i) { return _xdrChoices[i]; }
    static int getNumChoices() { return _xdrChoices.size(); }
    static int getClosestChoice(double val);
    static int getClosestChoiceIdx(double val);
    static int getCeilChoice(double val);
    static int getCeilChoiceIdx(double val);
    static int getFloorChoice(double val);
    static int getFloorChoiceIdx(double val);

    db::Group &getGroup(int instId);

    void checkFeasibility(int &nChoice);

    void report() const;
    void reportSol();
    void reportTdmAssignment() const;
    void reportTronconUsage() const;

    void saveSol(vector<double> &sol) const;
    void recoverSol(vector<double> &sol);
    void writeSol(string filename) const;
    void readSol(string filename);

private:
    void constructTimingGraph();
    bool isOptVar(int idx) const { return _isOptVar[idx]; };

    int _nDevice;
    int _nTroncon;
    vector<int> _instToDevice;
    vector<TdmNet *> _nets;
    vector<db::Group> *_groups;
    TimingGraph *_timingGraph;
    vector<XdrVar *> _xdrVars;
    vector<XdrVar *> _optXdrVars;
    vector<bool> _isOptVar;
    vector<vector<Troncon *>> _troncons;
    vector<pair<int, int>> _idxToTroncon;

    static vector<int> _xdrChoices;
    const int _tronconLimit = 20;
    const int _maxChoice = 1600;
};

class XdrVar {
public:
    XdrVar(TdmNet *net, bool _forward, int _id);

    double getVal() const { return _val; }
    void setVal(double val) { _val = val; }
    TdmNet *getNet() const { return _net; }
    bool isForward() const { return _forward; }
    bool isLegal() const { return _val == (int)_val && (_val == 1 || ((int)_val % 8 == 0)); }
    int _id;

private:
    double _val;
    TdmNet *_net;
    bool _forward;
};

class Troncon {
public:
    void addNet(TdmNet *net);
    TdmNet *getNet(int i) const { return _nets[i]; }
    int getNumNets() const { return _nets.size(); }
    vector<TdmNet *> &getNets() { return _nets; }

    Troncon(int device1, int device2, int limit) {
        _devices.first = device1;
        _devices.second = device2;
        _limit = limit;
    }

    double getContUsage() const;
    int getUsage() const;

    pair<int, int> _devices;
    int _limit;
    int _id;

private:
    vector<TdmNet *> _nets;
};

extern TdmDB tdmDatabase;