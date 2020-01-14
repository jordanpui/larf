#include "tdm_db.h"
#include "db/db.h"
#include "db/group.h"
#include "tdm_net.h"
#include "tdm_part.h"
#include "timing_graph.h"

vector<int> TdmDB::_xdrChoices;
TdmDB tdmDatabase;

void TdmDB::init(int nDevice, vector<db::Group>* groups) {
    _nDevice = nDevice;
    _groups = groups;

    // gen mapping from inst to device
    _instToDevice.resize(_groups->size());
    ifstream instDeviceFile("instance.device");
    for (unsigned i = 0; i < _groups->size(); i++) {
        string name;
        instDeviceFile >> name >> _instToDevice[i];
    }
    instDeviceFile.close();

    // gen all the tdm nets
    getTdmNets(_instToDevice, _nets);

    // gen troncon, xdr var
    _troncons.resize(_nDevice, vector<Troncon*>(_nDevice, NULL));
    for (int i = 0; i < _nDevice; i++) {
        for (int j = i + 1; j < _nDevice; j++) {
            _troncons[i][j] = new Troncon(i, j, _tronconLimit);
        }
    }
    for (auto net : _nets) {
        if (net->isInterNet()) {
            int fromDevice = net->getFromDevice();
            int toDevice = net->getToDevice();

            Troncon* troncon = getTroncon(fromDevice, toDevice);
            troncon->addNet(net);

            XdrVar* xdrVar = new XdrVar(net, fromDevice < toDevice, _xdrVars.size());
            _xdrVars.push_back(xdrVar);
        }
    }
    for (int i = 0; i < _nDevice; i++) {
        for (int j = i + 1; j < _nDevice; j++) {
            if (_troncons[i][j]->getNumNets() > 0) {
                _troncons[i][j]->_id = _idxToTroncon.size();
                _idxToTroncon.emplace_back(i, j);
            } else {
                delete _troncons[i][j];
                _troncons[i][j] = NULL;
            }
        }
    }
    _nTroncon = _idxToTroncon.size();

    // gen xdr choices
    _xdrChoices = {1};
    for (int i = 1; i * 8 <= _maxChoice; i++) _xdrChoices.push_back(i * 8);

    // gen xdrVar need to be optimized
    for (int i = 0; i < _nTroncon; i++) {
        Troncon* troncon = getTroncon(i);
        if (troncon->getNumNets() > troncon->_limit) {
            for (auto net : troncon->getNets()) _optXdrVars.push_back(net->getXdrVar());
        } else {
            for (auto net : troncon->getNets()) net->getXdrVar()->setVal(1);
        }
    }
    _isOptVar.assign(_xdrVars.size(), false);
    for (auto var : _optXdrVars) _isOptVar[var->_id] = true;

    // construct timing graph
    constructTimingGraph();

    // report();
}

void TdmDB::report() const {
    log() << "---- reprot TdmDB ----" << endl;

    int maxXdrNum = INT_MIN, minXdrNum = INT_MAX;
    for (auto pair : _idxToTroncon) {
        Troncon* troncon = getTroncon(pair.first, pair.second);
        maxXdrNum = max(maxXdrNum, troncon->getNumNets());
        minXdrNum = min(minXdrNum, troncon->getNumNets());
    }

    printlog(LOG_INFO,
             "#troncon=%lu, avgXdrNum=%f, maxXdrNum=%d, minXdrNum=%d",
             _idxToTroncon.size(),
             _xdrVars.size() * 1.0 / _idxToTroncon.size(),
             maxXdrNum,
             minXdrNum);

    int nOptimalTroncon = 0;
    for (int i = 0; i < _nTroncon; i++)
        if (getTroncon(i)->getNumNets() <= getTroncon(i)->_limit) nOptimalTroncon++;
    printlog(LOG_INFO, "active_troncon=%d, varNeedOpt=%lu", _idxToTroncon.size() - nOptimalTroncon, _optXdrVars.size());

    _timingGraph->report();
    reportTronconUsage();
}

void TdmDB::reportSol() {
    updateTiming();

    double choiceVio = getChoiceVio();

    if (choiceVio == 0) {
        printlog(LOG_INFO, "at=%f, LimitVio=%d, choiceVio=%f", getArrivalTime(), getLimitVio(), choiceVio);
    } else {
        printlog(LOG_INFO, "at=%f, contLimitVio=%f, choiceVio=%f", getArrivalTime(), getContLimitVio(), choiceVio);
    }

    // for (auto xdr : _optXdrVars) cout << xdr->getVal() << " ";
    // cout << endl;
}

Troncon* TdmDB::getTroncon(XdrVar* xdrVar) const { return getTroncon(xdrVar->getNet()); }

Troncon* TdmDB::getTroncon(TdmNet* tdmNet) const { return getTroncon(tdmNet->getFromDevice(), tdmNet->getToDevice()); }

void TdmDB::constructTimingGraph() {
    _timingGraph = new TimingGraph();
    for (auto& group : *_groups) {
        _timingGraph->addNode(group.instances[0]);
    }

    for (auto net : _nets) {
        db::Instance* driver = net->getPin(net->getDriverIdx())->instance;

        set<db::Instance*> instances;
        for (auto pin : net->getPins()) instances.insert(pin->instance);

        for (auto instance : instances) {
            if (instance == driver) continue;
            Edge* edge = _timingGraph->addEdge(driver->id, instance->id, net);

            if (net->isInterNet()) _timingGraph->addMapping(net->getXdrVar(), edge);
        }
    }

    _timingGraph->breakCycle();
    _timingGraph->setConstDelay();
    _timingGraph->setSrcSink();
    _timingGraph->levelize();
    _timingGraph->removeAbnEdges();
}

void TdmDB::updateTiming() {
    _timingGraph->updateArrivalTime();
    _timingGraph->updateRequireTime();
}

double TdmDB::getArrivalTime() const { return _timingGraph->getSinkAT(); }

XdrVar::XdrVar(TdmNet* net, bool forward, int id) {
    _net = net;
    if (net->_xdrVar) assert(false);
    net->_xdrVar = this;
    _val = 1;
    _forward = forward;
    _id = id;
}

double TdmDB::getContLimitVio() const {
    double vio = 0;
    for (int t = 0; t < getNumTroncon(); t++) vio += getContLimitVio(getTroncon(t));
    return vio;
}

int TdmDB::getLimitVio() const {
    int vio = 0;
    for (int t = 0; t < getNumTroncon(); t++) vio += getLimitVio(getTroncon(t));
    return vio;
}

double TdmDB::getChoiceVio() const {
    double vio = 0;
    for (int t = 0; t < getNumTroncon(); t++) vio += getChoiceVio(getTroncon(t));
    return vio;
}

double TdmDB::getContLimitVio(Troncon* troncon) const { return max(0.0, troncon->getContUsage() - troncon->_limit); }

int TdmDB::getLimitVio(Troncon* troncon) const { return max(0, troncon->getUsage() - troncon->_limit); }

double TdmDB::getChoiceVio(Troncon* troncon) const {
    double vio = 0;

    for (auto net : troncon->getNets()) {
        double val = net->getXdrVar()->getVal();
        int closestChoice = getClosestChoice(val);

        vio += abs(val - closestChoice);

        // if (abs(val - closestChoice) != 0) {
        //     cout << vio << " " << val << endl;
        //     getchar();
        // }
    }

    return vio;
}

int TdmDB::getClosestChoiceIdx(double val) {
    if (val <= _xdrChoices.front()) return 0;
    if (val >= _xdrChoices.back()) return _xdrChoices.size() - 1;

    int lb = floor(val / 8);
    int ub = ceil(val / 8);

    if (lb == ub)
        return lb;
    else if (abs(val - _xdrChoices[lb]) > abs(val - _xdrChoices[ub]))
        return ub;
    else if (abs(val - _xdrChoices[lb]) < abs(val - _xdrChoices[ub]))
        return lb;
    else
        return ub;
}

int TdmDB::getClosestChoice(double val) { return _xdrChoices[getClosestChoiceIdx(val)]; }

int TdmDB::getCeilChoiceIdx(double val) {
    val = max(val, (double)_xdrChoices.front());

    if (val <= _xdrChoices.front()) return 0;
    return ceil(val / 8);
}

int TdmDB::getCeilChoice(double val) { return _xdrChoices[getCeilChoiceIdx(val)]; }

int TdmDB::getFloorChoiceIdx(double val) {
    val = min(val, (double)_xdrChoices.back());

    if (val <= _xdrChoices.front()) return 0;
    return floor(val / 8);
}

int TdmDB::getFloorChoice(double val) { return _xdrChoices[getFloorChoiceIdx(val)]; }

void Troncon::addNet(TdmNet* net) {
    _nets.push_back(net);
    net->_troncon = this;
}

double Troncon::getContUsage() const {
    double usage = 0;

    for (auto net : _nets) {
        double val = net->getXdrVar()->getVal();
        usage += 1 / val;
    }

    return usage;
}

int Troncon::getUsage() const {
    unordered_map<double, int> forwardTdmCount, backwardTdmCount;
    int usage = 0;

    for (auto net : _nets) {
        double val = net->getXdrVar()->getVal();
        bool isLegal = net->getXdrVar()->isLegal();
        assert(isLegal);

        bool isForward = net->getXdrVar()->isForward();

        if (isForward)
            forwardTdmCount[val]++;
        else
            backwardTdmCount[val]++;
    }
    for (auto pair : forwardTdmCount) usage += ceil(pair.second / pair.first);
    for (auto pair : backwardTdmCount) usage += ceil(pair.second / pair.first);

    return usage;
}

void TdmDB::reportTdmAssignment() const {
    log() << "report xdr assignment of opt vars" << endl;
    for (auto var : _optXdrVars) cout << var->getVal() << " ";
    cout << endl << endl;
}

bool TdmDB::isOptVar(XdrVar* xdrVar) const { return isOptVar(xdrVar->_id); }

void TdmDB::checkFeasibility(int& nChoice) {
    map<Troncon*, vector<int>> tronconToXdrVar;
    for (int i = 0, sz = _optXdrVars.size(); i < sz; i++) {
        Troncon* troncon = _optXdrVars[i]->getNet()->_troncon;
        tronconToXdrVar[troncon].push_back(i);
    }
    for (auto& pair : tronconToXdrVar) {
        int numForwVars = 0, numBackVars = 0;
        for (auto i : pair.second) {
            if (_optXdrVars[i]->isForward())
                numForwVars++;
            else
                numBackVars++;
        }
        int limit = pair.first->_limit;
        int forwUsage = ceil(numForwVars * 1.0 / _xdrChoices[nChoice - 1]);
        int backUsage = ceil(numBackVars * 1.0 / _xdrChoices[nChoice - 1]);
        if (forwUsage + backUsage > limit) {
            printlog(LOG_WARN,
                     "troncon#%d: #forwVars=%d #backVars=%d, max_choice=%d, forwUsage/backUsage/Limit=%d/%d/%d",
                     pair.first->_id,
                     numForwVars,
                     numBackVars,
                     _xdrChoices[nChoice - 1],
                     forwUsage,
                     backUsage,
                     limit);
            int idx = nChoice - 1;
            for (; nChoice * 8 <= 1600; nChoice++) {
                _xdrChoices.push_back(nChoice * 8);
                forwUsage = ceil(numForwVars * 1.0 / _xdrChoices[nChoice]);
                backUsage = ceil(numBackVars * 1.0 / _xdrChoices[nChoice]);
                if (forwUsage + backUsage <= limit) {
                    printlog(LOG_WARN, "max_choice increase:%d->%d", idx * 8, _xdrChoices[nChoice]);
                    nChoice++;
                    break;
                }
            }
        }
    }
}

void TdmDB::reportTronconUsage() const {
    log() << "---- report troncon vars ----" << endl;
    map<Troncon*, vector<int>> tronconToXdrVar;
    for (int i = 0, sz = _optXdrVars.size(); i < sz; i++) {
        Troncon* troncon = _optXdrVars[i]->getNet()->_troncon;
        tronconToXdrVar[troncon].push_back(i);
    }

    for (auto& pair : tronconToXdrVar) {
        int numForwVars = 0, numBackVars = 0;
        for (auto i : pair.second) {
            if (_optXdrVars[i]->isForward())
                numForwVars++;
            else
                numBackVars++;
        }
        int limit = pair.first->_limit;

        printlog(LOG_INFO,
                 "troncon#%d: #forwVars=%d #backVars=%d, limit=%d",
                 pair.first->_id,
                 numForwVars,
                 numBackVars,
                 limit);
    }
}

db::Group& TdmDB::getGroup(int instId) { return (*_groups)[instId]; }

void TdmDB::saveSol(vector<double>& sol) const {
    sol.resize(_xdrVars.size());
    for (int i = 0, sz = _xdrVars.size(); i < sz; i++) sol[i] = _xdrVars[i]->getVal();
}

void TdmDB::recoverSol(vector<double>& sol) {
    for (int i = 0, sz = _xdrVars.size(); i < sz; i++) _xdrVars[i]->setVal(sol[i]);
}

void TdmDB::writeSol(string filename) const {
    ofstream file(filename);
    for (int i = 0, sz = _xdrVars.size(); i < sz; i++) file << _xdrVars[i]->getVal() << endl;
    file.close();
}

void TdmDB::readSol(string filename) {
    ifstream file(filename);
    for (int i = 0, sz = _xdrVars.size(); i < sz; i++) {
        double val;
        file >> val;
        _xdrVars[i]->setVal(val);
    }
    file.close();
}
