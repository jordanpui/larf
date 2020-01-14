#include "timing_graph.h"
#include "db/db.h"
#include "db/group.h"
#include "db/instance.h"
#include "tdm_db.h"
#include "tdm_net.h"
#include "db/site.h"

Edge::Edge(Node* fromNode, Node* toNode, TdmNet* net, int id) {
    _driver = fromNode;
    _fanout = toNode;
    _id = id;
    _net = net;
    _constDelay = fromNode->_delay;

    fromNode->addFanout(this);
    toNode->addDriver(this);
}

void Edge::updateDelay() {
    if (!_net) {
        _delay = 0;
    } else {
        _delay = _constDelay;
        if (_net->isInterNet()) {
            _delay += _tdmCoef * _net->getXdrVar()->getVal();
        }
    }
    _arrivalTimeAlongEdge = _delay + _driver->getArrivalTime();
}

double Edge::getDelay(int val) const {
    double delay;

    if (!_net) {
        delay = 0;
    } else {
        delay = _constDelay;
        if (_net->isInterNet()) delay += _tdmCoef * val;
    }
    return delay;
}

Node::Node(db::Instance* instance, int id) {
    _id = id;
    _instance = instance;
    _arrivalTime = -1;

    if (instance && instance->IsLUT())
        _delay = 2;
    else
        _delay = 0;
}

vector<Edge*>& TimingGraph::getEdges(XdrVar* var) { return _xdrToEdges[var]; }

vector<Edge*>& TimingGraph::getEdges(TdmNet* net) { return getEdges(net->getXdrVar()); }

Edge* TimingGraph::addEdge(int u, int v, TdmNet* net) { return addEdge(_nodes[u], _nodes[v], net); }

Node* TimingGraph::addNode(db::Instance* instance) {
    Node* node = new Node(instance, _nodes.size());
    _nodes.push_back(node);

    return node;
}
Edge* TimingGraph::addEdge(Node* u, Node* v, TdmNet* net) {
    Edge* edge = new Edge(u, v, net, _edges.size());
    _edges.push_back(edge);
    return edge;
}

void TimingGraph::addMapping(XdrVar* var, Edge* edge) { _xdrToEdges[var].push_back(edge); }

void TimingGraph::breakCycle() {
    int sz = _nodes.size();
    vector<pair<int, int>> edgeToAdd;

    for (int i = 0; i < sz; i++) {
        Node* curNode = _nodes[i];
        db::Instance* inst = curNode->_instance;

        if (inst->IsLUT() || inst->IsIO()) continue;

        Node* pseudoNode = addNode(inst);
        for (auto edge : curNode->_fanouts) {
            edge->_driver = pseudoNode;
            pseudoNode->_fanouts.push_back(edge);
        }
        curNode->_fanouts.clear();
    }

    if (isCyclic()) {
        cout << "error: cannot remove all cycle by breaking the current set of instances" << endl;
        exit(1);
    }

    // printlog(LOG_INFO, "after breaking loops: diff_nodes=%d", getNumNodes() - sz);
}

void TimingGraph::removeAbnEdges() {
    vector<double> intraNetDelay;
    for (auto edge : _edges) {
        if (edge->_net && edge->isConstEdge() && edge->_constDelay > 0) {
            intraNetDelay.push_back(edge->_constDelay);
        }
    }
    double sum = 0;
    for (auto val : intraNetDelay) sum += val;
    sort(intraNetDelay.begin(), intraNetDelay.end());
    double outlierRatio = 0.02;
    double outlierThreshold = intraNetDelay[intraNetDelay.size() * (1 - outlierRatio)];
    for (auto edge : _edges) {
        if (edge->_net && edge->isConstEdge()) {
            if (!edge->_driver->_instance->IsLUTFF() || !edge->_fanout->_instance->IsLUTFF()) edge->_constDelay = 0;
            if (false && edge->_constDelay > outlierThreshold) edge->_constDelay = 0;
        }
    }
}

void TimingGraph::setConstDelay() {
    // some stats
    double sumConstDelay = 0;
    double maxConstDelay = DBL_MIN, minConstDelay = DBL_MAX;
    int cnt = 0;

    for (auto edge : _edges) {
        if (!edge->_net) continue;

        if (edge->_net->isIntraNet()) {
            db::Group& driver = tdmDatabase.getGroup(edge->_driver->_instance->id);
            db::Group& fanout = tdmDatabase.getGroup(edge->_fanout->_instance->id);

            db::Site* driverSite = db::database.getSite(driver.x, driver.y);
            db::Site* fanoutSite = db::database.getSite(fanout.x, fanout.y);

            if (driverSite == fanoutSite && edge->_driver->_instance->IsLUT() && edge->_fanout->_instance->IsFF())
                continue;

            double wireDelay = wireDelayCoef * max(1.0, abs(driver.x - fanout.x) + abs(driver.y - fanout.y));
            edge->_constDelay += wireDelay;

            maxConstDelay = max(wireDelay, maxConstDelay);
            minConstDelay = min(wireDelay, minConstDelay);
            sumConstDelay += wireDelay;
            cnt++;
        } else {
            edge->_constDelay += 0;
        }
    }

    int gateCnt = 0;
    double gateDelay = 0;
    for (auto node : _nodes) {
        if (node->_instance->IsLUT()) {
            gateCnt++;
            gateDelay += node->_delay;
        }
    }

    // printlog(LOG_INFO,
    //          "logicDelay:%.3f(%.3f/%d), wireDelay:%.3f(%.3f/%d), min/max=%.3f/%.3f",
    //          gateDelay / gateCnt,
    //          gateDelay,
    //          gateCnt,
    //          sumConstDelay / cnt,
    //          sumConstDelay,
    //          cnt,
    //          minConstDelay,
    //          maxConstDelay);
}

void TimingGraph::setSrcSink() {
    vector<Node*> pis, pos;

    for (auto node : _nodes) {
        if (node->_drivers.empty()) pis.push_back(node);
        if (node->_fanouts.empty()) pos.push_back(node);
    }

    _source = addNode(NULL);
    _sink = addNode(NULL);

    for (auto node : pis) {
        addEdge(_source, node, NULL);
    }
    for (auto node : pos) {
        addEdge(node, _sink, NULL);
    }
}

void TimingGraph::forwardPropagateST() {
    queue<Node*> q;
    q.push(_source);

    vector<int> numDrivers(_nodes.size(), 0);
    for (auto node : _nodes) {
        numDrivers[node->_id] = node->_drivers.size();
    }

    while (!q.empty()) {
        Node* top = q.front();
        q.pop();

        for (auto edge : top->_fanouts) {
            edge->updateDelay();
            edge->_fanout->updateArrivalTime(edge->getArrivalTimeAlongEdge());
            numDrivers[edge->_fanout->_id]--;
            if (numDrivers[edge->_fanout->_id] == 0) q.push(edge->_fanout);
        }
    }
}

void TimingGraph::forwardPropagateMT() {
    for (unsigned l = 0; l < _levels.size(); l++) {
        std::mutex idx_mutex;
        int vIdx = 0;

        const int batchSize = 20;

        auto updateNodeArrivalTime = [&]() {
            while (true) {
                int idx;
                idx_mutex.lock();
                idx = vIdx;
                vIdx += batchSize;
                idx_mutex.unlock();

                if (idx >= (int)_levels[l].size()) break;

                for (int i = idx; i < min((int)_levels[l].size(), idx + batchSize); i++) {
                    Node* node = _levels[l][i];
                    for (auto driver : node->_drivers) node->updateArrivalTime(driver->getArrivalTimeAlongEdge());
                    for (auto fanout : node->_fanouts) fanout->updateDelay();
                }
            }
        };

        int nThreads = max(1, setting.nThreads);
        std::thread threads[nThreads];
        for (int i = 0; i < nThreads; i++) threads[i] = std::thread(updateNodeArrivalTime);
        for (int i = 0; i < nThreads; i++) threads[i].join();
    }
}

void TimingGraph::backwardPropagateST() {
    queue<Node*> q;
    q.push(_sink);

    vector<int> numFanouts(_nodes.size(), 0);
    for (auto node : _nodes) {
        numFanouts[node->_id] = node->_fanouts.size();
    }

    while (!q.empty()) {
        Node* top = q.front();
        q.pop();

        for (auto edge : top->_drivers) {
            edge->updateDelay();
            edge->_driver->updateRequireTime(top->getRequireTime() - edge->getDelay());
            numFanouts[edge->_driver->_id]--;
            if (numFanouts[edge->_driver->_id] == 0) q.push(edge->_driver);
        }
    }
}

void TimingGraph::backwardPropagateMT() {
    // TODO: implement multi-thread structure
    for (unsigned l = 0; l < _revLevels.size(); l++) {
        for (auto node : _revLevels[l]) {
            for (auto edge : node->_drivers) {
                edge->updateDelay();
                edge->_driver->updateRequireTime(node->getRequireTime() - edge->getDelay());
            }
        }
    }
}

void TimingGraph::updateArrivalTime() {
    resetArrivalTime();
    _source->updateArrivalTime(0);
    forwardPropagateMT();
}

void TimingGraph::updateRequireTime() {
    resetRequireTime();
    _sink->updateRequireTime(getSinkAT());
    backwardPropagateMT();
}

void TimingGraph::getSRCoef(XdrVar* var, double& k, double& b) {
    vector<Edge*> edges = getEdges(var);
    double sinkAT = getSinkAT();

    // Note: select the edge with the worst slack ratio in the current assignment
    double maxSR = DBL_MIN;

    for (auto edge : edges) {
        double curK =
            1 + (edge->_driver->getArrivalTime() - edge->_fanout->getRequireTime() + edge->_constDelay) / sinkAT;
        double curB = edge->_tdmCoef / sinkAT;

        if (curK + curB * var->getVal() > maxSR) {
            k = curK;
            b = curB;
            maxSR = curK + curB * var->getVal();
        }
    }
}

void TimingGraph::resetTiming() {
    resetArrivalTime();
    resetRequireTime();
}

void TimingGraph::resetArrivalTime() {
    for (auto node : _nodes) node->resetArrivalTime();
}

void TimingGraph::resetRequireTime() {
    for (auto node : _nodes) node->resetRequireTime();
}

void TimingGraph::levelize() {
    if (isCyclic()) {
        cout << "Error:cyclic graph" << endl;
        return;
    }

    forLevelize();
    revLevelize();
}

void TimingGraph::revLevelize() {
    queue<Node*> q;
    vector<int> numFanouts(_nodes.size(), 0);
    for (auto node : _nodes) {
        numFanouts[node->_id] = node->_fanouts.size();
    }

    vector<Node*>* level;
    q.push(_sink);

    int curLevelNode = 1;
    int prevLevelNode = 0;

    while (true) {
        if (prevLevelNode == 0) {
            _revLevels.resize(_revLevels.size() + 1);
            level = &_revLevels.back();
            prevLevelNode = curLevelNode;
            curLevelNode = 0;
        }
        Node* curNode = q.front();
        q.pop();

        prevLevelNode--;
        level->push_back(curNode);

        if (curNode == _source) break;

        for (auto edge : curNode->_drivers) {
            numFanouts[edge->_driver->_id]--;
            if (numFanouts[edge->_driver->_id] == 0) {
                q.push(edge->_driver);
                curLevelNode++;
            }
        }
    }
}

void TimingGraph::forLevelize() {
    queue<Node*> q;
    vector<int> numDrivers(_nodes.size(), 0);
    for (auto node : _nodes) {
        numDrivers[node->_id] = node->_drivers.size();
    }

    vector<Node*>* level;
    q.push(_source);

    int curLevelNode = 1;
    int prevLevelNode = 0;

    while (true) {
        if (prevLevelNode == 0) {
            _levels.resize(_levels.size() + 1);
            level = &_levels.back();
            prevLevelNode = curLevelNode;
            curLevelNode = 0;
        }
        Node* curNode = q.front();
        q.pop();

        prevLevelNode--;
        level->push_back(curNode);

        if (curNode == _sink) break;

        for (auto edge : curNode->_fanouts) {
            numDrivers[edge->_fanout->_id]--;
            if (numDrivers[edge->_fanout->_id] == 0) {
                q.push(edge->_fanout);
                curLevelNode++;
            }
        }
    }
}

void TimingGraph::report() const {
    log() << "---- report timing graph ----" << endl;
    printlog(LOG_INFO, "#level/#revLevel=%lu/%lu", _levels.size(), _revLevels.size());
    int cnt = 0;
    for (auto edge : _edges) {
        if (edge->_net && edge->_net->isInterNet()) cnt++;
    }
    printlog(LOG_INFO, "#nodes=%d, #edges=%d, #xdr-edges=%d", getNumNodes(), getNumEdges(), cnt);

    vector<double> intraNetDelay;
    for (auto edge : _edges)
        if (edge->_net && edge->isConstEdge()) intraNetDelay.push_back(edge->_constDelay);
    double sum = 0;
    for (auto val : intraNetDelay) sum += val;
    sort(intraNetDelay.begin(), intraNetDelay.end());
    printlog(LOG_INFO,
             "intra-net delay: min=%f, max=%f, avg=%f, 0.25/0.5/0.75/0.9/0.95/0.98=%.1f/%.1f/%.1f/%.1f/%.1f/%.1f",
             intraNetDelay.front(),
             intraNetDelay.back(),
             sum / intraNetDelay.size(),
             intraNetDelay[intraNetDelay.size() * 0.25],
             intraNetDelay[intraNetDelay.size() * 0.5],
             intraNetDelay[intraNetDelay.size() * 0.75],
             intraNetDelay[intraNetDelay.size() * 0.9],
             intraNetDelay[intraNetDelay.size() * 0.95],
             intraNetDelay[intraNetDelay.size() * 0.98]);

    int maxFanin = -1, maxFanout = -1;
    for (auto node : _nodes) {
        if (node != _source && node != _sink) {
            maxFanin = max((int)node->_drivers.size(), maxFanin);
            maxFanout = max((int)node->_fanouts.size(), maxFanout);
        }
    }
    printlog(LOG_INFO,
             "#src=%lu, #sink=%lu, maxFanin=%d, maxFanout=%d",
             _source->_fanouts.size(),
             _sink->_drivers.size(),
             maxFanin,
             maxFanout);

    vector<int> levelConstEdge(_levels.size(), 0), levelXdrEdge(_levels.size(), 0), levelOptXdrEdge(_levels.size(), 0);
    for (unsigned l = 0; l < _levels.size(); l++) {
        for (auto node : _levels[l]) {
            for (auto edge : node->_fanouts) {
                if (edge->isConstEdge()) {
                    levelConstEdge[l]++;
                } else {
                    levelXdrEdge[l]++;
                    if (tdmDatabase.isOptVar(edge->_net->_xdrVar)) levelOptXdrEdge[l]++;
                }
            }
        }
        printlog(LOG_INFO,
                 "level %lu: #node=%lu, #const_edge=%d, #xdr_edge/opt_xdr_edge=%d/%d",
                 l,
                 _levels[l].size(),
                 levelConstEdge[l],
                 levelXdrEdge[l],
                 levelOptXdrEdge[l]);
    }
}

bool TimingGraph::isCyclicUtil(int v, vector<bool>& visited, vector<bool>& recStack) {
    if (!visited[v]) {
        // Mark the current node as visited and part of recursion stack
        visited[v] = true;
        recStack[v] = true;

        // Recur for all the vertices adjacent to this vertex
        list<int>::iterator i;
        for (auto edge : _nodes[v]->_fanouts) {
            int i = edge->_fanout->_id;
            if (!visited[i] && isCyclicUtil(i, visited, recStack)) {
                return true;
            } else if (recStack[i]) {
                return true;
            }
        }
    }
    recStack[v] = false;  // remove the vertex from recursion stack
    return false;
}

// Returns true if the graph contains a cycle, else false.
bool TimingGraph::isCyclic() {
    // Mark all the vertices as not visited and not part of recursion
    // stack
    vector<bool> visited(_nodes.size(), false);
    vector<bool> recStack(_nodes.size(), false);

    // Call the recursive helper function to detect cycle in different
    // DFS trees
    for (unsigned i = 0; i < _nodes.size(); i++)
        if (isCyclicUtil(i, visited, recStack)) return true;

    return false;
}

bool TimingGraph::DFSUtil(int v, int dest, vector<bool>& visited) {
    visited[v] = true;

    if (dest == v) return true;

    for (auto edge : _nodes[v]->_fanouts) {
        int i = edge->_fanout->_id;
        if (!visited[i]) {
            if (DFSUtil(i, dest, visited)) {
                return true;
            }
        }
    }

    return false;
}

void TimingGraph::DFS(int v, int dest) {
    vector<bool> visited(_nodes.size(), false);

    DFSUtil(v, dest, visited);
}

bool Edge::isConstEdge() const { return !_net || _net->isIntraNet(); }

bool Edge::isCritical() const { return _fanout->getArrivalTime() == _driver->getArrivalTime() + getDelay(); }

bool TimingGraph::isOptEdge(const Edge* edge) const {
    return edge->_net && !edge->isConstEdge() && tdmDatabase.isOptVar(edge->_net->_xdrVar);
}

void Node::updateArrivalTime(double value) { _arrivalTime = max(_arrivalTime, value); }

vector<Edge*> TimingGraph::getCriticalPath() const {
    vector<Edge*> path;

    Node* sink = getSink();

    while (sink->_drivers.size() > 0) {
        for (auto fanin : sink->_drivers) {
            if (fanin->isCritical()) {
                sink = fanin->_driver;
                path.push_back(fanin);
                break;
            }
        }
    }

    return path;
}
