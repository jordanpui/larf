#pragma once

#include "global.h"

class Node;
class Edge;
class TimingGraph;
class TdmNet;
namespace db {
class Instance;
class Group;
}  // namespace db
class XdrVar;
class TdmDB;

class Node {
public:
    void updateArrivalTime(double value);
    void updateRequireTime(double value) { _requireTime = min(_requireTime, value); }
    void resetArrivalTime() { _arrivalTime = -1; }
    void resetRequireTime() { _requireTime = DBL_MAX; }
    double getArrivalTime() const { return _arrivalTime; }
    double getRequireTime() const { return _requireTime; }
    double getSlack() const { return _requireTime - _arrivalTime; }

    void addDriver(Edge* driver) { _drivers.push_back(driver); }
    void addFanout(Edge* fanout) { _fanouts.push_back(fanout); }

    Node(db::Instance* instance, int id);

    vector<Edge*> _drivers;
    vector<Edge*> _fanouts;
    int _id;

    db::Instance* _instance;
    double _delay;

private:
    double _arrivalTime;
    double _requireTime;
};

class Edge {
public:
    Edge(Node* fromNode, Node* toNode, TdmNet* net, int id);
    void updateDelay();
    double getDelay() const { return _delay; }
    double getDelay(int val) const;
    double getArrivalTimeAlongEdge() const { return _arrivalTimeAlongEdge; }
    bool isConstEdge() const;
    bool isCritical() const;

    Node* _driver;
    Node* _fanout;
    int _id;

    TdmNet* _net;
    double _constDelay;
    const double _tdmCoef = 5;

private:
    double _delay;
    double _arrivalTimeAlongEdge;
};

class TimingGraph {
public:
    void updateArrivalTime();
    void updateRequireTime();
    void resetTiming();
    void resetArrivalTime();
    void resetRequireTime();

    double getSinkAT() const { return _sink->getArrivalTime(); }

    Edge* addEdge(int u, int v, TdmNet* net);
    Edge* addEdge(Node* u, Node* v, TdmNet* net);
    Node* addNode(db::Instance* instance);

    void breakCycle();
    void setConstDelay();
    void setSrcSink();
    void removeAbnEdges();

    void levelize();
    vector<vector<Node*>>& getRevLevels() { return _revLevels; }
    vector<vector<Node*>>& getLevels() { return _levels; }

    int getNumNodes() const { return _nodes.size(); }
    int getNumEdges() const { return _edges.size(); }

    Node* getNode(int i) { return _nodes[i]; }
    Edge* getEdge(int i) { return _edges[i]; }
    vector<Edge*>& getEdges(XdrVar* var);
    vector<Edge*>& getEdges(TdmNet* net);
    void getSRCoef(XdrVar* var, double& k, double& b);

    bool isOptEdge(const Edge* edge) const;

    void addMapping(XdrVar* var, Edge* edge);

    void report() const;

    bool isCyclic();
    void DFS(int v, int dest);

    Node* getSink() const { return _sink; }
    Node* getSource() const { return _source; }

    vector<Edge*> getCriticalPath() const;

private:
    Node* _source;
    Node* _sink;

    map<XdrVar*, vector<Edge*>> _xdrToEdges;

    vector<Node*> _nodes;
    vector<Edge*> _edges;

    vector<vector<Node*>> _levels;     // from source to sink
    vector<vector<Node*>> _revLevels;  // frome sink to source

    const double outlierRatio = 0.02;
    const double wireDelayCoef = 1;

    void forLevelize();
    void revLevelize();

    bool isCyclicUtil(int v, vector<bool>& visited, vector<bool>& recStack);
    bool DFSUtil(int v, int dest, vector<bool>& visited);

    void forwardPropagateST();
    void forwardPropagateMT();
    void backwardPropagateST();
    void backwardPropagateMT();
};