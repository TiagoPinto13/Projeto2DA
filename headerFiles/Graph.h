#ifndef PROJETO1DA_GRAPH_H
#define PROJETO1DA_GRAPH_H

#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include "MutablePriorityQueue.h"
#include "DeliverySites.h"
#include "PumpingStations.h"
#include "WaterReservoir.h"


class Edge;

class Vertex {
public:
    Vertex(const std::string& in, const int type);
    Edge* addEdge(Vertex* dest, double w);
    bool removeEdge(const std::string& in);
    void removeOutgoingEdges();

    bool operator<(const Vertex& vertex) const;

    std::string getInfo() const;
    int getType() const; //water reservoirs-0; pumping stations-1 ; delivery sites-2
    std::vector<Edge*> getAdj() const;
    bool isVisited() const;
    bool isProcessing() const;
    unsigned int getIndegree() const;
    double getDist() const;
    Edge* getPath() const;
    std::vector<Edge*> getIncoming() const;

    void setInfo(const std::string& info);
    void setVisited(bool visited);
    void setProcessing(bool processing);
    void setIndegree(unsigned int indegree);
    void setDist(double dist);
    void setPath(Edge* path);

    friend class MutablePriorityQueue<Vertex>;

protected:
    std::string info;
    std::vector<Edge*> adj;
    int type_;
    bool visited;
    bool processing;
    unsigned int indegree;
    double dist;
    Edge* path;
    std::vector<Edge*> incoming;

    int queueIndex;
    void deleteEdge(Edge* edge);
};

class Edge {
public:
    Edge(Vertex* orig, Vertex* dest, double w);

    Vertex* getDest() const;
    double getWeight() const;
    bool isSelected() const;
    Vertex* getOrig() const;
    Edge* getReverse() const;
    double getFlow() const;

    void setSelected(bool selected);
    void setReverse(Edge* reverse);
    void setFlow(double flow);

protected:
    Vertex* dest;
    double weight;
    bool selected;
    Vertex* orig;
    Edge* reverse;
    double flow;
};

class Graph {
public:
    Vertex* findVertex(const std::string& in) const;
    bool addVertex(const std::string& in, const int type);
    bool removeVertex(const std::string& in);
    bool addEdge(const std::string& source, const std::string& dest, double w);
    bool removeEdge(const std::string& source, const std::string& dest);
    bool addBidirectionalEdge(const std::string& source, const std::string& dest, double w);

    int getNumVertex() const;
    std::vector<Vertex*> getVertexSet() const;
    std::unordered_map<std::string, Vertex*> getVertexMap() const;

    std::vector<std::string> dfs() const;
    std::vector<std::string> dfs(const std::string& source) const;
    void dfsVisit(Vertex* v, std::vector<std::string>& res) const;
    std::vector<std::string> bfs(const std::string& source) const;

    bool isDAG() const;
    bool dfsIsDAG(Vertex* v) const;
    std::vector<std::string> topsort() const;

protected:
    std::vector<Vertex*> vertexSet;
    std::unordered_map<std::string, Vertex*> vertexMap;

    double** distMatrix;
    int** pathMatrix;

    int findVertexIdx(const std::string& in) const;
    void deleteMatrix(int** m, int n);
    void deleteMatrix(double** m, int n);
};
#endif