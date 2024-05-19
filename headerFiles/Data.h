#ifndef PROJ2DA_DATA_H
#define PROJ2DA_DATA_H

#include "Graph.h"
#include <string>
#include <map>
#include <unordered_set>

class Data {
public:

    Data() = default;
    void readNodes(std::string nodeFilePath, int numberOfNodes);
    void readEdges(bool realWorldGraphs, std::string edgesFilePath);
    void parseTOY(bool tourismCSV, std::string edgesFilePath);
    Graph getNetwork();
    double getCost();
    bool isTourism();
    std::vector<Vertex*> getBestTour();
    std::map<std::string,std::string> getTourismLabels();

    void backtrackingTSP();
    void backtrack(std::vector<Vertex*>& currentTour, double currentCost);

    double calculateTourCost(const std::vector<Vertex*>& tour) const;
    double haversineDistance(double lat1, double lon1, double lat2, double lon2);
    void resetNodesVisitation();
    Vertex* findNearestNeighbor(Vertex* v);
    void triangularHeuristicAproximation(const std::string& startNodeId);
    void dfsMST(Vertex* v, const std::vector<Vertex*>& mst);
    std::vector<Vertex*> getAproximationTour();
    double getAproximationTourCost();
    void clusterApproximationTSP(const std::string& startNodeId);
    Vertex* findNearestNeighborCluster(Vertex* v, const std::unordered_set<Vertex*>& cluster);
    std::vector<Vertex*> getClusterTour();
    std::vector<Vertex *> prim(Graph * g);
    double getClusterTourCost();
    void preorderTraversalMST(Vertex* u);
    void mstApproximationTSP(const std::string& startNodeId);
    std::vector<Vertex*> getMSTTour();
    double getMSTTourCost();






    std::string bfs_farthest_node(const std::string& start);
    std::unordered_map<std::string, int> dijkstra(const std::string& start);
    std::vector<std::string> tsp_real_world1( std::string start);
    void removeVertex(std::string id);
    void removeEdge(std::string id1, std::string id2);

    std::vector<std::string> twoOpt(const std::vector<std::string>& tour);
    bool isConnected(const std::string& start);
    std::vector<std::string> tsp_real_world2( std::string start);


private:
        std::vector<Vertex*> bestTour;
        double bestCost;

        std::vector<Vertex*> aproximation_tour_;
        double aproximation_tourCost_;

        std::vector<Vertex*> cluster_tour_;
        double cluster_tourCost_;

        std::vector<Vertex*> mst_tour_;
        double mst_tourCost_;

        bool tourism=false;
        Graph network_;
        std::map<std::string,std::string> tourismLabels;

    };

#endif //PROJ2DA_DATA_H
