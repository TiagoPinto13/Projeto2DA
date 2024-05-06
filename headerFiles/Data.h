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
    Vertex* findNearestNeighbor(Vertex* v);
    void triangularHeuristicAproximation(const std::string& startNodeId);
    std::vector<Vertex*> getAproximationTour();
    double getAproximationTourCost();
    void clusterApproximationTSP(const std::string& startNodeId);
    Vertex* findNearestNeighborInCluster(Vertex* v, const std::unordered_set<Vertex*>& cluster);
    std::vector<Vertex*> getClusterTour();
    double getClusterTourCost();
    void preorderTraversalMST(Vertex* u);
    void mstApproximationTSP(const std::string& startNodeId);
    std::vector<Vertex*> getMSTTour();
    double getMSTTourCost();
    std::vector<std::string> tsp_subgraph(const Graph& graph, std::string start);
    std::vector<std::string> merge_tours(const std::vector<std::vector<std::string> >& tours);
    std::vector<std::string> tsp_real_world(const std::string& start_node);


    bool isInBestTour(Vertex* v);

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
