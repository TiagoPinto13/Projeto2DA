//
// Created by antero on 25-04-2024.
//

#ifndef PROJ2DA_DATA_H
#define PROJ2DA_DATA_H


#include "Graph.h"
#include <string>
#include <map>
#include <unordered_set>

class Data {
public:

    Data() = default;
    void readNodes(std::string nodeFilePath);
    void readEdges(bool realWorldGraphs, std::string edgesFilePath);
    void parseTOY(bool tourismCSV, std::string edgesFilePath);
    Graph getNetwork();
    double getCost();
    std::map<std::string,std::string> getTourismLabels();

    std::vector<Vertex*> backtrackingTSP();
    double calculateTourCost(const std::vector<Vertex*>& tour) const;
    void backtrack(std::vector<Vertex*>& currentTour, double currentCost);
    double haversineDistance(double lat1, double lon1, double lat2, double lon2);
    Vertex* findNearestNeighbor(Vertex* v);
    void triangularHeuristicAproximation(const std::string& startNodeId);
    std::vector<Vertex*> getAproximationTour();
    double getAproximationTourCost();
    void clusterApproximationTSP(const std::string& startNodeId);
    Vertex* findNearestNeighborInCluster(Vertex* v, const std::unordered_set<Vertex*>& cluster);
    std::vector<Vertex*> getClusterTour();
    double getClusterTourCost();

private:
    std::vector<Vertex*> bestTour;
    double bestCost;

    std::vector<Vertex*> aproximation_tour_;
    double aproximation_tourCost_;

    std::vector<Vertex*> cluster_tour_;
    double cluster_tourCost_;

    bool tourism=false;
    Graph network_;
    std::map<std::string,std::string> tourismLabels;
};

#endif //PROJ2DA_DATA_H
