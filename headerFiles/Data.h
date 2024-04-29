//
// Created by antero on 25-04-2024.
//

#ifndef PROJ2DA_DATA_H
#define PROJ2DA_DATA_H


#include "Graph.h"
#include <string>
#include <map>

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

        private:
    std::vector<Vertex*> bestTour;
    double bestCost;

    bool tourism=false;
    Graph network_;
    std::map<std::string,std::string> tourismLabels;
};

#endif //PROJ2DA_DATA_H
