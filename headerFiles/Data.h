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
    void readEdges(bool realWorldGraphs, std::string edgesFilePath); //bool to skip the first line, since in the realWorldGraphs there's a 1st line to skip
    void parseTOY(bool tourismCSV, std::string edgesFilePath); //bool to store the names of locals only present in the tourismCSV
    Graph getNetwork();
    std::map<std::string,std::string> getTourismLabels();

private:
    bool tourism=false;
    Graph network_;
    std::map<std::string,std::string> tourismLabels;
};

#endif //PROJ2DA_DATA_H
