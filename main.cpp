#include <iostream>
#include "headerFiles/Data.h"
#include "headerFiles/Graph.h"



int main() {
    Data data;

    data.readNodes("../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/nodes.csv");
    data.readEdges(false,"../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/edges_25.csv");


    data.backtrackingTSP();
    std::cout<<data.getCost();

    return 0;
}
