#include <iostream>
#include <chrono>
#include "headerFiles/Data.h"
#include "headerFiles/Graph.h"
#include "headerFiles/Menu.h"


int main() {
    /*Data data;
    Data data2;

    data.readNodes("../dataset/Real-world Graphs/Real-world Graphs/graph1/nodes.csv");
    data.readEdges(true,"../dataset/Real-world Graphs/Real-world Graphs/graph1/edges.csv");

    auto start1 = std::chrono::high_resolution_clock::now();
    data.triangularHeuristicAproximation("400");
    auto end1 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration1 = end1 - start1;

    std::cout << "Triangular Heuristic: " << data.getAproximationTourCost() << " " << data.getAproximationTour().size() << " Time taken: " << duration1.count() << " seconds" << std::endl;

    data2.readNodes("../dataset/Real-world Graphs/Real-world Graphs/graph1/nodes.csv");
    data2.readEdges(true,"../dataset/Real-world Graphs/Real-world Graphs/graph1/edges.csv");

    auto start2 = std::chrono::high_resolution_clock::now();
    data2.clusterApproximationTSP("400");
    auto end2 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration2 = end2 - start2;

    std::cout << "Cluster Approximation: " << data2.getClusterTourCost() << " " << data2.getClusterTour().size() << " Time taken: " << duration2.count() << " seconds" << std::endl;
    */
    Menu menu;
    menu.drawMenu();
    return 0;
}
