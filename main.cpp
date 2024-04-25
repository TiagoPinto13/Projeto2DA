#include <iostream>
#include "headerFiles/Data.h"
#include "headerFiles/Graph.h"



int main() {
    Data data;
    int k=0;
    data.parseTOY(true, "../dataset/Toy-Graphs/Toy-Graphs/tourism.csv");

    for(auto p:data.getTourismLabels()) {
        std::cout<<p.first<<" "<<p.second<<std::endl;
    }
    return 0;
}
