#include <iostream>
#include <chrono>
#include <iomanip>
#include "headerFiles/Data.h"
#include "headerFiles/Graph.h"
#include "headerFiles/Menu.h"
using namespace std;

int main(){
    Data d;
    d.parseTOY(true, "../dataset/Toy-Graphs/Toy-Graphs/tourism.csv");
    d.backtrackingTSP();
    std::cout<<d.getCost();
    for(auto v: d.getBestTour())
    return 0;
}
