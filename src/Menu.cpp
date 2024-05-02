#include "Menu.h"
#include <iomanip>
#include <iostream>
using namespace std;

Menu::Menu() {
    data_ = Data();
}

void Menu::drawTop(){
    cout << "┌──────────────────────────────────────────────────┐" << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << setw(50) << "│         ██████   █████  ██████  ██████           │" << endl;
    cout << setw(50) << "│         ██   ██ ██   ██ ██   ██      ██          │" << endl;
    cout << setw(50) << "│         ██   ██ ███████ ██████   █████           │" << endl;
    cout << setw(50) << "│         ██   ██ ██   ██ ██      ██               │" << endl;
    cout << setw(50) << "│         ██████  ██   ██ ██      ███████          │" << endl;
    cout << "│" << setw(53) << "│" << endl;
}

void Menu::drawBottom(){
    cout << "└──────────────────────────────────────────────────┘" << endl;
}

void Menu::drawMenu() {
    char key;
    bool flag = true;
    while (flag) {
        drawTop();
        cout << "│" << setw(53) << "│" << endl;
        cout << "│    Options:                                      │" << endl;
        cout << "│     [1] Backtracking Algorithm                   │" << endl;
        cout << "│     [2] Triangular Approximation Heuristic       │" << endl;
        cout << "│     [3] Cluster Approximation Heuristic          │" << endl;
        cout << "│     [4] MST Approximation Heuristic              │" << endl;
        cout << "│     [5] Approximation Heuristic Analysis         │" << endl;
        cout << "│     [Q] Exit                                     │" << endl;
        cout << "│" << setw(53) << "│" << endl;
        drawBottom();
        cout << "Choose an option: ";
        cin >> key;
        switch (key) {
            case '1':{
                string vertex_id;
                cout << "Enter the vertex id: ";
                cin >> vertex_id;
                drawBacktracking(vertex_id);
                break;
            }
            case '2': {
                string vertex_id;
                cout << "Enter the vertex id: ";
                cin >> vertex_id;
                drawTriangular(vertex_id);
                break;
            }
            case '3': {
                string vertex_id;
                cout << "Enter the vertex id: ";
                cin >> vertex_id;
                drawCluster(vertex_id);
                break;
            }
            case '4': {
                string vertex_id;
                cout << "Enter the vertex id: ";
                cin >> vertex_id;
                drawMST(vertex_id);
                break;
            }
            case '5': {
                string vertex_id;
                cout << "Enter the vertex id: ";
                cin >> vertex_id;
                drawApproximationAnalysis(vertex_id);
                break;
            }
            case 'Q':
            case 'q': {
                cout << "Exiting..." << endl;
                flag = false;
                break;
            }
            default:
                cout << "Invalid option" << endl;
                break;
        }
    }
}

void Menu::drawBacktracking(string vertex_id) {
    data_.readNodes("../dataset/Real-world Graphs/Real-world Graphs/graph1/nodes.csv");
    data_.readEdges(true,"../dataset/Real-world Graphs/Real-world Graphs/graph1/edges.csv");

    auto start = chrono::high_resolution_clock::now();
    data_.backtrackingTSP();
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> duration = end - start;

    cout << "┌─ Backtracking Algorithm  ────────────────────────┐" << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getCost() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Tour size: " << right << left << setw(36) << data_.getNetwork().getVertexSet().size() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration.count()) +  " seconds" << "│" << right << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << "└──────────────────────────────────────────────────┘" << endl;
}

void Menu::drawTriangular(string vertex_id) {
    data_.readNodes("../dataset/Real-world Graphs/Real-world Graphs/graph1/nodes.csv");
    data_.readEdges(true,"../dataset/Real-world Graphs/Real-world Graphs/graph1/edges.csv");

    auto start = chrono::high_resolution_clock::now();
    data_.triangularHeuristicAproximation(vertex_id);
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> duration = end - start;

    cout << "┌─ Triangular Approximation Heuristic ─────────────┐" << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getAproximationTourCost() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Tour size: " << right << left << setw(36) << data_.getAproximationTour().size() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration.count()) +  " seconds" << "│" << right << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << "└──────────────────────────────────────────────────┘" << endl;
}

void Menu::drawCluster(string vertex_id) {
    data_.readNodes("../dataset/Real-world Graphs/Real-world Graphs/graph1/nodes.csv");
    data_.readEdges(true,"../dataset/Real-world Graphs/Real-world Graphs/graph1/edges.csv");

    auto start = chrono::high_resolution_clock::now();
    data_.clusterApproximationTSP(vertex_id);
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> duration = end - start;

    cout << "┌─ Cluster Approximation Heuristic ────────────────┐" << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getClusterTourCost() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Tour size: " << right << left << setw(36) << data_.getClusterTour().size() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration.count()) +  " seconds" << "│" << right << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << "└──────────────────────────────────────────────────┘" << endl;
}

void Menu::drawMST(std::string vertex_id) {
    data_.readNodes("../dataset/Real-world Graphs/Real-world Graphs/graph1/nodes.csv");
    data_.readEdges(true,"../dataset/Real-world Graphs/Real-world Graphs/graph1/edges.csv");

    auto start = chrono::high_resolution_clock::now();
    data_.mstApproximationTSP(vertex_id);
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> duration = end - start;

    cout << "┌─ MST Approximation Heuristic ────────────────────┐" << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getMSTTourCost() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Tour size: " << right << left << setw(36) << data_.getMSTTour().size() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration.count()) +  " seconds" << "│" << right << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << "└──────────────────────────────────────────────────┘" << endl;
}

void Menu::drawApproximationAnalysis(std::string vertex_id) {
    /*data_.readNodes("../dataset/Real-world Graphs/Real-world Graphs/graph1/nodes.csv");
    data_.readEdges(true,"../dataset/Real-world Graphs/Real-world Graphs/graph1/edges.csv");

    auto start = chrono::high_resolution_clock::now();
    data_.backtrackingTSP();
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> duration = end - start;
*/
    cout << "┌─ Approximation Heuristic Analysis ───────────────┐" << endl;
    cout << "│" << setw(53) << "│" << endl;
/*
    cout << "│  Backtracking Algorithm :                        │" << endl;
    cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getCost() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Tour size: " << right << left << setw(36) << data_.getNetwork().getVertexSet().size() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration.count()) +  " seconds" << "│" << right << endl;
    cout << "│" << setw(53) << "│" << endl;
*/

    data_.readNodes("../dataset/Real-world Graphs/Real-world Graphs/graph1/nodes.csv");
    data_.readEdges(true,"../dataset/Real-world Graphs/Real-world Graphs/graph1/edges.csv");

    auto start = chrono::high_resolution_clock::now();
    data_.triangularHeuristicAproximation(vertex_id);
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> duration = end - start;

    cout << "│  Triangular Approximation Heuristic :            │" << endl;
    cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getAproximationTourCost() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Tour size: " << right << left << setw(36) << data_.getAproximationTour().size() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration.count()) +  " seconds" << "│" << right << endl;
    cout << "│" << setw(53) << "│" << endl;

    data_.readNodes("../dataset/Real-world Graphs/Real-world Graphs/graph1/nodes.csv");
    data_.readEdges(true,"../dataset/Real-world Graphs/Real-world Graphs/graph1/edges.csv");

    auto start2 = chrono::high_resolution_clock::now();
    data_.clusterApproximationTSP(vertex_id);
    auto end2 = chrono::high_resolution_clock::now();

    chrono::duration<double> duration2 = end2 - start2;

    cout << "│  Cluster Approximation Heuristic :               │" << endl;
    cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getClusterTourCost() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Tour size: " << right << left << setw(36) << data_.getClusterTour().size() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration2.count()) +  " seconds" << "│" << right << endl;
    cout << "│" << setw(53) << "│" << endl;

    data_.readNodes("../dataset/Real-world Graphs/Real-world Graphs/graph1/nodes.csv");
    data_.readEdges(true,"../dataset/Real-world Graphs/Real-world Graphs/graph1/edges.csv");

    auto start3 = chrono::high_resolution_clock::now();
    data_.mstApproximationTSP(vertex_id);
    auto end3 = chrono::high_resolution_clock::now();

    chrono::duration<double> duration3 = end3 - start3;

    cout << "│  MST Approximation Heuristic :                   │" << endl;
    cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getMSTTourCost() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Tour size: " << right << left << setw(36) << data_.getMSTTour().size() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration3.count()) +  " seconds" << "│" << right << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << "└──────────────────────────────────────────────────┘" << endl;
}