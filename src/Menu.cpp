#include "Menu.h"
#include <iomanip>
#include <iostream>
#include <chrono>

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

void Menu::firstMenu(){
    drawTop();
    cout << "│     Choose the dataset:                          │" << endl;
    cout << "│          [1] TOY Graphs                          │" << endl;
    cout << "│          [2] Real World Graphs                   │" << endl;
    cout << "│          [3] Fully Connected Graphs              │" << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << "└──────────────────────────────────────────────────┘" << endl;
    cout << "Enter your choice: ";
    int choice;
    cin >> choice;
    if(choice == 1){
        drawTop();
        cout << "│     Choose the file:                             │" << endl;
        cout << "│          [1] Shipping                            │" << endl;
        cout << "│          [2] Stadiums                            │" << endl;
        cout << "│          [3] Tourism                             │" << endl;
        cout << "│" << setw(53) << "│" << endl;
        cout << "└──────────────────────────────────────────────────┘" << endl;
        cout << "Enter your choice: ";
        int choice2;
        cin >> choice2;
        if(choice2 == 1){
            drawMenu(1);
        }
        else if(choice2 == 2){
            drawMenu(2);
        }
        else if(choice2 == 3){
            drawMenu(3);
        }
        else{
            cout << "Invalid choice!" << endl;
        }
    }
    else if(choice == 2){
        drawTop();
        cout << "│     Choose the Graph:                            │" << endl;
        cout << "│          [1] Graph 1                             │" << endl;
        cout << "│          [2] Graph 2                             │" << endl;
        cout << "│          [3] Graph 3                             │" << endl;
        cout << "│" << setw(53) << "│" << endl;
        cout << "└──────────────────────────────────────────────────┘" << endl;
        cout << "Enter your choice: ";
        int choice2;
        cin >> choice2;
        if(choice2 == 1){
            drawMenu(4);
        }
        else if(choice2 == 2){
            drawMenu(5);
        }
        else if(choice2 == 3){
            drawMenu(6);
        }
    }
    else if(choice == 3){
        drawTop();
        cout << "│     Choose the Number of Edges:                  │" << endl;
        cout << "│          [1] 25                                  │" << endl;
        cout << "│          [2] 50                                  │" << endl;
        cout << "│          [3] 75                                  │" << endl;
        cout << "│          [4] 100                                 │" << endl;
        cout << "│          [5] 200                                 │" << endl;
        cout << "│          [6] 300                                 │" << endl;
        cout << "│          [7] 400                                 │" << endl;
        cout << "│          [8] 500                                 │" << endl;
        cout << "│          [9] 600                                 │" << endl;
        cout << "│          [10] 700                                │" << endl;
        cout << "│          [11] 800                                │" << endl;
        cout << "│          [12] 900                                │" << endl;
        cout << "│" << setw(53) << "│" << endl;
        cout << "└──────────────────────────────────────────────────┘" << endl;
        cout << "Enter your choice: ";
        int choice2;
        cin >> choice2;
        if(choice2 == 1){
            drawMenu(7);
        }
        else if(choice2 == 2){
            drawMenu(8);
        }
        else if(choice2 == 3){
            drawMenu(9);
        }
        else if(choice2 == 4){
            drawMenu(10);
        }
        else if(choice2 == 5){
            drawMenu(11);
        }
        else if(choice2 == 6){
            drawMenu(12);
        }
        else if(choice2 == 7){
            drawMenu(13);
        }
        else if(choice2 == 8){
            drawMenu(14);
        }
        else if(choice2 == 9){
            drawMenu(15);
        }
        else if(choice2 == 10){
            drawMenu(16);
        }
        else if(choice2 == 11){
            drawMenu(17);
        }
        else if(choice2 == 12) {
            drawMenu(18);
        }
    }
    else{
        cout << "Invalid choice!" << endl;
    }
}

void Menu::drawBottom(){
    cout << "└──────────────────────────────────────────────────┘" << endl;
}

void Menu::drawMenu(int option) {
    if(option ==1){
        data_.parseTOY(false,"../dataset/Toy-Graphs/Toy-Graphs/shipping.csv");
    }
    else if(option == 2){
        data_.parseTOY(false,"../dataset/Toy-Graphs/Toy-Graphs/stadiums.csv");
    }
    else if(option == 3){
        data_.parseTOY(true,"../dataset/Toy-Graphs/Toy-Graphs/tourism.csv");
    }
    else if(option == 4){
        data_.readNodes("../dataset/Real-world Graphs/Real-world Graphs/graph1/nodes.csv", -1);
        data_.readEdges(true,"../dataset/Real-world Graphs/Real-world Graphs/graph1/edges.csv");
    }
    else if(option == 5){
        data_.readNodes("../dataset/Real-world Graphs/Real-world Graphs/graph2/nodes.csv", -1);
        data_.readEdges(true,"../dataset/Real-world Graphs/Real-world Graphs/graph2/edges.csv");
    }
    else if(option == 6){
        data_.readNodes("../dataset/Real-world Graphs/Real-world Graphs/graph3/nodes.csv", -1);
        data_.readEdges(true,"../dataset/Real-world Graphs/Real-world Graphs/graph3/edges.csv");
    }
    else if(option == 7){
        data_.readNodes("../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/nodes.csv", 25);
        data_.readEdges(false,"../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/edges_25.csv");
    }
    else if(option == 8){
        data_.readNodes("../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/nodes.csv", 50);
        data_.readEdges(false,"../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/edges_50.csv");
    }
    else if(option == 9){
        data_.readNodes("../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/nodes.csv", 75);
        data_.readEdges(false,"../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/edges_75.csv");
    }
    else if(option == 10){
        data_.readNodes("../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/nodes.csv", 100);
        data_.readEdges(false,"../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/edges_100.csv");
    }
    else if(option == 11){
        data_.readNodes("../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/nodes.csv", 200);
        data_.readEdges(false,"../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/edges_200.csv");
    }
    else if(option == 12){
        data_.readNodes("../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/nodes.csv", 300);
        data_.readEdges(false,"../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/edges_300.csv");
    }
    else if(option == 13){
        data_.readNodes("../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/nodes.csv", 400);
        data_.readEdges(false,"../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/edges_400.csv");
    }
    else if(option == 14){
        data_.readNodes("../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/nodes.csv", 500);
        data_.readEdges(false,"../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/edges_500.csv");
    }
    else if(option == 15){
        data_.readNodes("../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/nodes.csv", 600);
        data_.readEdges(false,"../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/edges_600.csv");
    }
    else if(option == 16){
        data_.readNodes("../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/nodes.csv", 700);
        data_.readEdges(false,"../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/edges_700.csv");
    }
    else if(option == 17){
        data_.readNodes("../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/nodes.csv", 800);
        data_.readEdges(false,"../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/edges_800.csv");
    }
    else if(option == 18){
        data_.readNodes("../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/nodes.csv", 900);
        data_.readEdges(false,"../dataset/Extra_Fully_Connected_Graphs/Extra_Fully_Connected_Graphs/edges_900.csv");
    }
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
        cout << "│     [6] TSP in Real World                        │" << endl;
        cout << "│     [7] Change current dataset                   │" << endl;
        cout << "│     [Q] Exit                                     │" << endl;
        cout << "│" << setw(53) << "│" << endl;
        drawBottom();
        cout << "Choose an option: ";
        cin >> key;
        switch (key) {
            case '1':{
                /*string vertex_id;
                cout << "Enter the vertex id: ";
                cin >> vertex_id;*/
                drawBacktracking(/*vertex_id*/);  //diz no enunciado que é sempre com o vertex 0
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
            case '6':{
                string vertex_id;
                cout << "Enter the vertex id: ";
                cin >> vertex_id;
                drawTspRealWorld(vertex_id);
                break;
            }
            case '7': {
                flag = false;
                firstMenu();
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

void Menu::drawBacktracking(/*string vertex_id*/) {

    auto start = chrono::high_resolution_clock::now();
    data_.backtrackingTSP();
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> duration = end - start;

    if(data_.isTourism()) {
        cout << "┌─ Backtracking Algorithm  ────────────────────────┐" << endl;
        cout << "│" << setw(53) << "│" << endl;
        cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getCost() << right << " │" << endl;
        cout << "│ " << left << setw(12) << "Tour size: " << right << left << setw(36) << data_.getBestTour().size() << right << " │" << endl;
        cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration.count()) +  " seconds" << "│" << right << endl;
        cout << "│" << setw(53) << "│" << endl;
        for(auto v: data_.getBestTour()) {
            cout << "│ " << left << "vertex: " << v->getInfo() << " - " << setw(10) << data_.getTourismLabels()[v->getInfo()] << right <<setw(30) << "│"  << endl;

        }
        cout << "│" << setw(53) << right <<"│" << endl;
        cout << "└──────────────────────────────────────────────────┘" << endl;
        waitForEnter();
    } else {

    }
    cout << "┌─ Backtracking Algorithm  ────────────────────────┐" << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getCost() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Tour size: " << right << left << setw(36) << data_.getBestTour().size() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration.count()) +  " seconds" << "│" << right << endl;
    cout << "│" << setw(53) << "│" << endl;
    for(auto v: data_.getBestTour()) {
        cout << "│ " << left << "vertex: " << setw(4) << v->getInfo() << right <<setw(40) << "│"  << endl;

    }
    cout << "│" << setw(53) << right <<"│" << endl;
    cout << "└──────────────────────────────────────────────────┘" << endl;
    waitForEnter();
}

void Menu::drawTriangular(string vertex_id) {
    auto start = chrono::high_resolution_clock::now();
    data_.triangularHeuristicAproximation(vertex_id);
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> duration = end - start;

    cout << "┌─ Triangular Approximation Heuristic ─────────────┐" << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getAproximationTourCost() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration.count()) +  " seconds" << "│" << right << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << "└──────────────────────────────────────────────────┘" << endl;
    waitForEnter();
}

void Menu::drawCluster(string vertex_id) {
    auto start = chrono::high_resolution_clock::now();
    data_.clusterApproximationTSP(vertex_id);
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> duration = end - start;

    cout << "┌─ Cluster Approximation Heuristic ────────────────┐" << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getClusterTourCost() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration.count()) +  " seconds" << "│" << right << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << "└──────────────────────────────────────────────────┘" << endl;
    waitForEnter();
}

void Menu::drawMST(std::string vertex_id) {
    auto start = chrono::high_resolution_clock::now();
    data_.mstApproximationTSP(vertex_id);
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> duration = end - start;

    cout << "┌─ MST Approximation Heuristic ────────────────────┐" << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getMSTTourCost() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration.count()) +  " seconds" << "│" << right << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << "└──────────────────────────────────────────────────┘" << endl;
    waitForEnter();
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
    cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration.count()) +  " seconds" << "│" << right << endl;
    cout << "│" << setw(53) << "│" << endl;
*/

    auto start = chrono::high_resolution_clock::now();
    data_.triangularHeuristicAproximation(vertex_id);
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> duration = end - start;

    cout << "│  Triangular Approximation Heuristic :            │" << endl;
    cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getAproximationTourCost() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration.count()) +  " seconds" << "│" << right << endl;
    cout << "│" << setw(53) << "│" << endl;

    auto start2 = chrono::high_resolution_clock::now();
    data_.clusterApproximationTSP(vertex_id);
    auto end2 = chrono::high_resolution_clock::now();

    chrono::duration<double> duration2 = end2 - start2;

    cout << "│  Cluster Approximation Heuristic :               │" << endl;
    cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getClusterTourCost() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration2.count()) +  " seconds" << "│" << right << endl;
    cout << "│" << setw(53) << "│" << endl;

    auto start3 = chrono::high_resolution_clock::now();
    data_.mstApproximationTSP(vertex_id);
    auto end3 = chrono::high_resolution_clock::now();

    chrono::duration<double> duration3 = end3 - start3;

    cout << "│  MST Approximation Heuristic :                   │" << endl;
    cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getMSTTourCost() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration3.count()) +  " seconds" << "│" << right << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << "└──────────────────────────────────────────────────┘" << endl;
    waitForEnter();
}

double Menu::calculate_tour_cost(const std::vector<std::string>& tour) {
    double total_cost = 0.0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        double edge_weight = data_.getNetwork().getEdgeWeight(tour[i], tour[i + 1]);
        total_cost += edge_weight;
    }
    return total_cost;
}

void Menu::drawTspRealWorld(std::string vertex_id) {
    auto start = chrono::high_resolution_clock::now();
    std::vector<std::string> tour = data_.tsp_real_world(vertex_id);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;

    cout << "┌─ TSP in Real World ──────────────────────────────┐" << endl;
    cout << "│                                                  │" << endl;
    cout << "│ Start Node:                                      │" << endl;

    if (!tour.empty()) {
        cout << "│ " << left << setw(12) << "Time taken:" << right << left << setw(37) << to_string(duration.count()) +  " seconds" << "│" << endl;

        cout << "│ " << left << setw(12) << "Tour Cost:" << left << setw(39) << calculate_tour_cost(tour) << "│" << endl;
        cout << "│ " << left << setw(12) << "Tour:" << left << setw(39) << "[";
        for (size_t i = 0; i < tour.size(); ++i) {
            cout << tour[i];
            if (i < tour.size() - 1) {
                cout << ", ";
            }
        }
        cout << "]" << "│" << endl;
    } else {
        cout << "│ No feasible tour exists.                         │" << endl;
    }
    cout << "│                                                  │" << endl;
    cout << "└──────────────────────────────────────────────────┘" << endl;
    waitForEnter();
}

void Menu::waitForEnter() {
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    cout << "Press ENTER to continue...";
    getchar();
}