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
        cout << "│     Choose the Number of Nodes:                  │" << endl;
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
    data_ = Data();
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
        cout << "│     [4] Approximation Heuristic Analysis         │" << endl;
        cout << "│     [5] TSP in Real World                        │" << endl;
        cout << "│     [6] Change current dataset                   │" << endl;
        cout << "│     [7] Remove Vertex/Edge                       │" << endl;
        cout << "│     [Q] Exit                                     │" << endl;
        cout << "│" << setw(53) << "│" << endl;
        drawBottom();
        cout << "Choose an option: ";
        cin >> key;
        switch (key) {
            case '1':{
                drawBacktracking();
                break;
            }
            case '2': {
                drawTriangular("0");
                break;
            }
            case '3': {
                drawCluster("0");
                break;
            }
            case '4': {

                drawApproximationAnalysis("0");
                break;
            }

            case '5': {
                drawTop();
                cout << "│" << setw(53) << "│" << endl;
                cout << "│    Options:                                      │" << endl;
                cout << "│     [1] Fast Method                              │" << endl;
                cout << "│     [2] Cost eficient Method                     │" << endl;
                drawBottom();
                char key10;
                cout << "Choose an option: ";
                cin >> key10;
                if (key10 == '1') {
                    string input;
                    cout << "Choose a start vertex: ";
                    cin >> input;
                    drawTspRealWorld(input);
                } else if (key10 == '2') {
                    string input;
                    cout << "Choose a start vertex: ";
                    cin >> input;
                    drawTspRealWorld2(input);
                } else {
                    cout << "Invalid option" << endl;
                }
            }
                
            case '6': {
                flag = false;
                firstMenu();
                break;
            }
            case '7': {
                drawRemoveVertexEdge();
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

/**
 * @brief Draw the menu for removing a vertex or an edge from the graph.
 *
 * This function displays a menu with options to remove either a vertex or an edge
 * from the graph. It waits for user input and performs the corresponding operation.
 *
 * The available options are:
 * - [1] Remove Vertex
 * - [2] Remove Edge
 * - [3] Back
 */

void Menu::drawRemoveVertexEdge(){
    char key;
    bool flag = true;
    while (flag) {
        drawTop();
        cout << "│" << setw(53) << "│" << endl;
        cout << "│    Options:                                      │" << endl;
        cout << "│     [1] Remove Vertex                            │" << endl;
        cout << "│     [2] Remove Edge                              │" << endl;
        cout << "│     [3] Back                                     │" << endl;
        cout << "│" << setw(53) << "│" << endl;
        drawBottom();
        cout << "Choose an option: ";
        cin >> key;
        switch (key) {
            case '1':{
                string vertex_id;
                cout << "Enter the vertex id: ";
                cin >> vertex_id;
                data_.removeVertex(vertex_id);
                break;
            }
            case '2': {
                string vertex_id1, vertex_id2;
                cout << "Enter the first vertex id: ";
                cin >> vertex_id1;
                cout << "Enter the second vertex id: ";
                cin >> vertex_id2;
                data_.removeEdge(vertex_id1, vertex_id2);
                break;
            }
            case '3': {
                flag = false;
                break;
            }
            default:
                cout << "Invalid option" << endl;
                break;
        }
    }

}


/**
 * @brief Draw the results of the backtracking algorithm for the Traveling Salesman Problem (TSP).
 *
 * This function executes the backtracking algorithm to find the solution for the TSP and displays
 * the results including the tour cost, tour size, time taken for computation, and the vertices
 * included in the best tour.
 *
 * If the problem involves tourism, it also displays additional information about each vertex,
 * such as tourism labels.
 */
void Menu::drawBacktracking() {

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


/**
 * @brief Draw the result of the Triangular Approximation Heuristic for the Traveling Salesman Problem (TSP).
 *
 * This function calculates the Triangular Approximation Heuristic solution for the TSP
 * starting from the specified vertex and draws the result along with the tour cost
 * and the time taken for the calculation.
 *
 * @param vertex_id The ID of the starting vertex for the TSP tour.
 */
void Menu::drawTriangular(string vertex_id) {
    auto start = chrono::high_resolution_clock::now();
    data_.triangularHeuristicAproximation(vertex_id);
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> duration = end - start;

    cout << "┌─ Triangular Approximation Heuristic ─────────────┐" << endl;
    cout << "│" << setw(53) << "│" << endl;
    if(data_.getAproximationTourCost()==-1) {
        cout << "│" << setw(27)<< "no tour"<< setw(23)<<" " << "│" << endl;
    }
    else {
        cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getAproximationTourCost() << right << " │" << endl;
        cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration.count()) +  " seconds" << "│" << right << endl;
        for(auto v: data_.getAproximationTour()) {
            cout << "│ " << left << "vertex: " << setw(4) << v->getInfo() << right <<setw(40) << "│"  << endl;
        }
    }

    cout << "│" << setw(53) << "│" << endl;
    cout << "└──────────────────────────────────────────────────┘" << endl;
    waitForEnter();
}


/**
 * @brief Draw the result of the Cluster Approximation Heuristic for the Traveling Salesman Problem (TSP).
 *
 * This function calculates the Cluster Approximation Heuristic solution for the TSP
 * starting from the specified vertex and draws the result along with the tour cost
 * and the time taken for the calculation.
 *
 * @param vertex_id The ID of the starting vertex for the TSP tour.
 */
void Menu::drawCluster(string vertex_id) {
    auto start = chrono::high_resolution_clock::now();
    data_.clusterApproximationTSP(vertex_id);
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> duration = end - start;

    cout << "┌─ Cluster Approximation Heuristic ────────────────┐" << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getClusterTourCost() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration.count()) +  " seconds" << "│" << right << endl;
    for(auto v: data_.getClusterTour()) {
        cout << "│ " << left << "vertex: " << setw(4) << v->getInfo() << right <<setw(40) << "│"  << endl;
    }
    cout << "│" << setw(53) << "│" << endl;
    cout << "└──────────────────────────────────────────────────┘" << endl;
    waitForEnter();
}






/**
 * @brief Draw the result of the Approximation Heuristic Analysis for the Traveling Salesman Problem (TSP).
 *
 * This function calculates and compares the results of three approximation heuristics
 * (Triangular Approximation, Cluster Approximation, and MST Approximation) for the TSP
 * starting from the specified vertex. It draws the results of each heuristic along with
 * the tour cost and the time taken for the calculation.
 *
 * @param vertex_id The ID of the starting vertex for the TSP tour.
 */
void Menu::drawApproximationAnalysis(std::string vertex_id) {
    cout << "┌─ Approximation Heuristic Analysis ───────────────┐" << endl;
    cout << "│" << setw(53) << "│" << endl;

    auto start = chrono::high_resolution_clock::now();
    data_.backtrackingTSP();
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> duration = end - start;

    cout << "│  Backtracking Algorithm:                         │" << endl;
    cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getCost() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Tour size: " << right << left << setw(36) << data_.getBestTour().size() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration.count()) +  " seconds" << "│" << right << endl;
    cout << "│" << setw(53) << "│" << endl;
    
    auto start1 = chrono::high_resolution_clock::now();
    data_.triangularHeuristicAproximation(vertex_id);
    auto end1 = chrono::high_resolution_clock::now();

    chrono::duration<double> duration1 = end1 - start1;

    cout << "│  Triangular Approximation Heuristic :            │" << endl;
    cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getAproximationTourCost() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Tour size: " << right << left << setw(36) << data_.getAproximationTour().size() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration1.count()) +  " seconds" << "│" << right << endl;
    cout << "│" << setw(53) << "│" << endl;

    auto start2 = chrono::high_resolution_clock::now();
    data_.clusterApproximationTSP(vertex_id);
    auto end2 = chrono::high_resolution_clock::now();

    chrono::duration<double> duration2 = end2 - start2;

    cout << "│  Cluster Approximation Heuristic :               │" << endl;
    cout << "│ " << left << setw(12) << "Tour cost: " << right << left << setw(36) << data_.getClusterTourCost() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Tour size: " << right << left << setw(36) << data_.getClusterTour().size() << right << " │" << endl;
    cout << "│ " << left << setw(12) << "Time taken: " << right << left <<  setw(37) << to_string(duration2.count()) +  " seconds" << "│" << right << endl;
    cout << "│" << setw(53) << "│" << endl;
    cout << "└──────────────────────────────────────────────────┘" << endl;
    waitForEnter();
}

double Menu::calculate_tour_cost(const std::vector<std::string>& tour) {
    int cost = 0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        Vertex* v1 = data_.getNetwork().findVertex(tour[i]);
        Vertex* v2 = data_.getNetwork().findVertex(tour[i + 1]);
        for (const auto& edge : v1->getAdj()) {
            if (edge->getDest()->getInfo() == v2->getInfo()) {
                cost += edge->getWeight();
                break;
            }
        }
    }
    return cost;
}

/**
 * @brief Executes the Traveling Salesman Problem (TSP) algorithm on real-world data.
 *
 * This function executes the TSP algorithm on real-world data using the fast method.
 * It prompts the user to input a start vertex ID and then computes the TSP tour
 * using the fast method. It displays the start node, the time taken to compute
 * the tour, the tour cost, and the tour itself. If no feasible tour exists,
 * it notifies the user.
 *
 * @param vertex_id The ID of the start vertex for the TSP algorithm.
 */
void Menu::drawTspRealWorld(std::string vertex_id) {
    auto start = chrono::high_resolution_clock::now();
    std::vector<std::string> tour = data_.tsp_real_world1(vertex_id);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;

    int count = 0;
    cout << "┌─ TSP in Real World ──────────────────────────────┐" << endl;
    cout << "│                                                  │" << endl;
    if(std::stoi(vertex_id) > 9)
        cout << "│ Start Node: " << left << setw(36) << vertex_id << "│" << endl;
    else
        cout << "│ Start Node: " << left << setw(37) << vertex_id << "│" << endl;

    if (!tour.empty()) {
        cout << "│ " << left << setw(12) << "Time taken:" << right << left << setw(37) << to_string(duration.count()) + " seconds" << "│" << endl;
        cout << "│ " << left << setw(12) << "Tour Cost:" << left << setw(37) << fixed << setprecision(2) << calculate_tour_cost(tour) << "│" << endl;
        cout << "│ " << left << setw(12) << "Tour:" << "[";
        for (size_t i = 0; i < tour.size(); ++i) {
            std::cout << "│ " << std::left << std::setw(48) << tour[i];
            if (i < tour.size() - 1) {
                std::cout << ", ";
            }
            std::cout << "│" << std::endl;
            count++;
        }
        std::cout << "│" << std::left << std::setw(48) << "]" << "│" << std::endl;
    } else {
        cout << "│ No feasible tour exists.                         │" << endl;
    }
    cout << "│                                                  │" << endl;
    cout << "└──────────────────────────────────────────────────┘" << endl;
    waitForEnter();
}

/**
 * @brief Executes the Traveling Salesman Problem (TSP) algorithm on real-world data.
 *
 * This function executes the TSP algorithm on real-world data using the cost-efficient method.
 * It prompts the user to input a start vertex ID and then computes the TSP tour
 * using the cost-efficient method. It displays the start node, the time taken to compute
 * the tour, the tour cost, and the tour itself. If no feasible tour exists,
 * it notifies the user.
 *
 * @param vertex_id The ID of the start vertex for the TSP algorithm.
 */
void Menu::drawTspRealWorld2(std::string vertex_id) {
    auto start = chrono::high_resolution_clock::now();
    std::vector<std::string> tour = data_.tsp_real_world2(vertex_id);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;
    int count2=0;
    cout << "┌─ TSP in Real World ──────────────────────────────┐" << endl;
    cout << "│                                                  │" << endl;
    if(std::stoi(vertex_id) > 9)
        cout << "│ Start Node: " << left << setw(36) << vertex_id << "│" << endl;
    else
        cout << "│ Start Node: " << left << setw(37) << vertex_id << "│" << endl;

    if (!tour.empty()) {
        cout << "│ " << left << setw(12) << "Time taken:" << right << left << setw(37) << to_string(duration.count()) + " seconds" << "│" << endl;
        cout << "│ " << left << setw(12) << "Tour Cost:" << left << setw(37) << fixed << setprecision(2) << calculate_tour_cost(tour) << "│" << endl;
        cout << "│ " << left << setw(12) << "Tour:" << "[";
        for (size_t i = 0; i < tour.size(); ++i) {
            cout << tour[i];
            if (i < tour.size() - 1) {
                cout << ", ";
            }
            if(count2 % 10==0)
                cout << endl;
            count2++;
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
