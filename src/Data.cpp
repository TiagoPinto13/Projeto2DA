//
// Created by antero on 25-04-2024.
//
#include "../headerFiles/Data.h"
#include <fstream>
#include <cmath>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <limits>
#include <climits>
#include <stack>
#include <climits>
#include <list>

using namespace std;


void Data::readNodes(string nodeFilePath, int numberOfNodes) {
    ifstream nodesFile(nodeFilePath);
    if(nodesFile.fail()){
        ostringstream error_message;
        error_message << "Could not open file \"" << nodeFilePath << '"';
        throw ios_base::failure(error_message.str());
    }

    string id_str, longitude_str, latitude_str;

    getline(nodesFile, id_str);
    if(numberOfNodes==-1) {
        while(getline(nodesFile, id_str, ','), getline(nodesFile, longitude_str, ','),
                getline(nodesFile, latitude_str, '\n')) {


            double latitude = stod(latitude_str);
            double longitude = stod(longitude_str);

            network_.addVertex(id_str, longitude, latitude);
        }
    }
    else {
        int nr = numberOfNodes;
        while(nr!=0) {
            nr--;
            getline(nodesFile, id_str, ',');
            getline(nodesFile, longitude_str, ',');
            getline(nodesFile, latitude_str, '\n');

            double latitude = stod(latitude_str);
            double longitude = stod(longitude_str);

            network_.addVertex(id_str, longitude, latitude);
        }
    }

}
void Data::readEdges(bool realWorldGraphs, string edgesFilePath) {  //bool to skip the first line, since in the realWorldGraphs there's a 1st line to skip
    ifstream edgesFile(edgesFilePath);
    if(edgesFile.fail()){
        ostringstream error_message;
        error_message << "Could not open file \"" << edgesFilePath << '"';
        throw ios_base::failure(error_message.str());
    }

    tourism = false;

    string c1, c2, c3;

    if(realWorldGraphs) {
        getline(edgesFile,c1);
    }
    while(getline(edgesFile, c1, ','), getline(edgesFile, c2, ','),
            getline(edgesFile, c3, '\n')) {

        double weight = stod(c3);

        network_.addEdge(c1, c2, weight);
        network_.addEdge(c2, c1, weight);
    }
}

void Data::parseTOY(bool tourismCSV, string edgesFilePath) {  //bool to store the names of locals only present in the tourismCSV
    ifstream edgesFile(edgesFilePath);
    if(edgesFile.fail()){
        ostringstream error_message;
        error_message << "Could not open file \"" << edgesFilePath << '"';
        throw ios_base::failure(error_message.str());
    }

    if(tourismCSV) {
        tourism = tourismCSV; //so we can know, after processing the tourism.csv, that we processed it (if needed)
    }
    else {
        tourism = false;
    }

    string origem, destino, distancia, labelOrigem, labelDestino;

    getline(edgesFile, origem);

    while(getline(edgesFile, origem, ','), getline(edgesFile, destino, ',')) {



        if(tourismCSV) {
            getline(edgesFile, distancia, ',');
            getline(edgesFile, labelOrigem, ',');
            getline(edgesFile, labelDestino, '\n');
            if (!(tourismLabels.find(origem) != tourismLabels.end() && tourismLabels.find(destino) != tourismLabels.end())) {
                tourismLabels[origem] = labelOrigem;

                labelDestino.erase(remove(labelDestino.begin(), labelDestino.end(), '\r'), labelDestino.end());
                tourismLabels[destino] = labelDestino;
            }
        }

        else {
            getline(edgesFile, distancia, '\n');
        }

        double weight = stod(distancia);

        if(!network_.addEdge(origem, destino, weight)) {
            network_.addVertex(origem,0,0);
            network_.addVertex(destino,0,0);
            network_.addEdge(origem, destino, weight);
        }
        network_.addEdge(destino, origem, weight);
    }
}

void Data::resetNodesVisitation() {
    for(auto vertex : network_.getVertexSet()) {
        vertex->setVisited(false);
    }
}

Graph Data::getNetwork() {
    return network_;
}

map<string, string> Data::getTourismLabels() {
    return tourismLabels;
}

void Data::backtrackingTSP() {
    bestTour.clear();
    bestCost = numeric_limits<double>::max();
    vector<Vertex*> currentTour;
    Vertex* v = network_.findVertex("0");


    for(auto vertex : network_.getVertexSet()) {
        vertex->setVisited(false);
    }

    currentTour.push_back(v);
    backtrack(currentTour, 0);
    resetNodesVisitation();
}

void Data::backtrack(vector<Vertex*>& currentTour, double currentCost) {
    if (currentTour.size() == network_.getVertexSet().size()+1 && currentTour.back() == currentTour.front()) {
        if (currentCost < bestCost) {
            bestTour = currentTour;
            bestCost = currentCost;
        }
        return;
    }

    Vertex* lastVertex = currentTour.back();
    for (auto edge : lastVertex->getAdj()) {
        Vertex* neighbor = edge->getDest();
        if (!neighbor->isVisited()) {
            currentTour.push_back(neighbor);
            neighbor->setVisited(true);
            backtrack(currentTour, currentCost + edge->getWeight());
            neighbor->setVisited(false);
            currentTour.pop_back();
        }
    }
}

double Data::calculateTourCost(const vector<Vertex*>& tour) const {
    double cost = 0;

    for (size_t i = 0; i < tour.size() - 1; ++i) {
        Vertex* v1 = tour[i];
        Vertex* v2 = tour[i + 1];

        for (Edge* edge : v1->getAdj()) {
            if (edge->getDest() == v2) {
                cost += edge->getWeight();
                break;
            }
        }
    }
    return cost;
}

double Data::getCost() {
    return bestCost;
}

double haversineDistance(double lat1, double lon1, double lat2, double lon2){
    lat1 *= M_PI / 180.0;
    lon1 *= M_PI / 180.0;
    lat2 *= M_PI / 180.0;
    lon2 *= M_PI / 180.0;

    const double EARTHRADIUS = 6371.0;

    double dLat = lat2 - lat1;
    double dLon = lon2 - lon1;
    double a = sin(dLat / 2) * sin(dLat / 2) +
               cos(lat1) * cos(lat2) *
               sin(dLon / 2) * sin(dLon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    return EARTHRADIUS * c;
}

Vertex* Data::findNearestNeighbor(Vertex* v) {
    Vertex* nearestNeighbor = nullptr;
    double minDistance = numeric_limits<double>::max();
    for (Edge* edge : v->getAdj()) {
        Vertex* neighbor = edge->getDest();
        if (!neighbor->isVisited()) {
            double distance = edge->getWeight();
            if (distance < minDistance) {
                minDistance = distance;
                nearestNeighbor = neighbor;
            }
        }
    }

    return nearestNeighbor;
}

std::vector<Vertex *> Data::prim(Graph * g) {
    if (g->getVertexSet().empty()) {
        return g->getVertexSet();
    }
    for(auto v : g->getVertexSet()) {
        v->setDist(INT_MAX);
        v->setPath(nullptr);
        v->setVisited(false);
    }
    Vertex* s = g->getVertexSet().front();
    s->setDist(0);
    MutablePriorityQueue<Vertex> q;
    q.insert(s);
    while( ! q.empty() ) {
        auto v = q.extractMin();
        v->setVisited(true);
        for(auto &e : v->getAdj()) {
            Vertex* w = e->getDest();
            if (!w->isVisited()) {
                auto oldDist = w->getDist();
                if(e->getWeight() < oldDist) {
                    w->setDist(e->getWeight());

                    w->setPath(e);
                    if (oldDist == INT_MAX) {
                        q.insert(w);
                    }
                    else {
                        q.decreaseKey(w);
                    }
                }
            }
        }
    }

    return g->getVertexSet();
}

void Data::dfsMST(Vertex* v, const std::vector<Vertex*>& mst) {
    v->setVisited(true);

    if (std::find(aproximation_tour_.begin(), aproximation_tour_.end(), v) == aproximation_tour_.end()) {
        aproximation_tour_.push_back(v);
    }

    for (auto& edge : v->getAdj()) {
        Vertex* neighbor = edge->getDest();
        if (!neighbor->isVisited()) {
            aproximation_tourCost_ += edge->getWeight();
            dfsMST(neighbor, mst);
        }
    }
}



void Data::triangularHeuristicAproximation(const string& startNodeId) {
    aproximation_tour_.clear();
    aproximation_tourCost_ = 0.0;

    Vertex* startVertex = network_.findVertex(startNodeId);
    if (!startVertex) {
        cerr << "Start node not found in the graph.\n";
        return;
    }

    std::vector<Vertex*> mst = prim(&network_);

    for(auto v : network_.getVertexSet()) {
        v->setVisited(false);
    }
    Graph mstGraph;
    for(auto v : mst) {
        mstGraph.addVertex(v->getInfo(),v->getLong(),v->getLat());
        auto ep = v->getPath();
        if (ep != nullptr) {
            if(!mstGraph.addBidirectionalEdge(ep->getOrig()->getInfo(),ep->getDest()->getInfo(),ep->getWeight())) {
                mstGraph.addVertex(ep->getOrig()->getInfo(),ep->getOrig()->getLong(),ep->getOrig()->getLat());
                mstGraph.addVertex(ep->getDest()->getInfo(),ep->getDest()->getLong(),ep->getDest()->getLat());
                mstGraph.addBidirectionalEdge(ep->getOrig()->getInfo(),ep->getDest()->getInfo(),ep->getWeight());
            }
        }
    }
    //dfsMST(startVertex, mstGraph.getVertexSet());
    auto vector1 = mstGraph.dfs();
    for(auto s: vector1) {
        aproximation_tour_.push_back(network_.findVertex(s));
    }
    aproximation_tour_.push_back(startVertex);

    aproximation_tourCost_ = calculateTourCost(aproximation_tour_);
    //resetNodesVisitation();
}

vector<Vertex*> Data::getAproximationTour() {
    return aproximation_tour_;
}

double Data::getAproximationTourCost() {
    return aproximation_tourCost_;
}

void Data::clusterApproximationTSP(const string& startNodeId){
    const auto& vertices = network_.getVertexSet();
    unordered_set<Vertex*> unvisited(vertices.begin(), vertices.end());

    cluster_tour_.clear();
    cluster_tourCost_ = 0.0;

    Vertex* startVertex = network_.findVertex(startNodeId);
    if (!startVertex) {
        cerr << "Start node not found in the graph.\n";
        return;
    }

    cluster_tour_.push_back(startVertex);
    startVertex->setVisited(true);
    unvisited.erase(startVertex);

    while (cluster_tour_.size() < vertices.size()) {
        Vertex* lastVertex = cluster_tour_.back();
        Vertex* nearestNeighbor = findNearestNeighborCluster(lastVertex, unvisited);

        if (nearestNeighbor) {
            nearestNeighbor->setVisited(true);
            cluster_tour_.push_back(nearestNeighbor);

            for (Edge* edge : lastVertex->getAdj()) {
                if (edge->getDest() == nearestNeighbor) {
                    cluster_tourCost_ += edge->getWeight();
                    break;
                }
            }

            unvisited.erase(nearestNeighbor);
        } else {
            break;
        }
    }

    for (Edge* edge : startVertex->getAdj()) {
        if (edge->getDest() == cluster_tour_.back()) {
            cluster_tourCost_ += edge->getWeight();
            break;
        }
    }

    cluster_tour_.push_back(startVertex);
    resetNodesVisitation();
}


Vertex* Data::findNearestNeighborCluster(Vertex* v, const unordered_set<Vertex*>& unvisited) {
    Vertex* nearestNeighbor = nullptr;
    double minDistance = numeric_limits<double>::max();

    for (Edge* edge : v->getAdj()) {
        Vertex* neighbor = edge->getDest();
        if (!neighbor->isVisited() && unvisited.count(neighbor)) {
            double distance = edge->getWeight();
            if (distance < minDistance) {
                minDistance = distance;
                nearestNeighbor = neighbor;
            }
        }
    }

    return nearestNeighbor;
}



vector<Vertex*> Data::getClusterTour() {
    return cluster_tour_;
}

double Data::getClusterTourCost() {
    return cluster_tourCost_;
}


void Data::mstApproximationTSP(const string& startNodeId) {
    mst_tour_.clear();
    mst_tourCost_ = 0.0;

    Vertex* startVertex = network_.findVertex(startNodeId);
    if (!startVertex) {
        cerr << "Start node not found in the graph.\n";
        return;
    }

    MutablePriorityQueue<Vertex> pq;
    pq.insert(startVertex);

    vector<bool> visited(network_.getVertexSet().size(), false);
    visited[stoi(startVertex->getInfo())] = true;

    while (!pq.empty()) {
        Vertex* u = pq.extractMin();

        for (Edge* edge : u->getAdj()) {
            Vertex* v = edge->getDest();
            double weight = edge->getWeight();
            int vIndex = stoi(v->getInfo());

            if (!visited[vIndex]) {
                v->setParent(u);
                pq.insert(v);
                visited[vIndex] = true;
            }
        }
    }

    preorderTraversalMST(startVertex);

    mst_tourCost_ = calculateTourCost(mst_tour_);
    resetNodesVisitation();
}

void Data::preorderTraversalMST(Vertex* u) {
    mst_tour_.push_back(u);

    for (Edge* edge : u->getAdj()) {
        Vertex* v = edge->getDest();

        if (v->getParent() == u) {
            preorderTraversalMST(v);
        }
    }
}


vector<Vertex*> Data::getMSTTour() {
    return mst_tour_;
}

double Data::getMSTTourCost() {
    return mst_tourCost_;
}




/*TSP REAL WORLD
vector<string> Data::merge_tours(const vector<vector<string> >& tours) {
    vector<string> merged_tour;
    for (const auto& tour : tours) {
        if (merged_tour.empty()) {
            merged_tour = tour;
        } else {
            for (size_t i = 0; i < merged_tour.size(); ++i) {
                if (merged_tour[i] == tour[0]) {
                    merged_tour.insert(merged_tour.begin() + i, tour.begin() + 1, tour.end());
                    break;
                }
            }
        }
    }
    return merged_tour;
}

vector<string> Data::tsp_subgraph(const Graph& subgraph, string start) {
    vector<string> tour;
    tour.push_back(start);
    unordered_set<string> unvisited;
    for (const auto& vertex : subgraph.getVertexSet()) {
        unvisited.insert(vertex->getInfo());
    }
    unvisited.erase(start);
    string current_node = start;
    while (!unvisited.empty()) {
        vector<Edge*> adj;
        for (auto temp : subgraph.getVertexSet()) {
            if (temp->getInfo() == current_node) {
                adj = temp->getAdj();
                break;
            }
        }
        string next_node;
        double min_weight = numeric_limits<double>::infinity();
        for (const auto& edge : adj) {
            string dest_info = edge->getDest()->getInfo();
            if (unvisited.count(dest_info) && subgraph.getEdgeWeight(current_node, dest_info) < min_weight) {
                min_weight = subgraph.getEdgeWeight(current_node, dest_info);
                next_node = dest_info;
            }
        }
        if (!next_node.empty()) {
            tour.push_back(next_node);
            unvisited.erase(next_node);
            current_node = next_node;
        } else {
            break;
        }
    }
    tour.push_back(start);
    resetNodesVisitation();
    return tour;
}


vector<string> Data::tsp_real_world(const string& start_node) {
    vector<vector<string> > subgraph_tours;
    for (const auto& vertex : network_.getVertexSet()) {
        Graph subgraph;
        subgraph.addVertex(vertex->getInfo(), vertex->getLong(), vertex->getLat());
        for (const Edge* edge : vertex->getAdj()) {
            string dest_node = edge->getDest()->getInfo();
            double weight = edge->getWeight();
            subgraph.addEdge(vertex->getInfo(), dest_node, weight);
        }
        if (subgraph.getVertexSet().size() > 1) {
            subgraph_tours.push_back(tsp_subgraph(subgraph, vertex->getInfo()));
        }
    }

    vector<string> tour = merge_tours(subgraph_tours);

    // Verifica se o tour é válido
    if (tour.size() == network_.getVertexSet().size() + 1 && tour.front() == tour.back() && tour.front() == start_node) {
        return tour;
    } else {
        // Se o tour não for válido, retorna um tour vazio
        return vector<string>();
    }
}

std::vector<Vertex *> Data::getBestTour() {
    return bestTour;
}

bool Data::isTourism() {
    return tourism;
}*/






const int INF = std::numeric_limits<int>::max();


std::string Data::bfs_farthest_node(const std::string& start) {
    std::unordered_set<std::string> visited;
    std::queue<std::pair<std::string, int>> q; // par (nó, distância)
    q.push({start, 0});
    visited.insert(start);

    std::string farthest_node = start;
    int max_distance = 0;

    while (!q.empty()) {
        auto [node, dist] = q.front();
        q.pop();

        if (dist > max_distance) {
            max_distance = dist;
            farthest_node = node;
        }

        Vertex* vertex = network_.findVertex(node);
        if (vertex != nullptr) {
            for (const auto& edge : vertex->getAdj()) {
                std::string neighbor = edge->getDest()->getInfo();
                if (visited.find(neighbor) == visited.end()) {
                    visited.insert(neighbor);
                    q.push({neighbor, dist + 1});
                }
            }
        }
    }

    return farthest_node;
}

std::unordered_map<std::string, int> Data::dijkstra(const std::string& start) {
    std::unordered_map<std::string, int> distances;
    for (const auto& vertex : network_.getVertexSet()) {
        distances[vertex->getInfo()] = INF;
    }
    distances[start] = 0;

    std::priority_queue<std::pair<int, std::string>, std::vector<std::pair<int, std::string>>, std::greater<>> pq;
    pq.push({0, start});

    while (!pq.empty()) {
        auto [current_distance, current_node] = pq.top();
        pq.pop();

        if (current_distance > distances[current_node]) continue;

        Vertex* node = network_.findVertex(current_node);
        for (auto& edge : node->getAdj()) {
            Vertex* neighbor = edge->getDest();
            int weight = edge->getWeight();

            int distance = current_distance + weight;
            if (distance < distances[neighbor->getInfo()]) {
                distances[neighbor->getInfo()] = distance;
                pq.push({distance, neighbor->getInfo()});
            }
        }
    }

    return distances;
}
bool Data::isConnected(const std::string& start) {
    if (network_.getVertexSet().empty()) return false;

    std::unordered_set<std::string> visited;
    std::queue<std::string> q;
    q.push(start);
    visited.insert(start);

    while (!q.empty()) {
        std::string node = q.front();
        q.pop();

        Vertex* vertex = network_.findVertex(node);
        if (vertex != nullptr) {
            for (const auto& edge : vertex->getAdj()) {
                std::string neighbor = edge->getDest()->getInfo();
                if (visited.find(neighbor) == visited.end()) {
                    visited.insert(neighbor);
                    q.push(neighbor);
                }
            }
        }
    }

    return visited.size() == network_.getVertexSet().size();
}




std::vector<std::string> Data::tsp_real_world1(std::string start) {
    if (network_.getVertexSet().empty() || !isConnected(start)) {
        std::cout << "No feasible tour exists" << std::endl;
        return {};
    }

    std::string farthest_node = bfs_farthest_node(start);
    auto distances_from_farthest = dijkstra(farthest_node);

    std::unordered_set<std::string> visited;
    std::vector<std::string> path;
    std::string current_node = start;
    visited.insert(current_node);
    path.push_back(current_node);

    while (visited.size() < network_.getVertexSet().size()) {
        std::string next_node;
        int min_distance = INF;

        Vertex* node = network_.findVertex(current_node);
        if (node != nullptr) {
            for (const auto& edge : node->getAdj()) {
                std::string neighbor = edge->getDest()->getInfo();
                int weight = edge->getWeight();

                if (visited.find(neighbor) == visited.end() && weight < min_distance) {
                    min_distance = weight;
                    next_node = neighbor;
                }
            }
        }

        // Verificar se encontrou o próximo nó
        if (next_node.empty()) {
            std::cout << "No feasible tour exists" << std::endl;
            return {};
        }

        visited.insert(next_node);
        path.push_back(next_node);
        current_node = next_node;
    }

    path.push_back(start);  // Retornar ao início para completar o tour
    return path;
}





void Data::removeVertex(string id) {
    if(network_.findVertex(id) != nullptr) {
        network_.removeVertex(id);
    }
    else {
        cerr << "Vertex not found in the graph.\n";
    }
}

void Data::removeEdge(string id1, string id2) {
    if(network_.findVertex(id1) != nullptr && network_.findVertex(id2) != nullptr) {
        network_.removeEdge(id1, id2);
    }
    else {
        cerr << "One or both vertices not found in the graph.\n";
    }
}







Graph Data::copyGraph() const {
    Graph copy;
    for (const auto& vertex : network_.getVertexSet()) {
        copy.addVertex(vertex->getInfo(), vertex->getLong(), vertex->getLat());
    }
    for (const auto& vertex : network_.getVertexSet()) {
        for (const auto& edge : vertex->getAdj()) {
            copy.addEdge(vertex->getInfo(), edge->getDest()->getInfo(), edge->getWeight());
        }
    }
    return copy;
}

Graph Data::primMST(const std::string& start) {
    Graph mst;
    std::unordered_set<std::string> visited;
    auto comp = [](const std::pair<int, std::pair<std::string, std::string>>& a, const std::pair<int, std::pair<std::string, std::string>>& b) {
        return a.first > b.first;
    };
    std::priority_queue<std::pair<int, std::pair<std::string, std::string>>, std::vector<std::pair<int, std::pair<std::string, std::string>>>, decltype(comp)> pq(comp);

    auto addEdges = [&](const std::string& node) {
        Vertex* vertex = network_.findVertex(node);
        for (const auto& edge : vertex->getAdj()) {
            if (visited.find(edge->getDest()->getInfo()) == visited.end()) {
                pq.push({edge->getWeight(), {node, edge->getDest()->getInfo()}});
            }
        }
    };

    Vertex* startVertex = network_.findVertex(start);
    mst.addVertex(startVertex->getInfo(), startVertex->getLong(), startVertex->getLat());
    visited.insert(start);
    addEdges(start);

    while (visited.size() < network_.getVertexSet().size() && !pq.empty()) {
        auto [weight, nodes] = pq.top();
        pq.pop();
        auto [u, v] = nodes;
        if (visited.find(v) == visited.end()) {
            Vertex* vertex = network_.findVertex(v);
            mst.addVertex(vertex->getInfo(), vertex->getLong(), vertex->getLat());
            mst.addEdge(u, v, weight);
            visited.insert(v);
            addEdges(v);
        }
    }

    return mst;
}

std::vector<std::string> Data::findOddDegreeVertices(const Graph& mst) {
    std::unordered_map<std::string, int> degree;
    for (const auto& vertex : mst.getVertexSet()) {
        degree[vertex->getInfo()] = 0;
    }
    for (const auto& vertex : mst.getVertexSet()) {
        for (const auto& edge : vertex->getAdj()) {
            degree[vertex->getInfo()]++;
        }
    }

    std::vector<std::string> oddVertices;
    for (const auto& [vertex, deg] : degree) {
        if (deg % 2 != 0) {
            oddVertices.push_back(vertex);
        }
    }

    return oddVertices;
}

std::vector<std::pair<std::string, std::string>> Data::minimumWeightMatching(const std::vector<std::string>& oddVertices) {
    std::vector<std::pair<std::string, std::string>> matching;
    std::unordered_set<std::string> unmatched(oddVertices.begin(), oddVertices.end());

    while (!unmatched.empty()) {
        std::string u = *unmatched.begin();
        unmatched.erase(u);

        std::string bestMatch;
        int bestWeight = INF;
        for (const std::string& v : unmatched) {
            Vertex* vertex = network_.findVertex(u);
            Vertex* dest = network_.findVertex(v);
            for (const auto& edge : vertex->getAdj()) {
                if (edge->getDest()->getInfo() == dest->getInfo()) {
                    int weight = edge->getWeight();
                    if (weight < bestWeight) {
                        bestWeight = weight;
                        bestMatch = v;
                    }
                    break;
                }
            }
        }

        matching.emplace_back(u, bestMatch);
        unmatched.erase(bestMatch);
    }

    return matching;
}

std::vector<std::string> Data::findEulerCircuit(Graph& graph) {
    std::vector<std::string> circuit;
    std::stack<std::string> stack;
    std::unordered_map<std::string, std::list<std::string>> adjList;

    for (const auto& vertex : graph.getVertexSet()) {
        for (const auto& edge : vertex->getAdj()) {
            adjList[vertex->getInfo()].push_back(edge->getDest()->getInfo());
        }
    }

    std::string start = graph.getVertexSet()[0]->getInfo();
    stack.push(start);

    while (!stack.empty()) {
        std::string current = stack.top();
        if (adjList[current].empty()) {
            circuit.push_back(current);
            stack.pop();
        } else {
            std::string next = adjList[current].front();
            adjList[current].pop_front();
            stack.push(next);
        }
    }

    std::reverse(circuit.begin(), circuit.end());
    return circuit;
}

std::vector<std::string> Data::eulerToHamiltonian(const std::vector<std::string>& eulerCircuit) {
    std::unordered_set<std::string> visited;
    std::vector<std::string> hamiltonianCircuit;

    for (const std::string& node : eulerCircuit) {
        if (visited.find(node) == visited.end()) {
            hamiltonianCircuit.push_back(node);
            visited.insert(node);
        }
    }

    hamiltonianCircuit.push_back(hamiltonianCircuit[0]); // Completar o circuito
    return hamiltonianCircuit;
}

std::vector<std::string> Data::twoOpt(const std::vector<std::string>& tour) {
    auto calculateTourCost = [&](const std::vector<std::string>& t) {
        int cost = 0;
        for (size_t i = 0; i < t.size() - 1; ++i) {
            Vertex* v1 = network_.findVertex(t[i]);
            Vertex* v2 = network_.findVertex(t[i + 1]);
            for (const auto& edge : v1->getAdj()) {
                if (edge->getDest()->getInfo() == v2->getInfo()) {
                    cost += edge->getWeight();
                    break;
                }
            }
        }
        return cost;
    };

    std::vector<std::string> bestTour = tour;
    int bestCost = calculateTourCost(bestTour);

    bool improvement = true;
    while (improvement) {
        improvement = false;
        for (size_t i = 1; i < bestTour.size() - 2; ++i) {
            for (size_t j = i + 1; j < bestTour.size() - 1; ++j) {
                std::vector<std::string> newTour = bestTour;
                std::reverse(newTour.begin() + i, newTour.begin() + j);

                int newCost = calculateTourCost(newTour);
                if (newCost < bestCost) {
                    bestTour = newTour;
                    bestCost = newCost;
                    improvement = true;
                    break; // Adicionar um break para sair do loop interno
                }
            }
            if (improvement) {
                break; // Adicionar um break para sair do loop externo
            }
        }
    }

    return bestTour;
}


std::vector<std::string> Data::tsp_real_world2(std::string start) {
    if (network_.getVertexSet().empty() || !isConnected(start)) {
        std::cout << "No feasible tour exists" << std::endl;
        return {};
    }

    Graph originalGraph = copyGraph();
    Graph mst = primMST(start);
    std::vector<std::string> oddVertices = findOddDegreeVertices(mst);
    std::vector<std::pair<std::string, std::string>> matching = minimumWeightMatching(oddVertices);

    for (const auto& pair : matching) {
        auto firstVertex = originalGraph.findVertex(pair.first);
        auto secondVertex = originalGraph.findVertex(pair.second);
        for (const auto& edge : firstVertex->getAdj()) {
            if (edge->getDest()->getInfo() == pair.second) {
                mst.addEdge(pair.first, pair.second, edge->getWeight());
                break; // Adicionando o break para evitar adições duplicadas
            }
        }
    }

    std::vector<std::string> eulerCircuit = findEulerCircuit(mst);
    std::vector<std::string> hamiltonianCircuit = eulerToHamiltonian(eulerCircuit);
    std::vector<std::string> optimizedTour = twoOpt(hamiltonianCircuit);

    return optimizedTour;
}




/*
 Advantages of this approach:

Handles real-world graphs that are not fully connected.
Allows arbitrary starting points specified by the user.
Employs scalable algorithms like DFS or BFS for subgraph identification.
Utilizes heuristic methods for TSP within subgraphs, balancing efficiency and optimality.

Disadvantages:

May not guarantee an optimal solution due to the heuristic nature of TSP algorithms.
Complexity increases with the number of disconnected components in the graph.
Requires careful handling of merging subgraph tours to minimize tour length.

The branch-and-bound method:
 The problem is broken down into sub-problems in this approach. The solution of those individual sub-problems would provide an optimal solution.
Complexity: O(N^2 * 2^N)
 */