//
// Created by antero on 25-04-2024.
//
#include "../headerFiles/Data.h"
#include <fstream>
#include <cmath>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <limits>
#include <climits>
#include <stack>
#include <climits>
#include <list>

using namespace std;

/**
 * @brief Reads node data from a file and adds vertices to the network.
 *
 * @param nodeFilePath Path to the file containing node data.
 * @param numberOfNodes Number of nodes to read, or -1 to read all nodes.
 * @throws ios_base::failure if the file cannot be opened.
 *
 * @complexity O(N) where N is the number of nodes being read.
 */
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

            network_.addVertex(id_str, longitude, latitude, true);

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

            network_.addVertex(id_str, longitude, latitude, true);
        }
    }

}

/**
 * @brief Reads edge data from a file and adds edges to the network.
 *
 * @param realWorldGraphs Indicates if the input file has a header line to skip.
 * @param edgesFilePath Path to the file containing edge data.
 * @throws ios_base::failure if the file cannot be opened.
 *
 * @complexity O(E) where E is the number of edges being read.
 */
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


/**
 * @brief Parses a file to process edges and optionally store tourism labels.
 *
 * @param tourismCSV Indicates if the input file contains tourism labels.
 * @param edgesFilePath Path to the file containing edge data.
 * @throws ios_base::failure if the file cannot be opened.
 *
 * @complexity O(E) where E is the number of edges being read.
 */
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
            network_.addVertex(origem,0,0, false);
            network_.addVertex(destino,0,0, false);
            network_.addEdge(origem, destino, weight);
        }
        network_.addEdge(destino, origem, weight);
    }
}


/**
 * @brief Resets the visitation state of all nodes in the network.
 *
 * @complexity O(V) where V is the number of vertices in the network.
 */
void Data::resetNodesVisitation() {
    for(auto vertex : network_.getVertexSet()) {
        vertex->setVisited(false);
    }
}


/**
 * @brief Returns the network graph.
 *
 * @return Graph The network graph.
 *
 * @complexity O(1)
 */
Graph Data::getNetwork() {
    return network_;
}

/**
 * @brief Returns the map of tourism labels.
 *
 * @return map<string, string> The map of tourism labels.
 *
 * @complexity O(1)
 */
map<string, string> Data::getTourismLabels() {
    return tourismLabels;
}



/**
 * @brief Solves the TSP using a backtracking approach.
 *
 * @complexity O(N!) where N is the number of nodes.
 */
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

/**
 * @brief Helper function for backtracking TSP to explore all tours.
 *
 * @param currentTour Current tour path.
 * @param currentCost Current cost of the tour.
 *
 * @complexity O(N!) where N is the number of nodes.
 */
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


/**
 * @brief Calculates the cost of a given tour.
 *
 * @param tour The tour path as a vector of vertices.
 * @return double The total cost of the tour.
 *
 * @complexity O(N^2) where N is the number of nodes in the tour.
 */
double Data::calculateTourCost(const vector<Vertex*>& tour) const {
    double cost = 0;
    int nodenr=0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        Vertex* v1 = tour[i];
        Vertex* v2 = tour[i + 1];
        bool haveEdge=false;
        for (Edge* edge : v1->getAdj()) {
            if (edge->getDest()->getInfo() == v2->getInfo()) {
                cost += edge->getWeight();
                nodenr++;
                haveEdge=true;
                break;
            }
        }
        if(!haveEdge && v1->hasCoord() && v2->hasCoord()) {
            cost += haversineDistance(v1->getLat(),v1->getLong(),v2->getLat(),v2->getLong());
            nodenr++;
        }
    }
    if(nodenr+1 < aproximation_tour_.size()) {
        return -1;
    }
    return cost;
}


/**
 * @brief Returns the best tour cost found by the backtracking algorithm.
 *
 * @return double The cost of the best tour.
 *
 * @complexity O(1)
 */
double Data::getCost() {
    return bestCost;
}



/**
 * @brief Calculates the Haversine distance between two geographical points.
 *
 * @param lat1 Latitude of the first point.
 * @param lon1 Longitude of the first point.
 * @param lat2 Latitude of the second point.
 * @param lon2 Longitude of the second point.
 * @return double The distance between the two points.
 *
 * @complexity O(1)
 */


double Data::haversineDistance(double lat1, double lon1, double lat2, double lon2) const {

    lat1 *= M_PI / 180.0;
    lon1 *= M_PI / 180.0;
    lat2 *= M_PI / 180.0;
    lon2 *= M_PI / 180.0;

    const double EARTHRADIUS = 6371000;

    double dLat = lat2 - lat1;
    double dLon = lon2 - lon1;
    double a = sin(dLat / 2) * sin(dLat / 2) +
               cos(lat1) * cos(lat2) *
               sin(dLon / 2) * sin(dLon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    return EARTHRADIUS * c;
}


/**
 * @brief Finds the nearest unvisited neighbor of a given vertex.
 *
 * @param v The vertex to find the nearest neighbor for.
 * @return Vertex* Pointer to the nearest neighbor vertex.
 *
 * @complexity O(E) where E is the number of edges adjacent to the vertex.
 */
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


/**
 * @brief Finds the Minimum Spanning Tree (MST) of a given graph using Prim's algorithm.
 *
 * @param g Pointer to the graph.
 * @return vector<Vertex*> The vertices of the MST.
 *
 * @complexity O((V + E) log V) where V is the number of vertices and E is the number of edges.
 */
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


/**
 * @brief Performs a depth-first search (DFS) on the MST starting from a given vertex.
 *
 * @param v The starting vertex for DFS.
 * @param mst The vertices of the MST.
 *
 * @complexity O(V + E) where V is the number of vertices and E is the number of edges in the MST.
 */
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






/**
 * @brief Returns the tour found by the approximation algorithm.
 *
 * @return vector<Vertex*> The tour as a vector of vertices.
 *
 * @complexity O(1)
 */
void Data::createMstGraph(Graph &mstGraph, std::vector<Vertex*>  mst) {
    for(auto v : mst) {
        mstGraph.addVertex(v->getInfo(),v->getLong(),v->getLat(), v->hasCoord());
        auto ep = v->getPath();
        if (ep != nullptr) {
            if(!mstGraph.addBidirectionalEdge(ep->getOrig()->getInfo(),ep->getDest()->getInfo(),ep->getWeight())) {
                mstGraph.addVertex(ep->getOrig()->getInfo(),ep->getOrig()->getLong(),ep->getOrig()->getLat(), ep->getOrig()->hasCoord());
                mstGraph.addVertex(ep->getDest()->getInfo(),ep->getDest()->getLong(),ep->getDest()->getLat(), ep->getDest()->hasCoord());
                mstGraph.addBidirectionalEdge(ep->getOrig()->getInfo(),ep->getDest()->getInfo(),ep->getWeight());
            }
        }
    }
}


/**
 * @brief Approximates the TSP solution using a triangular heuristic starting from a given node.
 *
 * @param startNodeId The ID of the starting node.
 *
 * @complexity O((V + E) log V) where V is the number of vertices and E is the number of edges.
 */

void Data::triangularHeuristicAproximation(const string& startNodeId) {
    aproximation_tour_.clear();
    aproximation_tourCost_ = 0.0;
    Vertex* startVertex = network_.findVertex("0");
    if (!startVertex) {
        cerr << "Start node not found in the graph.\n";
        return;
    }
    std::vector<Vertex*> mst = prim(&network_);
    resetNodesVisitation();
    Graph mstGraph;
    createMstGraph(mstGraph,mst);
    auto tour = mstGraph.dfs();
    for(auto s: tour) {
        aproximation_tour_.push_back(network_.findVertex(s));
    }
    aproximation_tour_.push_back(startVertex);
    aproximation_tourCost_ = calculateTourCost(aproximation_tour_);
}

/**
 * @brief Returns the tour found by the approximation algorithm.
 *
 * @return vector<Vertex*> The tour as a vector of vertices.
 *
 * @complexity O(1)
 */
vector<Vertex*> Data::getAproximationTour() {
    return aproximation_tour_;
}


/**
 * @brief Returns the cost of the tour found by the approximation algorithm.
 *
 * @return double The cost of the tour.
 *
 * @complexity O(1)
 */
double Data::getAproximationTourCost() {
    return aproximation_tourCost_;
}


/**
 * @brief Approximates the TSP solution using a clustering approach starting from a given node.
 *
 * @param startNodeId The ID of the starting node.
 *
 * @complexity O(V^2) where V is the number of vertices.
 */
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


/**
 * @brief Finds the nearest neighbor for clustering approach considering unvisited nodes.
 *
 * @param v The current vertex.
 * @param unvisited Set of unvisited vertices.
 * @return Vertex* Pointer to the nearest neighbor vertex.
 *
 * @complexity O(E) where E is the number of edges adjacent to the vertex.
 */
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



/**
 * @brief Returns the tour found by the clustering approximation algorithm.
 *
 * @return vector<Vertex*> The tour as a vector of vertices.
 *
 * @complexity O(1)
 */
vector<Vertex*> Data::getClusterTour() {
    return cluster_tour_;
}


/**
 * @brief Returns the cost of the tour found by the clustering approximation algorithm.
 *
 * @return double The cost of the tour.
 *
 * @complexity O(1)
 */
double Data::getClusterTourCost() {
    return cluster_tourCost_;
}

/**
 * @brief Approximates the TSP solution using MST starting from a given node.
 *
 * @param startNodeId The ID of the starting node.
 *
 * @complexity O((V + E) log V) where V is the number of vertices and E is the number of edges.
 */
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


/**
 * @brief Performs a preorder traversal on the MST starting from a given vertex.
 *
 * @param u The starting vertex.
 *
 * @complexity O(V + E) where V is the number of vertices and E is the number of edges in the MST.
 */
void Data::preorderTraversalMST(Vertex* u) {
    mst_tour_.push_back(u);

    for (Edge* edge : u->getAdj()) {
        Vertex* v = edge->getDest();

        if (v->getParent() == u) {
            preorderTraversalMST(v);
        }
    }
}

/**
 * @brief Returns the tour found by the MST approximation algorithm.
 *
 * @return vector<Vertex*> The tour as a vector of vertices.
 *
 * @complexity O(1)
 */
vector<Vertex*> Data::getMSTTour() {
    return mst_tour_;
}


/**
 * @brief Returns the cost of the tour found by the MST approximation algorithm.
 *
 * @return double The cost of the tour.
 *
 * @complexity O(1)
 */
double Data::getMSTTourCost() {
    return mst_tourCost_;
}




/**
 * @brief Returns the best tour found by the backtracking algorithm.
 *
 * @return vector<Vertex*> The tour as a vector of vertices.
 *
 * @complexity O(1)
 */
std::vector<Vertex *> Data::getBestTour() {
    return bestTour;
}

/**
 * @brief Returns the cost of the best tour found by the backtracking algorithm.
 *
 * @return double The cost of the best tour.
 *
 * @complexity O(1)
 */
bool Data::isTourism() {
    return tourism;
}

const int INF = std::numeric_limits<int>::max();

/**
 * @brief Performs BFS to find the farthest node from a given start node.
 *
 * @param start The starting node.
 * @return string The farthest node found.
 *
 * @complexity O(V + E) where V is the number of vertices and E is the number of edges.
 */
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


/**
 * @brief Performs Dijkstra's algorithm to find shortest paths from a start node.
 *
 * @param start The starting node.
 * @return unordered_map<string, int> The shortest path distances from the start node.
 *
 * @complexity O((V + E) log V) where V is the number of vertices and E is the number of edges.
 */
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



/**
 * @brief Solves the Traveling Salesman Problem (TSP) using a heuristic approach for real-world graphs.
 *
 * @param start The starting node ID.
 * @return std::vector<std::string> The tour path as a vector of node IDs.
 *
 * @complexity O(V^2) where V is the number of vertices. This includes the BFS, Dijkstra, and the main loop.
 */
std::vector<std::string> Data::tsp_real_world1(std::string start) {
    if (network_.getVertexSet().empty() || !isConnected(start)) {
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


        if (next_node.empty()) {
            return {};
        }

        visited.insert(next_node);
        path.push_back(next_node);
        current_node = next_node;
    }

    path.push_back(start);
    return path;
}




/**
 * @brief Removes a vertex from the network.
 *
 * @param id The ID of the vertex to remove.
 *
 * @complexity O(V + E) where V is the number of vertices and E is the number of edges.
 */
void Data::removeVertex(string id) {
    if(network_.findVertex(id) != nullptr) {
        network_.removeVertex(id);
    }
    else {
        cerr << "Vertex not found in the graph.\n";
    }
}


/**
 * @brief Removes an edge between two vertices from the network.
 *
 * @param id1 The ID of the first vertex.
 * @param id2 The ID of the second vertex.
 *
 * @complexity O(E) where E is the number of edges.
 */
void Data::removeEdge(string id1, string id2) {
    if(network_.findVertex(id1) != nullptr && network_.findVertex(id2) != nullptr) {
        network_.removeEdge(id1, id2);
    }
    else {
        cerr << "One or both vertices not found in the graph.\n";
    }
}




/**
 * @brief Checks if the network is connected starting from a given node.
 *
 * @param start The starting node ID.
 * @return bool True if the network is connected, false otherwise.
 *
 * @complexity O(V + E) where V is the number of vertices and E is the number of edges.
 */
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



/**
 * @brief Optimizes a given TSP tour using the two-opt algorithm.
 *
 * @param tour The initial tour path.
 * @return std::vector<std::string> The optimized tour path.
 *
 * @complexity O(N^2 * E) where N is the number of nodes and E is the number of edges.
 */

std::vector<std::string> Data::twoOpt(const std::vector<std::string>& tour) {
    auto calculateTourCost = [&](const std::vector<std::string>& t) {
        int cost = 0;
        for (size_t i = 0; i < t.size() - 1; ++i) {
            Vertex* v1 = network_.findVertex(t[i]);
            Vertex* v2 = network_.findVertex(t[i + 1]);
            if (v1 && v2) {
                for (const auto& edge : v1->getAdj()) {
                    if (edge->getDest()->getInfo() == v2->getInfo()) {
                        cost += edge->getWeight();
                        break;
                    }
                }
            } else {
                return std::numeric_limits<int>::max();
            }
        }
        return cost;
    };

    std::vector<std::string> bestTour = tour;
    int bestCost = calculateTourCost(bestTour);

    bool improvement = true;
    int iteration = 0;
    int noImprovementCount = 0;
    const int maxNoImprovement = 100;

    while (improvement && noImprovementCount < maxNoImprovement) {
        improvement = false;

        for (size_t i = 1; i < bestTour.size() - 1; ++i) {
            for (size_t j = i + 2; j < bestTour.size(); ++j) {
                std::vector<std::string> newTour = bestTour;
                std::reverse(newTour.begin() + i, newTour.begin() + j);

                int newCost = calculateTourCost(newTour);

                if (newCost < bestCost) {
                    bestTour = newTour;
                    bestCost = newCost;
                    improvement = true;
                    noImprovementCount = 0;
                }
            }
        }

        iteration++;
        noImprovementCount++;
    }

    return bestTour;
}



/**
 * @brief Solves the TSP using a hybrid approach for real-world graphs, combining heuristic approach and two-opt heuristic.
 *
 * @param start The starting node ID.
 * @return std::vector<std::string> The optimized tour path.
 *
 * @complexity O((V + E) log V + N^2 * E) where V is the number of vertices, E is the number of edges, and N is the number of nodes in the tour.
 */
std::vector<std::string> Data::tsp_real_world2(const std::string start) {

    if (network_.getVertexSet().empty() || !isConnected(start)) {
        return {};
    }
    std::vector<std::string> tour = tsp_real_world1(start);
    std::vector<std::string> optimizedTour = twoOpt(tour);

    return optimizedTour;
}





/*
Vantagens:

Simplicidade de implementação.
Rápido para grafos pequenos.
Pode ser usado como base para algoritmos mais complexos.

 
Desvantagens:

Complexidade exponencial.
Não garante a solução ótima.
Não é eficiente para grafos grandes.

 */


