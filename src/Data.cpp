//
// Created by antero on 25-04-2024.
//
#include "../headerFiles/Data.h"
#include <fstream>
#include <cmath>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <queue>
#include <limits>
#include "../headerFiles/MutablePriorityQueue.h"

using namespace std;


void Data::readNodes(std::string nodeFilePath) {
    ifstream nodesFile(nodeFilePath);
    if(nodesFile.fail()){
        ostringstream error_message;
        error_message << "Could not open file \"" << nodeFilePath << '"';
        throw ios_base::failure(error_message.str());
    }

    string id_str, longitude_str, latitude_str;

    getline(nodesFile, id_str);

    while(getline(nodesFile, id_str, ','), getline(nodesFile, longitude_str, ','),
    getline(nodesFile, latitude_str, '\n')) {


        double latitude = std::stod(latitude_str);
        double longitude = std::stod(longitude_str);

        network_.addVertex(id_str, longitude, latitude);
    }
}
void Data::readEdges(bool realWorldGraphs, std::string edgesFilePath) {  //bool to skip the first line, since in the realWorldGraphs there's a 1st line to skip
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

        double weight = std::stod(c3);

        network_.addEdge(c1, c2, weight);
    }
}

void Data::parseTOY(bool tourismCSV, std::string edgesFilePath) {  //bool to store the names of locals only present in the tourismCSV
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

                labelDestino.erase(std::remove(labelDestino.begin(), labelDestino.end(), '\r'), labelDestino.end());
                tourismLabels[destino] = labelDestino;
            }
        }

        else {
            getline(edgesFile, distancia, '\n');
        }

        double weight = std::stod(distancia);

        if(!network_.addEdge(origem, destino, weight)) {
            network_.addVertex(origem,0,0);
            network_.addVertex(destino,0,0);
            network_.addEdge(origem, destino, weight);
        }
    }
}

Graph Data::getNetwork() {
    return network_;
}

std::map<std::string, std::string> Data::getTourismLabels() {
    return tourismLabels;
}

std::vector<Vertex*> Data::backtrackingTSP() {
    bestTour.clear();
    bestCost = std::numeric_limits<double>::max();

    std::vector<Vertex*> currentTour;

    Vertex* v = network_.findVertex("0");
    currentTour.push_back(v);
    v->setVisited(true);

    backtrack(currentTour, 0);

    for (Vertex* vertex : bestTour) {
        vertex->setVisited(false);
    }

    return bestTour;
}

void Data::backtrack(std::vector<Vertex*>& currentTour, double currentCost) {
    if (currentTour.size() == network_.getVertexSet().size()) {
        double tourCost = calculateTourCost(currentTour);
        if (tourCost < bestCost) {
            bestTour = currentTour;
            bestCost = tourCost;
        }
        return;
    }


    Vertex* lastVertex = currentTour.back();
    for (Edge* edge : lastVertex->getAdj()) {
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
    double a = std::sin(dLat / 2) * std::sin(dLat / 2) +
               std::cos(lat1) * std::cos(lat2) *
               std::sin(dLon / 2) * std::sin(dLon / 2);
    double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));
    return EARTHRADIUS * c;
}

Vertex* Data::findNearestNeighbor(Vertex* v) {
    Vertex* nearestNeighbor = nullptr;
    double minDistance = std::numeric_limits<double>::max();
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

void Data::triangularHeuristicAproximation(const std::string& startNodeId) {
    aproximation_tour_.clear();
    aproximation_tourCost_ = 0.0;

    Vertex* startVertex = network_.findVertex(startNodeId);
    if (!startVertex) {
        std::cerr << "Start node not found in the graph.\n";
        return;
    }

    aproximation_tour_.push_back(startVertex);
    startVertex->setVisited(true);
    while (aproximation_tour_.size() < network_.getVertexSet().size()) {
        Vertex* lastVertex = aproximation_tour_.back();
        Vertex* nearestNeighbor = findNearestNeighbor(lastVertex);
        if (nearestNeighbor) {
            nearestNeighbor->setVisited(true);
            aproximation_tour_.push_back(nearestNeighbor);
            for (Edge* edge : lastVertex->getAdj()) {
                if (edge->getDest() == nearestNeighbor) {
                    aproximation_tourCost_ += edge->getWeight();
                    break;
                }
            }
        }else{
            break;
        }
    }
    for (Edge* edge : startVertex->getAdj()) {
        if (edge->getDest() == aproximation_tour_.back()) {
            aproximation_tourCost_ += edge->getWeight();
            break;
        }
    }
    aproximation_tour_.push_back(startVertex);
}

std::vector<Vertex*> Data::getAproximationTour() {
    return aproximation_tour_;
}

double Data::getAproximationTourCost() {
    return aproximation_tourCost_;
}

void Data::clusterApproximationTSP(const string& startNodeId){
    const auto& vertices = network_.getVertexSet();
    unordered_set<Vertex*> clusters;

    for (const auto& vertex : vertices) {
        clusters.insert(vertex);
    }

    cluster_tour_.clear();
    cluster_tourCost_ = 0.0;

    Vertex* startVertex = network_.findVertex(startNodeId);
    if (!startVertex) {
        cerr << "Start node not found in the graph.\n";
        return;
    }

    cluster_tour_.push_back(startVertex);
    startVertex->setVisited(true);

    while (!clusters.empty()) {
        Vertex* lastVertex = cluster_tour_.back();
        Vertex* nearestNeighbor = findNearestNeighborInCluster(lastVertex, clusters);

        if (nearestNeighbor) {
            nearestNeighbor->setVisited(true);
            cluster_tour_.push_back(nearestNeighbor);

            for (Edge* edge : lastVertex->getAdj()) {
                if (edge->getDest() == nearestNeighbor) {
                    cluster_tourCost_ += edge->getWeight();
                    break;
                }
            }

            clusters.erase(nearestNeighbor);
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
}

Vertex* Data::findNearestNeighborInCluster(Vertex* v, const unordered_set<Vertex*>& cluster) {
    Vertex* nearestNeighbor = nullptr;
    double minDistance = numeric_limits<double>::max();

    for (Edge* edge : v->getAdj()) {
        Vertex* neighbor = edge->getDest();
        if (!neighbor->isVisited() && cluster.count(neighbor)) {
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

    vector<bool> visited(network_.getVertexSet().size(), false); // Track visited vertices
    MutablePriorityQueue<Vertex> pq;

    pq.insert(startVertex);
    visited[stoi(startVertex->getInfo())] = true; // Mark start vertex as visited

    while (!pq.empty()) {
        Vertex* u = pq.extractMin();

        for (Edge* edge : u->getAdj()) {
            Vertex* v = edge->getDest();
            double weight = edge->getWeight();
            int vIndex = stoi(v->getInfo());
            if (!visited[vIndex]) {
                v->setParent(u); // Set parent pointer
                pq.insert(v);
                visited[vIndex] = true;
            }
        }
    }

    preorderTraversalMST(startVertex);
    mst_tourCost_ = calculateTourCost(mst_tour_);
}

void Data::preorderTraversalMST(Vertex* u) {
    mst_tour_.push_back(u);
    for (Edge* edge : u->getAdj()) {
        Vertex* v = edge->getDest();
        if (v->getParent() == u) { // Check if v is connected to u
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

