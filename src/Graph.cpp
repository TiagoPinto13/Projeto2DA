
#include "../includes/Graph.h"

Vertex::Vertex(const std::string& in, const int type) : type_(type), info(in), visited(false), processing(false), indegree(0), dist(0), path(nullptr), queueIndex(0) {}

/**
 * @brief Adds an edge from this vertex to the destination vertex with the given weight.
 *
 * @param d The destination vertex.
 * @param w The weight of the edge.
 *
 * @return The pointer to the newly created edge.
 *
 * @complexity Time Complexity: O(1)
 */
Edge* Vertex::addEdge(Vertex* d, double w) {
    auto newEdge = new Edge(this, d, w);
    adj.push_back(newEdge);
    d->incoming.push_back(newEdge);
    return newEdge;
}

/**
 * @brief Removes the edge with the given destination vertex info from this vertex's adjacency list.
 *
 * @param in The info of the destination vertex of the edge to be removed.
 *
 * @return True if an edge was removed, false otherwise.
 *
 * @complexity Time Complexity: O(E), where E is the number of edges adjacent to this vertex.
 */
bool Vertex::removeEdge(const std::string& in) {
    bool removedEdge = false;
    auto it = adj.begin();
    while (it != adj.end()) {
        Edge* edge = *it;
        Vertex* dest = edge->getDest();
        if (dest->getInfo() == in) {
            it = adj.erase(it);
            deleteEdge(edge);
            removedEdge = true;
        }
        else {
            it++;
        }
    }
    return removedEdge;
}

/**
 * @brief Removes all outgoing edges from this vertex.
 *
 * @complexity Time Complexity: O(E), where E is the number of edges adjacent to this vertex.
 */
void Vertex::removeOutgoingEdges() {
    auto it = adj.begin();
    while (it != adj.end()) {
        Edge* edge = *it;
        it = adj.erase(it);
        deleteEdge(edge);
    }
}

/**
 * @brief Overloaded less than operator for comparing vertices based on their distance.
 *
 * @param vertex The vertex to compare with.
 *
 * @return True if this vertex's distance is less than the other vertex's distance, false otherwise.
 *
 * @complexity Time Complexity: O(1)
 */
bool Vertex::operator<(const Vertex& vertex) const {
    return this->dist < vertex.dist;
}

/**
 * @brief Gets the information stored in this vertex.
 *
 * @return The information stored in this vertex.
 *
 * @complexity Time Complexity: O(1)
 */
std::string Vertex::getInfo() const {
    return this->info;
}

/**
 * @brief Gets the adjacency list of this vertex.
 *
 * @return The adjacency list of this vertex.
 *
 * @complexity Time Complexity: O(1)
 */
std::vector<Edge*> Vertex::getAdj() const {
    return this->adj;
}

/**
 * @brief Checks if this vertex has been visited.
 *
 * @return True if this vertex has been visited, false otherwise.
 *
 * @complexity Time Complexity: O(1)
 */
bool Vertex::isVisited() const {
    return this->visited;
}

/**
 * @brief Checks if this vertex is currently being processed.
 *
 * @return True if this vertex is currently being processed, false otherwise.
 *
 * @complexity Time Complexity: O(1)
 */
bool Vertex::isProcessing() const {
    return this->processing;
}

/**
 * @brief Gets the indegree of this vertex.
 *
 * @return The indegree of this vertex.
 *
 * @complexity Time Complexity: O(1)
 */
unsigned int Vertex::getIndegree() const {
    return this->indegree;
}

/**
 * @brief Gets the distance of this vertex.
 *
 * @return The distance of this vertex.
 *
 * @complexity Time Complexity: O(1)
 */
double Vertex::getDist() const {
    return this->dist;
}

/**
 * @brief Gets the path stored in this vertex.
 *
 * @return The path stored in this vertex.
 *
 * @complexity Time Complexity: O(1)
 */
Edge* Vertex::getPath() const {
    return this->path;
}

/**
 * @brief Gets the incoming edges of this vertex.
 *
 * @return The incoming edges of this vertex.
 *
 * @complexity Time Complexity: O(1)
 */
std::vector<Edge*> Vertex::getIncoming() const {
    return this->incoming;
}

/**
 * @brief Sets the information stored in this vertex.
 *
 * @param info The information to set.
 *
 * @complexity Time Complexity: O(1)
 */
void Vertex::setInfo(const std::string& info) {
    this->info = info;
}


/**
 * @brief Sets the visited status of this vertex.
 *
 * @param visited The visited status to set.
 *
 * @complexity Time Complexity: O(1)
 */
void Vertex::setVisited(bool visited) {
    this->visited = visited;
}

/**
 * @brief Sets the processing status of this vertex.
 *
 * @param processing The processing status to set.
 *
 * @complexity Time Complexity: O(1)
 */
void Vertex::setProcessing(bool processing) {
    this->processing = processing;
}

/**
 * @brief Sets the indegree of this vertex.
 *
 * @param indegree The indegree to set.
 *
 * @complexity Time Complexity: O(1)
 */
void Vertex::setIndegree(unsigned int indegree) {
    this->indegree = indegree;
}

/**
 * @brief Gets the type of this vertex.
 *
 * @return The type of this vertex.
 *
 * @complexity Time Complexity: O(1)
 */
int Vertex::getType() const {
    return type_;
}

/**
 * @brief Sets the distance of this vertex.
 *
 * @param dist The distance to set.
 *
 * @complexity Time Complexity: O(1)
 */
void Vertex::setDist(double dist) {
    this->dist = dist;
}

/**
 * @brief Sets the path stored in this vertex.
 *
 * @param path The path to set.
 *
 * @complexity Time Complexity: O(1)
 */
void Vertex::setPath(Edge* path) {
    this->path = path;
}

/**
 * @brief Deletes the given edge from the incoming edges of its destination vertex.
 *
 * @param edge The edge to delete.
 *
 * @complexity Time Complexity: O(E), where E is the number of incoming edges of the destination vertex.
 */
void Vertex::deleteEdge(Edge* edge) {
    Vertex* dest = edge->getDest();
    auto it = dest->incoming.begin();
    while (it != dest->incoming.end()) {
        if ((*it)->getOrig()->getInfo() == info) {
            it = dest->incoming.erase(it);
        }
        else {
            it++;
        }
    }
    delete edge;
}

/**
 * @brief Constructor for Edge class.
 *
 * @param orig Pointer to the origin vertex of the edge.
 * @param dest Pointer to the destination vertex of the edge.
 * @param w Weight of the edge.
 *
 * @complexity Time Complexity: O(1)
 */
Edge::Edge(Vertex* orig, Vertex* dest, double w) : orig(orig), dest(dest), weight(w), reverse(nullptr), flow(0) {}

/**
 * @brief Gets the destination vertex of this edge.
 *
 * @return Pointer to the destination vertex of this edge.
 *
 * @complexity Time Complexity: O(1)
 */
Vertex* Edge::getDest() const {
    return this->dest;
}

/**
 * @brief Gets the weight of this edge.
 *
 * @return The weight of this edge.
 *
 * @complexity Time Complexity: O(1)
 */
double Edge::getWeight() const {
    return this->weight;
}

/**
 * @brief Gets the origin vertex of this edge.
 *
 * @return Pointer to the origin vertex of this edge.
 *
 * @complexity Time Complexity: O(1)
 */
Vertex* Edge::getOrig() const {
    return this->orig;
}

/**
 * @brief Gets the reverse edge of this edge.
 *
 * @return Pointer to the reverse edge of this edge.
 *
 * @complexity Time Complexity: O(1)
 */
Edge* Edge::getReverse() const {
    return this->reverse;
}

/**
 * @brief Checks if this edge is selected.
 *
 * @return True if this edge is selected, false otherwise.
 *
 * @complexity Time Complexity: O(1)
 */
bool Edge::isSelected() const {
    return this->selected;
}

/**
 * @brief Gets the flow of this edge.
 *
 * @return The flow of this edge.
 *
 * @complexity Time Complexity: O(1)
 */
double Edge::getFlow() const {
    return flow;
}

/**
 * @brief Sets the selected status of this edge.
 *
 * @param selected The selected status to set.
 *
 * @complexity Time Complexity: O(1)
 */
void Edge::setSelected(bool selected) {
    this->selected = selected;
}

/**
 * @brief Sets the reverse edge of this edge.
 *
 * @param reverse Pointer to the reverse edge to set.
 *
 * @complexity Time Complexity: O(1)
 */
void Edge::setReverse(Edge* reverse) {
    this->reverse = reverse;
}

/**
 * @brief Sets the flow of this edge.
 *
 * @param flow The flow to set.
 *
 * @complexity Time Complexity: O(1)
 */
void Edge::setFlow(double flow) {
    this->flow = flow;
}


/**
 * @brief Finds a vertex in the graph given its information.
 *
 * @param in The information of the vertex to find.
 *
 * @return Pointer to the vertex if found, nullptr otherwise.
 *
 * @complexity Time Complexity: O(1) on average, O(n) in worst case, where n is the number of vertices in the graph.
 */
Vertex* Graph::findVertex(const std::string& in) const {
    auto it = vertexMap.find(in);
    if (it != vertexMap.end()) {
        return it->second;
    }
    return nullptr;
}

/**
 * @brief Adds a vertex to the graph.
 *
 * @param in The information of the vertex to add.
 * @param type The type of the vertex to add.
 *
 * @return True if the vertex was added successfully, false otherwise.
 *
 * @complexity Time Complexity: O(1) on average, O(n) in worst case, where n is the number of vertices in the graph.
 */
bool Graph::addVertex(const std::string& in, const int type) {
    if (findVertex(in) != nullptr) {
        return false;
    }
    Vertex* vertex = new Vertex(in,type);
    vertexSet.push_back(vertex);
    vertexMap[in] = vertex;
    return true;
}

/**
 * @brief Removes a vertex from the graph.
 *
 * @param in The information of the vertex to remove.
 *
 * @return True if the vertex was removed successfully, false otherwise.
 *
 * @complexity Time Complexity: O(n * m), where n is the number of vertices in the graph and m is the maximum number of edges per vertex.
 */
bool Graph::removeVertex(const std::string& in) {
    auto it = vertexMap.find(in);
    if (it == vertexMap.end()) {
        return false;
    }
    Vertex* v = it->second;
    for (auto u : vertexSet) {
        if (u != v) {
            u->removeEdge(v->getInfo());
        }
        else {
            u->removeOutgoingEdges();
        }
    }
    vertexSet.erase(std::remove(vertexSet.begin(), vertexSet.end(), v), vertexSet.end());
    vertexMap.erase(it);
    delete v;
    // Após remover o vértice, atualize o fluxo residual das arestas

    return true;
}

/**
 * @brief Adds an edge to the graph.
 *
 * @param source The information of the source vertex.
 * @param dest The information of the destination vertex.
 * @param w The weight of the edge.
 *
 * @return True if the edge was added successfully, false otherwise.
 *
 * @complexity Time Complexity: O(n), where n is the number of vertices in the graph.
 */
bool Graph::addEdge(const std::string& source, const std::string& dest, double w) {
    Vertex* v1 = findVertex(source);
    Vertex* v2 = findVertex(dest);
    if (v1 == nullptr || v2 == nullptr) {
        return false;
    }
    v1->addEdge(v2, w);
    return true;
}

/**
 * @brief Removes an edge from the graph.
 *
 * @param source The information of the source vertex.
 * @param dest The information of the destination vertex.
 *
 * @return True if the edge was removed successfully, false otherwise.
 *
 * @complexity Time Complexity: O(n), where n is the number of vertices in the graph.
 */
bool Graph::removeEdge(const std::string& source, const std::string& dest) {
    Vertex* srcVertex = findVertex(source);
    if (srcVertex == nullptr) {
        return false;
    }
    return srcVertex->removeEdge(dest);
}

/**
 * @brief Adds a bidirectional edge to the graph.
 *
 * @param source The information of the source vertex.
 * @param dest The information of the destination vertex.
 * @param w The weight of the edge.
 *
 * @return True if the bidirectional edge was added successfully, false otherwise.
 *
 * @complexity Time Complexity: O(n), where n is the number of vertices in the graph.
 */
bool Graph::addBidirectionalEdge(const std::string& source, const std::string& dest, double w) {
    Vertex* v1 = findVertex(source);
    Vertex* v2 = findVertex(dest);
    if (v1 == nullptr || v2 == nullptr) {
        return false;
    }
    Edge* e1 = v1->addEdge(v2, w);
    Edge* e2 = v2->addEdge(v1, w);
    e1->setReverse(e2);
    e2->setReverse(e1);
    return true;
}

/**
 * @brief Gets the number of vertices in the graph.
 *
 * @return The number of vertices in the graph.
 *
 * @complexity Time Complexity: O(1)
 */
int Graph::getNumVertex() const {
    return vertexSet.size();
}

/**
 * @brief Gets the vector containing all vertices in the graph.
 *
 * @return Vector containing all vertices in the graph.
 *
 * @complexity Time Complexity: O(1)
 */
std::vector<Vertex*> Graph::getVertexSet() const {
    return vertexSet;
}

/**
 * @brief Gets the map containing the information of each vertex in the graph.
 *
 * @return Map containing the information of each vertex in the graph.
 *
 * @complexity Time Complexity: O(1)
 */
std::unordered_map<std::string, Vertex*> Graph::getVertexMap() const {
    return vertexMap;
}

/**
 * @brief Performs a depth-first search (DFS) traversal of the graph.
 *
 * @return Vector containing the information of the vertices visited during the DFS traversal.
 *
 * @complexity Time Complexity: O(V + E), where V is the number of vertices and E is the number of edges in the graph.
 */
std::vector<std::string> Graph::dfs() const {
    std::vector<std::string> res;
    for (auto v : vertexSet) {
        v->setVisited(false);
    }
    for (auto v : vertexSet) {
        if (!v->isVisited()) {
            dfsVisit(v, res);
        }
    }
    return res;
}

/**
 * @brief Performs a depth-first search (DFS) traversal of the graph starting from a specific vertex.
 *
 * @param source The information of the source vertex to start the DFS traversal.
 *
 * @return Vector containing the information of the vertices visited during the DFS traversal.
 *
 * @complexity Time Complexity: O(V + E), where V is the number of vertices and E is the number of edges in the graph.
 */
std::vector<std::string> Graph::dfs(const std::string& source) const {
    std::vector<std::string> res;
    Vertex* s = findVertex(source);
    if (s == nullptr) {
        return res;
    }
    for (auto v : vertexSet) {
        v->setVisited(false);
    }
    dfsVisit(s, res);
    return res;
}

/**
 * @brief Helper function for depth-first search (DFS) traversal.
 *
 * @param v Pointer to the vertex to start the DFS traversal.
 * @param res Vector to store the information of the vertices visited during the DFS traversal.
 *
 * @complexity Time Complexity: O(V + E), where V is the number of vertices and E is the number of edges in the graph.
 */
void Graph::dfsVisit(Vertex* v, std::vector<std::string>& res) const {
    v->setVisited(true);
    res.push_back(v->getInfo());
    for (auto& e : v->getAdj()) {
        Vertex* w = e->getDest();
        if (!w->isVisited()) {
            dfsVisit(w, res);
        }
    }
}

/**
 * @brief Performs a breadth-first search (BFS) traversal of the graph starting from a specific vertex.
 *
 * @param source The information of the source vertex to start the BFS traversal.
 *
 * @return Vector containing the information of the vertices visited during the BFS traversal.
 *
 * @complexity Time Complexity: O(V + E), where V is the number of vertices and E is the number of edges in the graph.
 */
std::vector<std::string> Graph::bfs(const std::string& source) const {
    std::vector<std::string> res;
    Vertex* s = findVertex(source);
    if (s == nullptr) {
        return res;
    }
    for (auto v : vertexSet) {
        v->setVisited(false);
    }
    std::queue<Vertex*> q;
    q.push(s);
    s->setVisited(true);
    while (!q.empty()) {
        Vertex* v = q.front();
        q.pop();
        res.push_back(v->getInfo());
        for (auto& e : v->getAdj()) {
            Vertex* w = e->getDest();
            if (!w->isVisited()) {
                q.push(w);
                w->setVisited(true);
            }
        }
    }
    return res;
}

/**
 * @brief Checks if the graph is a directed acyclic graph (DAG).
 *
 * @return True if the graph is a DAG, false otherwise.
 *
 * @complexity Time Complexity: O(V + E), where V is the number of vertices and E is the number of edges in the graph.
 */
bool Graph::isDAG() const {
    for (auto v : vertexSet) {
        v->setVisited(false);
        v->setProcessing(false);
    }
    for (auto v : vertexSet) {
        if (!v->isVisited()) {
            if (!dfsIsDAG(v)) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief Helper function for checking if the graph is a directed acyclic graph (DAG) using depth-first search (DFS).
 *
 * @param v Pointer to the vertex to start the DFS traversal.
 *
 * @return True if the graph is a DAG, false otherwise.
 *
 * @complexity Time Complexity: O(V + E), where V is the number of vertices and E is the number of edges in the graph.
 */
bool Graph::dfsIsDAG(Vertex* v) const {
    v->setVisited(true);
    v->setProcessing(true);
    for (auto e : v->getAdj()) {
        Vertex* w = e->getDest();
        if (w->isProcessing()) {
            return false;
        }
        if (!w->isVisited()) {
            if (!dfsIsDAG(w)) {
                return false;
            }
        }
    }
    v->setProcessing(false);
    return true;
}

/**
 * @brief Performs a topological sort of the graph.
 *
 * @return Vector containing the information of the vertices in topological order.
 *
 * @complexity Time Complexity: O(V + E), where V is the number of vertices and E is the number of edges in the graph.
 */
std::vector<std::string> Graph::topsort() const {
    std::vector<std::string> res;
    for (auto v : vertexSet) {
        v->setIndegree(0);
    }
    for (auto v : vertexSet) {
        for (auto e : v->getAdj()) {
            unsigned int indegree = e->getDest()->getIndegree();
            e->getDest()->setIndegree(indegree + 1);
        }
    }
    std::queue<Vertex*> q;
    for (auto v : vertexSet) {
        if (v->getIndegree() == 0) {
            q.push(v);
        }
    }
    while (!q.empty()) {
        Vertex* v = q.front();
        q.pop();
        res.push_back(v->getInfo());
        for (auto e : v->getAdj()) {
            Vertex* w = e->getDest();
            w->setIndegree(w->getIndegree() - 1);
            if (w->getIndegree() == 0) {
                q.push(w);
            }
        }
    }
    if (res.size() != vertexSet.size()) {
        res.clear();
    }
    return res;
}

/**
 * @brief Finds the index of a vertex in the graph given its information.
 *
 * @param in The information of the vertex to find.
 *
 * @return The index of the vertex if found, -1 otherwise.
 *
 * @complexity Time Complexity: O(V), where V is the number of vertices in the graph.
 */
int Graph::findVertexIdx(const std::string& in) const {
    for (unsigned i = 0; i < vertexSet.size(); i++) {
        if (vertexSet[i]->getInfo() == in) {
            return i;
        }
    }
    return -1;
}

/**
 * @brief Deletes a 2D integer matrix.
 *
 * @param m The matrix to delete.
 * @param n The number of rows in the matrix.
 *
 * @complexity Time Complexity: O(n^2), where n is the number of rows in the matrix.
 */
void Graph::deleteMatrix(int** m, int n) {
    if (m != nullptr) {
        for (int i = 0; i < n; i++) {
            if (m[i] != nullptr) {
                delete[] m[i];
            }
        }
        delete[] m;
    }
}
    
/**
 * @brief Deletes a 2D double matrix.
 *
 * @param m The matrix to delete.
 * @param n The number of rows in the matrix.
 *
 * @complexity Time Complexity: O(n^2), where n is the number of rows in the matrix.
 */
void Graph::deleteMatrix(double** m, int n) {
    if (m != nullptr) {
        for (int i = 0; i < n; i++) {
            if (m[i] != nullptr) {
                delete[] m[i];
            }
        }
        delete[] m;
    }
}







