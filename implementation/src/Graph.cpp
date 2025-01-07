#include "Graph.hpp"
Graph::Graph(size_t order, bool isDirected, float probabilityOfEdge): 
	order(order), size(0), isDirected(isDirected), delta(0), Delta(0) {

    size_t connectedVertex {0};
    float probability {0.0};

    std::random_device randomNumber;
    std::mt19937 seed(randomNumber());
    std::uniform_int_distribution<int> gap(0, order - 1);
    std::uniform_real_distribution<float> probabilityGap(0.0, 1.0);

    for (size_t i {0}; i < order; ++i) {
        adjacency_list[i] = {};
    }

    for (size_t i {0}; i < order; ++i) {
        connectedVertex = gap(seed);

        while (i == connectedVertex) {
            connectedVertex = gap(seed);
		}
		
        if (!edgeExists(i, connectedVertex)) {
            addEdge(i, connectedVertex);
        }

        for (size_t j {i + 1}; j < order; ++j) {
            if (!edgeExists(i, j)) {
                probability = probabilityGap(seed);
                if (probabilityOfEdge >= probability) {
                    addEdge(i, j);
                }
           }
        }
    }
    
    delta = computeMinVertexDegree();
    Delta = computeMaxVertexDegree();
}

Graph::Graph(const std::string& filename, 
		  const std::vector<std::vector<std::string>>& cartesian_plane, 
		  const std::unordered_map<size_t, std::string>& map_vertices_to_coordenates,
		  const std::unordered_map<std::string, size_t>& map_coordenates_to_vertices,
		  bool isDirected): 
			order(0), size(0), isDirected(isDirected), delta(0), Delta(0),
			map_vertices_to_coordenates(map_vertices_to_coordenates),
			map_coordenates_to_vertices(map_coordenates_to_vertices),
			cartesian_plane(cartesian_plane) {
	
    std::ifstream file(filename);
    
    if (!file) {
        throw std::runtime_error("Error opening the file!");
    }
    
    size_t source, destination {0};
    std::string line {};
    
    while (std::getline(file, line)) {
        std::stringstream ssEdges(line);

        if (ssEdges >> source >> destination) {
            addVertex(source);
            addVertex(destination);
            addEdge(source, destination);
        }
        
    }
    
    delta = computeMinVertexDegree();
    Delta = computeMaxVertexDegree();
    
    file.close();
}

Graph::Graph(const Graph& graph): order(graph.order), size(graph.size),
	isDirected(graph.isDirected), delta(graph.delta),
	Delta(graph.Delta), adjacency_list(graph.adjacency_list), 
	map_vertices_to_coordenates(graph.map_vertices_to_coordenates),
	map_coordenates_to_vertices(graph.map_coordenates_to_vertices), 
	cartesian_plane(graph.cartesian_plane) {}

size_t Graph::computeMaxVertexDegree() {
	size_t maxDegree = getVertexDegree(getAdjacencyList().begin()->first);
    size_t temp {0};
    
    for (const auto& [i, _] : getAdjacencyList()) {
        temp = getVertexDegree(i);
        if (maxDegree < temp) {
            maxDegree = temp;
        }
    }
            
    return maxDegree;
}


size_t Graph::computeMinVertexDegree() {
 	size_t minDegree = getVertexDegree(getAdjacencyList().begin()->first);
    size_t temp {0};
    
    for (const auto& [i, _] : getAdjacencyList()) {
        temp = getVertexDegree(i);
        if (minDegree > temp) {
            minDegree = temp;
        }
    }
            
    return minDegree;
}

size_t Graph::getVertexDegree(size_t vertex) const {
    return !vertexExists(vertex) ? 0 : const_cast<const std::unordered_map<size_t, std::list<size_t>>&>(adjacency_list).at(vertex).size();
}

void Graph::addVertex(size_t source) {
	if (!vertexExists(source)) {
        adjacency_list[source] = {};
        ++order;
    }
}

void Graph::addEdge(size_t source, size_t destination) {
    if (source == destination) { return; }
    
    if (this->isDirected == false) {
        this->adjacency_list[source].emplace_back(destination);
        this->adjacency_list[destination].emplace_back(source);
    } 
    
    else {
        this->adjacency_list[source].emplace_back(destination);
    }

    this->size += 1;
}

size_t Graph::getMaxDegree() const { return Delta; }

size_t Graph::getMinDegree() const { return delta; }

size_t Graph::getSize() const { return this->size; }

size_t Graph::getOrder() const { return this->order; }

const std::unordered_map<size_t, std::list<size_t>>& Graph::getAdjacencyList() const { return this->adjacency_list; }

size_t Graph::coordenateToVertex(const std::string& coordenate) const {
    return this->map_coordenates_to_vertices.at(coordenate);
}

const std::string& Graph::vertexToCoordenate(size_t vertex) const {
    return this->map_vertices_to_coordenates.at(vertex);
}

const std::list<size_t>& Graph::getAdjacencyList(size_t vertex) const { return this->adjacency_list.at(vertex); }

const std::vector<std::vector<std::string>>& Graph::getCartesianPlane() const { return this->cartesian_plane; }

bool Graph::vertexExists(size_t vertex) const { return adjacency_list.find(vertex) != adjacency_list.end(); }

bool Graph::edgeExists(size_t u, size_t v) const {
    return std::find(adjacency_list.at(u).begin(), adjacency_list.at(u).end(), v) != adjacency_list.at(u).end();
}

size_t Graph::heuristic1(int x1, int y1, int x2, int y2) {
	return 10 * euclideanDistance(x1, y1, x2, y2);
}

size_t Graph::heuristic2(int x1, int y1, int x2, int y2) {
	return 10 * manhatannDistance(x1, y1, x2, y2);
}

size_t Graph::euclideanDistance(int x1, int y1, int x2, int y2) {
	return std::floor(std::sqrt(std::pow(std::abs(x1 - x2), 2) + std::pow(std::abs(y1 - y2), 2)));
}

size_t Graph::manhatannDistance(int x1, int y1, int x2, int y2) {
	return std::abs(x1 - x2) + std::abs(y1 - y2);
}

std::string Graph::f1(size_t x, size_t y) {
	return x == 0 ? "(" + std::to_string(x) + "," + std::to_string(y) + ")" :
					"(" + std::to_string(x - 1) + "," + std::to_string(y) + ")";
}

std::string Graph::f2(size_t x, size_t y) {
	return x == 30 ? "(" + std::to_string(x) + "," + std::to_string(y) + ")" :
					"(" + std::to_string(x + 1) + "," + std::to_string(y) + ")";
}

std::string Graph::f3(size_t x, size_t y) {
	return y == 0 ? "(" + std::to_string(x) + "," + std::to_string(y) + ")" :
					"(" + std::to_string(x) + "," + std::to_string(y - 1) + ")";
}

std::string Graph::f4(size_t x, size_t y) {
	return y == 30 ? "(" + std::to_string(x) + "," + std::to_string(y) + ")" :
					"(" + std::to_string(x) + "," + std::to_string(y + 1) + ")";
}

constexpr size_t Graph::c1() {
    return 10;
}

constexpr size_t Graph::c2() {
    return 15;
}

size_t Graph::c3(size_t steps) {
    return 10 + (std::abs(5 - static_cast<int>(steps)) % 6);
}

size_t Graph::c4(size_t steps) {
    return 5 + (std::abs(10 - static_cast<int>(steps)) % 11);
}

void Graph::printStats(const std::string& initial_state,
					   const std::string& search_objective,
	 				   const std::vector<std::string>& path,
	 			       size_t path_cost,
	 				   size_t generated_vertices_amount,
	 				   size_t visited_vertices_amount) {
	 			
	 			std::cout << "Initial state: " << initial_state << '\n'
	 					  << "Search objective: " << search_objective << '\n'
	 					  << "Path coordenates: " << search_objective << '\n';
	 					  
	 			for (const auto& coordenate: path) {
	 				std::cout << coordenate << ' ';
	 			}
	 			
 				std::cout << "\nPath cost: " << path_cost << '\n'
	 			 		  << "Venerated vertices amount: " << generated_vertices_amount << '\n'
	 			 		  << "Visited vertices amount: " << visited_vertices_amount << '\n';
}

void Graph::breadthFirstSearch(size_t u, size_t v) {
    std::vector<bool> visited(this->order, false);
    std::stack<std::string> last_visited;
    
    std::vector<std::string> path;
    path.reserve(this->order);
    
    std::queue<int> queue;
    
    size_t temp {0};
    size_t generated_vertices_amount {0};
    size_t visited_vertices_amount {0};
    size_t path_cost {0};
    size_t steps_root_to_objective_amount {0};

    visited[u] = true;
    queue.push(u);

    while (!queue.empty()) {
    	
        temp = queue.front();
        
        if (!last_visited.empty()) {
		    if ((last_visited.top().at(1) > vertexToCoordenate(temp).at(1)) || (last_visited.top().at(1) < vertexToCoordenate(temp).at(1))) {
		    	path_cost += c3(steps_root_to_objective_amount);
		    }
		    
		    else if ((last_visited.top().at(4) > vertexToCoordenate(temp).at(4)) || (last_visited.top().at(4) < vertexToCoordenate(temp).at(4))) {
		    	path_cost += c2();
		    }
	 	}
        
        last_visited.push(vertexToCoordenate(temp));
        queue.pop();
        
        ++visited_vertices_amount;
        path.push_back(vertexToCoordenate(temp));
      	
      	if (temp == v) {
  		    printStats(vertexToCoordenate(u), vertexToCoordenate(v), path, path_cost, generated_vertices_amount, visited_vertices_amount);
  		    return;
      	}

        for (const auto& neighbor : getAdjacencyList(temp)) {
        	if (!visited[neighbor]) {
            	visited[neighbor] = true;
                queue.push(neighbor);
                ++generated_vertices_amount;
            }
        }
    }
}

void Graph::depthFirstSearch(size_t u, size_t v) {
    std::vector<bool> visited(this->order, false);
    std::vector<std::string> path;
    std::stack<std::string> last_visited;
    
    path.reserve(this->order);
    
    std::stack<size_t> stack;
    
    size_t temp {0};
    size_t generated_vertices_amount {0};
    size_t visited_vertices_amount {0};
    size_t path_cost {0};
    size_t steps_root_to_objective_amount {0};
    
    visited[u] = true;
    stack.push(u);

    while (!stack.empty()) {
        temp = stack.top();
        
        if (!last_visited.empty()) {
		    if ((last_visited.top().at(1) > vertexToCoordenate(temp).at(1)) || (last_visited.top().at(1) < vertexToCoordenate(temp).at(1))) {
		    	path_cost += c3(steps_root_to_objective_amount);
		    }
		    
		    else if ((last_visited.top().at(4) > vertexToCoordenate(temp).at(4)) || (last_visited.top().at(4) < vertexToCoordenate(temp).at(4))) {
		    	path_cost += c2();
		    }
	 	}
	 
    	last_visited.push(vertexToCoordenate(temp));
        stack.pop();
        
        ++visited_vertices_amount;
        path.push_back(vertexToCoordenate(temp));
        
        if (temp == v) {
  		    printStats(vertexToCoordenate(u), vertexToCoordenate(v), path, path_cost, generated_vertices_amount, visited_vertices_amount);
  		    return;
      	}
      
        for (const auto& neighbor: getAdjacencyList(temp)) {
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                stack.push(neighbor);
                ++generated_vertices_amount;
            }
        }
    }
}

/*
void Graph::uniformCostSearch(size_t u, size_t v) {
    std::vector<bool> visited(this->order, false);

    std::priority_queue<std::pair<size_t, int>, std::vector<std::pair<size_t, int>>, std::greater<>> queue;

    queue.push({0, u});

    while (!queue.empty()) {
        auto [current_node, current_cost] = queue.top();
        queue.pop();

        if (visited[current_node]) {
        	continue;
        }

        visited[current_node] = true;

        if (current_node == v) {
            std::cout << "Destino alcançado com custo: " << current_cost << '\n';
            return;
        }

        for (const auto& neighbor : getAdjacencyList(current_node)) {
            if (!visited[neighbor]) {
                queue.push({neighbor, current_cost + 1});
            }
        }
    }

    std::cout << "Nenhum caminho encontrado de " << u << " para " << v << '\n';
}


void Graph::AStar(size_t u, size_t v) {
    std::vector<bool> visited(this->order, false);
    
    std::priority_queue<std::pair<size_t, int>, std::vector<std::pair<size_t, int>>, std::greater<>> queue;
    
    queue.push({u, 0}); 

    while (!queue.empty()) {
        auto [current_cost, current_node] = queue.top();
        queue.pop();

        visited[current_node] = true;

        if (current_node == v) {
            std::cout << "Destination reached with cost: " << current_cost << '\n';
            return;
        }
        
        for (const auto& [neighbor, edge_cost] : getAdjacencyList(current_node)) {
            if (!visited[neighbor]) {
                int newCost = current_cost + edge_cost; //+ heuristic1(u, v);
                queue.push({newCost, neighbor});
            }
        }
    }

    std::cout << "No path found from " << u << " to " << v << '\n';
}


void Graph::greedySearch(size_t u, size_t v) {
    std::vector<bool> visited(this->order, false);
    std::unordered_map<size_t, int> cost;         

    queue.push({0, u}); 
    cost[u] = 0;

    while (!queue.empty()) {
        auto [currentCost, currentNode] = queue.top();
        queue.pop();

        if (visited[currentNode]) { continue; };

        visited[currentNode] = true;

        if (currentNode == v) {
            std::cout << "Destination reached with cost: " << currentCost << '\n';
            return;
        }

        for (const auto& [neighbor, edgeCost] : getAdjacencyList(currentNode)) {
            if (!visited[neighbor]) {
                int newCost = currentCost + edgeCost;

                if (cost.find(neighbor) == cost.end() || newCost < cost[neighbor]) {
                    cost[neighbor] = newCost;
                    queue.push({newCost, neighbor});
                }
            }
        }
    }

    std::cout << "No path found from " << u << " to " << v << '\n';
}

*/
void Graph::deleteAdjacencyList(size_t vertex) {
    if (adjacency_list.find(vertex) == adjacency_list.end()) { return; }

    std::queue<size_t> toDelete;
    toDelete.push(vertex);
    
    int currentVertex { -1 };

    for (const auto& it: this->adjacency_list[vertex]) {
    	toDelete.push(it);
    }

    while (!toDelete.empty()) {
        currentVertex = toDelete.front();
        toDelete.pop();

        for (const auto& it: this->adjacency_list[currentVertex]) {
        	this->adjacency_list[it].remove(currentVertex);
       	}
       	
       	deleteVertex(currentVertex);
    }
}

void Graph::deleteVertex(size_t vertex) {
    this->size -= this->adjacency_list[vertex].size();
    this->adjacency_list.erase(vertex);
    --this->order;
}

std::ostream& operator<< (std::ostream& os, const Graph& graph) {
    for (const auto& [u, v] : graph.adjacency_list) {
        size_t vertex = u;  

        os << vertex << " ----> ";
        for (const auto& neighbor : v) { 
            os << neighbor << " ";  
    	}
    	
        os << '\n'; 	
 	}
 	   
    return os;
}
