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
		  const std::unordered_map<size_t, std::string>& map_vertices_to_coordinates,
		  const std::unordered_map<std::string, size_t>& map_coordinates_to_vertices,
		  bool isDirected): 
			order(0), size(0), isDirected(isDirected), delta(0), Delta(0),
			map_vertices_to_coordinates(map_vertices_to_coordinates),
			map_coordinates_to_vertices(map_coordinates_to_vertices),
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
	map_vertices_to_coordinates(graph.map_vertices_to_coordinates),
	map_coordinates_to_vertices(graph.map_coordinates_to_vertices), 
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

size_t Graph::coordinateToVertex(const std::string& coordinate) const {
    return this->map_coordinates_to_vertices.at(coordinate);
}

std::pair<size_t, size_t> Graph::vertexToCoordinate(size_t vertex) const {
    const std::string& coord_str = this->map_vertices_to_coordinates.at(vertex);
    
    std::string clean_str = coord_str.substr(1, coord_str.length() - 2);
    
    std::istringstream iss(clean_str);
    std::string x_str, y_str;
    std::getline(iss, x_str, ',');
    std::getline(iss, y_str);
    
    int x = std::stoi(x_str);
    int y = std::stoi(y_str);
    
    return {x, y};
}


const std::list<size_t>& Graph::getAdjacencyList(size_t vertex) const { return this->adjacency_list.at(vertex); }

const std::vector<std::vector<std::string>>& Graph::getCartesianPlane() const { return this->cartesian_plane; }

bool Graph::vertexExists(size_t vertex) const { return adjacency_list.find(vertex) != adjacency_list.end(); }

bool Graph::edgeExists(size_t u, size_t v) const {
    return std::find(adjacency_list.at(u).begin(), adjacency_list.at(u).end(), v) != adjacency_list.at(u).end();
}

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












// Apartir daqui

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
	 			 		  << "Generated vertices amount: " << generated_vertices_amount << '\n'
	 			 		  << "Visited vertices amount: " << visited_vertices_amount << '\n';
}

std::string Graph::coordinatesToString(const std::pair<size_t, size_t>& coordinates) {
	auto x = coordinates.first;
	auto y = coordinates.second;
	
	return "(" + std::to_string(x) + "," + std::to_string(y) + ")";
}

void Graph::cenaryChoose(std::stack<std::pair<size_t, size_t>>& last_visited,
                         size_t& path_cost, size_t steps_root_to_objective_amount, size_t cenary, size_t temp) {
                
    auto last_coords = last_visited.top(); 
    auto current_coords = vertexToCoordinate(temp);
    
    if (cenary == 1) {
        path_cost += c1();
    }
    
    else {      
        // Compara se "x" (eixo horizontal) mudou
        if ((current_coords.first > last_coords.first) || (current_coords.first < last_coords.first)) {
            if (cenary == 2) {
                path_cost += c2();
            }
            else if (cenary == 3) {
                path_cost += c3(steps_root_to_objective_amount);
            }
            else if (cenary == 4) {
                path_cost += c4(steps_root_to_objective_amount);
            }
        }
        // Compara se "y" (eixo vertical) mudou
        else if ((current_coords.second > last_coords.second) || (current_coords.second < last_coords.second)) {
            path_cost += c1();
        }		        
    }         
}

size_t Graph::neighborCost(size_t u, size_t v, size_t steps_root_to_objective_amount, size_t cenary) {
    auto current_coords = vertexToCoordinate(u);
    auto destination_coords = vertexToCoordinate(v);
    
    size_t path_cost {0};
    
    if (cenary == 1) {
        path_cost += c1();
    }
    
    else {      
        // Compara se "x" (eixo horizontal) mudou
        if ((current_coords.first > destination_coords.first) || (current_coords.first < destination_coords.first)) {
            if (cenary == 2) {
                path_cost += c2();
            }
            else if (cenary == 3) {
                path_cost += c3(steps_root_to_objective_amount);
            }
            else if (cenary == 4) {
                path_cost += c4(steps_root_to_objective_amount);
            }
        }
        // Compara se "y" (eixo vertical) mudou
        else if ((current_coords.second > destination_coords.second) || (current_coords.second < destination_coords.second)) {
            path_cost += c1();
        }		        
    }
    
    return path_cost;         
}

void Graph::breadthFirstSearch(const std::string& u, const std::string& v, size_t cenary) {
    std::vector<bool> visited(this->order, false);
    std::stack<std::pair<size_t, size_t>> last_visited;
    std::vector<std::string> path;
    path.reserve(this->order);
    
    std::queue<int> queue;
    
    size_t source { coordinateToVertex(u) };
    size_t destination { coordinateToVertex(v) };
    size_t temp {0};
    size_t generated_vertices_amount {0};
    size_t visited_vertices_amount {0};
    size_t path_cost {0};
    int steps_root_to_objective_amount {-1};
   
    visited[source] = true;
    queue.push(source);
    
    while (!queue.empty()) {
        temp = queue.front();
        
        ++steps_root_to_objective_amount;
        
        if (!last_visited.empty()) {
            cenaryChoose(last_visited, path_cost, steps_root_to_objective_amount, cenary, temp);
        }
        
        last_visited.push(vertexToCoordinate(temp));
        queue.pop();
        
        ++visited_vertices_amount;
        path.push_back(coordinatesToString(vertexToCoordinate(temp))); 
        
        if (temp == destination) {
            printStats(coordinatesToString(vertexToCoordinate(source)),
             			coordinatesToString(vertexToCoordinate(destination)),
             			path, path_cost, generated_vertices_amount,
             			visited_vertices_amount);
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

void Graph::depthFirstSearch(const std::string& u, const std::string& v, size_t cenary) {
    std::vector<bool> visited(this->order, false);
    std::vector<std::string> path;
    std::stack<std::pair<size_t, size_t>> last_visited;
    
    path.reserve(this->order);
    
    std::stack<size_t> stack;
    
    size_t source { coordinateToVertex(u) };
    size_t destination { coordinateToVertex(v) };
   	size_t temp {0};
    size_t generated_vertices_amount {0};
    size_t visited_vertices_amount {0};
    size_t path_cost {0};
    int steps_root_to_objective_amount {-1};
    
    stack.push(source);
    
    while (!stack.empty()) {
        temp = stack.top();
        stack.pop();
        
        if (visited[temp]) {
            continue;
        }
        
        visited[temp] = true;
        
        ++steps_root_to_objective_amount;
        
       	if (!last_visited.empty()) {
	        cenaryChoose(last_visited, path_cost, steps_root_to_objective_amount, cenary, temp);
        }
		
        last_visited.push(vertexToCoordinate(temp));

        
        ++visited_vertices_amount;
        path.push_back(coordinatesToString(vertexToCoordinate(temp)));

        if (temp == destination) {
            printStats(coordinatesToString(vertexToCoordinate(source)),
             			coordinatesToString(vertexToCoordinate(destination)),
             			path, path_cost, generated_vertices_amount,
             			visited_vertices_amount);
            return;
        }
        
        for (const auto& neighbor : getAdjacencyList(temp)) {
            if (!visited[neighbor]) {
                stack.push(neighbor);
                ++generated_vertices_amount;
            }
        }
    }
}

void Graph::uniformCostSearch(const std::string& u, const std::string& v, size_t cenary) {
    std::vector<bool> visited(this->order, false);
    std::vector<int> predecessor(this->order, -1);
    
    std::vector<std::string> path;

    using Node = std::pair<size_t, size_t>;
    std::priority_queue<Node, std::vector<Node>, std::greater<>> priority_queue;

    size_t source { coordinateToVertex(u) };
    size_t destination { coordinateToVertex(v) };

    size_t generated_vertices_amount {0};
    size_t visited_vertices_amount {0};
    size_t path_cost {0};
    int steps_root_to_objective_amount {-1};
    size_t new_cost {0};

    priority_queue.push({source, 0});

    while (!priority_queue.empty()) {
        auto [current_vertex, current_cost] = priority_queue.top();
        priority_queue.pop();

        if (visited[current_vertex]) {
            continue;
        }

        ++steps_root_to_objective_amount;
        
        visited[current_vertex] = true;
        ++visited_vertices_amount;
		
		if (current_vertex == destination) {
            path_cost = current_cost;
            
            size_t current = destination;
            while (current != source) {
                path.push_back(coordinatesToString(vertexToCoordinate(current)));
                current = predecessor[current];
            }
            
            path.push_back(coordinatesToString(vertexToCoordinate(source)));

            // Inverter o caminho para que ele fique da origem até o destino
            std::reverse(path.begin(), path.end());

            printStats(coordinatesToString(vertexToCoordinate(source)),
                       coordinatesToString(vertexToCoordinate(destination)),
                       path, path_cost, generated_vertices_amount,
                       visited_vertices_amount);
            return;
        }
        
        for (const auto& neighbor : getAdjacencyList(current_vertex)) {
            if (!visited[neighbor]) {
                new_cost = current_cost + neighborCost(current_vertex, neighbor, steps_root_to_objective_amount, cenary);
                priority_queue.push({neighbor, new_cost});
                ++generated_vertices_amount;
                
                // Armazenar o predecessor para reconstrução do caminho.
                if (predecessor[neighbor] == -1) {
                    predecessor[neighbor] = current_vertex;
                }
            }
        }
    }
}

void Graph::AStar(const std::string& u, const std::string& v, size_t cenary) {
    std::vector<bool> visited(this->order, false);
    std::vector<int> predecessor(this->order, -1);
    std::vector<std::string> path;

    // Node: {vertex, {cost_with_heuristic, cost_without_heuristic}}
    using Node = std::pair<size_t, std::pair<size_t, size_t>>;
    std::priority_queue<Node, std::vector<Node>, std::greater<>> priority_queue;

    size_t source { coordinateToVertex(u) };
    size_t destination { coordinateToVertex(v) };

    size_t generated_vertices_amount {0};
    size_t visited_vertices_amount {0};
    int steps_root_to_objective_amount {-1};
    size_t path_cost {0};
    size_t new_cost_with_heuristic {0};
    size_t acumulated_cost {0};

    priority_queue.push({source, {0, 0}});

    while (!priority_queue.empty()) {
        auto [current_vertex, costs] = priority_queue.top();
        auto [cost_with_heuristic, cost_without_heuristic] = costs;
      
        priority_queue.pop();

        if (visited[current_vertex]) {
            continue;
        }
		
        ++steps_root_to_objective_amount;
        
		visited[current_vertex] = true;
        ++visited_vertices_amount;
        
        if (current_vertex == destination) {
            path_cost = cost_without_heuristic;
            
            // Reconstrução do caminho a partir do destino
            size_t current = destination;
            
            while (current != source) {
                path.push_back(coordinatesToString(vertexToCoordinate(current)));
                current = predecessor[current];
            }
            
            path.push_back(coordinatesToString(vertexToCoordinate(source)));

            // Inverter o caminho para que ele fique da origem até o destino
            std::reverse(path.begin(), path.end());

            printStats(
                coordinatesToString(vertexToCoordinate(source)),
                coordinatesToString(vertexToCoordinate(destination)),
                path,
                path_cost,
                generated_vertices_amount,
                visited_vertices_amount
            );
            return;
        }

        for (const auto& neighbor : getAdjacencyList(current_vertex)) {
            if (!visited[neighbor]) {
                acumulated_cost = cost_without_heuristic + neighborCost(current_vertex, neighbor, steps_root_to_objective_amount, cenary);

                new_cost_with_heuristic = 
                    acumulated_cost + euclideanDistance(
                        vertexToCoordinate(neighbor).first, vertexToCoordinate(neighbor).second,
                        vertexToCoordinate(destination).first, vertexToCoordinate(destination).second
                    );

                priority_queue.push({neighbor, {new_cost_with_heuristic, acumulated_cost}});
                ++generated_vertices_amount;

                // Armazenar o predecessor para reconstrução do caminho.
                if (predecessor[neighbor] == -1) {
                    predecessor[neighbor] = current_vertex;
                }
            }
        }
    }
}


void Graph::greedySearch(const std::string& u, const std::string& v, size_t cenary) {
    std::vector<bool> visited(this->order, false);
    std::vector<int> predecessor(this->order, -1);
    std::vector<std::string> path;

    // Node: {vertex, {cost_with_heuristic, cost_without_heuristic}}
    using Node = std::pair<size_t, std::pair<size_t, size_t>>;
    std::priority_queue<Node, std::vector<Node>, std::greater<>> priority_queue;

    size_t source { coordinateToVertex(u) };
    size_t destination { coordinateToVertex(v) };

    size_t generated_vertices_amount {0};
    size_t visited_vertices_amount {0};
    int steps_root_to_objective_amount {-1};
    size_t path_cost {0};
    size_t new_cost_with_heuristic {0};
    size_t acumulated_cost {0};

    priority_queue.push({source, {0, 0}});

    while (!priority_queue.empty()) {
        auto [current_vertex, costs] = priority_queue.top();
        auto [cost_with_heuristic, cost_without_heuristic] = costs;
      
        priority_queue.pop();

        if (visited[current_vertex]) {
            continue;
        }
		
        ++steps_root_to_objective_amount;
        
		visited[current_vertex] = true;
        ++visited_vertices_amount;
        
        if (current_vertex == destination) {
            path_cost = cost_without_heuristic;

            // Reconstrução do caminho a partir do destino
            size_t current = destination;
            
            while (current != source) {
                path.push_back(coordinatesToString(vertexToCoordinate(current)));
                current = predecessor[current];
            }
            
            path.push_back(coordinatesToString(vertexToCoordinate(source)));

            // Inverter o caminho para que ele fique da origem até o destino
            std::reverse(path.begin(), path.end());

            printStats(
                coordinatesToString(vertexToCoordinate(source)),
                coordinatesToString(vertexToCoordinate(destination)),
                path,
                path_cost,
                generated_vertices_amount,
                visited_vertices_amount
            );
            
            return;
        }

        for (const auto& neighbor : getAdjacencyList(current_vertex)) {
            if (!visited[neighbor]) {
                acumulated_cost = cost_without_heuristic + neighborCost(current_vertex, neighbor, steps_root_to_objective_amount, cenary);

                new_cost_with_heuristic = 
                     	euclideanDistance(
                        vertexToCoordinate(neighbor).first, vertexToCoordinate(neighbor).second,
                        vertexToCoordinate(destination).first, vertexToCoordinate(destination).second
                    );

                priority_queue.push({neighbor, {new_cost_with_heuristic, acumulated_cost}});
                ++generated_vertices_amount;

                // Armazenar o predecessor para reconstrução do caminho.
                if (predecessor[neighbor] == -1) {
                    predecessor[neighbor] = current_vertex;
                }
            }
        }
    }
}
