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

double Graph::heuristic1(int x1, int y1, int x2, int y2) {
	return 10 * euclideanDistance(x1, y1, x2, y2);
}

size_t Graph::heuristic2(int x1, int y1, int x2, int y2) {
	return 10 * manhatannDistance(x1, y1, x2, y2);
}

double Graph::euclideanDistance(int x1, int y1, int x2, int y2) {
	return std::sqrt(std::pow(std::abs(x1 - x2), 2) + std::pow(std::abs(y1 - y2), 2));
}

size_t Graph::manhatannDistance(int x1, int y1, int x2, int y2) {
	return std::abs(x1 - x2) + std::abs(y1 - y2);
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

std::string Graph::getStats(const std::string& initial_state,
					   const std::string& search_objective,
	 				   const std::vector<std::string>& path,
	 			       size_t path_cost,
	 				   size_t generated_vertices_amount,
	 				   size_t visited_vertices_amount) {
	 			
	 			std::cout << "Initial state: " << initial_state << '\n'
	 					  << "Search objective: " << search_objective << '\n'
	 					  << "Path coordenates: " << '\n';
	 					  
	 			for (const auto& coordenate: path) {
	 				std::cout << coordenate << ' ';
	 			}
	 			
 				std::cout << "\nPath cost: " << path_cost << '\n'
	 			 		  << "Generated vertices amount: " << generated_vertices_amount << '\n'
	 			 		  << "Visited vertices amount: " << visited_vertices_amount << '\n';
	return "string";
}

std::string Graph::coordinatesToString(const std::pair<size_t, size_t>& coordinates) {
	auto x = coordinates.first;
	auto y = coordinates.second;
	
	return "(" + std::to_string(x) + "," + std::to_string(y) + ")";
}

size_t Graph::neighborCost(size_t u, size_t v, size_t steps_root_to_objective_amount, size_t cenary) {
    auto current_coords { vertexToCoordinate(u) };
    auto destination_coords { vertexToCoordinate(v) };
    
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

std::vector<std::string> Graph::buildPath(size_t source, size_t destination, const std::vector<int>& predecessor) {
    std::vector<std::string> path;
    size_t current = destination;
    
    while (current != source) {
        path.push_back(coordinatesToString(vertexToCoordinate(current)));
        current = predecessor[current];
    }
    
    path.push_back(coordinatesToString(vertexToCoordinate(source)));
    std::reverse(path.begin(), path.end());
    
    return path;
}

std::string Graph::breadthFirstSearch(const std::string& u, const std::string& v, size_t cenary) {
    std::unordered_set<size_t> visited; 
    std::vector<int> predecessor(this->order, -1);
    std::queue<Node> queue;
    std::vector<std::string> path;
    path.reserve(this->order);
    
    size_t source { coordinateToVertex(u) };
    size_t destination { coordinateToVertex(v) };
    size_t generated_vertices_amount {0};
    size_t visited_vertices_amount {0};
    size_t new_path_cost {0};

    Node current_node;
    
    visited.emplace(source);
    
    queue.push(Node(source));
    
    while (!queue.empty()) {
        current_node = queue.front();
        queue.pop();

        ++visited_vertices_amount;

        if (current_node.vertex == destination) {
            path = buildPath(source, destination, predecessor);

            return getStats(
                coordinatesToString(vertexToCoordinate(source)),
                coordinatesToString(vertexToCoordinate(destination)),
                path,
                current_node.path_cost,
                generated_vertices_amount,
                visited_vertices_amount
            );
        }

        // Explora os vizinhos do nó atual
        for (const auto& neighbor : getAdjacencyList(current_node.vertex)) {
            if (!visited.count(neighbor)) {
                visited.emplace(neighbor);
                predecessor[neighbor] = current_node.vertex;

                // Calcula o custo acumulado do caminho
                new_path_cost = current_node.path_cost + 
                                       neighborCost(current_node.vertex, neighbor, current_node.height, cenary);

                // Cria o próximo nó e adiciona à fila
                queue.push(Node(neighbor, 0, 0, new_path_cost, current_node.height + 1));
                ++generated_vertices_amount;
            }
        }
    }

    // Caso o destino não seja alcançável
    return "Destination not reachable from source.\n";
}

std::string Graph::depthFirstSearch(const std::string& u, const std::string& v, size_t cenary) {
    std::unordered_set<size_t> visited;  
    std::vector<int> predecessor(this->order, -1);  
    std::stack<Node> stack;
    std::vector<std::string> path;
    path.reserve(this->order);

    size_t source { coordinateToVertex(u) };
    size_t destination { coordinateToVertex(v) };
    size_t generated_vertices_amount {0};
    size_t visited_vertices_amount {0};
    size_t new_path_cost {0};
	Node current_node;
	
    // Criação do nó inicial
    stack.push(Node(source));
    
    visited.emplace(source);

    while (!stack.empty()) {
        current_node = stack.top();
        stack.pop();
		
		if (!visited.count(current_node.vertex)) {
            visited.emplace(current_node.vertex);
		}
		
        ++visited_vertices_amount;

        // Verifica se o destino foi alcançado
        if (current_node.vertex == destination) {
            path = buildPath(source, destination, predecessor);

            return getStats(
                coordinatesToString(vertexToCoordinate(source)),
                coordinatesToString(vertexToCoordinate(destination)),
                path,
                current_node.path_cost,
                generated_vertices_amount,
                visited_vertices_amount
            );
        }

        // Explora os vizinhos do nó atual
        for (const auto& neighbor : getAdjacencyList(current_node.vertex)) {
            if (!visited.count(neighbor)) {

                predecessor[neighbor] = current_node.vertex;

                // Calcula o custo acumulado do caminho
                new_path_cost = current_node.path_cost + neighborCost(current_node.vertex, neighbor, current_node.height, cenary);

                // Cria o próximo nó e adiciona à pilha
                stack.push(Node(neighbor, 0, 0, new_path_cost, current_node.height + 1));
                ++generated_vertices_amount;
            }
        }
    }

    // Caso o destino não seja alcançável
    return "Destination not reachable from source.\n";
}


std::string Graph::AStar(const std::string& u, const std::string& v, size_t cenary) {
    std::unordered_map<size_t, size_t> cost_map; // Armazena f(n) = g(n) + h(n)
    std::priority_queue<Node, std::vector<Node>, std::greater<>> priority_queue;
    std::vector<std::string> path;

    std::vector<int> predecessor(this->order, -1);
    
    size_t source { coordinateToVertex(u) };
    size_t destination { coordinateToVertex(v) };

    size_t generated_vertices_amount {0};
    size_t visited_vertices_amount {0};

    Node current_node;

    // Insere o nó inicial na fila de prioridade
    priority_queue.push(Node(source));
    cost_map[source] = 0;

    while (!priority_queue.empty()) {
        current_node = priority_queue.top();
        priority_queue.pop();

        ++visited_vertices_amount;

        // Verifica se o destino foi alcançado
        if (current_node.vertex == destination) {
            path = buildPath(source, destination, predecessor);
            
            // Exibe as estatísticas
            return getStats(
                coordinatesToString(vertexToCoordinate(source)),
                coordinatesToString(vertexToCoordinate(destination)),
                path,
                current_node.path_cost,
                generated_vertices_amount,
                visited_vertices_amount
            );
        }

        // Verifica os vizinhos do vértice atual
        for (const auto& neighbor : getAdjacencyList(current_node.vertex)) {
            // Calcula o custo acumulado g(n)
            size_t new_acumulated_path_cost = current_node.path_cost +
                                              neighborCost(current_node.vertex, neighbor, current_node.height, cenary);

            // Calcula o custo heurístico h(n)
            size_t heuristic_cost = heuristic2(
                vertexToCoordinate(neighbor).first, vertexToCoordinate(neighbor).second,
                vertexToCoordinate(destination).first, vertexToCoordinate(destination).second
            );

            // f(n) = g(n) + h(n)
            size_t new_cost_with_heuristic = new_acumulated_path_cost + heuristic_cost;

            // Adiciona o nó à fila se:
            // - Não foi visitado ainda
            // - Ou encontrou um custo menor
            if ((cost_map.find(neighbor) == cost_map.end()) || (new_cost_with_heuristic < cost_map[neighbor])) {
				cost_map[neighbor] = new_cost_with_heuristic;
				predecessor[neighbor] = current_node.vertex;
				priority_queue.push(Node(neighbor, new_cost_with_heuristic, heuristic_cost, 
						                 new_acumulated_path_cost, current_node.height + 1));
				++generated_vertices_amount;
			}
        }
    }

    // Caso o destino não seja alcançável
    return "Destination not reachable from source.\n";
}


std::string Graph::uniformCostSearch(const std::string& u, const std::string& v, size_t cenary) {
    std::vector<int> predecessor(this->order, -1);
    std::vector<std::string> path;
    std::unordered_map<size_t, size_t> cost_map; // Armazena g(n)

    std::priority_queue<Node, std::vector<Node>, std::greater<>> priority_queue;

    size_t source { coordinateToVertex(u) };
    size_t destination { coordinateToVertex(v) };

    size_t generated_vertices_amount {0};
    size_t visited_vertices_amount {0};

    Node current_node;

    // Insere o nó inicial na fila de prioridade com custo acumulado 0
    priority_queue.push( { Node(source) } );
    cost_map[source] = 0;

    while (!priority_queue.empty()) {
        current_node = priority_queue.top();
        priority_queue.pop();

        ++visited_vertices_amount;

       	if (current_node.vertex == destination) {
        	path = buildPath(source, destination, predecessor);

            return getStats(
                coordinatesToString(vertexToCoordinate(source)),
                coordinatesToString(vertexToCoordinate(destination)),
                path,
                current_node.path_cost,
                generated_vertices_amount,
                visited_vertices_amount
            );
        }

        // Verifica os vizinhos do vértice atual
        for (const auto& neighbor : getAdjacencyList(current_node.vertex)) {
            // Calcula o custo acumulado g(n)
            size_t new_acumulated_path_cost = current_node.path_cost + 
                                         neighborCost(current_node.vertex, neighbor, current_node.height, cenary);

            // Adiciona o nó à fila se:
            // - Não foi visitado ainda
            // - Ou encontrou um custo menor
            if (cost_map.find(neighbor) == cost_map.end() || new_acumulated_path_cost < cost_map[neighbor]) {
                cost_map[neighbor] = new_acumulated_path_cost;
                predecessor[neighbor] = current_node.vertex;

                // Cria e insere o nó na fila de prioridade com custo acumulado atualizado
                priority_queue.push( { Node(neighbor, new_acumulated_path_cost, 0, new_acumulated_path_cost, current_node.height + 1) } );
                 					 
                ++generated_vertices_amount;
            }
        }
    }

    // Caso o destino não seja alcançável
    return "Destination not reachable from source.\n";
}

std::string Graph::greedySearch(const std::string& u, const std::string& v, size_t cenary) {
    std::vector<int> predecessor(this->order, -1);
    std::vector<std::string> path;
    std::unordered_map<size_t, size_t> cost_map; // Armazena f(n) = g(n) + h(n)

	std::priority_queue<Node, std::vector<Node>, std::greater<>> priority_queue;


    size_t source { coordinateToVertex(u) };
    size_t destination { coordinateToVertex(v) };

    size_t generated_vertices_amount {0};
    size_t visited_vertices_amount {0};

    Node current_node;

    // Insere o nó inicial na fila de prioridade
    priority_queue.push( {Node(source)} );
    cost_map[source] = 0;

    while (!priority_queue.empty()) {
        current_node = priority_queue.top();
        priority_queue.pop();

        ++visited_vertices_amount;

        if (current_node.vertex == destination) {
            path = buildPath(source, destination, predecessor);

            return getStats(
                coordinatesToString(vertexToCoordinate(source)),
                coordinatesToString(vertexToCoordinate(destination)),
                path,
                current_node.path_cost,
                generated_vertices_amount,
                visited_vertices_amount
            );
        }

        // Verifica os vizinhos do vértice atual
        for (const auto& neighbor : getAdjacencyList(current_node.vertex)) {
            // Calcula o custo acumulado g(n)
            size_t new_acumulated_path_cost = current_node.path_cost + 
                                         neighborCost(current_node.vertex, neighbor, current_node.height, cenary);

            // Calcula o custo heurístico h(n)
            size_t heuristic_cost = heuristic1(
                vertexToCoordinate(neighbor).first, vertexToCoordinate(neighbor).second,
                vertexToCoordinate(destination).first, vertexToCoordinate(destination).second
            );

            // Adiciona o nó à fila se:
            // - Não foi visitado ainda
            // - Ou encontrou um custo menor
            if (cost_map.find(neighbor) == cost_map.end()) {
                cost_map[neighbor] = heuristic_cost;
                predecessor[neighbor] = current_node.vertex;

                // Cria e insere o nó na fila de prioridade
                priority_queue.push( {Node(neighbor, heuristic_cost, 0, new_acumulated_path_cost, current_node.height + 1)} );
                ++generated_vertices_amount;
            }
        }
    }
    // Caso o destino não seja alcançável
    return "Destination not reachable from source.\n";
}
