#include "Graph.hpp"

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
    std::getline(iss, x_str, ' ');
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
	
	std::string coordinates {};
	
	for (const auto& it : path) {
		coordinates += it + ' ';
	}  
	
	return initial_state + "," + search_objective + "," + 
	 	   std::to_string(generated_vertices_amount) + "," +
	 	   std::to_string(visited_vertices_amount) + "," + std::to_string(path_cost) + "," + coordinates;
}

std::string Graph::coordinatesToString(const std::pair<size_t, size_t>& coordinates) {
	auto x = coordinates.first;
	auto y = coordinates.second;
	
	return "(" + std::to_string(x) + " " + std::to_string(y) + ")";
}

/**
 * @brief Calculates the cost of moving from one vertex to another based on the provided scenario.
 * 
 * This function calculates the cost associated with moving between two vertices in a graph, considering the provided scenario. 
 * The costs may vary depending on whether the movement is horizontal or vertical, and different scenarios apply 
 * specific cost functions.
 * 
 * @param u Source vertex.
 * @param v Destination vertex.
 * @param steps_root_to_objective_amount Number of steps from the root node to the objective.
 * @param cenary The scenario used to calculate the cost. The possible scenarios are:
 *               - Scenario 1: A fixed cost of 10 for all movements.
 *               - Scenario 2: A fixed cost of 15 for horizontal movements.
 *               - Scenario 3: A cost based on the distance to the objective. The cost is calculated as 
 *                 10 + (|5 - steps| % 6), where `steps` is the number of steps from the root to the objective.
 *               - Scenario 4: A cost based on the distance to the objective. The cost is calculated as 
 *                 5 + (|10 - steps| % 11), where `steps` is the number of steps from the root to the objective.
 * @return size_t The calculated cost for the path between the vertices.
 * 
 * @details The calculation distinguishes between horizontal and vertical movements. If the movement is horizontal 
 * (change in the x-axis), specific costs are applied depending on the scenario. Otherwise, for vertical movements 
 * (change in the y-axis), the cost from scenario 1 is used.
 */

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

/**
 * @brief Performs a Breadth-First Search (BFS) on the graph.
 * 
 * @param u Coordinates of the source vertex in string format.
 * @param v Coordinates of the destination vertex in string format.
 * @param cenary Scenario used for cost calculation.
 * @param randomize Indicates whether the neighbors should be explored in random order.
 * @return std::string Search statistics: initial state, search objective, number of generated and visited vertices,
 *         cost, and coordinates of the found path.
 * 
 * @details
 * Breadth-First Search (BFS) is an algorithm for traversing or searching through a graph, starting at the initial
 * vertex and exploring all neighboring vertices at the present depth level before moving on to vertices at the next depth level.
 * 
 * In BFS, the algorithm explores all possible neighbors of a vertex before proceeding to the next level of neighbors.
 * This guarantees that the shortest path to the destination (if it exists) will be found.
 * 
 * The search expands by visiting vertices level by level. BFS does not backtrack and guarantees finding the shortest path
 * (if one exists) in an unweighted graph or graph with equal weights on edges.
 * 
 * If the `randomize` parameter is set to `true`, the neighbors of each node will be explored in random order, which can
 * introduce variability into the search process. If `randomize` is `false`, the neighbors will be explored in a fixed order.
 * 
 * The search statistics include details such as the initial state, goal state, the number of vertices generated and visited
 * during the search, the total cost of the found path, and the coordinates of the path.
 *
 * @note Nodes can't be repeated.
 *
 */

std::string Graph::breadthFirstSearch(const std::string& u, const std::string& v, size_t cenary, bool randomize) {
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
        
          // Armazena os vizinhos em um vetor
        

        // Embaralha os vizinhos se randomize estiver ativo
		if (randomize) {
			std::random_device rd;
			std::mt19937 gen(rd());
			
			// Extrai os elementos da fila para um vetor temporário
			std::vector<Node> neighbors;
			
			while (!queue.empty()) {
				neighbors.push_back(queue.front());
				queue.pop();
			}
			
			// Embaralha os elementos no vetor
			std::shuffle(neighbors.begin(), neighbors.end(), gen);
			
			// Reconstrói a fila com a ordem embaralhada
			for (const auto& neighbor : neighbors) {
				queue.push(neighbor);
			}
		}
    }

    // Caso o destino não seja alcançável
    return "Destination not reachable from source.\n";
}

/**
 * @brief Performs a depth-first search (DFS) on the graph.
 * 
 * @param u Coordinates of the source vertex in string format.
 * @param v Coordinates of the destination vertex in string format.
 * @param cenary Scenario used for cost calculation.
 * @param randomize Indicates whether the neighbors should be explored in random order.
 * @return std::string Search statistics: initial state, search objective, number of generated and visited vertices,
 *         cost, and coordinates of the found path.
 * 
 * @details
 * Depth-first search (DFS) is an algorithm that explores as far as possible along each branch before backtracking. 
 * Starting at the source vertex, the algorithm recursively visits the next unexplored neighbor until the destination 
 * vertex is found or all paths have been explored. If the `randomize` parameter is set to true, the order in which 
 * neighbors are explored is randomized, introducing variability in the traversal.
 * 
 * This approach is useful for searching deep paths in the graph, but it does not guarantee the shortest path is found. 
 * DFS may explore paths unnecessarily, especially in graphs with cycles or large branching factors. The algorithm 
 * can be adapted with cycle detection to avoid infinite loops.
 *
 * @note Nodes can't be repeated.
 *
 */

std::string Graph::depthFirstSearch(const std::string& u, const std::string& v, size_t cenary, bool randomize) {
    std::unordered_set<size_t> visited;  
    std::vector<int> predecessor(this->order, -1);  
    std::unordered_set<Node, NodeHash> tree;
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
				
                Node node(neighbor, 0, 0, new_path_cost, current_node.height + 1);
                // Cria o próximo nó e adiciona à pilha
                stack.push(node);
                
                tree.emplace(node);

                ++generated_vertices_amount;
            }
        }
        
        if (randomize) {
			std::random_device rd;
			std::mt19937 gen(rd());
			
			// Extrai os elementos da fila para um vetor temporário
			std::vector<Node> neighbors;
			
			while (!stack.empty()) {
				neighbors.push_back(stack.top());
				stack.pop();
			}
			
			// Embaralha os elementos no vetor
			std::shuffle(neighbors.begin(), neighbors.end(), gen);
			
			// Reconstrói a fila com a ordem embaralhada
			for (const auto& neighbor : neighbors) {
				stack.push(neighbor);
			}
		}
    }

    // Caso o destino não seja alcançável
    return "Destination not reachable from source.\n";
}

/**
 * @brief Performs A* pathfinding algorithm to find the optimal path between two points, using search trees.
 *
 * @param u The starting point (initial state) as a string representation of coordinates.
 * @param v The destination point (goal state) as a string representation of coordinates.
 * @param cenary The scenario or environment configuration that influences the search.
 * @param heuristic The heuristic function used to guide the search (e.g., Manhattan, Euclidean).
 * 
 * @details
 * This function implements the A* algorithm for finding the shortest path from the starting node (u) 
 * to the destination node (v). The algorithm combines the actual cost to reach a node (g(n)) with an estimate
 * of the cost from that node to the goal (h(n)), where the total cost function is defined as f(n) = g(n) + h(n).
 * 
 * The search expands nodes by always choosing the one with the lowest total cost (f(n)), ensuring that 
 * the most promising nodes are explored first. At each step, the algorithm evaluates neighboring nodes 
 * and adds them to the open set if they are not already visited or if a cheaper path to them is found.
 * 
 * The A* search continues until the goal is reached or all possible paths have been explored. The function 
 * also takes into account the given scenario (cenary) and heuristic function to guide the search efficiently.
 *
 * @note Nodes can be repeated because the same vertex may appear with a lower cost than its previous value.
 *
 *
 */
 
std::pair<std::string, size_t> Graph::AStar(const std::string& u, const std::string& v, size_t cenary, size_t heuristic, const std::unordered_set<std::string>& drugstores) {
    std::unordered_map<size_t, size_t> cost_map; // Armazena f(n) = g(n) + h(n)
    std::unordered_set<Node, NodeHash> tree;
    std::priority_queue<Node, std::vector<Node>, std::greater<>> priority_queue;
    std::vector<std::string> path;

    std::vector<int> predecessor(this->order, -1);

    size_t source { coordinateToVertex(u) };
    size_t destination { coordinateToVertex(v) };

    size_t generated_vertices_amount {0};
    size_t visited_vertices_amount {0};
    size_t heuristic_cost {0};
    size_t new_acumulated_path_cost {0};
    bool drugstore_visited {false};
    bool destination_reached {false};
    Node current_node(source);

    // Insere o nó inicial na fila de prioridade
    priority_queue.push(current_node);
    cost_map[source] = 0;
    tree.emplace(current_node);

    while (!priority_queue.empty()) {
        current_node = priority_queue.top();
        priority_queue.pop();

        ++visited_vertices_amount;
		
		if (!drugstores.empty()) {
		    // Atualiza se uma farmácia foi visitada
		    drugstore_visited = drugstore_visited || (drugstores.count(coordinatesToString(vertexToCoordinate(current_node.vertex))) > 0);

		    // Verifica se o destino foi alcançado
		    if (current_node.vertex == destination) {
		        destination_reached = true;
		    }

		    // Verifica se o destino foi alcançado e a farmácia foi visitada
		    if (destination_reached && drugstore_visited) {
		        path = buildPath(source, destination, predecessor);

		        // Exibe as estatísticas
		        return std::make_pair(
		            getStats(
		                coordinatesToString(vertexToCoordinate(source)),
		                coordinatesToString(vertexToCoordinate(destination)),
		                path,
		                current_node.path_cost,
		                generated_vertices_amount,
		                visited_vertices_amount),
		            current_node.path_cost);
		    }
		}
		
		else {
		    if (current_node.vertex == destination) {
		        path = buildPath(source, destination, predecessor);

		        // Exibe as estatísticas
		        return std::make_pair(
		            getStats(
		                coordinatesToString(vertexToCoordinate(source)),
		                coordinatesToString(vertexToCoordinate(destination)),
		                path,
		                current_node.path_cost,
		                generated_vertices_amount,
		                visited_vertices_amount),
		            current_node.path_cost);
		    }
		}

        // Verifica os vizinhos do vértice atual
        for (const auto& neighbor : getAdjacencyList(current_node.vertex)) {
            // Calcula o custo acumulado g(n)
            new_acumulated_path_cost = current_node.path_cost + 
                                       neighborCost(current_node.vertex, neighbor, current_node.height, cenary);

            // Calcula o custo heurístico h(n)
            heuristic_cost = (heuristic == 1 ? heuristic1(
                vertexToCoordinate(neighbor).first, vertexToCoordinate(neighbor).second,
                vertexToCoordinate(destination).first, vertexToCoordinate(destination).second
            ) : heuristic2(
                vertexToCoordinate(neighbor).first, vertexToCoordinate(neighbor).second,
                vertexToCoordinate(destination).first, vertexToCoordinate(destination).second
            ));

            // Adiciona o nó à fila de prioridade se houver melhoria de custo ou nova visita a uma farmácia
            if (cost_map.find(neighbor) == cost_map.end() || (new_acumulated_path_cost < cost_map[neighbor])) {
                cost_map[neighbor] = new_acumulated_path_cost;
                predecessor[neighbor] = current_node.vertex;

                Node node(neighbor, 
                    cost_map[neighbor], 
                    heuristic_cost, 
                    new_acumulated_path_cost, 
                    current_node.height + 1, 
                    false); // Atualiza o status da farmácia
                    
                priority_queue.push(node);

                tree.emplace(node);

                ++generated_vertices_amount;
            }
        }
    }

    // Verifica se uma farmácia foi visitada, mesmo sem alcançar o destino
    if (!destination_reached && drugstore_visited) {
        return std::make_pair("Farmácia visitada, mas destino não alcançado.\n", 0);
    }

    return std::make_pair<std::string, size_t>("Destino inalcançável a partir da origem.\n", 0);
}


/**
 * @brief Performs uniform-cost search on the graph, using search trees
 * 
 * @param u Coordinates of the source vertex in string format.
 * @param v Coordinates of the destination vertex in string format.
 * @param cenary Scenario used for cost calculation.
 * @return std::string Search statistics: initial state, search objective, number of generated and visited vertices,
 *         cost, and coordinates of the found path.
 * 
 * @details
 * Uniform-cost search is a variant of the A* algorithm where the heuristic function h(n) is set to 0. 
 * This means the algorithm does not consider an estimate of the cost to the destination, and is guided solely 
 * by the accumulated cost of each node (g(n)).
 * 
 * At each expanded node, the algorithm checks its neighbors and adds them to the priority queue based on 
 * the total accumulated cost to reach that point. The node with the lowest accumulated cost is expanded first. 
 * The search continues until the objective is reached, or until all possible nodes are explored.
 * 
 * This algorithm is guaranteed to find the least-cost path in graphs with non-negative edge costs.
 *
 * @note Nodes can be repeated because the same vertex may appear with a lower cost than its previous value.
 *
 */


std::string Graph::uniformCostSearch(const std::string& u, const std::string& v, size_t cenary) {
    std::vector<int> predecessor(this->order, -1);
    std::vector<std::string> path;
    std::unordered_map<size_t, size_t> cost_map; // Armazena g(n)
    std::unordered_set<Node, NodeHash> tree;

    std::priority_queue<Node, std::vector<Node>, std::greater<>> priority_queue;

    size_t source { coordinateToVertex(u) };
    size_t destination { coordinateToVertex(v) };

    size_t generated_vertices_amount {0};
    size_t visited_vertices_amount {0};
    size_t new_acumulated_path_cost {0};

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
            new_acumulated_path_cost = current_node.path_cost + 
                                         neighborCost(current_node.vertex, neighbor, current_node.height, cenary);

            // Adiciona o nó à fila se:
            // - Não foi visitado ainda
            // - Ou encontrou um custo menor
            if (cost_map.find(neighbor) == cost_map.end() || new_acumulated_path_cost < cost_map[neighbor]) {
                cost_map[neighbor] = new_acumulated_path_cost;
                predecessor[neighbor] = current_node.vertex;

                // Cria e insere o nó na fila de prioridade, onde o valor de comparação do nó é somente seu custo acumulado desde o nó raiz.

                
                Node node(neighbor, new_acumulated_path_cost, 0, new_acumulated_path_cost, current_node.height + 1);
                    
                priority_queue.push(node);

                tree.emplace(node);
                 					 
                ++generated_vertices_amount;
            }
        }
    }

    // Caso o destino não seja alcançável
    return "Destination not reachable from source.\n";
}

/**
 * @brief Performs greedy search on the graph, using search trees.
 * 
 * @param u Coordinates of the source vertex in string format.
 * @param v Coordinates of the destination vertex in string format.
 * @param cenary Scenario used for cost calculation.
 * @return std::string Search statistics: initial state, search objective, number of generated and visited vertices,
 *         cost, and coordinates of the found path.
 * 
 * @details
 * In the greedy search algorithm, the initial node is visited, and its adjacent vertices are generated. The search is 
 * guided solely by the heuristic estimate h(n), which approximates the cost from the current node to the goal. 
 * The cost function f(n) is given by f(n) = 0 + h(n), making it a variant of A* where the accumulated cost g(n) is ignored.
 * 
 * The algorithm expands nodes based on the smallest estimated cost to reach the destination (h(n)), 
 * without considering the cost already incurred to reach the current node.
 * 
 * This approach may not always find the optimal path as it doesn't account for the actual cost of the path taken.
 *
 * @note Nodes can be repeated because the same vertex may appear with a lower cost than its previous value.
 *
 */

 
std::string Graph::greedySearch(const std::string& u, const std::string& v, size_t cenary, size_t heuristic) {
    std::vector<int> predecessor(this->order, -1);
    std::vector<std::string> path;
    std::unordered_map<size_t, size_t> cost_map; // Armazena f(n) = g(n) + h(n)

	std::priority_queue<Node, std::vector<Node>, std::greater<>> priority_queue;


    size_t source { coordinateToVertex(u) };
    size_t destination { coordinateToVertex(v) };

    size_t generated_vertices_amount {0};
    size_t visited_vertices_amount {0};
    
    size_t new_acumulated_path_cost {0};
    size_t heuristic_cost {0};

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
            new_acumulated_path_cost = current_node.path_cost + 
                                         neighborCost(current_node.vertex, neighbor, current_node.height, cenary);

            // Calcula o custo heurístico h(n)
            heuristic_cost = (heuristic == 1 ? heuristic1(
				vertexToCoordinate(neighbor).first, vertexToCoordinate(neighbor).second,
				vertexToCoordinate(destination).first, vertexToCoordinate(destination).second
      		) : heuristic2(
				vertexToCoordinate(neighbor).first, vertexToCoordinate(neighbor).second,
				vertexToCoordinate(destination).first, vertexToCoordinate(destination).second
			  ));

            // Adiciona o nó à fila se:
            // - Não foi visitado ainda
            if (cost_map.find(neighbor) == cost_map.end()) {
                cost_map[neighbor] = heuristic_cost;
                predecessor[neighbor] = current_node.vertex;

                // Cria e insere o nó na fila de prioridade, onde o valor de comparação do nó é somente seu custo heurístico
                priority_queue.push( { Node(neighbor, heuristic_cost, 0, new_acumulated_path_cost, current_node.height + 1) } );
                ++generated_vertices_amount;
            }
        }
    }
    // Caso o destino não seja alcançável
    return "Destination not reachable from source.\n";
}
