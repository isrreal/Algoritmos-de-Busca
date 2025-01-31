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

std::vector<std::string> Graph::buildPath(NodePtr node_ptr) {
    std::vector<std::string> path;
    
    while (node_ptr) {
        path.push_back(coordinatesToString(vertexToCoordinate(node_ptr->vertex)));
        node_ptr = node_ptr->parent;
    }

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
    std::queue<NodePtr> queue;
    std::vector<std::string> path;
    std::vector<NodePtr> tree;
    
    path.reserve(this->order);
    
    size_t source { coordinateToVertex(u) };
    size_t destination { coordinateToVertex(v) };
    size_t generated_vertices_amount {0};
    size_t visited_vertices_amount {0};
    size_t new_path_cost {0};
	
    NodePtr current_node { std::make_shared<Node>(source) };
    
    visited.emplace(source);
    
    queue.push(current_node);
    
    while (!queue.empty()) {
        current_node = queue.front();
        queue.pop();

        ++visited_vertices_amount;

        if (current_node->vertex == destination) {
            path = buildPath(current_node);

            return getStats(
                coordinatesToString(vertexToCoordinate(source)),
                coordinatesToString(vertexToCoordinate(destination)),
                path,
                current_node->path_cost,
                generated_vertices_amount,
                visited_vertices_amount
            );
        }

        // Explora os vizinhos do nó atual
        for (const auto& neighbor : getAdjacencyList(current_node->vertex)) {
            if (!visited.count(neighbor)) {
                visited.emplace(neighbor);

                // Calcula o custo acumulado do caminho
                new_path_cost = current_node->path_cost + 
                                       neighborCost(current_node->vertex, neighbor, current_node->height, cenary);

                // Cria o próximo nó e adiciona à fila
                queue.push(std::make_shared<Node>(neighbor, 0, 0, new_path_cost, current_node, current_node->height + 1));
                ++generated_vertices_amount;
            }
        }
        
          // Armazena os vizinhos em um vetor
        

        // Embaralha os vizinhos se randomize estiver ativo
		if (randomize) {
			std::random_device rd;
			std::mt19937 gen(rd());
			
			// Extrai os elementos da fila para um vetor temporário
			std::vector<NodePtr> neighbors;
			
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
 *         cost, and coordinates of the found path.s
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
    std::stack<NodePtr> stack;
    std::vector<std::string> path;
    std::vector<NodePtr> tree;
    
    path.reserve(this->order);
    
    size_t source { coordinateToVertex(u) };
    size_t destination { coordinateToVertex(v) };
    size_t generated_vertices_amount {0};
    size_t visited_vertices_amount {0};
    size_t new_path_cost {0};
	
    NodePtr current_node { std::make_shared<Node>(source) };
    
    visited.emplace(source);
    
    stack.push(current_node);
    
    while (!stack.empty()) {
        current_node = stack.top();
        stack.pop();

        ++visited_vertices_amount;

        if (current_node->vertex == destination) {
            path = buildPath(current_node);

            return getStats(
                coordinatesToString(vertexToCoordinate(source)),
                coordinatesToString(vertexToCoordinate(destination)),
                path,
                current_node->path_cost,
                generated_vertices_amount,
                visited_vertices_amount
            );
        }

        // Explora os vizinhos do nó atual
        for (const auto& neighbor : getAdjacencyList(current_node->vertex)) {
            if (!visited.count(neighbor)) {
                visited.emplace(neighbor);

                // Calcula o custo acumulado do caminho
                new_path_cost = current_node->path_cost + 
                                       neighborCost(current_node->vertex, neighbor, current_node->height, cenary);

                // Cria o próximo nó e adiciona à fila
                stack.push(std::make_shared<Node>(neighbor, 0, 0, new_path_cost, current_node, current_node->height + 1));
                ++generated_vertices_amount;
            }
        }
        
          // Armazena os vizinhos em um vetor
        

        // Embaralha os vizinhos se randomize estiver ativo
		if (randomize) {
			std::random_device rd;
			std::mt19937 gen(rd());
			
			// Extrai os elementos da fila para um vetor temporário
			std::vector<NodePtr> neighbors;
			
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
 * @brief Implements the A* pathfinding algorithm to find the optimal path between two points in a graph.
 *
 * @param u The starting point as a string representation of coordinates (e.g., "x,y").
 * @param v The destination point as a string representation of coordinates (e.g., "x,y").
 * @param cenary An integer representing the environmental configuration that influences path costs.
 * @param heuristic An integer representing the heuristic function to be used (e.g., 1 for Manhattan, 2 for Euclidean).
 * @param drugstores A set of string coordinates representing drugstore locations. If provided, the algorithm prioritizes visiting a drugstore before continuing to the destination.
 * 
 * @return A string containing details of the path found, including the reconstructed path, total cost, 
 *         number of generated vertices, and number of visited vertices. If the destination is unreachable, 
 *         returns an appropriate error message.
 * 
 * @details
 * This function uses the A* algorithm to compute the shortest path from a starting node (u) to a destination 
 * node (v) in a weighted graph. The A* algorithm evaluates paths by minimizing the total cost function 
 * - f(n) = g(n) + h(n), where:
 * - g(n): The accumulated cost to reach node (n).
 * - h(n): The heuristic estimate of the cost to reach the destination from node (n).
 * 
 * The algorithm processes nodes by expanding the one with the lowest f(n) value first, ensuring an 
 * efficient and guided search toward the destination. 
 * 
 * If drugstores are provided, the algorithm prioritizes visiting the nearest drugstore before resuming 
 * the search for the destination. Once a drugstore is visited, the search restarts from the drugstore.
 * 
 * The algorithm maintains a priority queue for open nodes, a cost map to track f(n) values, and a 
 * search tree for backtracking the final path. Nodes are revisited if a cheaper path to them is found during exploration.
 * 
 * @note
 * - The heuristic function must be consistent (admissible and monotonic) to guarantee optimality.
 * - Drugstores are optional; if provided, the algorithm ensures at least one is visited before reaching the destination.
 * - Nodes are uniquely identified by their coordinates and mapped to vertices in the graph.
 * 
 * @exception
 * Returns an error message if the destination is unreachable from the source.
 */


std::string Graph::AStar(const std::string& u, const std::string& v, size_t cenary, size_t heuristic, const std::unordered_set<std::string>& drugstores) {
    std::unordered_map<size_t, size_t> cost_map; // guarda f(n) = g(n) + h(n)
    std::vector<NodePtr> tree;
    std::priority_queue<NodePtr, std::vector<NodePtr>, std::greater<>> priority_queue;
    std::vector<std::string> path;
    size_t source { coordinateToVertex(u) };
    size_t destination { coordinateToVertex(v) };

    size_t generated_vertices_amount {0};
    size_t visited_vertices_amount {0};
    size_t heuristic_cost {0};
    size_t new_acumulated_path_cost {0};
    
    NodePtr current_node { std::make_shared<Node>(source) };
    NodePtr drugstore_node;

    priority_queue.push(current_node);
    cost_map[source] = 0;

    while (!priority_queue.empty()) {
        current_node = priority_queue.top();
        priority_queue.pop();
        tree.push_back(current_node);
        ++visited_vertices_amount;

        // Verifica se a farmácia foi visitada
        if (!drugstore_node && !drugstores.empty()) {
            if (drugstores.count(coordinatesToString(vertexToCoordinate(current_node->vertex))) > 0) {
                drugstore_node = std::make_shared<Node>(*current_node);

                // Limpa a fila de prioridade para reiniciar a busca a partir da farmácia
                while (!priority_queue.empty()) { priority_queue.pop(); }
                
                priority_queue.push(drugstore_node);
                cost_map.clear();
                cost_map[drugstore_node->vertex] = drugstore_node->path_cost;
                continue;
            }
        }

        // Verifica se o destino foi alcançado
        if (current_node->vertex == destination) {
        
            // Reconstrução do caminho    
            path = buildPath(current_node);

            return getStats(coordinatesToString(vertexToCoordinate(source)),
                    coordinatesToString(vertexToCoordinate(destination)),
                    path,
                    current_node->path_cost,
                    generated_vertices_amount,
                    visited_vertices_amount);
        }

        for (const auto& neighbor : getAdjacencyList(current_node->vertex)) {
            // Caso haja farmácias e nenhuma não foi visitada, ignora o destino.
            if ((!drugstore_node && !drugstores.empty()) && (neighbor == destination)) { continue; }

            new_acumulated_path_cost = current_node->path_cost +
                                       neighborCost(current_node->vertex, neighbor, current_node->height, cenary);
                                       
            heuristic_cost = (heuristic == 1 ? heuristic1(
                                  vertexToCoordinate(neighbor).first, vertexToCoordinate(neighbor).second,
                                  vertexToCoordinate(destination).first, vertexToCoordinate(destination).second
                              ) : heuristic2(
                                  vertexToCoordinate(neighbor).first, vertexToCoordinate(neighbor).second,
                                  vertexToCoordinate(destination).first, vertexToCoordinate(destination).second
                              ));
			
            // Adiciona o nó à fila se:
            // - Não foi visitado ainda
            // - Ou encontrou um custo menor
            // f(n) = g(n) + h(n)
            if (cost_map.find(neighbor) == cost_map.end() || (new_acumulated_path_cost < cost_map[neighbor])) {
                cost_map[neighbor] = new_acumulated_path_cost;

                NodePtr node { std::make_shared<Node>(neighbor, cost_map[neighbor], heuristic_cost, new_acumulated_path_cost, current_node, current_node->height + 1) };

                priority_queue.push(node);
                tree.push_back(node);

                ++generated_vertices_amount;
            }
        }
    }

    return "Destination not reachable from source.\n";
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
 */

std::string Graph::uniformCostSearch(const std::string& u, const std::string& v, size_t cenary) {
    std::unordered_map<size_t, size_t> cost_map; // guarda f(n) = g(n) + 0
    std::vector<NodePtr> tree;
    std::priority_queue<NodePtr, std::vector<NodePtr>, std::greater<>> priority_queue;
    std::vector<std::string> path;
    size_t source { coordinateToVertex(u) };
    size_t destination { coordinateToVertex(v) };

    size_t generated_vertices_amount {0};
    size_t visited_vertices_amount {0};
    size_t new_acumulated_path_cost {0};
    NodePtr current_node { std::make_shared<Node>(source) };

    // Insere o nó inicial na fila de prioridade
    priority_queue.push(current_node);
    cost_map[source] = 0;

    while (!priority_queue.empty()) {
        current_node = priority_queue.top();
        priority_queue.pop();
        tree.push_back(current_node);
        ++visited_vertices_amount;

       	if (current_node->vertex == destination) {
        	path = buildPath(current_node);

            return getStats(
                coordinatesToString(vertexToCoordinate(source)),
                coordinatesToString(vertexToCoordinate(destination)),
                path,
                current_node->path_cost,
                generated_vertices_amount,
                visited_vertices_amount
            );
        }

        // Verifica os vizinhos do vértice atual
        for (const auto& neighbor : getAdjacencyList(current_node->vertex)) {
            // Calcula o custo acumulado g(n)
            new_acumulated_path_cost = current_node->path_cost + 
                                         neighborCost(current_node->vertex, neighbor, current_node->height, cenary);

            // Adiciona o nó à fila se:
            // - Não foi visitado ainda
            // - Ou encontrou um custo menor
            // f(n) = g(n) + 0
            if (cost_map.find(neighbor) == cost_map.end() || new_acumulated_path_cost < cost_map[neighbor]) {
                cost_map[neighbor] = new_acumulated_path_cost;

                // Cria e insere o nó na fila de prioridade, onde o valor de comparação do nó é somente seu custo acumulado desde o nó raiz.
         
                NodePtr node { std::make_shared<Node>(neighbor, new_acumulated_path_cost, 0, new_acumulated_path_cost, current_node, current_node->height + 1) };
                    
                priority_queue.push(node);

                tree.push_back(node);
                 					 
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
    std::unordered_map<size_t, size_t> cost_map; // guarda f(n) = 0 + h(n)
    std::vector<NodePtr> tree;
    std::priority_queue<NodePtr, std::vector<NodePtr>, std::greater<>> priority_queue;
    std::vector<std::string> path;
    
    size_t source { coordinateToVertex(u) };
    size_t destination { coordinateToVertex(v) };

    size_t generated_vertices_amount {0};
    size_t visited_vertices_amount {0};
    size_t heuristic_cost {0};
    size_t new_acumulated_path_cost {0};
    NodePtr current_node { std::make_shared<Node>(source) };

    // Insere o nó inicial na fila de prioridade
    priority_queue.push(current_node);
    cost_map[source] = 0;

    while (!priority_queue.empty()) {
        current_node = priority_queue.top();
        priority_queue.pop();
        
        tree.push_back(current_node);
        ++visited_vertices_amount;

       	if (current_node->vertex == destination) {
        	path = buildPath(current_node);

            return getStats(
                coordinatesToString(vertexToCoordinate(source)),
                coordinatesToString(vertexToCoordinate(destination)),
                path,
                current_node->path_cost,
                generated_vertices_amount,
                visited_vertices_amount
            );
        }

        // Verifica os vizinhos do vértice atual
        for (const auto& neighbor : getAdjacencyList(current_node->vertex)) {
            // Calcula o custo acumulado g(n)
            new_acumulated_path_cost = current_node->path_cost + 
                                         neighborCost(current_node->vertex, neighbor, current_node->height, cenary);
                                         
			heuristic_cost = (heuristic == 1 ? heuristic1(
                                  vertexToCoordinate(neighbor).first, vertexToCoordinate(neighbor).second,
                                  vertexToCoordinate(destination).first, vertexToCoordinate(destination).second
                              ) : heuristic2(
                                  vertexToCoordinate(neighbor).first, vertexToCoordinate(neighbor).second,
                                  vertexToCoordinate(destination).first, vertexToCoordinate(destination).second
                              ));
            // Adiciona o nó à fila se:
            // - Não foi visitado ainda
            // - Ou encontrou um custo menor
            // f(n) = 0 + h(n)
            		
            if (cost_map.find(neighbor) == cost_map.end() || (heuristic_cost < cost_map[neighbor])) {
                cost_map[neighbor] = heuristic_cost;

               	NodePtr node { std::make_shared<Node>(neighbor, heuristic_cost, heuristic_cost, new_acumulated_path_cost, current_node, current_node->height + 1) };

                priority_queue.push(node);
                tree.push_back(node);

                ++generated_vertices_amount;
            }
        }
    }

    // Caso o destino não seja alcançável
    return "Destination not reachable from source.\n";
} 
