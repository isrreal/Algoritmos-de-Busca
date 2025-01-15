#include <iostream>
#include <unordered_map>
#include <fstream>
#include <string>
#include <vector>
#include "Graph.hpp"
#include <random>

std::unordered_map<std::string, size_t> map_coordinates_to_vertices;
std::unordered_map<size_t, std::string> map_vertices_to_coordinates;
std::vector<std::vector<std::string>> cartesian_plane;
constexpr size_t COORDINATES_AMOUNT = 31 * 31;


void buildMatrix() {
    size_t vertex {0};
    cartesian_plane.resize(31, std::vector<std::string>(31));
    std::string coordenate {};
    
    map_vertices_to_coordinates.reserve(COORDINATES_AMOUNT);
    map_coordinates_to_vertices.reserve(COORDINATES_AMOUNT);

    for (size_t x {0}; x < 31; ++x) {
        for (size_t y {0}; y < 31; ++y) {
            coordenate = "(" + std::to_string(x) + " " + std::to_string(y) + ")";
            map_vertices_to_coordinates[vertex] = coordenate;
            map_coordinates_to_vertices[coordenate] = vertex++;
            cartesian_plane[x][y] = coordenate;
        }
    }
}

std::string coordinateToString(size_t x, size_t y) {
    return "(" + std::to_string(x) + " " + std::to_string(y) + ")";
}

void buildGraphFile() {
    buildMatrix();
    std::ofstream graph_edges("graph.txt");

    if (!graph_edges.is_open()) {
        std::cerr << "Erro ao abrir o arquivo 'graph.txt' para escrita." << std::endl;
        return;
    }
    
    std::string current_coordenate;
    std::string neighbor_down;
    std::string neighbor_right;
    
    for (size_t x {0}; x < 31; ++x) {
        for (size_t y {0}; y < 31; ++y) {

            current_coordenate = coordinateToString(x, y);
            
            if ((x + 1) < 31) {
                neighbor_down = coordinateToString(x + 1, y);
                graph_edges << map_coordinates_to_vertices[current_coordenate] << " " 
                            << map_coordinates_to_vertices[neighbor_down] << '\n';
            }
            
            if ((y + 1) < 31) {
                neighbor_right = coordinateToString(x, y + 1);
                graph_edges << map_coordinates_to_vertices[current_coordenate] << " " 
                            << map_coordinates_to_vertices[neighbor_right] << '\n';
            }
        }
    }

    graph_edges.close();
}

void generateColumns(std::ofstream& file) {
    file << "iteration,algorithm,initial state,search objective,generated vertices amount,visited vertices amount,path cost,path coordinates\n";
}

void part1(const std::string& source, const std::string& destination, Graph& graph, size_t cenary) {
    std::ofstream generalCSV("part1-results-cost" + std::to_string(cenary) + ".csv");
    
    if (generalCSV.is_open()) {
        generateColumns(generalCSV);

        for (size_t i {0}; i < 50; ++i) {
 
            generalCSV << i << ",BFS," << graph.breadthFirstSearch(source, destination, cenary) << '\n';

            generalCSV << i << ",DFS," << graph.depthFirstSearch(source, destination, cenary) << '\n';

            generalCSV << i << ",UCS," << graph.uniformCostSearch(source, destination, cenary) << '\n';
        }
    } 
    
    else {
        std::cerr << "Error opening general CSV file!" << std::endl;
    }
}

void part2(const std::string& source, const std::string& destination, Graph& graph, size_t cenary) {
    std::ofstream generalCSV("part2-results-cost" + std::to_string(cenary) + ".csv");
    
    if (generalCSV.is_open()) {
        generateColumns(generalCSV);

        for (size_t i {0}; i < 50; ++i) {
            generalCSV << i << ",UCS," << graph.uniformCostSearch(source, destination, cenary) << '\n';
            generalCSV << i << ",A*," << graph.AStar(source, destination, cenary).first << '\n';
        }
    } 
    
    else {
        std::cerr << "Error opening general CSV file!" << std::endl;
    }
}

void part3(const std::string& source, const std::string& destination, Graph& graph, size_t cenary) {
    std::ofstream generalCSV("part3-results-cost" + std::to_string(cenary) + ".csv");
    
    if (generalCSV.is_open()) {
        generateColumns(generalCSV);

        for (size_t i {0}; i < 50; ++i) {
            generalCSV << i << ",GS," << graph.greedySearch(source, destination, cenary) << '\n';
            generalCSV << i << ",A*," << graph.AStar(source, destination, cenary).first << '\n';
        }
    } 
    
    else {
        std::cerr << "Error opening general CSV file!" << std::endl;
    }
}

void part4(const std::string& source, const std::string& destination, Graph& graph, size_t cenary) {
    std::ofstream generalCSV("part4-results-cost" + std::to_string(cenary) + ".csv");
    
    if (generalCSV.is_open()) {
        generateColumns(generalCSV);

        for (size_t i {0}; i < 20; ++i) {
            generalCSV << i << ",BFS," << graph.breadthFirstSearch(source, destination, cenary, true) << '\n';

            generalCSV << i << ",DFS," << graph.depthFirstSearch(source, destination, cenary, true) << '\n';            
        }
    } 
    
    else {
        std::cerr << "Error opening general CSV file!" << std::endl;
    }
}

void part5(const std::string& source, const std::string& destination, const std::vector<std::string>& drugstores, Graph& graph, size_t cenary) {
	std::string nearest_drugstore;
	size_t shortest_path = std::numeric_limits<size_t>::max(); 

	for (const auto& drugstore : drugstores) {
		auto result = graph.AStar(source, drugstore, cenary);
		
		if (result.second < shortest_path) {
		    shortest_path = result.second;
		    nearest_drugstore = drugstore;
		}
	}
	
    std::ofstream generalCSV("part5-results-cost" + std::to_string(cenary) + ".csv");
    
    if (generalCSV.is_open()) {
        generateColumns(generalCSV);

        for (size_t i {0}; i < 25; ++i) {
            generalCSV << i << ",A*," << graph.AStar(nearest_drugstore, destination, cenary).first << '\n';          
        }
    } 
    
    else {
        std::cerr << "Error opening general CSV file!" << std::endl;
    }
}


int main(void) {

	// instanciações necessárias para gerar o número pseudo-aleatório
	
	std::random_device random_number;
	std::mt19937 seed(random_number());
	std::uniform_int_distribution<size_t> gap(0, COORDINATES_AMOUNT - 1);
	
	buildGraphFile();
	
    Graph graph("graph.txt", cartesian_plane, map_vertices_to_coordinates, map_coordinates_to_vertices, false);
	
	std::string source {};
	std::string destination {};
	std::unordered_set<std::string> set;
	std::vector<std::string> drugstores;
	std::string drugstore {};
	
	
    for (size_t i {0}; i < 4; ++i) {
    	drugstore = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
    	
 		while (set.count(drugstore)) {
			drugstore = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
			drugstores.push_back(drugstore);
		}
		
    	set.emplace(drugstore);
    	drugstores.push_back(drugstore);
    }
    
    source = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
	destination = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
	
	for (size_t i {1}; i < 5; ++i) {
		part1(source, destination, graph, i);
		part2(source, destination, graph, i);
		part3(source, destination, graph, i);
		part4(source, destination, graph, i);
		part5(source, destination, drugstores, graph, i);
	}
	
    return 0;
}

