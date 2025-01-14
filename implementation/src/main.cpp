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
            coordenate = "(" + std::to_string(x) + "," + std::to_string(y) + ")";
            map_vertices_to_coordinates[vertex] = coordenate;
            map_coordinates_to_vertices[coordenate] = vertex++;
            cartesian_plane[x][y] = coordenate;
        }
    }
}

std::string coordinateToString(size_t x, size_t y) {
    return "(" + std::to_string(x) + "," + std::to_string(y) + ")";
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

int main(void) {

	// instanciações necessárias para gerar o número pseudo-aleatório
	
	std::random_device random_number;
	std::mt19937 seed(random_number());
	std::uniform_int_distribution<size_t> gap(0, COORDINATES_AMOUNT - 1);
	
	std::string x_source {};
	std::string y_source {};
	std::string x_destination {};
	std::string y_destination {};
		
    buildGraphFile();
    
    Graph graph("graph.txt", cartesian_plane, map_vertices_to_coordinates, map_coordinates_to_vertices, false);
    
    x_source = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
	y_source = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
	
	x_destination = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
	y_destination = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
    
    std::cout << "********* BFS ********* \n";
    graph.breadthFirstSearch(x_source, y_source, 3);
    
    std::cout << "\n\n ********* DFS ********* \n";
    graph.depthFirstSearch(x_source, y_source, 3);
    
    std::cout << "\n\n ********* Uniform Search (Dijkstra) ********* \n";
    graph.uniformCostSearch(x_source, y_source, 3);
    
    std::cout << "\n\n ********* A* ********* \n";
    graph.AStar(x_source, y_source, 3);
    
    std::cout << "\n\n ********* Greedy Search ********* \n";
    graph.greedySearch(x_source, y_source, 3);
   
    return 0;
}

