#include <iostream>
#include <unordered_map>
#include <fstream>
#include <string>
#include <vector>
#include "Graph.hpp"

std::unordered_map<std::string, size_t> map_coordenates_to_vertices;
std::unordered_map<size_t, std::string> map_vertices_to_coordenates;
std::vector<std::vector<std::string>> cartesian_plane;

void buildMatrix() {
    size_t vertex {0};
    cartesian_plane.resize(31, std::vector<std::string>(31));
	std::string coordenate {};
	
    for (size_t x {0}; x < 31; ++x) {
        for (size_t y {0}; y < 31; ++y) {
            coordenate = "(" + std::to_string(x) + "," + std::to_string(y) + ")";
            map_vertices_to_coordenates[vertex] = coordenate;
            map_coordenates_to_vertices[coordenate] = vertex++;
            cartesian_plane[x][y] = coordenate;
        }
    }
}

void buildGraphFile() {
    buildMatrix();
    std::ofstream graph_edges("graph.txt");

    if (!graph_edges.is_open()) {
        std::cerr << "Erro ao abrir o arquivo 'graph.txt' para escrita." << std::endl;
        return;
    }
    
    std::string current_coordenate {};
    std::string neighbor_down {};
    std::string neighbor_right {};
    
    for (size_t x {0}; x < 31; ++x) {
        for (size_t y {0}; y < 31; ++y) {

            current_coordenate = "(" + std::to_string(x) + "," + std::to_string(y) + ")";
            
            if (x + 1 < 31) {
                neighbor_down = "(" + std::to_string(x + 1) + "," + std::to_string(y) + ")";
                
                graph_edges << map_coordenates_to_vertices[current_coordenate] << " " 
                			<< map_coordenates_to_vertices[neighbor_down] << '\n';
            }
            
            if (y + 1 < 31) {
                neighbor_right = "(" + std::to_string(x) + "," + std::to_string(y + 1) + ")";
                
                graph_edges << map_coordenates_to_vertices[current_coordenate] << " " 
                			<< map_coordenates_to_vertices[neighbor_right] << '\n';
            }
        }
    }

    graph_edges.close();
}

int main(void) {
    buildGraphFile();
    Graph graph("graph.txt", cartesian_plane,map_vertices_to_coordenates, map_coordenates_to_vertices, false);
 	std::cout << graph.f1(0, 6);
    return 0;
}

