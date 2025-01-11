#include <iostream>
#include <unordered_map>
#include <fstream>
#include <string>
#include <vector>
#include "Graph.hpp"

std::unordered_map<std::string, size_t> map_coordinates_to_vertices;
std::unordered_map<size_t, std::string> map_vertices_to_coordinates;
std::vector<std::vector<std::string>> cartesian_plane;

void buildMatrix() {
    size_t vertex {0};
    cartesian_plane.resize(31, std::vector<std::string>(31));
    std::string coordenate {};
    
    map_vertices_to_coordinates.reserve(961);
    map_coordinates_to_vertices.reserve(961);

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
            
            if (x + 1 < 31) {
                neighbor_down = coordinateToString(x + 1, y);
                graph_edges << map_coordinates_to_vertices[current_coordenate] << " " 
                            << map_coordinates_to_vertices[neighbor_down] << '\n';
            }
            
            if (y + 1 < 31) {
                neighbor_right = coordinateToString(x, y + 1);
                graph_edges << map_coordinates_to_vertices[current_coordenate] << " " 
                            << map_coordinates_to_vertices[neighbor_right] << '\n';
            }
        }
    }

    graph_edges.close();
}

int main(void) {
    buildGraphFile();
    Graph graph("graph.txt", cartesian_plane, map_vertices_to_coordinates, map_coordinates_to_vertices, false);
    
    std::cout << "********* BFS ********* \n";
    graph.breadthFirstSearch("(0,0)", "(5,30)", 4);
    
    std::cout << "\n\n ********* DFS ********* \n";
    graph.depthFirstSearch("(0,0)", "(5,30)", 4);
    
    std::cout << "\n\n ********* Uniform Search (Dikjstra) ********* \n";
    graph.uniformCostSearch("(0,0)", "(5,30)", 4);
    
    std::cout << "\n\n ********* A* ********* \n";
    graph.AStar("(0,0)", "(5,30)", 4);
    return 0;
}

