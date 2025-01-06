#include <iostream>
#include <unordered_map>
#include <fstream>
#include <string>
#include <vector>
#include "Graph.hpp"

std::unordered_map<std::string, size_t> matrix;
std::vector<std::vector<std::string>> cartesian_plane;

void buildMatrix() {
    size_t counter {0};
    cartesian_plane.resize(31, std::vector<std::string>(31));

    for (size_t x {0}; x < 31; ++x) {
        for (size_t y {0}; y < 31; ++y) {
            std::string coord = "(" + std::to_string(x) + "," + std::to_string(y) + ")";
            matrix[coord] = counter++;
            cartesian_plane[x][y] = coord;
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

    for (size_t x {0}; x < 31; ++x) {
        for (size_t y {0}; y < 31; ++y) {
            std::string current = "(" + std::to_string(x) + "," + std::to_string(y) + ")";
            
            if (x + 1 < 31) {
                std::string neighbor_down = "(" + std::to_string(x + 1) + "," + std::to_string(y) + ")";
                graph_edges << matrix[current] << " " << matrix[neighbor_down] << '\n';
            }
            
            if (y + 1 < 31) {
                std::string neighbor_right = "(" + std::to_string(x) + "," + std::to_string(y + 1) + ")";
                graph_edges << matrix[current] << " " << matrix[neighbor_right] << '\n';
            }
        }
    }

    graph_edges.close();
}

int main(void) {
    buildGraphFile();
    Graph graph("graph.txt", cartesian_plane, false);

    std::cout << graph.getCartesianPlane()[0][5] << '\n';
 
    return 0;
}

