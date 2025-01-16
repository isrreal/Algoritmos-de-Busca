#include <iostream>
#include <unordered_map>
#include <fstream>
#include <string>
#include <vector>
#include "Graph.hpp"
#include "util_functions.hpp"

void part1(Graph& graph) {
    std::string source {};
    std::string destination {};

    // Abre os arquivos CSV fora dos loops
    std::ofstream generalCSV("part1-results.csv");

    if (generalCSV.is_open()) {
        // Escreve os cabeçalhos no arquivo CSV
        generalCSV << "Iteration,Algorithm,Cenary,Initial State,Objective,Generated vertices amount,Visited vertices amount,Cost,Path\n";

        for (size_t i {0}; i < 50; ++i) { 
            source = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
            destination = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
            
            for (size_t cenary {1}; cenary <= 4; ++cenary) {  
                generalCSV << i << ",BFS," << cenary << "," 
                           << graph.breadthFirstSearch(source, destination, cenary) << '\n';

                generalCSV << i << ",DFS," << cenary << "," 
                           << graph.depthFirstSearch(source, destination, cenary) << '\n';

                generalCSV << i << ",UCS," << cenary << "," 
                           << graph.uniformCostSearch(source, destination, cenary) << '\n';
            }
        }
    } 
    else {
        std::cerr << "Error opening general CSV file!" << std::endl;
    }
}

void part2(Graph& graph) {
    std::string source {};
    std::string destination {};

    // Abre os arquivos CSV fora dos loops
    std::ofstream generalCSV1("part2-results-ucs.csv");
    std::ofstream generalCSV2("part2-results-astar.csv");

    if (generalCSV1.is_open() && generalCSV2.is_open()) {
        // Escreve os cabeçalhos no arquivo CSV
        generalCSV1 << "Iteration,Algorithm,Cenary,Initial State,Objective,Generated vertices amount,Visited vertices amount,Cost,Path\n";
        generalCSV2 << "Iteration,Algorithm,Cenary,Initial State,Objective,Generated vertices amount,Visited vertices amount,Cost,Path\n";

        for (size_t i {0}; i < 50; ++i) {
            source = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
            destination = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
            
            for (size_t cenary {1}; cenary <= 4; ++cenary) {  
                generalCSV1 << i << "," << cenary << "," 
                            << graph.uniformCostSearch(source, destination, cenary) << '\n';
            }

            for (size_t cenary {1}; cenary <= 4; ++cenary) {  
                for (size_t heuristic {1}; heuristic <= 2; ++heuristic) {  
                    generalCSV2 << i << "," << cenary << "," << heuristic << ","
                                << graph.AStar(source, destination, cenary, heuristic).first << '\n';
                }
            }
        }
    } 
    else {
        std::cerr << "Error opening general CSV files!" << std::endl;
    }
}

void part3(Graph& graph) {
    std::string source {};
    std::string destination {};

    std::ofstream generalCSV("part3-results.csv");

    if (generalCSV.is_open()) {
        generalCSV << "Iteration,Algorithm,Cenary,Heuristic,Initial State,Objective,Generated vertices amount,Visited vertices amount,Cost,Path\n";

        for (size_t i {0}; i < 50; ++i) {
            source = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
            destination = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
            
            for (size_t cenary {1}; cenary <= 4; ++cenary) {  
                for (size_t heuristic {1}; heuristic <= 2; ++heuristic) {  
                    generalCSV << i << ",GS," << cenary << "," << heuristic << ","
                               << graph.greedySearch(source, destination, cenary, heuristic) << '\n';

                    generalCSV << i << ",A*," << cenary << "," << heuristic << ","
                               << graph.AStar(source, destination, cenary, heuristic).first << '\n';
                }
            }
        }
    } 
    
    else {
        std::cerr << "Error opening general CSV file!" << std::endl;
    }
}

void part4(Graph& graph) {
    std::string source {};
    std::string destination {};

    std::ofstream generalCSV("part4-results.csv");

    if (generalCSV.is_open()) {
        generalCSV << "Iteration,Algorithm,Cenary,Heuristic,Initial State,Objective,Generated vertices amount,Visited vertices amount,Cost,Path\n";

        for (size_t i {0}; i < 20; ++i) {
            source = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
            destination = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));

            for (size_t cenary {1}; cenary <= 4; ++cenary) {
                for (size_t j {0}; j < 10; ++j) {
                    generalCSV << i << ",BFS," << cenary << ","
                               << graph.breadthFirstSearch(source, destination, cenary, true) << '\n';

                    generalCSV << i << ",DFS," << cenary << ","
                               << graph.depthFirstSearch(source, destination, cenary, true) << '\n';
                }
            }
        }
    } 
    
    else {
        std::cerr << "Error opening general CSV file!" << std::endl;
    }
}

void part5(Graph& graph) {
    std::string source {};
    std::string destination {};
    std::vector<std::string> drugstores;
    std::string nearest_drugstore {};
    
    std::ofstream generalCSV("part5-results.csv");

    if (generalCSV.is_open()) {
        generalCSV << "Iteration,Algorithm,Cenary,Heuristic,Initial State,Objective,Generated vertices amount,Visited vertices amount,Cost,Path\n";

        for (size_t i {0}; i < 25; ++i) {
            source = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
            destination = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
            
            for (size_t cenary {1}; cenary <= 4; ++cenary) {
                for (size_t heuristic {1}; heuristic <= 2; ++heuristic) {
                    drugstores = generateDrugstores(graph);
                    nearest_drugstore = getNearestDrugstore(source, drugstores, graph, cenary, heuristic);

                    generalCSV << i << ",A*," << cenary << "," << heuristic << ","
                               << graph.AStar(nearest_drugstore, destination, cenary, heuristic).first << '\n';
                }
            }
        }
    } 
    
    else {
        std::cerr << "Error opening general CSV file!" << std::endl;
    }
}


int main(void) {	
	buildGraphFile();
	
    Graph graph("graph.txt", cartesian_plane, map_vertices_to_coordinates, map_coordinates_to_vertices, false);
	part1(graph);
	part2(graph);
	part3(graph);
	part4(graph);
    part5(graph);
    
    return 0;
}

