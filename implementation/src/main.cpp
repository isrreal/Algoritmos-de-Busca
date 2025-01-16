#include <iostream>
#include <unordered_map>
#include <fstream>
#include <string>
#include <vector>
#include "Graph.hpp"
#include "util_functions.hpp"

/**
 * @brief Executes experiments comparing BFS, DFS, and UCS algorithms under various scenarios.
 * 
 * This function performs an experimental evaluation of three search algorithms: 
 * Breadth-First Search (BFS), Depth-First Search (DFS), and Uniform Cost Search (UCS). 
 * The goal is to analyze the performance and behavior of these algorithms under different 
 * cost functions and scenarios, storing the results in a CSV file.
 * 
 * @param graph Reference to the graph object representing the search space (city map).
 * 
 * @details 
 * For each of the 50 iterations, two random coordinates are generated: one as the 
 * initial state (source) and the other as the objective (destination). Each algorithm 
 * is executed using four different cost functions (cenary 1 to 4). The results, including 
 * the algorithm used, the cost function, the generated vertices amount, visited vertices 
 * amount, path cost, and the path itself, are recorded in a CSV file.
 * 
 * The steps for each iteration are:
 * 1. Generate random coordinates `(x1, y1)` for the initial state and `(x2, y2)` for the objective.
 * 2. For each cost function (cenary 1 to 4):
 *    - Execute BFS and log its results.
 *    - Execute DFS and log its results.
 *    - Execute UCS and log its results.
 * 
 * Results are saved in the CSV file `part1-results.csv` with the following headers:
 * Iteration,Algorithm,Cenary,Initial State,Objective,Generated vertices amount,Visited vertices amount,Cost,Path.
 */
 
void part1(Graph& graph) {
    std::string source {};
    std::string destination {};

    std::ofstream generalCSV("part1-results.csv");

    if (generalCSV.is_open()) {
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

/**
 * @brief Executes experiments comparing UCS and A* algorithms under various scenarios and heuristics.
 * 
 * This function performs an experimental evaluation of two search algorithms: 
 * Uniform Cost Search (UCS) and A* (A-Star). The goal is to analyze their 
 * performance and behavior under different cost functions and heuristics, 
 * storing the results in two separate CSV files.
 * 
 * @param graph Reference to the graph object representing the search space (city map).
 * 
 * @details 
 * For each of the 50 iterations, two random coordinates are generated: one as the 
 * initial state (source) and the other as the objective (destination). UCS is executed 
 * using four different cost functions (cenary 1 to 4), while A* is executed using 
 * all combinations of the four cost functions and two heuristic functions (H1 and H2). 
 * The results are logged in the following CSV files:
 * 
 * - `part2-results-ucs.csv`: Contains results from the UCS algorithm.
 * - `part2-results-astar.csv`: Contains results from the A* algorithm.
 * 
 * Each CSV file contains the following headers:
 * `Iteration,Algorithm,Cenary,Initial State,Objective,Generated vertices amount,Visited vertices amount,Cost,Path`.
 * 
 */
 
void part2(Graph& graph) {
    std::string source {};
    std::string destination {};

    std::ofstream generalCSV1("part2-results-ucs.csv");
    std::ofstream generalCSV2("part2-results-astar.csv");

    if (generalCSV1.is_open() && generalCSV2.is_open()) {
        generalCSV1 << "Iteration,Cenary,Initial State,Objective,Generated vertices amount,Visited vertices amount,Cost,Path\n";
        generalCSV2 << "Iteration,Heuristic,Cenary,Initial State,Objective,Generated vertices amount,Visited vertices amount,Cost,Path\n";

        for (size_t i {0}; i < 50; ++i) {
            source = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
            destination = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
            
            for (size_t cenary {1}; cenary <= 4; ++cenary) {  
                generalCSV1 << i << ',' << cenary << ',' 
                            << graph.uniformCostSearch(source, destination, cenary) << '\n';
            }

            for (size_t cenary {1}; cenary <= 4; ++cenary) {  
                for (size_t heuristic {1}; heuristic <= 2; ++heuristic) {  
                    generalCSV2 << i << ',' << heuristic << ',' << cenary << ','
                                << graph.AStar(source, destination, cenary, heuristic, {}).first << '\n';
                }
            }
        }
    } 
    else {
        std::cerr << "Error opening general CSV files!" << std::endl;
    }
}

/**
 * @brief Executes experiments comparing Greedy Search and A* algorithms under various scenarios and heuristics.
 * 
 * This function performs an experimental evaluation of two search algorithms: 
 * Greedy Search (GS) and A* (A-Star). The goal is to compare their performance and 
 * behavior using different cost functions and heuristics, logging the results in a CSV file.
 * 
 * @param graph Reference to the graph object representing the search space (city map).
 * 
 * @details
 * For each of the 50 iterations, two random coordinates are generated: one as the 
 * initial state (source) and the other as the objective (destination). 
 * Both GS and A* algorithms are executed using all combinations of:
 * - Four cost functions (`cenary` 1 to 4)
 * - Two heuristic functions (`H1` and `H2`)
 * 
 * The results are logged in the CSV file `part3-results.csv` with the following format:
 * `Iteration,Algorithm,Cenary,Heuristic,Initial State,Objective,Generated vertices amount,Visited vertices amount,Cost,Path`.
 */

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
                               << graph.AStar(source, destination, cenary, heuristic, {}).first << '\n';
                }
            }
        }
    } 
    
    else {
        std::cerr << "Error opening general CSV file!" << std::endl;
    }
}

/**
 * @brief Executes experiments comparing Breadth-First Search (BFS) and Depth-First Search (DFS) 
 * with randomization of neighborhood generation.
 * 
 * This function evaluates the behavior of BFS and DFS under randomized neighbor generation. 
 * It introduces random shuffling of the order in which neighbors are generated during the search, 
 * allowing for variations in the paths found by these algorithms across multiple executions.
 * 
 * @param graph Reference to the graph object representing the search space.
 * 
 * @details
 * For each of the 20 iterations, two random coordinates are generated: one as the 
 * initial state (source) and the other as the objective (destination). 
 * BFS and DFS algorithms are executed 10 times for each of the four cost functions (`cenary` 1 to 4).
 * 
 * The results are saved in the CSV file `part4-results.csv` with the following format:
 * `Iteration,Algorithm,Cenary,Heuristic,Initial State,Objective,Generated vertices amount,Visited vertices amount,Cost,Path`.
 */
 
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

/**
 * @brief Executes experiments using the A* algorithm to solve the pathfinding problem 
 * with constraints: passing through a drugstore before reaching the destination.
 * 
 * This function simulates a scenario where an agent must return home from work but needs 
 * to stop at one of several drugstores along the way. The goal is to find the least-cost path 
 * that satisfies these constraints, using randomized coordinates and cost functions for variation.
 * 
 * @param graph Reference to the graph object representing the city or search space.
 * 
 * @details
 * The experiment evaluates the behavior of the A* algorithm in a constrained search problem, 
 * where randomization and multiple cost functions (`f1`, `f2`, `f3`, `f4`) combined with 
 * heuristics (`H1` and `H2`) are applied to simulate real-world unpredictability.
 * 
 * The results are saved in the CSV file `part5-results.csv` with the following format:
 * `Iteration,Algorithm,Cenary,Heuristic,Initial State,Objective,Generated vertices amount,Visited vertices amount,Cost,Path`.
 */
 
void part5(Graph& graph) {
    std::string source {};
    std::string destination {};
    std::unordered_set<std::string> drugstores;
    
    std::ofstream generalCSV("part5-results.csv");

    if (generalCSV.is_open()) {
        generalCSV << "Iteration,Algorithm,Cenary,Heuristic,Initial State,Objective,Generated vertices amount,Visited vertices amount,Cost,Path\n";

        for (size_t i {0}; i < 25; ++i) {
            source = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
            destination = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
            
            for (size_t cenary {1}; cenary <= 4; ++cenary) {
                for (size_t heuristic {1}; heuristic <= 2; ++heuristic) {
                    drugstores = generateDrugstores(graph);

                    generalCSV << i << ",A*," << cenary << "," << heuristic << ","
                               << graph.AStar(source, destination, cenary, heuristic, drugstores).first << '\n';
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

