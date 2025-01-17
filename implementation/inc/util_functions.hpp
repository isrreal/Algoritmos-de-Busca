#ifndef UTIL_FUNCTIONS_HPP
#define UTIL_FUNCTIONS_HPP

#include <vector>
#include <random>
#include <unordered_set>

std::unordered_map<std::string, size_t> map_coordinates_to_vertices;
std::unordered_map<size_t, std::string> map_vertices_to_coordinates;
std::vector<std::vector<std::string>> cartesian_plane;
constexpr size_t COORDINATES_AMOUNT = 31 * 31;

std::random_device random_number;
std::mt19937 seed(random_number());
std::uniform_int_distribution<size_t> gap(0, COORDINATES_AMOUNT - 1);

inline void buildMatrix() {
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

inline std::string coordinatesToString(size_t x, size_t y) {
    return "(" + std::to_string(x) + " " + std::to_string(y) + ")";
}

inline void buildGraphFile() {
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

            current_coordenate = coordinatesToString(x, y);
            
            if ((x + 1) < 31) {
                neighbor_down = coordinatesToString(x + 1, y);
                graph_edges << map_coordinates_to_vertices[current_coordenate] << " " 
                            << map_coordinates_to_vertices[neighbor_down] << '\n';
            }
            
            if ((y + 1) < 31) {
                neighbor_right = coordinatesToString(x, y + 1);
                graph_edges << map_coordinates_to_vertices[current_coordenate] << " " 
                            << map_coordinates_to_vertices[neighbor_right] << '\n';
            }
        }
    }

	graph_edges.close();
}

/**
 * @brief Generates a list of random drugstore locations in the graph.
 * 
 * This function is responsible for creating a list of four random and unique 
 * coordinates within the graph to represent the locations of drugstores. Each 
 * coordinate is converted to a string format and added to a set to ensure no duplicates.
 * 
 * @param graph Reference to the graph object representing the city.
 * @return A vector containing string representations of the coordinates of the drugstores.
 * 
 * @note Drugstore locations are chosen using a random seed and are converted to 
 * a string format for easier processing.
 */
 
inline std::unordered_set<std::string> generateDrugstores(Graph& graph) {
	std::unordered_set<std::string> drugstores;
	std::string drugstore {};
	
	for (size_t i {0}; i < 4; ++i) {
    	drugstore = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
    	
 		while (drugstores.count(drugstore)) {
			drugstore = graph.coordinatesToString(graph.vertexToCoordinate(gap(seed)));
		}
		
    	drugstores.emplace(drugstore);
    }
    
    return drugstores;
}

inline size_t getRandomNumber() { return gap(seed); }

#endif
