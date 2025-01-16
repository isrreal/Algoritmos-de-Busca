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

inline std::string coordinateToString(size_t x, size_t y) {
    return "(" + std::to_string(x) + " " + std::to_string(y) + ")";
}

inline void buildGraphFile() {
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

inline std::vector<std::string> generateDrugstores(Graph& graph) {
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
    
    return drugstores;
}

inline std::string getNearestDrugstore(const std::string& source, std::vector<std::string> drugstores, Graph& graph, size_t cenary, size_t heuristic) {
	std::string nearest_drugstore {};
	size_t shortest_path { std::numeric_limits<size_t>::max() }; 
	
	for (const auto& drugstore : drugstores) {
		auto result = graph.AStar(source, drugstore, cenary, heuristic);
		
		if (result.second < shortest_path) {
			shortest_path = result.second;
			nearest_drugstore = drugstore;
		}
	}
	
	return nearest_drugstore;
}

#endif
