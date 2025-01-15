#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <iostream>
#include <random>
#include <algorithm>
#include <list>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <sstream>
#include <fstream>
#include "Node.hpp"
#include <memory>

class Graph {
private:
    size_t order;
    size_t size;
    bool isDirected; 
    size_t delta; 
    size_t Delta;
    
    std::unordered_map<size_t, std::list<size_t>> adjacency_list;
	std::unordered_map<size_t, std::string> map_vertices_to_coordinates;
	std::unordered_map<std::string, size_t> map_coordinates_to_vertices;
	
    std::vector<std::vector<std::string>> cartesian_plane;
    
    void addVertex(size_t source);
    void addEdge(size_t source, size_t destination);
    Graph readGraph(const std::string& filename);
    size_t computeMaxVertexDegree();
    size_t computeMinVertexDegree();
    
public:	
	Graph(size_t order, bool isDirected, float probabilityOfEdge);	
    Graph(const std::string& filename, 
		  const std::vector<std::vector<std::string>>& cartesian_plane, 
		  const std::unordered_map<size_t, std::string>& map_vertices_to_coordinates,
		  const std::unordered_map<std::string, size_t>& map_coordinates_to_vertices,
		  bool isDirected);

    Graph() = default;
    ~Graph() = default;
  
    size_t getSize() const;
    size_t getOrder() const;
    size_t getVertexDegree(size_t vertex) const;
    size_t getMaxDegree() const;
    size_t getMinDegree() const;
   	const std::vector<std::vector<std::string>>& getCartesianPlane() const;
   	
	size_t coordinateToVertex(const std::string& coordinate) const;
	
	std::pair<size_t, size_t> vertexToCoordinate(size_t vertex) const;
    
    double heuristic1(int x1, int y1, int x2, int y2);
	size_t heuristic2(int x1, int y1, int x2, int y2);

	double euclideanDistance(int x1, int y1, int x2, int y2);
	size_t manhatannDistance(int x1, int y1, int x2, int y2);
	
	constexpr size_t c1();
	constexpr size_t c2();
	size_t c3(size_t steps);
	size_t c4(size_t steps);
	
	std::string getStats(const std::string& initial_state,
				   const std::string& search_objective,
 				   const std::vector<std::string>& path,
 			       size_t path_cost,
 				   size_t generated_vertices_amount,
 				   size_t visited_vertices_amount);
 				   
   	size_t pathCompute(std::stack<std::pair<size_t, size_t>>& last_visited,
                         size_t path_cost, size_t steps_root_to_objective_amount, size_t cenary, size_t temp);
                     
 	size_t neighborCost(size_t u, size_t v, size_t steps_root_to_objective_amount, size_t cenary);
 	
	std::vector<std::string> buildPath(size_t source, size_t destination, const std::vector<int>& predecessor);
 	
   	std::string coordinatesToString(const std::pair<size_t, size_t>& coordinates);

	const std::unordered_map<size_t, std::list<size_t>>& getAdjacencyList() const;
    
    const std::list<size_t>& getAdjacencyList(size_t vertex) const;
    
    bool edgeExists(size_t u, size_t v) const;
    
    bool vertexExists(size_t vertex) const;	

    std::string breadthFirstSearch(const std::string& u, const std::string& v, size_t cenary);  
    
    std::string depthFirstSearch(const std::string& u, const std::string& v, size_t cenary);  
    
    std::string uniformCostSearch(const std::string& u, const std::string& v, size_t cenary);
    
    std::string AStar(const std::string& u, const std::string& v, size_t cenary);
    
	std::string greedySearch(const std::string& u, const std::string& v, size_t cenary);
    
    friend std::ostream& operator<< (std::ostream& os, const Graph& graph);
};

#endif
