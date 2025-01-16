#ifndef NODE_HPP
#define NODE_HPP

#include <functional>

struct Node {
    size_t vertex;
    size_t value;
    size_t heuristic_cost;
    size_t path_cost;
    size_t height;
    bool drugstore_visited;

    // Construtor genérico
    Node(size_t vertex, size_t value = 0, size_t heuristic_cost = 0, size_t path_cost = 0, size_t height = 0, bool drugstore_visited = false) :
         vertex(vertex), value(value), heuristic_cost(heuristic_cost),
         path_cost(path_cost), height(height), drugstore_visited(drugstore_visited) {}

    Node() = default;
    
    ~Node() = default;

    // Sobrecarga do operador =
    Node& operator=(const Node& other) {
        if (this == &other) {
            return *this; 
        }

        this->vertex = other.vertex;
        this->value = other.value;
        this->heuristic_cost = other.heuristic_cost;
        this->path_cost = other.path_cost;
        this->height = other.height;
        this->drugstore_visited = other.drugstore_visited;

        return *this;
    }

    bool operator == (const Node& node) const {
        return (vertex == node.vertex) && (value == node.value) && (drugstore_visited == node.drugstore_visited);
    }

    bool operator < (const Node& node) const {
        return (value < node.value);
    }

    bool operator > (const Node& node) const {
        return (value > node.value);
    }
};

// Estrutura de hash personalizada para o Node
struct NodeHash {
    size_t operator()(const Node& node) const {
        size_t h1 = std::hash<size_t>{}(node.vertex);
        size_t h2 = std::hash<size_t>{}(node.value);
        size_t h3 = std::hash<bool>{}(node.drugstore_visited);

        // Combine os hashes de vertex e value (podemos usar XOR ou soma)
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

#endif

