#ifndef NODE_HPP
#define NODE_HPP

#include <iostream>
#include <memory>
#include <functional>

struct Node {
    size_t vertex;
    size_t value;
    size_t heuristic_cost;
    size_t path_cost;
    size_t height;
    std::shared_ptr<Node> parent; // Alterado para shared_ptr

    // Construtor genérico
    Node(size_t vertex, size_t value = 0, size_t heuristic_cost = 0, size_t path_cost = 0, std::shared_ptr<Node> parent = nullptr, size_t height = 0) :
         vertex(vertex), value(value), heuristic_cost(heuristic_cost),
         path_cost(path_cost), height(height), parent(parent) {}

    Node() = default;
    
    // Destruidor padrão, sem necessidade de implementação manual para shared_ptr
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
        this->parent = other.parent; // Copiando shared_ptr
        this->height = other.height;

        return *this;
    }

    bool operator == (const Node& node) const {
        return (vertex == node.vertex) && (value == node.value);
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

        // Combine os hashes de vertex, value e drugstore_visited de forma mais robusta
        return h1 ^ (h2 << 1);
    }
};

#endif

