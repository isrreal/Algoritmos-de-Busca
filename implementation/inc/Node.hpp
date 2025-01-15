#ifndef NODE_HPP
#define NODE_HPP

#include <cstddef> // Para size_t

// Definições para os tipos de busca
struct AStarSearch {};
struct UniformCostSearch {};
struct GreedySearch {};

// Template base para Node
template <typename T>
struct Node {
    size_t vertex;
    size_t acumulated_path_with_heuristic;
    size_t acumulated_path_cost;
    Node<T>* parent;
    size_t height;

    // Construtor genérico
    Node(size_t vertex, size_t acumulated_path_with_heuristic, size_t acumulated_path_cost, Node<T>* parent = nullptr, size_t height = 1)
        : vertex(vertex), acumulated_path_with_heuristic(acumulated_path_with_heuristic), acumulated_path_cost(acumulated_path_cost), parent(parent), height(height) {}

    Node() = default;

    ~Node() = default;
};

// Especialização para AStarSearch
template <>
struct Node<AStarSearch> {
    size_t vertex;
    size_t acumulated_path_with_heuristic;
    size_t acumulated_path_cost;
    Node<AStarSearch>* parent;
    size_t height;

    Node(size_t vertex, size_t acumulated_path_with_heuristic, size_t acumulated_path_cost, Node<AStarSearch>* parent = nullptr, size_t height = 1)
        : vertex(vertex), acumulated_path_with_heuristic(acumulated_path_with_heuristic), acumulated_path_cost(acumulated_path_cost), parent(parent), height(height) {}

    Node() = default;

    bool operator == (const Node<AStarSearch>& node) const {
        return (vertex == node.vertex) && (acumulated_path_with_heuristic == node.acumulated_path_with_heuristic);
    }

    bool operator < (const Node<AStarSearch>& node) const {
        return (acumulated_path_cost + acumulated_path_with_heuristic) < (node.acumulated_path_cost + node.acumulated_path_with_heuristic); // g(n) + h(n)
    }

    bool operator > (const Node<AStarSearch>& node) const {
        return (acumulated_path_cost + acumulated_path_with_heuristic) > (node.acumulated_path_cost + node.acumulated_path_with_heuristic); // g(n) + h(n)
    }

    ~Node() = default;
};

// Especialização para Busca Uniforme
template <>
struct Node<UniformCostSearch> {
    size_t vertex;
    size_t acumulated_path_with_heuristic;
    size_t acumulated_path_cost;
    Node<UniformCostSearch>* parent;
    size_t height;

    Node(size_t vertex, size_t acumulated_path_with_heuristic, size_t acumulated_path_cost, Node<UniformCostSearch>* parent = nullptr, size_t height = 1)
        : vertex(vertex), acumulated_path_with_heuristic(acumulated_path_with_heuristic), acumulated_path_cost(acumulated_path_cost), parent(parent), height(height) {}

    Node() = default;

    bool operator==(const Node<UniformCostSearch>& node) const {
        return (vertex == node.acumulated_path_cost) && (acumulated_path_with_heuristic == node.acumulated_path_cost);
    }

    bool operator<(const Node<UniformCostSearch>& node) const {
        return acumulated_path_cost < node.acumulated_path_cost; // Apenas g(n)
    }

    bool operator>(const Node<UniformCostSearch>& node) const {
        return acumulated_path_cost > node.acumulated_path_cost; // Apenas g(n)
    }

    ~Node() = default;
};

// Especialização para Busca Gulosa
template <>
struct Node<GreedySearch> {
    size_t vertex;
    size_t acumulated_path_with_heuristic;
    size_t acumulated_path_cost;
    Node<GreedySearch>* parent;
    size_t height;

    Node(size_t vertex, size_t acumulated_path_with_heuristic, size_t acumulated_path_cost, Node<GreedySearch>* parent = nullptr, size_t height = 1)
        : vertex(vertex), acumulated_path_with_heuristic(acumulated_path_with_heuristic), acumulated_path_cost(acumulated_path_cost), parent(parent), height(height) {}

    Node() = default;

    bool operator==(const Node<GreedySearch>& node) const {
        return (vertex == node.vertex) && (acumulated_path_with_heuristic == node.acumulated_path_with_heuristic);
    }

    bool operator<(const Node<GreedySearch>& node) const {
        return acumulated_path_with_heuristic < node.acumulated_path_with_heuristic; // Apenas h(n)
    }

    bool operator>(const Node<GreedySearch>& node) const {
        return acumulated_path_with_heuristic > node.acumulated_path_with_heuristic; // Apenas h(n)
    }

    ~Node() = default;
};

#endif

