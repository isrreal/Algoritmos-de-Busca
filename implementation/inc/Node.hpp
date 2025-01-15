#ifndef NODE_HPP
#define NODE_HPP

struct Node {
    size_t vertex;
    size_t value;
    size_t heuristic_cost;
    size_t path_cost;
    size_t height;

    // Construtor genérico
    Node(size_t vertex, size_t value = 0, size_t heuristic_cost = 0, size_t path_cost = 0, size_t height = 0) :
         vertex(vertex), value(value), heuristic_cost(heuristic_cost),
         path_cost(path_cost), height(height) {}

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

#endif
