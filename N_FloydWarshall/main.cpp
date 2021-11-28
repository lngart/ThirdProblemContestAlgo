#include <iostream>
#include <limits>
#include <vector>

template <typename T>
class IGraph {
public:
    using Vertex = T;

    struct Edge {
        Edge(Vertex end, T weight_edge) : end_edge{end}, weight{weight_edge} {
        }

        Vertex end_edge;
        T weight;

        bool operator>(const Edge &other) const {
            return weight > other.weight;
        }
    };

    virtual const std::vector<Edge> &GetNeighbors(Vertex v) const = 0;
    virtual void AddEdge(Vertex from, Vertex to, T weight) = 0;
    virtual size_t Size() const = 0;
    virtual ~IGraph() = default;
};

template <typename T>
class ListGraph : public IGraph<T> {
    using typename IGraph<T>::Vertex;
    using typename IGraph<T>::Edge;

    std::vector<std::vector<Edge>> vertices_;

public:
    explicit ListGraph(size_t num_vertices) : vertices_(num_vertices) {
    }

    const std::vector<Edge> &GetNeighbors(Vertex v) const override {
        return vertices_[v];
    }

    void AddEdge(Vertex from, Vertex to, T weight) override {
        vertices_[from].emplace_back(to, weight);
    }

    size_t Size() const override {
        return vertices_.size();
    }
};

enum class Color { WHITE, BLACK };

template <typename T>
void FindMinDistancesBetweenAllVertices(IGraph<T> &graph, std::vector<std::vector<T>> &distance) {
    for (size_t vertex{0}; vertex < graph.Size(); ++vertex) {
        for (const auto &neighbor : graph.GetNeighbors(vertex)) {
            auto adjacent_vertex = neighbor.end_edge, current_weight = neighbor.weight;
            distance[vertex][adjacent_vertex] = current_weight;
        }
    }

    for (size_t i{0}; i < graph.Size(); ++i) {
        for (size_t j{0}; j < graph.Size(); ++j) {
            for (size_t k{0}; k < graph.Size(); ++k) {
                distance[j][k] = std::min(distance[j][k], distance[j][i] + distance[i][k]);
            }
        }
    }
}

template <typename T>
void PrintMinDistancesBetweenAllVertices(IGraph<T> &graph) {
    std::vector<std::vector<T>> distance(graph.Size(), std::vector<T>(graph.Size()));

    FindMinDistancesBetweenAllVertices(graph, distance);

    for (size_t i{0}; i < graph.Size(); ++i) {
        for (size_t j{0}; j < graph.Size(); ++j) {
            std::cout << distance[i][j] << " ";
        }
        std::cout << '\n';
    }
}

void Initialization(size_t &num_vertices) {
    std::cin >> num_vertices;
}

template <typename T>
void AddingEdge(IGraph<T> &graph, size_t num_vertices) {
    for (size_t from{0}; from < num_vertices; ++from) {
        for (size_t to{0}; to < num_vertices; ++to) {
            T weight{0};
            std::cin >> weight;
            if (from != to) {
                graph.AddEdge(from, to, weight);  // We use numbering from 0
            }
        }
    }
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);

    size_t num_vertices{0};

    Initialization(num_vertices);

    ListGraph<int> graph(num_vertices);

    AddingEdge(graph, num_vertices);

    PrintMinDistancesBetweenAllVertices(graph);

    return 0;
}
