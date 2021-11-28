#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>
#include <numeric>

template <typename T>
class IGraph {
public:
    using Vertex = T;

    struct Edge {
        Edge(Vertex start, Vertex end, T weight) : start_edge{start}, end_edge{end}, weight_edge{weight} {
        }

        Vertex start_edge;
        Vertex end_edge;
        T weight_edge;
    };

    virtual const std::vector<Vertex> &GetNeighbors(Vertex v) const = 0;
    virtual void AddEdge(Vertex from, Vertex to, T weight) = 0;
    virtual const std::vector<Edge> &GetEdges() const = 0;
    virtual size_t Size() const = 0;
    virtual ~IGraph() = default;
};

template <typename T>
class ListGraph : public IGraph<T> {
    using typename IGraph<T>::Vertex;
    using typename IGraph<T>::Edge;

    std::vector<std::vector<Vertex>> vertices_;
    std::vector<Edge> edges_;

public:
    explicit ListGraph(size_t num_vertices) : vertices_(num_vertices) {
    }

    const std::vector<Vertex> &GetNeighbors(Vertex v) const override {
        return vertices_[v];
    }

    void AddEdge(Vertex from, Vertex to, T weight) override {
        vertices_[from].emplace_back(to);
        edges_.emplace_back(from, to, weight);
    }

    const std::vector<Edge> &GetEdges() const override {
        return edges_;
    }

    size_t Size() const override {
        return vertices_.size();
    }
};

template <typename T>
void FindtMinDistanceFromVertex(std::unique_ptr<IGraph<T>> &graph, std::vector<T> &distance, T start_vertex) {
    distance[start_vertex] = 0;

    for (size_t i{0}; i < graph->Size() - 1; ++i) {
        for (const auto &edge : graph->GetEdges()) {
            if (distance[edge.start_edge] < std::numeric_limits<T>::max()) {
                distance[edge.end_edge] =
                    std::min(distance[edge.end_edge], distance[edge.start_edge] + edge.weight_edge);
            }
        }
    }
}

template <typename T>
void PrintMinDistanceFromVertex(std::unique_ptr<IGraph<T>> &graph) {
    T start_vertex{0};

    std::vector<T> distance(graph->Size(), std::numeric_limits<T>::max());

    FindtMinDistanceFromVertex(graph, distance, start_vertex);

    for (const auto &el : distance) {
        std::cout << (el == std::numeric_limits<T>::max() ? 30000 : el) << " ";
    }
}

template <typename T>
void Initialization(std::unique_ptr<IGraph<T>> &graph, size_t &num_vertices, size_t &num_edges) {
    std::cin >> num_vertices >> num_edges;
    graph = std::make_unique<ListGraph<T>>(num_vertices);
}

template <typename T>
void AddingEdge(std::unique_ptr<IGraph<T>> &graph, size_t num_edges) {
    for (size_t i{0}; i < num_edges; ++i) {
        T from{0}, to{0}, weight{0};
        std::cin >> from >> to >> weight;
        graph->AddEdge(from - 1, to - 1, weight);  // We use numbering from 0
    }
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);

    size_t num_vertices{0}, num_edges{0};

    std::unique_ptr<IGraph<int>> graph;

    Initialization(graph, num_vertices, num_edges);

    AddingEdge(graph, num_edges);

    PrintMinDistanceFromVertex(graph);

    return 0;
}
