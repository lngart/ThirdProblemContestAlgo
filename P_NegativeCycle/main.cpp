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
void GetNegativeCycle(std::unique_ptr<IGraph<T>> &graph) {
    std::vector<T> distance(graph->Size()), parents(graph->Size(), std::numeric_limits<T>::max()), negative_cycle;

    T last_visited_vertex = std::numeric_limits<T>::max();
    for (size_t i{0}; i < graph->Size(); ++i) {
        last_visited_vertex = std::numeric_limits<T>::max();
        for (const auto &edge : graph->GetEdges()) {
            if (distance[edge.start_edge] + edge.weight_edge < distance[edge.end_edge]) {
                distance[edge.end_edge] = distance[edge.start_edge] + edge.weight_edge;
                parents[edge.end_edge] = edge.start_edge;
                last_visited_vertex = edge.end_edge;
            }
        }
    }

    if (last_visited_vertex == std::numeric_limits<T>::max()) {
        std::cout << "NO";
    } else {
        T start_vertex_in_negative_cycle = last_visited_vertex;
        for (size_t i{0}; i < graph->Size(); ++i) {
            start_vertex_in_negative_cycle = parents[start_vertex_in_negative_cycle];
        }

        for (T current_vertex = start_vertex_in_negative_cycle;; current_vertex = parents[current_vertex]) {
            negative_cycle.emplace_back(current_vertex);
            if (current_vertex == start_vertex_in_negative_cycle && negative_cycle.size() > 1) {
                break;
            }
        }

        std::cout << "YES\n" << negative_cycle.size() << '\n';

        for (int i = negative_cycle.size() - 1; i >= 0; --i) {
            std::cout << negative_cycle[i] + 1 << " ";
        }
    }
}

template <typename T>
void Initialization(std::unique_ptr<IGraph<T>> &graph, size_t &num_vertices) {
    std::cin >> num_vertices;
    graph = std::make_unique<ListGraph<T>>(num_vertices);
}

template <typename T>
void AddingEdge(std::unique_ptr<IGraph<T>> &graph, size_t num_vertices) {
    const uint16_t max_possible_edge_weight = 100000;

    for (size_t from{0}; from < num_vertices; ++from) {
        for (size_t to{0}; to < num_vertices; ++to) {
            T weight{0};
            std::cin >> weight;
            if (weight < max_possible_edge_weight) {
                graph->AddEdge(from, to, weight);  // We use numbering from 0
            }
        }
    }
}
int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);

    size_t num_vertices{0};

    std::unique_ptr<IGraph<int>> graph;

    Initialization(graph, num_vertices);

    AddingEdge(graph, num_vertices);

    GetNegativeCycle(graph);

    return 0;
}
