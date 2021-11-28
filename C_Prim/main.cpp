#include <iostream>
#include <limits>
#include <memory>
#include <vector>
#include <queue>

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
    virtual void AddEdge(Vertex from, Vertex to, T length) = 0;
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

    void AddEdge(Vertex from, Vertex to, T length) override {
        vertices_[from].push_back({to, length});
        vertices_[to].push_back({from, length});
    }

    size_t Size() const override {
        return vertices_.size();
    }
};

enum class Color { WHITE, BLACK };

template <typename T>
void GetWeightMinSpanning(std::unique_ptr<IGraph<T>> &graph) {
    using Edge = typename IGraph<T>::Edge;

    std::vector<T> min_weight_out_edge(graph->Size(), std::numeric_limits<T>::max());

    std::vector<Color> colors(graph->Size(), Color::WHITE);

    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> vertices_queue;
    vertices_queue.push({0, 0});

    T weight_tree{0};

    while (!vertices_queue.empty()) {
        auto lightest_edge = vertices_queue.top();
        vertices_queue.pop();

        auto current_vertex = lightest_edge.end_edge;
        if (colors[current_vertex] == Color::WHITE) {
            colors[current_vertex] = Color::BLACK;

            weight_tree += lightest_edge.weight;

            for (const auto &neighbor : graph->GetNeighbors(current_vertex)) {
                const auto ending_vertex = neighbor.end_edge, current_weight = neighbor.weight;

                if (current_weight < min_weight_out_edge[ending_vertex]) {
                    min_weight_out_edge[ending_vertex] = current_weight;

                    if (colors[ending_vertex] == Color::WHITE) {
                        vertices_queue.push({ending_vertex, min_weight_out_edge[ending_vertex]});
                    };
                }
            }
        }
    }

    std::cout << weight_tree;
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
    std::unique_ptr<IGraph<size_t>> graph;

    Initialization(graph, num_vertices, num_edges);

    AddingEdge(graph, num_edges);

    GetWeightMinSpanning(graph);

    return 0;
}
