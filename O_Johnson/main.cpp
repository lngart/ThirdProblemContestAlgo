#include <iostream>
#include <ios>
#include <limits>
#include <memory>
#include <utility>
#include <vector>
#include <queue>
#include <map>

template <typename T>
class IGraph {
public:
    using Vertex = T;

    struct Edge {
        Edge(Vertex start, Vertex end, T weight) : start_edge{start}, end_edge{end}, weight_edge{weight} {
        }

        Edge(Vertex end, T weight) : end_edge{end}, weight_edge{weight} {
        }

        Vertex start_edge;
        Vertex end_edge;
        T weight_edge;

        bool operator>(const Edge &other) const {
            return weight_edge > other.weight_edge;
        }

        bool operator<(const Edge &other) const {
            return weight_edge < other.weight_edge;
        }
    };

    struct EdgeVertices {
        EdgeVertices(Vertex start, Vertex end) : start_edge{start}, end_edge{end} {
        }

        Vertex start_edge;
        Vertex end_edge;

        bool operator<(const EdgeVertices &other) const {
            return (start_edge != other.start_edge ? start_edge < other.start_edge : end_edge < other.end_edge);
        }
    };

    virtual const std::vector<Edge> &GetNeighbors(Vertex v) const = 0;
    virtual void AddEdge(Vertex from, Vertex to, T weight) = 0;
    virtual const std::vector<Edge> &GetEdges() const = 0;
    virtual size_t Size() const = 0;
    virtual ~IGraph() = default;
};

template <typename T>
class ListGraph : public IGraph<T> {
    using typename IGraph<T>::Vertex;
    using typename IGraph<T>::Edge;

    std::vector<std::vector<Edge>> vertices_;
    std::vector<Edge> edges_;

public:
    explicit ListGraph(size_t num_vertices) : vertices_(num_vertices) {
    }

    const std::vector<Edge> &GetNeighbors(Vertex v) const override {
        return vertices_[v];
    }

    void AddEdge(Vertex from, Vertex to, T weight) override {
        vertices_[from].emplace_back(from, to, weight);
        edges_.emplace_back(Edge(from, to, weight));
    }

    const std::vector<Edge> &GetEdges() const override {
        return edges_;
    }

    size_t Size() const override {
        return vertices_.size();
    }
};

template <typename T>
void FindMinDistancesFromDummyVertex(std::unique_ptr<IGraph<T>> &graph, std::vector<T> &distance, T start_vertex) {
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

enum class Color { WHITE, BLACK };

template <typename T>
void GetMinDistancesBetweenTwoVertices(std::unique_ptr<IGraph<T>> &graph,
                                       std::map<typename IGraph<T>::EdgeVertices, T> &weights, std::vector<T> &distance,
                                       T start_vertex) {
    using EdgeVertices = typename IGraph<T>::EdgeVertices;
    using Edge = typename IGraph<T>::Edge;

    distance[start_vertex] = 0;

    std::vector<Color> colors(graph->Size() - 1, Color::WHITE);

    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> vertices_queue;
    vertices_queue.emplace(start_vertex, 0);

    while (!vertices_queue.empty()) {
        auto current_edge = vertices_queue.top();
        vertices_queue.pop();

        auto current_vertex = current_edge.end_edge;
        if (colors[current_vertex] == Color::WHITE) {
            colors[current_vertex] = Color::BLACK;

            for (const auto &neighbor : graph->GetNeighbors(current_vertex)) {
                auto adjacent_vertex = neighbor.end_edge,
                     current_weight = weights[EdgeVertices(current_vertex, adjacent_vertex)];

                if (distance[current_vertex] + current_weight < distance[adjacent_vertex]) {
                    distance[adjacent_vertex] = distance[current_vertex] + current_weight;
                    vertices_queue.emplace(adjacent_vertex, distance[adjacent_vertex]);
                }
            }
        }
    }
}

template <typename T>
void GetMaxShortestPath(std::unique_ptr<IGraph<T>> &graph) {
    using EdgeVertices = typename ListGraph<T>::EdgeVertices;

    T graph_size = graph->Size(), dummy_vertex = graph_size - 1;

    std::vector<T> distance_from_dummy_vertex(graph_size, std::numeric_limits<T>::max());

    FindMinDistancesFromDummyVertex(graph, distance_from_dummy_vertex, dummy_vertex);

    std::map<EdgeVertices, T> updated_weights;
    for (T vertex{0}; vertex < graph_size - 1; ++vertex) {
        for (const auto &neighbor : graph->GetNeighbors(vertex)) {
            updated_weights[(EdgeVertices(neighbor.start_edge, neighbor.end_edge))] +=
                (neighbor.weight_edge + distance_from_dummy_vertex[neighbor.start_edge] -
                 distance_from_dummy_vertex[neighbor.end_edge]);
        }
    }

    T max_shortest_path{0};
    for (T vertex{0}; vertex < graph_size - 1; ++vertex) {
        std::vector<T> distances_in_reweighting_graph(graph_size - 1, std::numeric_limits<T>::max());
        GetMinDistancesBetweenTwoVertices(graph, updated_weights, distances_in_reweighting_graph, vertex);

        for (T current_vertex{0}; current_vertex < graph_size - 1; ++current_vertex) {
            auto current_length_path =
                distances_in_reweighting_graph[current_vertex] -
                (distance_from_dummy_vertex[vertex] - distance_from_dummy_vertex[current_vertex]);

            if (distances_in_reweighting_graph[current_vertex] != std::numeric_limits<T>::max() &&
                max_shortest_path < current_length_path) {
                max_shortest_path = current_length_path;
            }
        }
    }

    std::cout << max_shortest_path;
}

template <typename T>
void Initialization(std::unique_ptr<IGraph<T>> &graph, size_t &num_vertices, size_t &num_edges) {
    std::cin >> num_vertices >> num_edges;
    ++num_vertices;
    graph = std::make_unique<ListGraph<T>>(num_vertices);
}

template <typename T>
void AddingEdge(std::unique_ptr<IGraph<T>> &graph, size_t &num_edges) {
    for (size_t i{0}; i < num_edges; ++i) {
        T from{0}, to{0}, weight{0};
        std::cin >> from >> to >> weight;
        graph->AddEdge(from, to, weight);
    }

    auto dummy_vertex = graph->Size() - 1;
    for (size_t vertex{0}; vertex < graph->Size() - 1; ++vertex) {
        graph->AddEdge(dummy_vertex, vertex, 0);
        ++num_edges;
    }
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);

    std::unique_ptr<IGraph<int>> graph;

    size_t num_vertices{0}, num_edges{0};

    Initialization(graph, num_vertices, num_edges);

    AddingEdge(graph, num_edges);

    GetMaxShortestPath(graph);

    return 0;
}
