#include <algorithm>
#include <ios>
#include <iostream>
#include <vector>
#include <numeric>
#include <map>

template <typename T, typename U>
class IGraph {
public:
    using Vertex = T;

    struct Edge {
        Edge(Vertex start, Vertex end, U new_weight = 1) : start_edge{start}, end_edge{end}, weight{new_weight} {
        }

        Vertex start_edge;
        Vertex end_edge;

        U weight;
    };

    virtual const std::vector<Vertex> &GetNeighbors(Vertex v) const = 0;
    virtual void AddEdge(Vertex from, Vertex to, U weight = 1) = 0;
    virtual U GetEdgeWeight(Vertex from, Vertex to) = 0;
    virtual size_t Size() const = 0;
    virtual ~IGraph() = default;
};

template <typename T, typename U>
struct EdgeLess {
    using Edge = typename IGraph<T, U>::Edge;

    bool operator()(const Edge &lhs, const Edge &rhs) const {
        return (lhs.start_edge != rhs.start_edge ? lhs.start_edge < rhs.start_edge : lhs.end_edge < rhs.end_edge);
    }
};

template <typename T, typename U>
class ListGraph : public IGraph<T, U> {
    using typename IGraph<T, U>::Vertex;
    using typename IGraph<T, U>::Edge;

    std::vector<std::vector<Vertex>> vertices_;
    std::map<Edge, U, EdgeLess<T, U>> weights_;

public:
    explicit ListGraph(size_t num_vertices) : vertices_(num_vertices) {
    }

    const std::vector<T> &GetNeighbors(Vertex v) const override {
        return vertices_[v];
    }

    void AddEdge(Vertex from, Vertex to, U weight = 1) override {
        vertices_[from].emplace_back(to);
        weights_[Edge(from, to)] += weight;
    }

    U GetEdgeWeight(Vertex from, Vertex to) override {
        return weights_[Edge(from, to)];
    }

    size_t Size() const override {
        return vertices_.size();
    }
};

enum class Color { WHITE, BLACK };

template <typename T, typename U>
U DFS(IGraph<T, U> &graph, T start_vertex, T finish_vertex, U current_flow, std::vector<Color> &colors,
      std::map<typename IGraph<T, U>::Edge, U, EdgeLess<T, U>> &edges_capacity) {
    using Edge = typename IGraph<T, U>::Edge;

    if (start_vertex == finish_vertex) {
        return current_flow;
    }

    colors[start_vertex] = Color::BLACK;

    for (T vertex{0}; vertex < graph.Size(); ++vertex) {
        if (colors[vertex] == Color::WHITE && edges_capacity[Edge(start_vertex, vertex)] > 0) {
            auto flow_change =
                DFS(graph, vertex, finish_vertex, std::min(current_flow, edges_capacity[Edge(start_vertex, vertex)]),
                    colors, edges_capacity);

            if (flow_change > 0) {
                edges_capacity[Edge(start_vertex, vertex)] -= flow_change;
                edges_capacity[Edge(vertex, start_vertex)] += flow_change;
                return flow_change;
            }
        }
    }

    return 0;
}

template <typename T, typename U>
U GetMaxFlow(IGraph<T, U> &graph, T start_vertex, T finish_vertex) {
    using Edge = typename IGraph<T, U>::Edge;

    std::map<Edge, U, EdgeLess<T, U>> edges_capacity;

    for (T vertex{0}; vertex < graph.Size(); ++vertex) {
        for (const auto &neighbor : graph.GetNeighbors(vertex)) {
            edges_capacity[Edge(vertex, neighbor)] = graph.GetEdgeWeight(vertex, neighbor);
        }
    }

    U flow_change = std::numeric_limits<U>::max();
    U max_flow{0};
    while (flow_change > 0) {
        std::vector<Color> colors(graph.Size(), Color::WHITE);
        flow_change = DFS(graph, start_vertex, finish_vertex, std::numeric_limits<U>::max(), colors, edges_capacity);
        max_flow += flow_change;
    }

    return max_flow;
}

template <typename T, typename U>
void PrintMaxFlow(IGraph<T, U> &graph) {
    T start_vertex = 0, finish_vertex = graph.Size() - 1;
    std::cout << GetMaxFlow(graph, start_vertex, finish_vertex);  // We use numbering from 0
}

void Initialization(size_t &num_vertices, size_t &num_edges) {
    std::cin >> num_vertices >> num_edges;
}

template <typename T, typename U>
void AddingEdge(IGraph<T, U> &graph, size_t num_edges) {
    for (size_t i{0}; i < num_edges; ++i) {
        T from{0}, to{0};
        U weight{0};
        std::cin >> from >> to >> weight;
        graph.AddEdge(from - 1, to - 1, weight);  // We use numbering from 0
    }
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);

    size_t num_vertices{0}, num_edges{0};

    Initialization(num_vertices, num_edges);

    ListGraph<size_t, size_t> graph(num_vertices);

    AddingEdge(graph, num_edges);

    PrintMaxFlow(graph);

    return 0;
}
