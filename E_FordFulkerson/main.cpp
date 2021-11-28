#include <algorithm>
#include <ios>
#include <iostream>
#include <vector>
#include <numeric>
#include <map>

template <typename T>
class IGraph {
public:
    using Vertex = T;

    struct Edge {
        Edge(Vertex start, Vertex end) : start_edge{start}, end_edge{end} {
        }

        Vertex start_edge;
        Vertex end_edge;

        bool operator<(const Edge &other) const {
            return (start_edge != other.start_edge ? start_edge < other.start_edge : end_edge < other.end_edge);
        }
    };

    virtual const std::vector<T> &GetNeighbors(Vertex v) const = 0;
    virtual void AddEdge(Vertex from, Vertex to, T weight) = 0;
    virtual T GetEdgeWeight(Vertex from, Vertex to) = 0;
    virtual size_t Size() const = 0;
    virtual ~IGraph() = default;
};

template <typename T>
class ListGraph : public IGraph<T> {
    using typename IGraph<T>::Vertex;
    using typename IGraph<T>::Edge;

    std::vector<std::vector<Vertex>> vertices_;
    std::map<Edge, T> flows_;

public:
    explicit ListGraph(size_t num_vertices) : vertices_(num_vertices) {
    }

    const std::vector<T> &GetNeighbors(Vertex v) const override {
        return vertices_[v];
    }

    void AddEdge(Vertex from, Vertex to, T weight) override {
        vertices_[from].push_back(to);
        flows_[Edge(from, to)] += weight;
    }

    T GetEdgeWeight(Vertex from, Vertex to) override {
        return flows_[Edge(from, to)];
    }

    size_t Size() const override {
        return vertices_.size();
    }
};

enum class Color { WHITE, BLACK };

template <typename T>
T DFS(IGraph<T> &graph, T start_vertex, T finish_vertex, T current_flow, std::vector<Color> &colors,
      std::map<typename IGraph<T>::Edge, T> &edges_capacity) {
    using Edge = typename IGraph<T>::Edge;

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

template <typename T>
T FindMaxFlow(IGraph<T> &graph, T start_vertex, T finish_vertex) {
    using Edge = typename IGraph<T>::Edge;

    std::map<Edge, T> edges_capacity;

    for (T vertex{0}; vertex < graph.Size(); ++vertex) {
        for (const auto &neighbor : graph.GetNeighbors(vertex)) {
            edges_capacity[Edge(vertex, neighbor)] = graph.GetEdgeWeight(vertex, neighbor);
        }
    }

    T flow_change = std::numeric_limits<T>::max();
    T max_flow{0};
    while (flow_change > 0) {
        std::vector<Color> colors(graph.Size(), Color::WHITE);
        flow_change = DFS(graph, start_vertex, finish_vertex, std::numeric_limits<T>::max(), colors, edges_capacity);
        max_flow += flow_change;
    }

    return max_flow;
}

template <typename T>
void GetMaxFlow(IGraph<T> &graph) {
    T start_vertex = 0, finish_vertex = graph.Size() - 1;
    std::cout << FindMaxFlow(graph, start_vertex, finish_vertex);  // We use numbering from 0
}

void Initialization(size_t &num_vertices, size_t &num_edges) {
    std::cin >> num_vertices >> num_edges;
}

template <typename T>
void AddingEdge(IGraph<T> &graph, size_t num_edges) {
    for (size_t i{0}; i < num_edges; ++i) {
        T from{0}, to{0}, weight{0};
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

    ListGraph<size_t> graph(num_vertices);

    AddingEdge(graph, num_edges);

    GetMaxFlow(graph);

    return 0;
}
