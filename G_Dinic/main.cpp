#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>
#include <numeric>
#include <queue>

template <typename T>
class IGraph {
public:
    using Vertex = T;

    virtual const std::vector<T> &GetNeighbors(Vertex v) const = 0;
    virtual void AddEdge(Vertex from, Vertex to, T weight) = 0;
    virtual T GetEdgeWeight(Vertex from, Vertex to) = 0;
    virtual size_t Size() const = 0;
    virtual ~IGraph() = default;
};

template <typename T>
class ListGraph : public IGraph<T> {
    using typename IGraph<T>::Vertex;

    std::vector<std::vector<T>> vertices_;
    std::vector<std::vector<T>> flows_;

public:
    explicit ListGraph(size_t num_vertices)
        : vertices_(num_vertices), flows_(num_vertices, std::vector<T>(num_vertices)) {
    }

    const std::vector<T> &GetNeighbors(Vertex v) const override {
        return vertices_[v];
    }

    void AddEdge(Vertex from, Vertex to, T weight) override {
        vertices_[from].emplace_back(to);
        flows_[from][to] += weight;
    }

    T GetEdgeWeight(Vertex from, Vertex to) override {
        return flows_[from][to];
    }

    size_t Size() const override {
        return vertices_.size();
    }
};

template <typename T>
bool BFS(std::unique_ptr<IGraph<T>> &graph, T start_vertex, T finish_vertex, std::vector<T> &layers,
         std::vector<std::vector<T>> &edges_capacity) {
    std::queue<T> vertices_queue;
    vertices_queue.emplace(start_vertex);

    layers[start_vertex] = 0;

    while (!vertices_queue.empty()) {
        auto current_vertex = vertices_queue.front();
        vertices_queue.pop();

        for (T vertex{0}; vertex < graph->Size(); ++vertex) {
            if (layers[vertex] == std::numeric_limits<T>::max() && edges_capacity[current_vertex][vertex] > 0) {
                layers[vertex] = layers[current_vertex] + 1;
                vertices_queue.emplace(vertex);
            }
        }
    }

    return (layers[finish_vertex] != std::numeric_limits<T>::max());
}

template <typename T>
T DFS(std::unique_ptr<IGraph<T>> &graph, T current_vertex, T finish_vertex, T current_flow, std::vector<T> &layers,
      std::vector<T> &ends_available_edges, std::vector<std::vector<T>> &edges_capacity) {
    if (current_flow == 0) {
        return 0;
    }

    if (current_vertex == finish_vertex) {
        return current_flow;
    }

    for (T &next_vertex = ends_available_edges[current_vertex]; next_vertex < graph->Size(); ++next_vertex) {
        if (layers[next_vertex] == layers[current_vertex] + 1) {
            T flow_change = DFS(graph, next_vertex, finish_vertex,
                                std::min(current_flow, edges_capacity[current_vertex][next_vertex]), layers,
                                ends_available_edges, edges_capacity);

            if (flow_change > 0) {
                edges_capacity[current_vertex][next_vertex] -= flow_change;
                edges_capacity[next_vertex][current_vertex] += flow_change;
                return flow_change;
            }
        }
    }

    return 0;
}

template <typename T>
T FindMaxFlow(std::unique_ptr<IGraph<T>> &graph, T start_vertex, T finish_vertex) {
    std::vector<std::vector<T>> edges_capacity(graph->Size(), std::vector<T>(graph->Size()));

    for (T vertex{0}; vertex < graph->Size(); ++vertex) {
        for (const auto &neighbor : graph->GetNeighbors(vertex)) {
            edges_capacity[vertex][neighbor] = graph->GetEdgeWeight(vertex, neighbor);
        }
    }

    T max_flow{0};

    std::vector<T> layers(graph->Size(), std::numeric_limits<T>::max());

    while (BFS(graph, start_vertex, finish_vertex, layers, edges_capacity)) {
        std::vector<T> ends_available_edges(graph->Size());

        T flow_change = DFS(graph, start_vertex, finish_vertex, std::numeric_limits<T>::max(), layers,
                            ends_available_edges, edges_capacity);

        while (flow_change > 0) {
            max_flow += flow_change;
            flow_change = DFS(graph, start_vertex, finish_vertex, std::numeric_limits<T>::max(), layers,
                              ends_available_edges, edges_capacity);
        }

        layers.assign(layers.size(), std::numeric_limits<T>::max());
    }
    return max_flow;
}

template <typename T>
void GetMaxFlow(std::unique_ptr<IGraph<T>> &graph) {
    T start_vertex = 0, finish_vertex = graph->Size() - 1;
    std::cout << FindMaxFlow(graph, start_vertex, finish_vertex);  // We use numbering from 0
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

    GetMaxFlow(graph);

    return 0;
}
