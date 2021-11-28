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

enum class Color { WHITE, BLACK };

template <typename T>
bool BFS(std::unique_ptr<IGraph<T>> &graph, T start_vertex, T finish_vertex,
         std::vector<std::vector<T>> &edges_capacity, std::vector<T> &parents) {
    std::vector<Color> colors(graph->Size(), Color::WHITE);

    std::queue<T> vertices_queue;
    vertices_queue.emplace(start_vertex);
    colors[start_vertex] = Color::BLACK;

    while (!vertices_queue.empty() && colors[finish_vertex] == Color::WHITE) {
        auto current_vertex = vertices_queue.front();
        vertices_queue.pop();

        for (T vertex{0}; vertex < graph->Size(); ++vertex) {
            if (colors[vertex] == Color::WHITE && edges_capacity[current_vertex][vertex] > 0) {
                parents[vertex] = current_vertex;
                colors[vertex] = Color::BLACK;
                vertices_queue.emplace(vertex);
            }
        }
    }

    return static_cast<bool>(colors[finish_vertex]);
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
    std::vector<T> parents(graph->Size(), std::numeric_limits<T>::max());

    while (BFS(graph, start_vertex, finish_vertex, edges_capacity, parents)) {
        T flow_change = std::numeric_limits<T>::max();

        for (T current_vertex{finish_vertex}; current_vertex != start_vertex;
             current_vertex = parents[current_vertex]) {
            auto current_parent = parents[current_vertex];
            flow_change = std::min(flow_change, edges_capacity[current_parent][current_vertex]);
        }

        max_flow += flow_change;

        for (T current_vertex{finish_vertex}; current_vertex != start_vertex;
             current_vertex = parents[current_vertex]) {
            auto current_parent = parents[current_vertex];
            edges_capacity[current_vertex][current_parent] += flow_change;
            edges_capacity[current_parent][current_vertex] -= flow_change;
        }
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
