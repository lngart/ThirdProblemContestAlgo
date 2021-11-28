#include <iostream>
#include <limits>
#include <memory>
#include <utility>
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
        vertices_[to].emplace_back(from, weight);
    }

    size_t Size() const override {
        return vertices_.size();
    }
};

enum class Color { WHITE, BLACK };

template <typename T>
void FindMinDistanceBetweenTwoVertices(std::unique_ptr<IGraph<T>> &graph, std::vector<T> &distance, T start_vertex) {
    using Edge = typename IGraph<T>::Edge;

    distance[start_vertex] = 0;

    std::vector<Color> colors(graph->Size(), Color::WHITE);

    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> vertices_queue;
    vertices_queue.emplace(start_vertex, 0);

    while (!vertices_queue.empty()) {
        auto current_edge = vertices_queue.top();
        vertices_queue.pop();

        auto current_vertex = current_edge.end_edge;
        if (colors[current_vertex] == Color::WHITE) {
            colors[current_vertex] = Color::BLACK;

            for (const auto &neighbor : graph->GetNeighbors(current_vertex)) {
                auto adjacent_vertex = neighbor.end_edge, current_weight = neighbor.weight;

                if (distance[current_vertex] + current_weight < distance[adjacent_vertex]) {
                    distance[adjacent_vertex] = distance[current_vertex] + current_weight;
                    vertices_queue.emplace(adjacent_vertex, distance[adjacent_vertex]);
                }
            }
        }
    }
}

template <typename T>
void PrintMinDistanceBetweenTwoVertices(std::unique_ptr<IGraph<T>> &graph, T start_vertex) {
    std::vector<T> distance(graph->Size(), std::numeric_limits<T>::max());

    FindMinDistanceBetweenTwoVertices(graph, distance, start_vertex);

    const unsigned no_path_between_vertices = 2009000999;
    for (const auto &el : distance) {
        std::cout << (el == std::numeric_limits<T>::max() ? no_path_between_vertices : el) << " ";
    }

    std::cout << '\n';
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
        graph->AddEdge(from, to, weight);
    }
}

void RunningAlgorithm() {
    size_t num_algorithm_runs{0};
    std::cin >> num_algorithm_runs;

    std::unique_ptr<IGraph<int>> graph;

    do {
        size_t num_vertices{0}, num_edges{0};
        Initialization(graph, num_vertices, num_edges);

        AddingEdge(graph, num_edges);

        int start_vertex{0};
        std::cin >> start_vertex;

        PrintMinDistanceBetweenTwoVertices(graph, start_vertex);
    } while (--num_algorithm_runs);
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);

    RunningAlgorithm();

    return 0;
}
