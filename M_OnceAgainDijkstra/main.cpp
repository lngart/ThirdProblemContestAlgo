#include <iostream>
#include <limits>
#include <memory>
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
void FindMinDistanceBetweenTwoVertices(std::unique_ptr<IGraph<T>> &graph, std::vector<T> &distance, T start_vertex) {
    distance[start_vertex] = 0;

    std::vector<Color> colors(graph->Size(), Color::WHITE);

    for (size_t i{0}; i < graph->Size(); ++i) {
        T nearest_undeveloped_vertex = std::numeric_limits<T>::max();

        for (size_t vertex{0}; vertex < graph->Size(); ++vertex) {
            if (colors[vertex] == Color::WHITE && (nearest_undeveloped_vertex == std::numeric_limits<T>::max() ||
                                                   distance[vertex] < distance[nearest_undeveloped_vertex])) {
                nearest_undeveloped_vertex = vertex;
            }
        }

        if (distance[nearest_undeveloped_vertex] == std::numeric_limits<T>::max()) {
            break;
        }

        colors[nearest_undeveloped_vertex] = Color::BLACK;

        for (const auto &neighbor : graph->GetNeighbors(nearest_undeveloped_vertex)) {
            auto adjacent_vertex = neighbor.end_edge, current_weight = neighbor.weight;

            if (distance[nearest_undeveloped_vertex] + current_weight < distance[adjacent_vertex]) {
                distance[adjacent_vertex] = distance[nearest_undeveloped_vertex] + current_weight;
            }
        }
    }
}

template <typename T>
void PrintMinDistanceBetweenTwoVertices(std::unique_ptr<IGraph<T>> &graph, T start_vertex, T end_vertex) {
    std::vector<T> distance(graph->Size(), std::numeric_limits<T>::max());

    FindMinDistanceBetweenTwoVertices(graph, distance, start_vertex);

    std::cout << (distance[end_vertex] == std::numeric_limits<T>::max() ? -1 : distance[end_vertex]);
}

template <typename T>
void Initialization(std::unique_ptr<IGraph<T>> &graph, size_t &num_vertices, T &start_vertex, T &end_vertex) {
    std::cin >> num_vertices >> start_vertex >> end_vertex;
    --start_vertex, --end_vertex;  // We use numbering from 0
    graph = std::make_unique<ListGraph<T>>(num_vertices);
}

template <typename T>
void AddingEdge(std::unique_ptr<IGraph<T>> &graph, size_t num_vertices) {
    for (size_t from{0}; from < num_vertices; ++from) {
        for (size_t to{0}; to < num_vertices; ++to) {
            T weight{0};
            std::cin >> weight;
            if (from != to && weight != -1) {
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
    int start_vertex{0}, end_vertex{0};
    std::unique_ptr<IGraph<int>> graph;

    Initialization(graph, num_vertices, start_vertex, end_vertex);

    AddingEdge(graph, num_vertices);

    PrintMinDistanceBetweenTwoVertices(graph, start_vertex, end_vertex);

    return 0;
}
