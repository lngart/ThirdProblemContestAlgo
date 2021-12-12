#include <iostream>
#include <limits>
#include <vector>
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
void FindMinDistancesBetweenAllVertices(IGraph<T, U> &graph, std::vector<std::vector<U>> &distance) {
    for (size_t current_vertex{0}; current_vertex < graph.Size(); ++current_vertex) {
        for (size_t next_vertex{0}; next_vertex < graph.Size(); ++next_vertex) {
            distance[current_vertex][next_vertex] = graph.GetEdgeWeight(current_vertex, next_vertex);
        }
    }

    for (size_t intermediate_vertex{0}; intermediate_vertex < graph.Size(); ++intermediate_vertex) {
        for (size_t current_vertex{0}; current_vertex < graph.Size(); ++current_vertex) {
            for (size_t next_vertex{0}; next_vertex < graph.Size(); ++next_vertex) {
                distance[current_vertex][next_vertex] =
                    std::min(distance[current_vertex][next_vertex], distance[current_vertex][intermediate_vertex] +
                                                                        distance[intermediate_vertex][next_vertex]);
            }
        }
    }
}

template <typename T, typename U>
void PrintMinDistancesBetweenAllVertices(IGraph<T, U> &graph) {
    std::vector<std::vector<U>> distance(graph.Size(), std::vector<U>(graph.Size()));

    FindMinDistancesBetweenAllVertices(graph, distance);

    for (size_t i{0}; i < graph.Size(); ++i) {
        for (size_t j{0}; j < graph.Size(); ++j) {
            std::cout << distance[i][j] << " ";
        }
        std::cout << '\n';
    }
}

void Initialization(size_t &num_vertices) {
    std::cin >> num_vertices;
}

template <typename T, typename U>
void AddingEdge(IGraph<T, U> &graph, size_t num_vertices) {
    for (size_t from{0}; from < num_vertices; ++from) {
        for (size_t to{0}; to < num_vertices; ++to) {
            U weight{0};
            std::cin >> weight;
            if (from != to) {
                graph.AddEdge(from, to, weight);  // We use numbering from 0
            }
        }
    }
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);

    size_t num_vertices{0};

    Initialization(num_vertices);

    ListGraph<int, int> graph(num_vertices);

    AddingEdge(graph, num_vertices);

    PrintMinDistancesBetweenAllVertices(graph);

    return 0;
}
