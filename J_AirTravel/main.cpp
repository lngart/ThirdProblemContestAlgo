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
void FindMinDistanceFromVertex(std::unique_ptr<IGraph<T>> &graph, std::vector<T> &distance,
                               std::vector<T> &new_distance, size_t way_length, size_t start_vertex) {
    distance[start_vertex] = 0;
    new_distance[start_vertex] = 0;

    for (size_t i{0}; i < way_length; ++i) {
        for (const auto &edge : graph->GetEdges()) {
            if (distance[edge.start_edge] < std::numeric_limits<T>::max()) {
                new_distance[edge.end_edge] =
                    std::min(new_distance[edge.end_edge], distance[edge.start_edge] + edge.weight_edge);
            }
        }

        distance = new_distance;
    }
}

template <typename T>
void PrintMinDistanceFromVertex(std::unique_ptr<IGraph<T>> &graph, size_t way_length, size_t start_vertex,
                                size_t end_vertex) {
    std::vector<T> distance(graph->Size(), std::numeric_limits<T>::max()),
        new_distance(graph->Size(), std::numeric_limits<T>::max());

    FindMinDistanceFromVertex(graph, distance, new_distance, way_length, start_vertex);

    std::cout << (new_distance[end_vertex] == std::numeric_limits<T>::max() ? -1 : new_distance[end_vertex]);
}

template <typename T>
void Initialization(std::unique_ptr<IGraph<T>> &graph, size_t &num_cities, size_t &num_flights, size_t &num_nights,
                    size_t &starting_city, size_t &ending_city) {
    std::cin >> num_cities >> num_flights >> num_nights >> starting_city >> ending_city;
    --starting_city, --ending_city;  // We use numbering from 0
    graph = std::make_unique<ListGraph<T>>(num_cities);
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

    size_t num_cities{0}, num_flights{0}, num_nights{0}, starting_city{0}, ending_city{0};

    std::unique_ptr<IGraph<int>> flight_route;

    Initialization(flight_route, num_cities, num_flights, num_nights, starting_city, ending_city);

    AddingEdge(flight_route, num_flights);

    PrintMinDistanceFromVertex(flight_route, num_nights, starting_city, ending_city);

    return 0;
}
