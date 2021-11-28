#include <algorithm>
#include <ios>
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
void GetMinDistancesBetweenAllVertices(std::unique_ptr<IGraph<T>> &graph, std::vector<std::vector<T>> &distance,
                                       std::vector<std::vector<T>> &parents) {

    for (size_t vertex{0}; vertex < graph->Size(); ++vertex) {
        for (const auto &neighbor : graph->GetNeighbors(vertex)) {
            auto adjacent_vertex = neighbor.end_edge, current_weight = neighbor.weight;
            distance[vertex][adjacent_vertex] = current_weight;
        }
    }

    for (size_t i{0}; i < graph->Size(); ++i) {
        distance[i][i] = 0;
    }

    for (size_t i{0}; i < graph->Size(); ++i) {
        for (size_t j{0}; j < graph->Size(); ++j) {
            for (size_t k{0}; k < graph->Size(); ++k) {
                if (distance[j][i] != std::numeric_limits<T>::max() &&
                    distance[i][k] != std::numeric_limits<T>::max() &&
                    distance[j][k] > distance[j][i] + distance[i][k]) {
                    distance[j][k] = distance[j][i] + distance[i][k];
                    parents[j][k] = i;
                }
            }
        }
    }
}

template <typename T>
void FindRoute(std::unique_ptr<IGraph<T>> &flights, std::vector<std::vector<T>> &parents, std::vector<T> &way, T from,
               T to) {
    if (parents[from][to] == std::numeric_limits<T>::max()) {
        way.emplace_back(from);
        return;
    }

    FindRoute(flights, parents, way, from, parents[from][to]);
    FindRoute(flights, parents, way, parents[from][to], to);
}

template <typename T>
void FindMostProfitableRoute(std::unique_ptr<IGraph<T>> &flights, std::vector<T> &way, std::vector<T> &concerts,
                             bool &is_route_with_infinite_weight) {
    std::vector<std::vector<T>> distance(flights->Size(),
                                         std::vector<T>(flights->Size(), std::numeric_limits<T>::max())),
        parents(flights->Size(), std::vector<T>(flights->Size(), std::numeric_limits<T>::max()));

    GetMinDistancesBetweenAllVertices(flights, distance, parents);

    for (size_t i{0}; i < concerts.size() - 1; ++i) {
        auto first_city = concerts[i], second_city = concerts[i + 1];

        for (size_t j{0}; j < flights->Size(); ++j) {
            if (distance[j][j] < 0 && distance[first_city][j] != std::numeric_limits<T>::max() &&
                distance[j][second_city] != std::numeric_limits<T>::max()) {
                is_route_with_infinite_weight = true;
                return;
            }
        }

        FindRoute(flights, parents, way, first_city, second_city);
    }

    way.emplace_back(concerts.back());
}

template <typename T>
void GetMostProfitableRoute(std::unique_ptr<IGraph<T>> &flights, std::vector<T> &concerts,
                            std::vector<std::vector<T>> &numbering_of_flights) {
    bool is_route_with_infinite_weight{false};
    std::vector<T> way;

    FindMostProfitableRoute(flights, way, concerts, is_route_with_infinite_weight);

    if (is_route_with_infinite_weight) {
        std::cout << "infinitely kind";
        return;
    }

    std::vector<T> most_profitable_flights;
    for (size_t i{0}; i < way.size() - 1; ++i) {
        auto current_city = way[i], next_city = way[i + 1];
        if (current_city != next_city) {
            most_profitable_flights.emplace_back(numbering_of_flights[current_city][next_city] + 1);
        }
    }

    std::cout << most_profitable_flights.size() << '\n';
    for (const auto &current_flight_number : most_profitable_flights) {
        std::cout << current_flight_number << " ";
    }
}

template <typename T>
void Initialization(std::unique_ptr<IGraph<T>> &flights, size_t &num_cities, size_t &num_flights,
                    size_t &num_concerts) {
    std::cin >> num_cities >> num_flights >> num_concerts;
    flights = std::make_unique<ListGraph<T>>(num_cities);
}

template <typename T>
void InitFlights(std::unique_ptr<IGraph<T>> &graph, std::vector<std::vector<T>> &numbering_of_flights,
                 size_t num_flights) {
    for (size_t i{0}; i < num_flights; ++i) {
        T from{0}, to{0}, weight{0};
        std::cin >> from >> to >> weight;
        --from, --to;  // We use numbering from 0

        graph->AddEdge(from, to, -weight);

        numbering_of_flights[from][to] = i;
    }
}

template <typename T>
void InitConcerts(std::vector<T> &concerts) {
    for (size_t i{0}; i < concerts.size(); ++i) {
        std::cin >> concerts[i];
        --concerts[i];  // We use numbering from 0
    }
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);

    size_t num_cities{0}, num_flights{0}, num_concerts{0};
    std::unique_ptr<IGraph<int>> flights;
    Initialization(flights, num_cities, num_flights, num_concerts);

    std::vector<std::vector<int>> numbering_of_flights(num_cities, std::vector<int>(num_cities));
    InitFlights(flights, numbering_of_flights, num_flights);

    std::vector<int> concerts(num_concerts);
    InitConcerts(concerts);

    GetMostProfitableRoute(flights, concerts, numbering_of_flights);

    return 0;
}
