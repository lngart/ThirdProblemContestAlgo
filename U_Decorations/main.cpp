#include <iostream>
#include <ios>
#include <limits>
#include <memory>
#include <utility>
#include <vector>
#include <numeric>
#include <queue>

template <typename T>
class IGraph {
public:
    using Vertex = T;

    struct EdgeEndWeight {
        EdgeEndWeight(Vertex end, T weight_edge) : end_edge{end}, weight{weight_edge} {
        }

        Vertex end_edge;
        T weight;

        bool operator>(const EdgeEndWeight &other) const {
            return weight > other.weight;
        }
    };

    struct Edge {
        Edge(Vertex start, Vertex end, T weight_edge) : start_edge{start}, end_edge{end}, weight{weight_edge} {
        }

        Vertex start_edge;
        Vertex end_edge;
        T weight;

        bool operator<(const Edge &other) const {
            return weight < other.weight || (weight == other.weight && start_edge < other.start_edge) ||
                   (weight == other.weight && start_edge == other.start_edge && end_edge < other.end_edge);
        }

        bool operator==(const Edge &other) {
            return weight == other.weight && end_edge == other.end_edge && start_edge == other.start_edge;
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
        vertices_[from].emplace_back(from, to, weight);
        vertices_[to].emplace_back(to, from, weight);
    }

    size_t Size() const override {
        return vertices_.size();
    }
};

template <typename T>
class DisjointSet {
public:
    explicit DisjointSet(size_t num_sets);

    void Union(T first, T second);
    bool ElementsInSameSet(T first, T second);

private:
    size_t Find(T element);

    std::vector<T> parent_;
    std::vector<size_t> size_;
};

template <typename T>
DisjointSet<T>::DisjointSet(size_t num_sets) : parent_(num_sets), size_(num_sets) {
    std::iota(parent_.begin(), parent_.end(), 0);
}

template <typename T>
size_t DisjointSet<T>::Find(T element) {
    if (element != parent_[element]) {
        parent_[element] = Find(parent_[element]);
    }
    return parent_[element];
}

template <typename T>
void DisjointSet<T>::Union(T first, T second) {
    size_t first_id = Find(first), second_id = Find(second);

    if (first_id == second_id) {
        return;
    }
    if (size_[first_id] > size_[second_id]) {
        parent_[second_id] = first_id;
        ++size_[first_id];
    } else {
        parent_[first_id] = second_id;
        ++size_[second_id];
    }
}

template <typename T>
bool DisjointSet<T>::ElementsInSameSet(T first, T second) {
    return Find(first) == Find(second);
}

enum class Color { WHITE, BLACK };

template <typename T>
void FindMinDistanceBetweenTwoVertices(std::unique_ptr<IGraph<T>> &graph, std::vector<T> &distance,
                                       std::vector<T> &start_vertices) {
    using EdgeEndWeight = typename IGraph<T>::EdgeEndWeight;

    std::vector<Color> colors(graph->Size(), Color::WHITE);

    std::priority_queue<EdgeEndWeight, std::vector<EdgeEndWeight>, std::greater<EdgeEndWeight>> vertices_queue;

    for (const auto &vertex : start_vertices) {
        distance[vertex] = 0;
        vertices_queue.emplace(EdgeEndWeight(vertex, 0));
    }

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
struct Route {
    using Edge = typename IGraph<T>::Edge;

    Route(Edge edge, int number, bool is_edge) : edge(edge), number(number), is_edge(is_edge) {
    }

    Edge edge;
    int number;
    bool is_edge;

    bool operator>(const Route &other) const {
        return edge.weight > other.edge.weight || (edge.weight == other.edge.weight && number > other.number);
    }
};

template <typename T>
void ItPossibleDrive(std::unique_ptr<IGraph<T>> &routes, std::vector<T> &intersections_with_stations,
                     std::vector<typename IGraph<T>::Edge> &roads,
                     std::priority_queue<Route<T>, std::vector<Route<T>>, std::greater<Route<T>>> &current_routes) {
    using Edge = typename IGraph<T>::Edge;

    auto num_requests = current_routes.size();

    std::vector<T> distance(routes->Size(), std::numeric_limits<T>::max());

    FindMinDistanceBetweenTwoVertices(routes, distance, intersections_with_stations);

    for (const auto &current_road : roads) {
        auto current_start = current_road.start_edge, current_end = current_road.end_edge,
             current_weight = current_road.weight + distance[current_start] + distance[current_end];
        current_routes.emplace(Edge(current_start, current_end, current_weight), std::numeric_limits<T>::min(), true);
    }

    DisjointSet<T> routes_sets(routes->Size());

    std::vector<bool> it_possible_drive(num_requests);

    while (!current_routes.empty()) {
        auto current_route = current_routes.top();
        current_routes.pop();

        bool is_edge = current_route.is_edge;
        auto number = current_route.number, current_end = current_route.edge.end_edge,
             current_start = current_route.edge.start_edge;

        if (!is_edge) {
            it_possible_drive[number] = routes_sets.ElementsInSameSet(current_end, current_start);
        } else {
            routes_sets.Union(current_end, current_start);
        }
    }

    for (const auto &current_requests : it_possible_drive) {
        std::cout << (current_requests ? "YES" : "NO") << '\n';
    }
}

template <typename T>
void Initialization(std::unique_ptr<IGraph<T>> &routes, size_t &num_intersections, size_t &num_stations,
                    size_t &num_roads) {
    std::cin >> num_intersections >> num_stations >> num_roads;
    routes = std::make_unique<ListGraph<T>>(num_intersections);
}

template <typename T>
void InitStations(std::vector<T> &intersections_with_stations) {
    for (size_t i{0}; i < intersections_with_stations.size(); ++i) {
        std::cin >> intersections_with_stations[i];
        --intersections_with_stations[i];  // We use numbering from 0
    }
}

template <typename T>
void InitRoads(std::unique_ptr<IGraph<T>> &routes, size_t num_roads, std::vector<typename IGraph<T>::Edge> &roads) {
    for (size_t i{0}; i < num_roads; ++i) {
        T from{0}, to{0}, length{0};
        std::cin >> from >> to >> length;
        --from, --to;  // We use numbering from 0

        routes->AddEdge(from, to, length);

        roads.emplace_back(from, to, length);
    }
}

template <typename T>
void GetRequests(size_t &num_requests,
                 std::priority_queue<Route<T>, std::vector<Route<T>>, std::greater<Route<T>>> &current_routes) {
    using Edge = typename IGraph<T>::Edge;

    std::cin >> num_requests;
    for (size_t i{0}; i < num_requests; ++i) {
        T from{0}, to{0}, capacity{0};
        std::cin >> from >> to >> capacity;
        --from, --to;  // We use numbering from 0

        current_routes.emplace(Edge(from, to, capacity), static_cast<int>(i), false);
    }
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);

    std::unique_ptr<IGraph<int>> routes;

    size_t num_intersections{0}, num_stations{0}, num_roads{0};
    Initialization(routes, num_intersections, num_stations, num_roads);

    std::vector<int> intersections_with_stations(num_stations);
    InitStations(intersections_with_stations);

    std::vector<typename IGraph<int>::Edge> roads;
    InitRoads(routes, num_roads, roads);

    size_t num_requests{0};
    std::priority_queue<Route<int>, std::vector<Route<int>>, std::greater<Route<int>>> current_routes;
    GetRequests(num_requests, current_routes);

    ItPossibleDrive(routes, intersections_with_stations, roads, current_routes);

    return 0;
}
