#include <iostream>
#include <limits>
#include <map>
#include <vector>
#include <queue>

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

    virtual const std::vector<Edge> &GetNeighbors(Vertex v) const = 0;
    virtual void AddEdge(Vertex from, Vertex to, U weight = 1) = 0;
    virtual size_t Size() const = 0;
    virtual ~IGraph() = default;
};

template <typename T, typename U>
class ListGraph : public IGraph<T, U> {
    using typename IGraph<T, U>::Vertex;
    using typename IGraph<T, U>::Edge;

    std::vector<std::vector<Edge>> vertices_;

public:
    explicit ListGraph(size_t num_vertices) : vertices_(num_vertices) {
    }

    const std::vector<Edge> &GetNeighbors(Vertex v) const override {
        return vertices_[v];
    }

    void AddEdge(Vertex from, Vertex to, U weight = 1) override {
        vertices_[from].emplace_back(from, to, weight);
    }

    size_t Size() const override {
        return vertices_.size();
    }
};

enum class Color { WHITE, BLACK };

template <typename T, typename U>
void FindMinDistanceBetweenTwoVertices(IGraph<T, U> &graph, std::vector<U> &distance, T start_vertex) {
    using Edge = typename IGraph<T, U>::Edge;

    distance[start_vertex] = 0;

    std::vector<Color> colors(graph.Size(), Color::WHITE);

    auto edges_weights_greater = [](const Edge &lhs, const Edge &rhs) { return lhs.weight > rhs.weight; };

    std::priority_queue<Edge, std::vector<Edge>, decltype(edges_weights_greater)> vertices_queue(edges_weights_greater);
    vertices_queue.emplace(start_vertex, start_vertex, 0);

    while (!vertices_queue.empty()) {
        auto current_edge = vertices_queue.top();
        vertices_queue.pop();

        auto current_vertex = current_edge.end_edge;
        if (colors[current_vertex] == Color::WHITE) {
            colors[current_vertex] = Color::BLACK;

            for (const auto &neighbor : graph.GetNeighbors(current_vertex)) {
                auto adjacent_vertex = neighbor.end_edge, current_weight = neighbor.weight;

                if (distance[current_vertex] + current_weight < distance[adjacent_vertex]) {
                    distance[adjacent_vertex] = distance[current_vertex] + current_weight;
                    vertices_queue.emplace(current_vertex, adjacent_vertex, distance[adjacent_vertex]);
                }
            }
        }
    }
}

template <typename T, typename U>
U GetMinDistanceBetweenTwoVertices(IGraph<T, U> &graph, T start_vertex, T end_vertex) {
    std::vector<U> distance(graph.Size(), std::numeric_limits<U>::max());

    FindMinDistanceBetweenTwoVertices(graph, distance, start_vertex);

    return distance[end_vertex];
}

template <typename T, typename U>
T GetMinCostLiftingSafe(IGraph<T, U> &skyscraper, size_t skyscraper_floor_number) {
    size_t start_floor{0};  // We use numbering from 0
    auto min_cost_lifting_safe = GetMinDistanceBetweenTwoVertices(skyscraper, start_floor, skyscraper_floor_number);

    return min_cost_lifting_safe;
}

template <typename T, typename U>
void PrintMinCostLiftingSafe(IGraph<T, U> &skyscraper, size_t skyscraper_floor_number) {
    auto min_cost_lifting_safe = GetMinCostLiftingSafe(skyscraper, skyscraper_floor_number);

    std::cout << min_cost_lifting_safe;
}

void Initialization(size_t &skyscraper_floor_number, size_t &cost_climbing_stairs, size_t &cost_descending_stairs,
                    size_t &cost_entry_into_elevator, size_t &cost_removal_from_elevator, size_t &num_elevators) {
    std::cin >> skyscraper_floor_number >> cost_climbing_stairs >> cost_descending_stairs >> cost_entry_into_elevator >>
        cost_removal_from_elevator >> num_elevators;

    --skyscraper_floor_number;  // We use numbering from 0
}

template <typename T, typename U>
void BuildMovementsElevators(IGraph<T, U> &skyscraper, size_t skyscraper_floor_number, size_t num_elevators,
                             size_t cost_entry_into_elevator, size_t cost_removal_from_elevator,
                             size_t cost_climbing_stairs, size_t cost_descending_stairs) {
    const size_t max_possible_floor = 1'000'000;

    size_t total_num_floors = skyscraper_floor_number;

    for (size_t i{0}; i < num_elevators; ++i) {
        size_t num_floors_with_elevator_stops{0};
        std::cin >> num_floors_with_elevator_stops;

        for (size_t j{0}; j < num_floors_with_elevator_stops; ++j) {
            T current_stop_floor{0};
            std::cin >> current_stop_floor;
            --current_stop_floor;  // We use numbering from 0

            total_num_floors = std::max(total_num_floors, current_stop_floor);

            skyscraper.AddEdge(current_stop_floor, max_possible_floor + i, cost_entry_into_elevator);
            skyscraper.AddEdge(max_possible_floor + i, current_stop_floor, cost_removal_from_elevator);
        }
    }

    for (size_t floor{0}; floor < total_num_floors; ++floor) {
        if (floor + 1 < max_possible_floor) {
            skyscraper.AddEdge(floor, floor + 1, cost_climbing_stairs);
        }
        if (floor > 0) {
            skyscraper.AddEdge(floor, floor - 1, cost_descending_stairs);
        }
    }
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);

    const size_t max_possible_floor = 1'000'000;

    size_t skyscraper_floor_number{0}, cost_climbing_stairs{0}, cost_descending_stairs{0}, cost_entry_into_elevator{0},
        cost_removal_from_elevator{0}, num_elevators{0};

    Initialization(skyscraper_floor_number, cost_climbing_stairs, cost_descending_stairs, cost_entry_into_elevator,
                   cost_removal_from_elevator, num_elevators);

    ListGraph<size_t, size_t> skyscraper(max_possible_floor + num_elevators);

    BuildMovementsElevators(skyscraper, skyscraper_floor_number, num_elevators, cost_entry_into_elevator,
                            cost_removal_from_elevator, cost_climbing_stairs, cost_descending_stairs);

    PrintMinCostLiftingSafe(skyscraper, skyscraper_floor_number);

    return 0;
}
