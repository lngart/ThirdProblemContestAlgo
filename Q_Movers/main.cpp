#include <iostream>
#include <limits>
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
    }

    size_t Size() const override {
        return vertices_.size();
    }
};

enum class Color { WHITE, BLACK };

template <typename T>
void FindMinDistanceBetweenTwoVertices(IGraph<T> &graph, std::vector<T> &distance, T start_vertex,
                                       size_t weight_moving_to_next_vertex, size_t weight_moving_to_prev_vertex) {
    using Edge = typename IGraph<T>::Edge;

    const size_t max_possible_floor = 1'000'000;

    distance[start_vertex] = 0;

    std::vector<Color> colors(graph.Size(), Color::WHITE);

    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> vertices_queue;
    vertices_queue.emplace(start_vertex, 0);

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
                    vertices_queue.emplace(adjacent_vertex, distance[adjacent_vertex]);
                }
            }

            if (current_vertex + 1 < max_possible_floor &&
                distance[current_vertex] + weight_moving_to_next_vertex < distance[current_vertex + 1]) {
                distance[current_vertex + 1] = distance[current_vertex] + weight_moving_to_next_vertex;
                vertices_queue.emplace(current_vertex + 1, distance[current_vertex + 1]);
            }

            if (current_vertex > 0 && current_vertex < max_possible_floor &&
                distance[current_vertex] + weight_moving_to_prev_vertex < distance[current_vertex - 1]) {
                distance[current_vertex - 1] = distance[current_vertex] + weight_moving_to_prev_vertex;
                vertices_queue.emplace(current_vertex - 1, distance[current_vertex - 1]);
            }
        }
    }
}

template <typename T>
void GetMinDistanceBetweenTwoVertices(IGraph<T> &graph, T start_vertex, T end_vertex,
                                      size_t weight_moving_to_next_vertex, size_t weight_moving_to_prev_vertex) {
    std::vector<T> distance(graph.Size(), std::numeric_limits<T>::max());

    FindMinDistanceBetweenTwoVertices(graph, distance, start_vertex, weight_moving_to_next_vertex,
                                      weight_moving_to_prev_vertex);

    std::cout << distance[end_vertex];
}

void Initialization(size_t &skyscraper_floor_number, size_t &cost_climbing_stairs, size_t &cost_descending_stairs,
                    size_t &cost_entry_into_elevator, size_t &cost_removal_from_elevator, size_t &num_elevators) {
    std::cin >> skyscraper_floor_number >> cost_climbing_stairs >> cost_descending_stairs >> cost_entry_into_elevator >>
        cost_removal_from_elevator >> num_elevators;

    --skyscraper_floor_number;  // We use numbering from 0
}

template <typename T>
void BuildMovementsElevators(IGraph<T> &skyscraper, size_t num_elevators, size_t cost_entry_into_elevator,
                             size_t cost_removal_from_elevator) {
    const size_t max_possible_floor = 1'000'000;

    for (size_t i{0}; i < num_elevators; ++i) {
        size_t num_floors_with_elevator_stops{0};
        std::cin >> num_floors_with_elevator_stops;

        for (size_t j{0}; j < num_floors_with_elevator_stops; ++j) {
            T current_stop_floor{0};
            std::cin >> current_stop_floor;
            --current_stop_floor;  // We use numbering from 0

            skyscraper.AddEdge(current_stop_floor, max_possible_floor + i, cost_entry_into_elevator);
            skyscraper.AddEdge(max_possible_floor + i, current_stop_floor, cost_removal_from_elevator);
        }
    }
}

template <typename T>
void GetMinCostLiftingSafe(IGraph<T> &skyscraper, size_t skyscraper_floor_number, size_t cost_climbing_stairs,
                           size_t cost_descending_stairs) {
    size_t start_floor{0};  // We use numbering from 0
    GetMinDistanceBetweenTwoVertices(skyscraper, start_floor, skyscraper_floor_number, cost_climbing_stairs,
                                     cost_descending_stairs);
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

    ListGraph<size_t> skyscraper(max_possible_floor + num_elevators);

    BuildMovementsElevators(skyscraper, num_elevators, cost_entry_into_elevator, cost_removal_from_elevator);

    GetMinCostLiftingSafe(skyscraper, skyscraper_floor_number, cost_climbing_stairs, cost_descending_stairs);

    return 0;
}
