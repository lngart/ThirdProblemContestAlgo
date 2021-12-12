#include <algorithm>
#include <ios>
#include <iostream>
#include <vector>
#include <numeric>

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

    virtual const std::vector<Vertex> &GetNeighbors(Vertex v) const = 0;
    virtual void AddEdge(Vertex from, Vertex to, U weight = 1) = 0;
    virtual const std::vector<Edge> &GetEdges() const = 0;
    virtual size_t Size() const = 0;
    virtual ~IGraph() = default;
};

template <typename T, typename U>
class ListGraph : public IGraph<T, U> {
    using typename IGraph<T, U>::Vertex;
    using typename IGraph<T, U>::Edge;

    std::vector<std::vector<Vertex>> vertices_;
    std::vector<Edge> edges_;

public:
    explicit ListGraph(size_t num_vertices) : vertices_(num_vertices) {
    }

    const std::vector<Vertex> &GetNeighbors(Vertex v) const override {
        return vertices_[v];
    }

    void AddEdge(Vertex from, Vertex to, U weight = 1) override {
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
struct EdgeWeights {
    T old_weight_edge;
    T new_weight_edge;
};

template <typename T>
void FindMinArrivalTime(IGraph<T, EdgeWeights<T>> &teleportation_stations, std::vector<T> &arrival_times_at_station,
                        T start_station) {
    arrival_times_at_station[start_station] = 0;

    const auto &teleportations = teleportation_stations.GetEdges();

    for (size_t i{0}; i < teleportations.size(); ++i) {
        for (const auto &current_teleportation : teleportations) {
            if (arrival_times_at_station[current_teleportation.start_edge] <=
                current_teleportation.weight.old_weight_edge) {
                arrival_times_at_station[current_teleportation.end_edge] =
                    std::min(arrival_times_at_station[current_teleportation.end_edge],
                             current_teleportation.weight.new_weight_edge);
            }
        }
    }
}

template <typename T>
void PrintMinArrivalTime(IGraph<T, EdgeWeights<T>> &teleportation_stations, T start_station, T end_station) {
    std::vector<T> arrival_times_at_station(teleportation_stations.Size(), std::numeric_limits<T>::max());

    FindMinArrivalTime(teleportation_stations, arrival_times_at_station, start_station);

    std::cout << arrival_times_at_station[end_station];
}

template <typename T>
void Initialization(size_t &num_stations, size_t &num_teleportations, T &start_station, T &end_station) {
    std::cin >> num_stations >> start_station >> end_station >> num_teleportations;
    --start_station, --end_station;  // We use numbering from 0
}

template <typename T>
void BuildTeleportationGraph(IGraph<T, EdgeWeights<T>> &graph, size_t num_edges) {
    for (size_t i{0}; i < num_edges; ++i) {
        T from{0}, old_weight{0}, to{0}, new_weight{0};
        std::cin >> from >> old_weight >> to >> new_weight;
        graph.AddEdge(from - 1, to - 1, {old_weight, new_weight});  // We use numbering from 0
    }
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);

    size_t num_stations{0}, num_teleportations{0};
    int start_station{0}, end_station{0};

    Initialization(num_stations, num_teleportations, start_station, end_station);

    ListGraph<int, EdgeWeights<int>> teleportation_stations(num_stations);

    BuildTeleportationGraph(teleportation_stations, num_teleportations);

    PrintMinArrivalTime(teleportation_stations, start_station, end_station);

    return 0;
}
