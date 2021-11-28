#include <iostream>
#include <ios>
#include <limits>
#include <memory>
#include <utility>
#include <vector>
#include <queue>
#include <map>

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

        bool operator<(const Edge &other) const {
            return weight < other.weight;
        }
    };

    virtual const std::vector<Edge> &GetNeighbors(Vertex v) const = 0;
    virtual bool EdgeExists(Vertex from, Vertex to) const = 0;
    virtual void AddEdge(Vertex from, Vertex to, T weight) = 0;
    virtual T GetEdgeWeight(Vertex from, Vertex to) = 0;
    virtual size_t Size() const = 0;
    virtual ~IGraph() = default;
};

template <typename T>
class ListGraph : public IGraph<T> {
    using typename IGraph<T>::Vertex;
    using typename IGraph<T>::Edge;

    std::vector<std::vector<Edge>> vertices_;
    std::map<Edge, T> weights_;

public:
    explicit ListGraph(size_t num_vertices) : vertices_(num_vertices) {
    }

    const std::vector<Edge> &GetNeighbors(Vertex v) const override {
        return vertices_[v];
    }

    bool EdgeExists(Vertex from, Vertex to) const override {
        if (vertices_[from].empty()) {
            return false;
        }

        for (const auto &neighbor : vertices_[from]) {
            if (neighbor.end_edge == to) {
                return true;
            }
        }

        return false;
    }

    void AddEdge(Vertex from, Vertex to, T weight) override {
        vertices_[from].emplace_back(to, weight);
        weights_[Edge(from, to)] += weight;

        vertices_[to].emplace_back(from, weight);
        weights_[Edge(to, from)] += weight;
    }

    T GetEdgeWeight(Vertex from, Vertex to) override {
        return weights_[Edge(from, to)];
    }

    size_t Size() const override {
        return vertices_.size();
    }
};

enum class Color { WHITE, BLACK };

template <typename T>
T GetMinDistanceBetweenTwoVertices(std::unique_ptr<IGraph<T>> &graph, const std::vector<T> &start_vertices,
                                   T end_vertex) {
    using Edge = typename IGraph<T>::Edge;

    std::vector<T> distance(graph->Size(), std::numeric_limits<T>::max());

    std::vector<Color> colors(graph->Size(), Color::WHITE);

    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> vertices_queue;
    for (const auto &vertex : start_vertices) {
        distance[vertex] = 0;
        vertices_queue.emplace(vertex, 0);
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

    return distance[end_vertex];
}

template <typename T>
void Initialization(std::unique_ptr<IGraph<T>> &countries, size_t &num_countries, size_t &num_transport_lines,
                    size_t &num_infected_countries) {
    std::cin >> num_countries >> num_transport_lines >> num_infected_countries;
    countries = std::make_unique<ListGraph<T>>(num_countries);
}

template <typename T>
void InitInfectedCountries(std::vector<T> &infected_countries) {
    for (size_t i{0}; i < infected_countries.size(); ++i) {
        std::cin >> infected_countries[i];
        --infected_countries[i];  // We use numbering from 0
    }
}

template <typename T>
void AddingEdge(std::unique_ptr<IGraph<T>> &graph, size_t num_edges) {
    for (size_t i{0}; i < num_edges; ++i) {
        T from{0}, to{0}, weight{0};
        std::cin >> from >> to >> weight;
        --from, --to;  // We use numbering from 0

        bool existing_edge = graph->EdgeExists(from, to);
        if (!existing_edge || weight < graph->GetEdgeWeight(from, to)) {
            graph->AddEdge(from, to, weight);
        }
    }
}

template <typename T>
void GetMinimumRequiredTime(std::unique_ptr<IGraph<T>> &graph, std::vector<T> &infected_countries) {
    T start_country{0}, finish_country{0};
    std::cin >> start_country >> finish_country;
    --start_country, --finish_country;  // We use numbering from 0

    size_t time_virus_spread = GetMinDistanceBetweenTwoVertices(graph, infected_countries, finish_country),
           arrival_time_people = GetMinDistanceBetweenTwoVertices(graph, {start_country}, finish_country);

    if (arrival_time_people >= time_virus_spread) {
        std::cout << -1;
    } else {
        std::cout << arrival_time_people;
    }
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);

    size_t num_countries{0}, num_transport_lines{0}, num_infected_countries{0};

    std::unique_ptr<IGraph<size_t>> countries;

    Initialization(countries, num_countries, num_transport_lines, num_infected_countries);

    std::vector<size_t> infected_countries(num_infected_countries);

    InitInfectedCountries(infected_countries);

    AddingEdge(countries, num_transport_lines);

    GetMinimumRequiredTime(countries, infected_countries);

    return 0;
}
