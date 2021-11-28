#include <algorithm>
#include <ios>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>
#include <numeric>

template <typename T>
class IGraph {
public:
    using Vertex = T;

    struct Edge {
        Edge(Vertex start, Vertex end, T weight_edge)
            : start_edge{std::min(start, end)}, end_edge{std::max(start, end)}, weight{weight_edge} {
        }

        Vertex start_edge;
        Vertex end_edge;
        T weight;

        bool operator<(const Edge &other) const {
            return weight < other.weight;
        }
    };

    virtual const std::vector<T> &GetNeighbors(Vertex v) const = 0;
    virtual const std::vector<Edge> &GetEdges() const = 0;
    virtual void AddEdge(Vertex from, Vertex to, T weight) = 0;
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

    const std::vector<T> &GetNeighbors(Vertex v) const override {
        return vertices_[v];
    }

    const std::vector<Edge> &GetEdges() const override {
        return edges_;
    }

    void AddEdge(Vertex from, Vertex to, T weight) override {
        vertices_[from].push_back(to);
        vertices_[to].push_back(from);
        edges_.push_back({from, to, weight});
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

template <typename T>
T FindMinCostSavingWorld(IGraph<T> &spies_costs) {
    using Edge = typename IGraph<T>::Edge;

    T min_cost{0};

    DisjointSet<T> spies_sets(spies_costs.Size());

    std::vector<Edge> costs = spies_costs.GetEdges();

    std::sort(costs.begin(), costs.end());

    for (const auto &cost : costs) {
        if (!spies_sets.ElementsInSameSet(cost.start_edge, cost.end_edge)) {
            spies_sets.Union(cost.end_edge, cost.start_edge);
            min_cost += cost.weight;
        }
    }

    return min_cost;
}

template <typename T>
void PrintMinCostSavingWorld(IGraph<T> &spies_costs) {
    T min_cost = FindMinCostSavingWorld(spies_costs);

    std::cout << min_cost;
}

void Initialization(size_t &num_vertices) {
    std::cin >> num_vertices;
}

template <typename T>
void AddingEdge(IGraph<T> &graph) {  // We use numbering from 0
    for (size_t from{0}; from < graph.Size() - 1; ++from) {
        for (size_t to{0}; to < graph.Size() - 1; ++to) {
            T weight{0};
            std::cin >> weight;
            if (from != to) {
                graph.AddEdge(from, to, weight);
            }
        }
    }

    for (size_t current_vertex{0}; current_vertex < graph.Size() - 1; ++current_vertex) {
        T weight{0};
        std::cin >> weight;
        graph.AddEdge(current_vertex, graph.Size() - 1, weight);
    }
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);

    size_t num_spies{0};

    Initialization(num_spies);

    ListGraph<size_t> spies_costs(
        num_spies + 1);  // We additionally store information about the cost of sending each spy on a mission

    AddingEdge(spies_costs);

    PrintMinCostSavingWorld(spies_costs);

    return 0;
}
