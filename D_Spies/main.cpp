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

    virtual const std::vector<T> &GetNeighbors(Vertex v) const = 0;
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

    const std::vector<T> &GetNeighbors(Vertex v) const override {
        return vertices_[v];
    }

    void AddEdge(Vertex from, Vertex to, U weight = 1) override {
        vertices_[from].emplace_back(to);
        vertices_[to].emplace_back(from);
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
    for (size_t i{0}; i < num_sets; ++i) {
        parent_[i] = i;
    }
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
    if (size_[first_id] < size_[second_id]) {
        std::swap(first_id, second_id);
    } else if (size_[first_id] == size_[second_id]) {
        ++size_[first_id];
    }

    parent_[second_id] = first_id;
}

template <typename T>
bool DisjointSet<T>::ElementsInSameSet(T first, T second) {
    return Find(first) == Find(second);
}

template <typename T, typename U>
std::vector<typename IGraph<T, U>::Edge> GetEdgesMinSpanningTree(IGraph<T, U> &graph) {
    using Edge = typename IGraph<T, U>::Edge;

    std::vector<Edge> edges_min_spannig_tree;

    DisjointSet<T> set_edges(graph.Size());

    std::vector<Edge> graph_edges = graph.GetEdges();

    auto edges_less_cmp = [](const Edge &lhs, const Edge &rhs) { return lhs.weight < rhs.weight; };

    std::sort(graph_edges.begin(), graph_edges.end(), edges_less_cmp);

    for (const auto &current_edge : graph_edges) {
        if (!set_edges.ElementsInSameSet(current_edge.start_edge, current_edge.end_edge)) {
            set_edges.Union(current_edge.end_edge, current_edge.start_edge);
            edges_min_spannig_tree.emplace_back(current_edge);
        }
    }

    return edges_min_spannig_tree;
}

template <typename T, typename U>
U GetWeightMinSpanningTree(IGraph<T, U> &graph) {
    using Edge = typename IGraph<T, U>::Edge;

    std::vector<Edge> edges_min_spannig_tree = GetEdgesMinSpanningTree(graph);

    U weight_min_spanning_tree{0};
    for (const auto &current_edge : edges_min_spannig_tree) {
        weight_min_spanning_tree += current_edge.weight;
    }

    return weight_min_spanning_tree;
}

template <typename T, typename U>
void PrintMinCostSavingWorld(IGraph<T, U> &spies_costs) {
    U min_cost = GetWeightMinSpanningTree(spies_costs);

    std::cout << min_cost;
}

void Initialization(size_t &num_vertices) {
    std::cin >> num_vertices;
}

template <typename T, typename U>
void SetMeetingCosts(IGraph<T, U> &spies_costs) {  // We use numbering from 0
    for (size_t first_spy{0}; first_spy < spies_costs.Size() - 1; ++first_spy) {
        for (size_t second_spy{0}; second_spy < spies_costs.Size() - 1; ++second_spy) {
            U cost{0};
            std::cin >> cost;
            if (first_spy != second_spy) {
                spies_costs.AddEdge(first_spy, second_spy, cost);
            }
        }
    }

    for (size_t current_spy{0}; current_spy < spies_costs.Size() - 1; ++current_spy) {
        U cost{0};
        std::cin >> cost;
        spies_costs.AddEdge(current_spy, spies_costs.Size() - 1, cost);
    }
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);

    size_t num_spies{0};

    Initialization(num_spies);

    ListGraph<size_t, size_t> spies_costs(
        num_spies + 1);  // We additionally store information about the cost of sending each spy on a mission

    SetMeetingCosts(spies_costs);

    PrintMinCostSavingWorld(spies_costs);

    return 0;
}
