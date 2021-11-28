#include <algorithm>
#include <iostream>
#include <memory>
#include <numeric>
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

        Edge &operator=(const Edge &right) {
            if (this != &right) {
                end_edge = std::move(right.end_edge);
                weight = std::move(right.weight);
            }

            return *this;
        }
    };

    virtual const std::vector<Edge> &GetNeighbors(Vertex v) const = 0;
    virtual void AddEdge(Vertex from, Vertex to, T length) = 0;
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

    void AddEdge(Vertex from, Vertex to, T length) override {
        vertices_[from].push_back({to, length});
        vertices_[to].push_back({from, length});
    }

    size_t Size() const override {
        return vertices_.size();
    }
};

template <typename T>
class DisjointSet {
public:
    explicit DisjointSet(size_t num_sets);

    void Union(T first, T second, size_t weight);
    size_t GetCurrentWeight(T element);
    size_t GetNumSets();
    size_t Find(T element);

private:
    std::vector<T> parent_;
    std::vector<size_t> size_, weight_;
};

template <typename T>
DisjointSet<T>::DisjointSet(size_t num_sets) : parent_(num_sets), size_(num_sets), weight_(num_sets) {
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
void DisjointSet<T>::Union(T first, T second, size_t weight) {
    size_t first_id = Find(first), second_id = Find(second);

    if (first_id == second_id) {
        return;
    }

    if (size_[first_id] >= size_[second_id]) {
        parent_[second_id] = first_id;
        weight_[first_id] += (weight_[second_id] + weight);

        if (size_[first_id] == size_[second_id]) {
            ++size_[first_id];
        }
    } else {
        parent_[first_id] = second_id;
        weight_[second_id] += (weight_[first_id] + weight);
    }
}

template <typename T>
size_t DisjointSet<T>::GetCurrentWeight(T element) {
    size_t element_id = Find(element);
    return weight_[element_id];
}

template <typename T>
size_t DisjointSet<T>::GetNumSets() {
    size_t cnt_sets{0};

    for (size_t i{0}; i < parent_.size(); ++i) {
        if (parent_[i] == i) {
            ++cnt_sets;
        }
    }

    return cnt_sets;
}

template <typename T>
void GetWeightMinSpanning(std::unique_ptr<IGraph<T>> &graph) {
    using Edge = typename IGraph<T>::Edge;

    DisjointSet<T> separate_trees(graph->Size());

    while (separate_trees.GetNumSets() > 1) {
        std::vector<Edge> lightest_edges(graph->Size(), {std::numeric_limits<T>::max(), std::numeric_limits<T>::max()});
        for (T vertex = 0; vertex < graph->Size(); ++vertex) {
            for (const auto &neighbor : graph->GetNeighbors(vertex)) {
                auto vertex_id = separate_trees.Find(vertex), neighbor_id = separate_trees.Find(neighbor.end_edge);

                if (vertex_id != neighbor_id && lightest_edges[vertex_id].weight >= neighbor.weight) {
                    lightest_edges[vertex_id] = neighbor;
                }
            }
        }

        for (T vertex = 0; vertex < lightest_edges.size(); ++vertex) {
            if (lightest_edges[vertex].end_edge != std::numeric_limits<T>::max()) {
                auto end_edge_id = separate_trees.Find(lightest_edges[vertex].end_edge);

                if (end_edge_id != vertex) {
                    separate_trees.Union(lightest_edges[vertex].end_edge, vertex, lightest_edges[vertex].weight);
                }
            }
        }
    }

    std::cout << separate_trees.GetCurrentWeight(0);
}

template <typename T>
void Initialization(std::unique_ptr<IGraph<T>> &graph, size_t &num_vertices, size_t &num_edges) {
    std::cin >> num_vertices >> num_edges;
    graph = std::make_unique<ListGraph<T>>(num_vertices);
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

    size_t num_vertices{0}, num_edges{0};
    std::unique_ptr<IGraph<size_t>> graph;

    Initialization(graph, num_vertices, num_edges);

    AddingEdge(graph, num_edges);

    GetWeightMinSpanning(graph);

    return 0;
}
