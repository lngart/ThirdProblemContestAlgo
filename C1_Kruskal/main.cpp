#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>

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
DisjointSet<T>::DisjointSet(size_t num_sets) : parent_(num_sets), size_(num_sets, 1) {
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
void GetWeightMinSpanning(DisjointSet<T> &graph, size_t num_edges) {
    T weight_tree{0};

    for (size_t i{0}; i < num_edges; ++i) {
        T starting_vertex{0}, ending_vertex{0}, weight_edge{0};
        std::cin >> starting_vertex >> ending_vertex >> weight_edge;

        if (!graph.ElementsInSameSet(starting_vertex - 1, ending_vertex - 1)) {
            graph.Union(ending_vertex - 1, starting_vertex - 1);
            weight_tree += weight_edge;
        }
    }
    std::cout << weight_tree;
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);

    size_t num_vertices{0}, num_edges{0};
    std::cin >> num_vertices >> num_edges;

    DisjointSet<size_t> graph(num_vertices);

    GetWeightMinSpanning(graph, num_edges);

    return 0;
}
