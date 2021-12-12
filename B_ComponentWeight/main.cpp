#include <algorithm>
#include <ios>
#include <iostream>
#include <numeric>
#include <vector>
#include <functional>

template <typename T, class Function>
class DisjointSet {
public:
    explicit DisjointSet(size_t num_sets);

    void Union(T first, T second, size_t weight = 1);
    size_t GetCurrentWeight(T element);

private:
    size_t Find(T element);

    std::vector<T> parent_;
    std::vector<size_t> size_, weight_;

    Function operation_;
};

template <typename T, class Function>
DisjointSet<T, Function>::DisjointSet(size_t num_sets) : parent_(num_sets), size_(num_sets), weight_(num_sets) {
    for (size_t i{0}; i < num_sets; ++i) {
        parent_[i] = i;
    }
}

template <typename T, class Function>
size_t DisjointSet<T, Function>::Find(T element) {
    if (element != parent_[element]) {
        parent_[element] = Find(parent_[element]);
    }
    return parent_[element];
}

template <typename T, class Function>
void DisjointSet<T, Function>::Union(T first, T second, size_t weight) {
    size_t first_id = Find(first), second_id = Find(second);

    if (first_id == second_id) {
        operation_(weight_[first_id], weight);
        return;
    }

    if (size_[first_id] < size_[second_id]) {
        std::swap(first_id, second_id);
    } else if (size_[first_id] == size_[second_id]) {
        ++size_[first_id];
    }

    parent_[second_id] = first_id;
    auto new_weight = weight_[second_id] + weight;
    operation_(weight_[first_id], new_weight);
}

template <typename T, class Function>
size_t DisjointSet<T, Function>::GetCurrentWeight(T element) {
    size_t element_id = Find(element);
    return weight_[element_id];
}

enum class TypeAction { ADD_EDGE = 1, GET_WEIGHT = 2 };

template <typename T, class Function>
void AddingEdge(DisjointSet<T, Function> &connectivity_components) {
    T start_vertex{0}, end_vertex{0};
    size_t weight_edge{0};
    std::cin >> start_vertex >> end_vertex >> weight_edge;
    connectivity_components.Union(start_vertex - 1, end_vertex - 1, weight_edge);
}

template <typename T, class Function>
void FindWeightEdgesInComponent(DisjointSet<T, Function> &connectivity_components, std::vector<T> &required_weights,
                                size_t num_requests) {
    for (size_t i{0}; i < num_requests; ++i) {
        uint16_t type_action{0};
        std::cin >> type_action;

        if (type_action == static_cast<uint16_t>(TypeAction::ADD_EDGE)) {
            AddingEdge(connectivity_components);
        } else if (type_action == static_cast<uint16_t>(TypeAction::GET_WEIGHT)) {
            T current_vertex{0};
            std::cin >> current_vertex;
            auto current_weight = connectivity_components.GetCurrentWeight(current_vertex - 1);
            required_weights.emplace_back(current_weight);
        }
    }
}

template <typename T, class Function>
void PrintRequiredWeights(DisjointSet<T, Function> &connectivity_components, size_t num_requests) {
    std::vector<T> required_weights;

    FindWeightEdgesInComponent(connectivity_components, required_weights, num_requests);

    for (const auto &weight : required_weights) {
        std::cout << weight << '\n';
    }
}

template <typename T>
class Sum {
public:
    void operator()(T &lhs, T &rhs) {
        lhs += rhs;
    }
};

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);

    size_t num_vertices{0}, num_requests{0};
    std::cin >> num_vertices >> num_requests;

    DisjointSet<size_t, Sum<size_t>> connectivity_components(num_vertices);

    PrintRequiredWeights(connectivity_components, num_requests);

    return 0;
}
