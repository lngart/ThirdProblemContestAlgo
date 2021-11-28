#include <algorithm>
#include <ios>
#include <iostream>
#include <numeric>
#include <vector>

template <typename T>
class DisjointSet {
public:
    explicit DisjointSet(size_t num_sets);

    void Union(T first, T second, size_t weight);
    size_t GetCurrentWeight(T element);

private:
    size_t Find(T element);

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
        weight_[first_id] += weight;
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

enum class TypeAction { ADD_EDGE = 1, GET_WEIGHT = 2 };

template <typename T>
void AddingEdge(DisjointSet<T> &connectivity_components) {
    T start_vertex{0}, end_vertex{0};
    size_t weight_edge{0};
    std::cin >> start_vertex >> end_vertex >> weight_edge;
    connectivity_components.Union(start_vertex - 1, end_vertex - 1, weight_edge);
}

template <typename T>
void FindWeightEdgesInComponent(DisjointSet<T> &connectivity_components, std::vector<T> &required_weights,
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

template <typename T>
void PrintRequiredWeights(DisjointSet<T> &connectivity_components, size_t num_requests) {
    std::vector<T> required_weights;

    FindWeightEdgesInComponent(connectivity_components, required_weights, num_requests);

    for (const auto &weight : required_weights) {
        std::cout << weight << '\n';
    }
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);

    size_t num_vertices{0}, num_requests{0};
    std::cin >> num_vertices >> num_requests;

    DisjointSet<int> connectivity_components(num_vertices);

    PrintRequiredWeights(connectivity_components, num_requests);

    return 0;
}
