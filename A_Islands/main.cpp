#include <algorithm>
#include <ios>
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
        size_[first_id] += size_[second_id];
    } else {
        parent_[first_id] = second_id;
        size_[second_id] += size_[first_id];
    }
}

template <typename T>
bool DisjointSet<T>::ElementsInSameSet(T first, T second) {
    return Find(first) == Find(second);
}

template <typename T>
void ProcessingIslandsBridges(DisjointSet<T> &islands_and_bridges, size_t num_islands, size_t num_bridges) {
    size_t cnt_unexplored_islands{num_islands}, cnt_required_num_bridges{0};

    do {
        ++cnt_required_num_bridges;

        T first_island{0}, second_island{0};
        std::cin >> first_island >> second_island;

        if (!islands_and_bridges.ElementsInSameSet(first_island, second_island)) {
            islands_and_bridges.Union(first_island, second_island);
            --cnt_unexplored_islands;
        }

        if (cnt_unexplored_islands == 1) {
            std::cout << cnt_required_num_bridges;
            return;
        }
    } while (--num_bridges);
}

int main() {
    std::ios_base::sync_with_stdio(false);

    size_t num_islands{0}, num_bridges{0};
    std::cin >> num_islands >> num_bridges;

    DisjointSet<size_t> islands_and_bridges(num_islands);

    ProcessingIslandsBridges(islands_and_bridges, num_islands, num_bridges);

    return 0;
}
