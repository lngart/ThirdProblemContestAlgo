#include <algorithm>
#include <iostream>
#include <vector>
#include <memory>
#include <set>
#include <queue>
#include <stack>
#include <cmath>
#include <numeric>

template <typename T>
class Field {
public:
    explicit Field(std::vector<T>& arrangement)
        : arrangement_(arrangement)
        , arrangement_size_{arrangement_.size()}
        , line_size_{static_cast<size_t>(sqrt(arrangement_size_))}
        , empty_cell_position_{0}
        , path_length_function_{0}
        , distance_traveled_function_{0}
        , remaining_distance_function_{0} {
        UpdateDate();
    }

    size_t GetLineSize() const;

    size_t GetNumInversions() const;

    const std::shared_ptr<Field>& GetParent() const;
    void UpdateParent(std::shared_ptr<Field>& new_parent);

    size_t GetEmptyCellPosition() const;

    void UpdatePathLengthFunction();

    T GetDistanceTraveledFunction() const;
    void UpdateDistanceTraveledFunction(T new_distance_traveled_function);

    const std::vector<T>& GetArrangement() const;
    void GetNeighbors(std::vector<Field<T>>& neighbors) const;

    bool operator>=(const Field& other) const {
        return path_length_function_ >= other.path_length_function_;
    }

    bool operator==(const Field& other) const {
        return arrangement_ == other.arrangement_;
    }

private:
    std::vector<T> arrangement_;
    size_t arrangement_size_;
    size_t line_size_;

    std::shared_ptr<Field> parent_;

    size_t empty_cell_position_;

    T path_length_function_;
    T distance_traveled_function_;
    T remaining_distance_function_;

    void UpdateDate();
    size_t FindEmptyCellPosition() const;
    T CalculateRemainingDistanceFunction() const;
};

template <typename T>
size_t Field<T>::GetLineSize() const {
    return line_size_;
}

template <typename T>
size_t Field<T>::GetNumInversions() const {
    size_t num_inversions{0};

    for (size_t i{0}; i < arrangement_size_; ++i) {
        for (size_t j{i + 1}; j < arrangement_size_; ++j) {
            if (arrangement_[i] > arrangement_[j] && arrangement_[i] != 0 && arrangement_[j] != 0) {
                ++num_inversions;
            }
        }
    }

    return num_inversions;
}

template <typename T>
const std::shared_ptr<Field<T>>& Field<T>::GetParent() const {
    return parent_;
}

template <typename T>
void Field<T>::UpdateParent(std::shared_ptr<Field<T>>& new_parent) {
    parent_ = new_parent;
}

template <typename T>
size_t Field<T>::GetEmptyCellPosition() const {
    return empty_cell_position_;
}

template <typename T>
void Field<T>::UpdatePathLengthFunction() {
    path_length_function_ = distance_traveled_function_ + remaining_distance_function_;
}

template <typename T>
T Field<T>::GetDistanceTraveledFunction() const {
    return distance_traveled_function_;
}

template <typename T>
void Field<T>::UpdateDistanceTraveledFunction(T new_distance_traveled_function) {
    distance_traveled_function_ = new_distance_traveled_function;
}

template <typename T>
const std::vector<T>& Field<T>::GetArrangement() const {
    return arrangement_;
}

template <typename T>
void Field<T>::UpdateDate() {
    empty_cell_position_ = FindEmptyCellPosition();
    remaining_distance_function_ = CalculateRemainingDistanceFunction();
    path_length_function_ = distance_traveled_function_ + remaining_distance_function_;
}

template <typename T>
void Field<T>::GetNeighbors(std::vector<Field<T>>& neighbors) const {
    Field<T> left_neighbour(*this), right_neighbour(*this), upper_neighbour(*this), lower_neighbour(*this);

    if (empty_cell_position_ % line_size_ >= 1) {
        std::swap(left_neighbour.arrangement_[empty_cell_position_],
                  left_neighbour.arrangement_[empty_cell_position_ - 1]);
        left_neighbour.UpdateDate();
        neighbors.emplace_back(left_neighbour);
    }

    if (empty_cell_position_ % line_size_ + 1 < line_size_) {
        std::swap(right_neighbour.arrangement_[empty_cell_position_],
                  right_neighbour.arrangement_[empty_cell_position_ + 1]);
        right_neighbour.UpdateDate();
        neighbors.emplace_back(right_neighbour);
    }

    if (empty_cell_position_ >= line_size_) {
        std::swap(upper_neighbour.arrangement_[empty_cell_position_],
                  upper_neighbour.arrangement_[empty_cell_position_ - line_size_]);
        upper_neighbour.UpdateDate();
        neighbors.emplace_back(upper_neighbour);
    }

    if (empty_cell_position_ + line_size_ < arrangement_size_) {
        std::swap(lower_neighbour.arrangement_[empty_cell_position_],
                  lower_neighbour.arrangement_[empty_cell_position_ + line_size_]);
        lower_neighbour.UpdateDate();
        neighbors.emplace_back(lower_neighbour);
    }
}

template <typename T>
size_t Field<T>::FindEmptyCellPosition() const {
    auto position_zero = std::find(arrangement_.begin(), arrangement_.end(), 0);
    return std::distance(arrangement_.begin(), position_zero);
}

template <typename T>
T Field<T>::CalculateRemainingDistanceFunction() const {
    T current_remaining_distance_function{0};

    for (size_t shift{0}; shift < arrangement_size_; ++shift) {
        if (arrangement_[shift] != 0) {
            current_remaining_distance_function += abs((shift / line_size_ - (arrangement_[shift] - 1) / line_size_)) +
                                                   abs((shift % line_size_ - (arrangement_[shift] - 1) % line_size_));
        }
    }

    return current_remaining_distance_function;
}

template <typename T>
struct Less {
    bool operator()(const Field<T>& lhs, const Field<T>& rhs) const {
        return lhs.GetArrangement() < rhs.GetArrangement();
    }
};

template <typename T>
void FindMinDistance(Field<T>& start_state, Field<T>& end_state, std::set<Field<T>, Less<T>>& treated_cells) {
    std::priority_queue<Field<T>, std::vector<Field<T>>, std::greater_equal<Field<T>>> cells_queue;
    cells_queue.push(start_state);

    while (!cells_queue.empty()) {
        auto current_cells = cells_queue.top();
        cells_queue.pop();

        if (current_cells == end_state) {
            end_state = current_cells;
            return;
        }

        treated_cells.insert(current_cells);

        std::vector<Field<T>> neighbors;
        current_cells.GetNeighbors(neighbors);

        for (auto& neighbor : neighbors) {
            T new_distance = current_cells.GetDistanceTraveledFunction() + 1;

            auto parent = std::make_shared<Field<T>>(*treated_cells.find(current_cells));

            if ((new_distance < neighbor.GetDistanceTraveledFunction()) ||
                treated_cells.find(neighbor) == treated_cells.end()) {
                neighbor.UpdateParent(parent);
                neighbor.UpdateDistanceTraveledFunction(new_distance);
                neighbor.UpdatePathLengthFunction();

                cells_queue.emplace(neighbor);
            }
        }
    }
}

template <typename T>
void FindPath(Field<T>& end_state) {
    std::stack<char> desired_path;

    T start_position{0}, end_position = end_state.GetEmptyCellPosition();
    std::shared_ptr<Field<T>> parent = end_state.GetParent();

    while (parent) {
        start_position = parent->GetEmptyCellPosition();

        if (start_position >= end_position) {
            if (start_position > end_position + 1) {
                desired_path.push('U');
            } else {
                desired_path.push('L');
            }
        } else {
            if (end_position > start_position + 1) {
                desired_path.push('D');
            } else {
                desired_path.push('R');
            }
        }

        end_position = start_position;
        parent = parent->GetParent();
    }

    std::cout << desired_path.size() << '\n';

    while (!desired_path.empty()) {
        std::cout << desired_path.top();
        desired_path.pop();
    }
}

template <typename T>
bool IsAchievable(Field<T>& field) {
    return field.GetNumInversions() % 2 != field.GetLineSize() % 2;
}

template <typename T>
void GetPath(std::vector<T>& start_arrangement, std::vector<T>& end_arrangement) {
    std::set<Field<T>, Less<T>> treated_cells;
    Field<T> start_state(start_arrangement), end_state(end_arrangement);

    if (!IsAchievable(start_state)) {
        std::cout << -1;
    } else {
        FindMinDistance(start_state, end_state, treated_cells);
        FindPath(end_state);
    }
}

template <typename T>
void InitArrangement(std::vector<T>& start_arrangement, std::vector<T>& end_arrangement) {
    for (size_t i{0}; i < start_arrangement.size(); ++i) {
        std::cin >> start_arrangement[i];
    }

    std::iota(end_arrangement.begin(), end_arrangement.end(), 1);
    end_arrangement.emplace_back(0);
}

int main() {
    const uint16_t field_size = 9;

    std::vector<int> start_arrangement(field_size), end_arrangement(field_size - 1);
    InitArrangement(start_arrangement, end_arrangement);

    GetPath(start_arrangement, end_arrangement);

    return 0;
}
