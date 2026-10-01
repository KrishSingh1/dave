///////////////////////////////////////////////////////////////////////////////////////
/// \file KDTree.h
/// \brief A k-d tree for coordinates and their associated values.
///
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_KDTree_H
#define LPJ_GUESS_KDTree_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>

template<typename Coordinate, std::size_t Dimensions>
class Point {
    static_assert(Dimensions > 0, "a point must have at least one dimension");

public:
    explicit Point(const std::array<Coordinate, Dimensions>& coordinates)
        : coordinates(coordinates) {
        validate();
    }

    Point(std::initializer_list<Coordinate> coordinates) {
        if (coordinates.size() != Dimensions)
            throw std::invalid_argument("incorrect number of point coordinates");
        std::copy(coordinates.begin(), coordinates.end(), this->coordinates.begin());
        validate();
    }

    Coordinate get(std::size_t dimension) const { return coordinates[dimension]; }

    // Convert before subtraction so integer coordinates cannot overflow in
    // their original type. The result is used for tree comparisons.
    double squared_distance(const Point& other) const {
        double squared = 0.0;
        for (std::size_t i = 0; i < Dimensions; ++i) {
            const double difference = static_cast<double>(coordinates[i]) -
                                      static_cast<double>(other.coordinates[i]);
            squared += difference * difference;
        }
        return squared;
    }

    double distance(const Point& other) const {
        return std::sqrt(squared_distance(other));
    }

private:
    void validate() const {
        for (const auto& coordinate : coordinates) {
            if (!std::isfinite(static_cast<double>(coordinate)))
                throw std::invalid_argument("point coordinates must be finite");
        }
    }

    std::array<Coordinate, Dimensions> coordinates;
};

template<typename Coordinate, std::size_t Dimensions>
std::ostream& operator<<(std::ostream& out, const Point<Coordinate, Dimensions>& point) {
    out << '(';
    for (std::size_t i = 0; i < Dimensions; ++i) {
        if (i != 0) out << ", ";
        out << point.get(i);
    }
    return out << ')';
}

template<typename Value, typename Coordinate = double, std::size_t Dimensions = 2>
class KDTree {
public:
    using PointType = Point<Coordinate, Dimensions>;

    struct Entry {
        PointType point;
        Value value;
    };

    // References remain valid while the index exists and is not moved or assigned.
    struct NearestResult {
        const PointType& point;
        const Value& value;
        double distance;
        std::size_t visited;
    };

    // Duplicate coordinates retain the last value in input order.
    explicit KDTree(std::vector<Entry> input_entries)
        : entries(std::move(input_entries)) {
        std::stable_sort(entries.begin(), entries.end(),
                         [](const Entry& a, const Entry& b) {
                             return lexicographically_less(a.point, b.point);
                         });
        std::vector<Entry> unique_entries;
        unique_entries.reserve(entries.size());
        for (auto& entry : entries) {
            if (!unique_entries.empty() &&
                same_point(unique_entries.back().point, entry.point)) {
                unique_entries.back() = std::move(entry);
            } else {
                unique_entries.push_back(std::move(entry));
            }
        }
        entries = std::move(unique_entries);
        build(0, entries.size(), 0);
    }

    bool empty() const { return entries.empty(); }
    std::size_t size() const { return entries.size(); }

    NearestResult nearest(const PointType& query) const {
        if (empty())
            throw std::logic_error("spatial index is empty");

        std::size_t best = entries.size();
        std::size_t visited = 0;
        double best_squared_distance = std::numeric_limits<double>::infinity();
        search(0, entries.size(), 0, query, best, best_squared_distance, visited);
        const Entry& entry = entries[best];
        return {entry.point, entry.value, std::sqrt(best_squared_distance), visited};
    }

    // Returns null if no coordinate matches, using exact coordinate equality.
    const Entry* find_exact(const PointType& query) const {
        std::size_t begin = 0;
        std::size_t end = entries.size();
        std::size_t axis = 0;
        while (begin != end) {
            const std::size_t middle = begin + (end - begin) / 2;
            const PointType& point = entries[middle].point;
            if (same_point(query, point))
                return &entries[middle];
            if (axis_less(query, point, axis))
                end = middle;
            else
                begin = middle + 1;
            axis = (axis + 1) % Dimensions;
        }
        return nullptr;
    }

private:
    static bool lexicographically_less(const PointType& a, const PointType& b) {
        for (std::size_t i = 0; i < Dimensions; ++i) {
            if (a.get(i) < b.get(i)) return true;
            if (b.get(i) < a.get(i)) return false;
        }
        return false;
    }

    static bool same_point(const PointType& a, const PointType& b) {
        return !lexicographically_less(a, b) && !lexicographically_less(b, a);
    }

    static bool axis_less(const PointType& a, const PointType& b, std::size_t axis) {
        if (a.get(axis) < b.get(axis)) return true;
        if (b.get(axis) < a.get(axis)) return false;
        return lexicographically_less(a, b);
    }

    // Each range is a subtree with its root at the midpoint.
    void build(std::size_t begin, std::size_t end, std::size_t axis) {
        if (begin == end) return;
        const std::size_t middle = begin + (end - begin) / 2;
        std::nth_element(entries.begin() + begin, entries.begin() + middle,
                         entries.begin() + end,
                         [axis](const Entry& a, const Entry& b) {
                             return axis_less(a.point, b.point, axis);
                         });
        const std::size_t next_axis = (axis + 1) % Dimensions;
        build(begin, middle, next_axis);
        build(middle + 1, end, next_axis);
    }

    void search(std::size_t begin, std::size_t end, std::size_t axis,
                const PointType& query, std::size_t& best,
                double& best_squared_distance, std::size_t& visited) const {
        if (begin == end) return;

        const std::size_t middle = begin + (end - begin) / 2;
        const Entry& candidate = entries[middle];
        ++visited;
        const double squared_distance = candidate.point.squared_distance(query);
        if (best == entries.size() || squared_distance < best_squared_distance ||
            (squared_distance == best_squared_distance &&
             lexicographically_less(candidate.point, entries[best].point))) {
            best = middle;
            best_squared_distance = squared_distance;
        }

        const double offset = static_cast<double>(query.get(axis)) -
                              static_cast<double>(candidate.point.get(axis));
        const std::size_t next_axis = (axis + 1) % Dimensions;
        if (offset < 0) {
            search(begin, middle, next_axis, query, best, best_squared_distance, visited);
            if (offset * offset <= best_squared_distance)
                search(middle + 1, end, next_axis, query, best, best_squared_distance, visited);
        } else {
            search(middle + 1, end, next_axis, query, best, best_squared_distance, visited);
            if (offset * offset <= best_squared_distance)
                search(begin, middle, next_axis, query, best, best_squared_distance, visited);
        }
    }

    std::vector<Entry> entries;
};

#endif // LPJ_GUESS_KDTree_H
