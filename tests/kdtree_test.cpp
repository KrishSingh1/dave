///////////////////////////////////////////////////////////////////////////////////////
/// \file KDTree_test.cpp
/// \brief Tests for coordinate and value lookup.
///
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
///////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "kdtree.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

TEST_CASE("KDTree returns the nearest coordinate, value and real distance", "[KDTree]") {
    using Index = KDTree<double>;
    std::vector<Index::Entry> entries = {
        {{0, 0}, 10}, {{3, 0}, 20}, {{5, 5}, 30}, {{6, 0}, 40}, {{9, 9}, 50}
    };
    const Index index(std::move(entries));

    const auto nearest = index.nearest({4, 0});
    REQUIRE(nearest.point.get(0) == 3);
    REQUIRE(nearest.point.get(1) == 0);
    REQUIRE(nearest.value == 20);
    REQUIRE(nearest.distance == Approx(1));
    REQUIRE(nearest.visited > 0);

    REQUIRE(index.nearest({5, 5}).distance == Approx(0));
    REQUIRE(index.nearest({0, 4}).distance == Approx(4));

    Point<double, 2> origin{0, 0};
    REQUIRE(origin.squared_distance({3, 4}) == Approx(25));
    REQUIRE(origin.distance({3, 4}) == Approx(5));
}

TEST_CASE("KDTree agrees with linear search in three dimensions", "[KDTree]") {
    using Index = KDTree<int, int, 3>;
    std::vector<Index::Entry> entries;
    for (int i = 0; i < 80; ++i)
        entries.push_back({{(i * 17) % 37 - 18, (i * 11) % 41 - 20, (i * 7) % 43 - 21}, i});

    const Index index(entries);
    for (const auto& entry : entries)
        REQUIRE(index.find_exact(entry.point)->value == entry.value);
    REQUIRE(index.find_exact({100, 100, 100}) == nullptr);

    for (int x = -22; x <= 22; x += 7) {
        for (int y = -24; y <= 24; y += 8) {
            for (int z = -25; z <= 25; z += 9) {
                Point<int, 3> query{x, y, z};
                const double expected = std::min_element(entries.begin(), entries.end(),
                    [&query](const Index::Entry& a, const Index::Entry& b) {
                        return a.point.squared_distance(query) < b.point.squared_distance(query);
                    })->point.distance(query);
                const auto result = index.nearest(query);
                REQUIRE(result.point.distance(query) == Approx(expected));
                REQUIRE(result.distance == Approx(expected));
            }
        }
    }
}

TEST_CASE("KDTree resolves ties and repeated coordinates", "[KDTree]") {
    using Index = KDTree<int, int, 1>;
    const Index index({{{1}, 10}, {{-1}, 20}, {{1}, 30}});
    REQUIRE(index.size() == 2);
    REQUIRE(index.nearest({0}).point.get(0) == -1);
    REQUIRE(index.nearest({0}).value == 20);
    REQUIRE(index.nearest({1}).value == 30);
    REQUIRE(index.find_exact({1})->value == 30);
    REQUIRE(index.find_exact({2}) == nullptr);
}

TEST_CASE("KDTree handles empty data and malformed points", "[KDTree]") {
    const KDTree<double> index({});
    REQUIRE(index.empty());
    REQUIRE(index.find_exact({0, 0}) == nullptr);
    REQUIRE_THROWS_AS(index.nearest({0, 0}), std::logic_error);
    REQUIRE_THROWS_AS((Point<double, 2>{1}), std::invalid_argument);
    REQUIRE_THROWS_AS((Point<double, 2>{1, 2, 3}), std::invalid_argument);
    REQUIRE_THROWS_AS((Point<double, 2>{std::numeric_limits<double>::quiet_NaN(), 0}),
                      std::invalid_argument);
}
