/***************************************************************************
 *            tests/test_polytope.cpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Set.
 *
 *  Ariadne Set is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Set is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Set. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file tests/test_polytope.cpp
 *  @brief Test exact polytopes and two-dimensional vertex projections.
 */

#include <initializer_list>
#include <vector>

#include "algebra/vector.hpp"
#include "logging/test_runner.hpp"
#include "set/polytope.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {

RationalPoint make_point(std::initializer_list<Int> values) {
    Vector<Rational> coordinates(values.size(), Rational(0));
    SizeType i = 0u;
    for (auto value: values) coordinates[i++] = Rational(value);
    return RationalPoint(coordinates);
}

bool has_vertex(const Polytope& p, const RationalPoint& vertex) {
    for (const auto& candidate: p.vertices()) {
        if (candidate.dimension() != vertex.dimension()) continue;
        bool same = true;
        for (SizeType i = 0u; i != vertex.dimension(); ++i) {
            same = same && (candidate[i] == vertex[i]);
        }
        if (same) return true;
    }
    return false;
}

} // namespace

class TestPolytope : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    Void test_dimension_and_empty() {
        Polytope p(3u);
        ARIADNE_TESTING_EXPECT_EQUAL(p.dimension(), 3u);
        ARIADNE_TESTING_PROVE(p.empty());
        ARIADNE_TESTING_EXPECT_EQUAL(p.number_of_vertices(), 0u);
        p.new_vertex(make_point({1, 2, 3}));
        ARIADNE_TESTING_PROVE(!p.empty());
        ARIADNE_TESTING_EXPECT_EQUAL(p.number_of_vertices(), 1u);
        ARIADNE_TESTING_EXPECT_EXCEPTION(p.new_vertex(make_point({1, 2})));
        ARIADNE_TESTING_EXPECT_EXCEPTION(p.reduce2d());
    }

    Void test_reduction_2d() {
        Polytope p(2u, {
            make_point({0, 0}), make_point({2, 0}), make_point({2, 2}),
            make_point({0, 2}), make_point({1, 1}), make_point({2, 0}),
            make_point({1, 0}), make_point({0, 1})
        });
        p.reduce2d();
        ARIADNE_TESTING_EXPECT_EQUAL(p.number_of_vertices(), 4u);
        ARIADNE_TESTING_PROVE(has_vertex(p, make_point({0, 0})));
        ARIADNE_TESTING_PROVE(has_vertex(p, make_point({2, 0})));
        ARIADNE_TESTING_PROVE(has_vertex(p, make_point({2, 2})));
        ARIADNE_TESTING_PROVE(has_vertex(p, make_point({0, 2})));
        ARIADNE_TESTING_EXPECT_EQUAL(p.vertex(0u)[0], Rational(0));
        ARIADNE_TESTING_EXPECT_EQUAL(p.vertex(0u)[1], Rational(0));
        // Repeating reduction must preserve the vertices and their coordinates.
        p.reduce2d();
        ARIADNE_TESTING_EXPECT_EQUAL(p.number_of_vertices(), 4u);
        ARIADNE_TESTING_PROVE(has_vertex(p, make_point({0, 0})));
        ARIADNE_TESTING_PROVE(has_vertex(p, make_point({2, 2})));
    }

    Void test_degenerate_2d() {
        Polytope empty(2u);
        empty.reduce2d();
        ARIADNE_TESTING_PROVE(empty.empty());
        Polytope segment(2u, {
            make_point({3, 3}), make_point({0, 0}),
            make_point({2, 2}), make_point({1, 1}), make_point({3, 3})
        });
        segment.reduce2d();
        ARIADNE_TESTING_EXPECT_EQUAL(segment.number_of_vertices(), 2u);
        ARIADNE_TESTING_PROVE(has_vertex(segment, make_point({0, 0})));
        ARIADNE_TESTING_PROVE(has_vertex(segment, make_point({3, 3})));
    }

    Void test_projection_and_hull() {
        Polytope p(3u, {
            make_point({0, 3, 2}), make_point({0, 4, 0}),
            make_point({1, 3, 0}), make_point({1, 4, 2}),
            make_point({0, 3, 1})
        });
        const Polytope xy = project2d(p, 2u, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(xy.dimension(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(xy.number_of_vertices(), 4u);
        ARIADNE_TESTING_PROVE(has_vertex(xy, make_point({0, 3})));
        ARIADNE_TESTING_PROVE(has_vertex(xy, make_point({2, 4})));
        const Polytope union_hull = convex_hull(xy, Polytope(2u, {make_point({3, 5})}));
        ARIADNE_TESTING_PROVE(has_vertex(union_hull, make_point({3, 5})));
        ARIADNE_TESTING_EXPECT_EXCEPTION(project2d(p, 0u, 0u));
        ARIADNE_TESTING_EXPECT_EXCEPTION(project2d(p, 0u, 3u));
        ARIADNE_TESTING_EXPECT_EXCEPTION(convex_hull(p, xy));
    }

    Void test_generator_mean() {
        Polytope p(2u, {make_point({0, 0}), make_point({1, 0}), make_point({1, 1})});
        const RationalPoint mean = mean_of_generators(p);
        ARIADNE_TESTING_PROVE(mean[0] == Rational(2) / Rational(3));
        ARIADNE_TESTING_PROVE(mean[1] == Rational(1) / Rational(3));
        ARIADNE_TESTING_EXPECT_EXCEPTION(mean_of_generators(Polytope(2u)));
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_polytope");
    ARIADNE_TESTING_ADD_TEST(TestPolytope, test_dimension_and_empty);
    ARIADNE_TESTING_ADD_TEST(TestPolytope, test_reduction_2d);
    ARIADNE_TESTING_ADD_TEST(TestPolytope, test_degenerate_2d);
    ARIADNE_TESTING_ADD_TEST(TestPolytope, test_projection_and_hull);
    ARIADNE_TESTING_ADD_TEST(TestPolytope, test_generator_mean);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
