/***************************************************************************
 *            tests/test_support_polygon.cpp
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

/** @file tests/test_support_polygon.cpp
 *  @brief Exact certified halfplane clipping and degenerate intersections.
 */

#include <array>
#include <limits>
#include <stdexcept>
#include <vector>

#include "logging/test_runner.hpp"
#include "set/exact_support_polygon.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {

Rational exact(double value) {
    return Rational(ExactDouble(value));
}

bool includes(const Polytope& polygon, const Rational& x, const Rational& y) {
    const auto& vertices = polygon.vertices();
    if (vertices.empty()) return false;
    if (vertices.size() == 1u) {
        if (vertices[0][0] != x) return false;
        return decide(vertices[0][1] == y);
    }
    if (vertices.size() == 2u) {
        const auto& a = vertices[0];
        const auto& b = vertices[1];
        const Rational cross = (b[0] - a[0]) * (y - a[1]) -
                               (b[1] - a[1]) * (x - a[0]);
        if (cross != Rational(0)) return false;
        const Rational dot = (x - a[0]) * (x - b[0]) +
                             (y - a[1]) * (y - b[1]);
        return dot <= Rational(0);
    }
    for (std::size_t i = 0u; i < vertices.size(); ++i) {
        const auto& a = vertices[i];
        const auto& b = vertices[(i + 1u) % vertices.size()];
        const Rational cross = (b[0] - a[0]) * (y - a[1]) -
                               (b[1] - a[1]) * (x - a[0]);
        if (cross < Rational(0)) return false;
    }
    return true;
}

}

class TestSupportPolygon : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    Void test_exact_clipping_and_inclusion() {
        const std::vector<CertifiedSupportHalfplane2d> cuts = {
            {3.0, 2.0, 0.3},
            {2.0, -3.0, 0.4},
            {-1.0, 0.0, 0.75},
            {0.0, 1.0, 0.6}
        };
        const auto polygon = certified_support_polygon_2d({-1.0, 1.0, -1.0, 1.0}, cuts);
        ARIADNE_TESTING_PROVE(polygon.number_of_vertices() >= 3u);

        // Non-dyadic intersection of two dyadic lines, checked exactly.
        const Rational x = (Rational(3) * exact(0.3) + Rational(2) * exact(0.4)) / Rational(13);
        const Rational y = (Rational(2) * exact(0.3) - Rational(3) * exact(0.4)) / Rational(13);
        ARIADNE_TESTING_PROVE(includes(polygon, x, y));

        bool all_contained = true;
        for (int ix = -16; ix <= 16; ++ix) {
            for (int iy = -16; iy <= 16; ++iy) {
                const Rational px = exact(static_cast<double>(ix) / 16.0);
                const Rational py = exact(static_cast<double>(iy) / 16.0);
                bool feasible = true;
                for (const auto& cut : cuts) {
                    if (exact(cut.nx) * px + exact(cut.ny) * py > exact(cut.upper)) {
                        feasible = false;
                        break;
                    }
                }
                if (feasible && !includes(polygon, px, py)) all_contained = false;
            }
        }
        ARIADNE_TESTING_PROVE(all_contained);
        // The result must be genuinely tighter than the initial box.
        bool has_cut_vertex = false;
        for (const auto& p : polygon.vertices()) {
            if (p[0] != exact(-1.0) && p[0] != exact(1.0) &&
                p[1] != exact(-1.0) && p[1] != exact(1.0)) {
                has_cut_vertex = true;
            }
        }
        ARIADNE_TESTING_PROVE(has_cut_vertex);
    }

    Void test_degenerate_intersections() {
        const auto empty = certified_support_polygon_2d(
            {-1.0, 1.0, -1.0, 1.0},
            {{1.0, 0.0, -0.5}, {-1.0, 0.0, -0.5}});
        ARIADNE_TESTING_PROVE(empty.empty());

        const auto singleton = certified_support_polygon_2d(
            {-1.0, 1.0, -1.0, 1.0},
            {{1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0},
             {0.0, 1.0, 0.0}, {0.0, -1.0, 0.0}});
        ARIADNE_TESTING_EXPECT_EQUAL(singleton.number_of_vertices(), 1u);
        ARIADNE_TESTING_PROVE(includes(singleton, Rational(0), Rational(0)));

        const auto segment = certified_support_polygon_2d(
            {-1.0, 1.0, -1.0, 1.0},
            {{1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}});
        ARIADNE_TESTING_EXPECT_EQUAL(segment.number_of_vertices(), 2u);
        ARIADNE_TESTING_PROVE(includes(segment, Rational(0), Rational(0)));
        ARIADNE_TESTING_PROVE(includes(segment, Rational(0), Rational(1)));
    }

    Void test_invalid_inputs_and_zero_cuts() {
        const std::vector<CertifiedSupportHalfplane2d> none;
        const auto box = certified_support_polygon_2d({-2.0, 2.0, -3.0, 3.0}, none);
        ARIADNE_TESTING_EXPECT_EQUAL(box.number_of_vertices(), 4u);
        ARIADNE_TESTING_PROVE(includes(box, Rational(0), Rational(0)));
        ARIADNE_TESTING_EXPECT_EXCEPTION(
            certified_support_polygon_2d({1.0, -1.0, -1.0, 1.0}, none));
        ARIADNE_TESTING_EXPECT_EXCEPTION(
            certified_support_polygon_2d({-1.0, 1.0, -1.0, 1.0}, {{1.0, 0.0, std::numeric_limits<double>::infinity()}}));
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_support_polygon");
    ARIADNE_TESTING_ADD_TEST(TestSupportPolygon, test_exact_clipping_and_inclusion);
    ARIADNE_TESTING_ADD_TEST(TestSupportPolygon, test_degenerate_intersections);
    ARIADNE_TESTING_ADD_TEST(TestSupportPolygon, test_invalid_inputs_and_zero_cuts);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
