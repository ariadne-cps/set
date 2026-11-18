/***************************************************************************
 *            tests/test_spre22_critical.cpp
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

/** @file tests/test_spre22_critical.cpp
 *  @brief Replay a genuine SPRE22 plotting input at the Set level.
 *
 *  The drawing replay was reported to stall at original index 11591.
 *  This does not call drawing's legacy boundary algorithm. It tests that
 *  the independent certified Set projection terminates and remains sound.
 */

#include <chrono>
#include <cstdint>

#include "ariadne-testing.hpp"
#include "fixtures/spre22_critical_11591.hpp"
#include "fixtures/spre22_snapshot_reader.hpp"
#include "logging/test_runner.hpp"
#include "set/affine_set.hpp"
#include "set/function_set.hpp"

using namespace Ariadne;

namespace {

// Exact convex-polygon predicate for checking a *relative* inclusion.
// The proof of soundness of each projection comes from its support bounds.
bool contains_vertex(const OuterPolygon2d& outer,
                     const Rational& x, const Rational& y) {
    const auto& vertices = outer.vertices();
    if (vertices.empty()) return false;
    if (vertices.size() == 1u)
        return vertices[0].first == x && vertices[0].second == y;
    if (vertices.size() == 2u) {
        const auto& a = vertices[0];
        const auto& b = vertices[1];
        const Rational cross = (b.first - a.first) * (y - a.second)
                             - (b.second - a.second) * (x - a.first);
        const Rational dot = (x - a.first) * (x - b.first)
                           + (y - a.second) * (y - b.second);
        return cross == Rational(0) && dot <= Rational(0);
    }
    for (std::size_t i = 0u; i != vertices.size(); ++i) {
        const auto& a = vertices[i];
        const auto& b = vertices[(i + 1u) % vertices.size()];
        if ((b.first - a.first) * (y - a.second)
          - (b.second - a.second) * (x - a.first) < Rational(0))
            return false;
    }
    return true;
}

}

class TestSPRE22Critical : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    Void test_11591_critical_projection() {
        const auto capture = SPRE22Fixture::decode_critical_11591(
            SPRE22Fixture::encoded_record_11591);
        ARIADNE_TESTING_EXPECT_EQUAL(capture.set.number_of_parameters(), 26u);
        ARIADNE_TESTING_EXPECT_EQUAL(capture.set.number_of_constraints(), 4u);
        ARIADNE_TESTING_EXPECT_EQUAL(capture.set.dimension(), 9u);
        ARIADNE_TESTING_EXPECT_EQUAL(capture.x, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(capture.y, 2u);

        const ValidatedAffineConstrainedImageSet affine =
            capture.set.affine_over_approximation();
        const SizeType ndirections = 24u;
        const auto begin = std::chrono::steady_clock::now();
        const auto optimised = outer_project2d_with_report(
            affine, capture.x, capture.y, ndirections, true);
        const double elapsed_ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - begin).count();
        const auto relaxed = outer_project2d_with_report(
            affine, capture.x, capture.y, ndirections, false);

        ARIADNE_TESTING_EXPECT_EQUAL(optimised.requested_directions, ndirections);
        ARIADNE_TESTING_EXPECT_EQUAL(relaxed.requested_directions, ndirections);
        ARIADNE_TESTING_PROVE(optimised.certified_cuts <= ndirections);
        ARIADNE_TESTING_PROVE(relaxed.certified_cuts <= ndirections);
        ARIADNE_TESTING_EXPECT_EQUAL(relaxed.simplex_cuts, 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(relaxed.simplex_failures, 0u);
        ARIADNE_TESTING_PROVE(!optimised.polygon.empty());
        ARIADNE_TESTING_PROVE(!relaxed.polygon.empty());

        // Both certified polygons must be inside the original validated
        // bounding box. This is a geometric check, not a claim that the
        // optimiser always improves over zero dual multipliers.
        const auto box = affine.bounding_box();
        const auto initial = certified_support_polygon_2d(
            {box[capture.x].lower_bound().raw().get_d(),
             box[capture.x].upper_bound().raw().get_d(),
             box[capture.y].lower_bound().raw().get_d(),
             box[capture.y].upper_bound().raw().get_d()}, {});
        bool inside_box = true;
        for (const auto& vertex : optimised.polygon.vertices())
            if (!contains_vertex(initial, vertex.first, vertex.second))
                inside_box = false;
        for (const auto& vertex : relaxed.polygon.vertices())
            if (!contains_vertex(initial, vertex.first, vertex.second))
                inside_box = false;
        ARIADNE_TESTING_PROVE(inside_box);
        ARIADNE_TESTING_NOTE("SPRE22 index=11591 x=1 y=2 parameters=26 constraints=4 "
                             << "projection_24_directions_ms=" << elapsed_ms
                             << " simplex_cuts=" << optimised.simplex_cuts
                             << " relaxed_cuts=" << optimised.relaxed_cuts
                             << " simplex_failures=" << optimised.simplex_failures);
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_spre22_critical");
    ARIADNE_TESTING_ADD_TEST(TestSPRE22Critical, test_11591_critical_projection);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
