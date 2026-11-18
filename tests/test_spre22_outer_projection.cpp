/***************************************************************************
 *            tests/test_spre22_outer_projection.cpp
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

/** @file tests/test_spre22_outer_projection.cpp
 *  @brief Replay two real SPRE22 enclosures through the affine outer projection.
 *
 *  This is an authentic Taylor-enclosure regression, not the original
 *  nc=3/nc=4 affine boundary cycling reproduction.
 */

#include <chrono>

#include "fixtures/spre22_affine_subset.hpp"
#include "logging/test_runner.hpp"
#include "set/affine_set.hpp"
#include "set/function_set.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {

bool inside_outer_polygon(const OuterPolygon2d& polygon, const Rational& x,
                          const Rational& y) {
    const auto& vertices = polygon.vertices();
    if (vertices.size() < 3u) return false;
    for (SizeType i = 0u; i != vertices.size(); ++i) {
        const auto& a = vertices[i];
        const auto& b = vertices[(i + 1u) % vertices.size()];
        const Rational cross =
            (b.first - a.first) * (y - a.second) -
            (b.second - a.second) * (x - a.first);
        if (cross < Rational(0)) return false;
    }
    return true;
}

} // namespace

class TestSPRE22OuterProjection : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    Void test_real_enclosure_replay() {
        const auto captured = SPRE22Fixture::affine_projection_samples();
        ARIADNE_TESTING_EXPECT_EQUAL(captured.size(), 2u);

        for (SizeType i = 0u; i != captured.size(); ++i) {
            const ValidatedAffineConstrainedImageSet affine =
                captured[i].affine_over_approximation();
            ARIADNE_TESTING_EXPECT_EQUAL(affine.dimension(), 5u);
            ARIADNE_TESTING_EXPECT_EQUAL(affine.number_of_constraints(), 1u);

            // The two real fixtures are snapshots 1 and 17 of the 48
            // captured SPRE22 enclosures. Cover both coordinate planes.
            for (const auto& axes : {std::pair<SizeType, SizeType>{1u, 2u},
                                    std::pair<SizeType, SizeType>{3u, 4u}}) {
                const auto start = std::chrono::steady_clock::now();
                const auto refined = outer_project2d_with_report(
                    affine, axes.first, axes.second, 24u, true);
                const auto elapsed = std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - start).count();

                const auto fallback = outer_project2d_with_report(
                    affine, axes.first, axes.second, 24u, false);
                ARIADNE_TESTING_EXPECT_EQUAL(refined.requested_directions, 24u);
                ARIADNE_TESTING_EXPECT_EQUAL(fallback.requested_directions, 24u);
                ARIADNE_TESTING_EXPECT_EQUAL(fallback.simplex_cuts, 0u);
                ARIADNE_TESTING_EXPECT_EQUAL(fallback.simplex_failures, 0u);
                ARIADNE_TESTING_PROVE(!refined.polygon.empty());
                ARIADNE_TESTING_PROVE(!fallback.polygon.empty());
                ARIADNE_TESTING_PROVE(refined.certified_cuts <= 24u);
                ARIADNE_TESTING_PROVE(fallback.certified_cuts <= 24u);

                // Both are mathematically conservative, and the simplex
                // certificate should not extend beyond the initial box.
                // Geometry is independently checked using rational predicates.
                bool refined_inside_fallback = true;
                for (const auto& vertex : refined.polygon.vertices()) {
                    if (!inside_outer_polygon(
                            fallback.polygon, vertex.first, vertex.second)) {
                        refined_inside_fallback = false;
                    }
                }
                ARIADNE_TESTING_PROVE(refined_inside_fallback);
                ARIADNE_TESTING_NOTE(
                    "SPRE22 real capture " << (i == 0u ? 1u : 17u)
                    << " projected axes=(" << axes.first << "," << axes.second
                    << ") simplex_24_directions_ms=" << elapsed
                    << " support_cuts=" << refined.certified_cuts
                    << " simplex_cuts=" << refined.simplex_cuts
                    << " relaxed_cuts=" << refined.relaxed_cuts
                    << " simplex_failures=" << refined.simplex_failures);
            }
        }
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_spre22_outer_projection");
    ARIADNE_TESTING_ADD_TEST(TestSPRE22OuterProjection, test_real_enclosure_replay);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
