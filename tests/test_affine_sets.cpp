/***************************************************************************
 *            tests/test_affine_sets.cpp
 *
 *  Copyright  2009-20  Pieter Collins
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

/** @file tests/test_affine_sets.cpp
 *  @brief Test affine constrained image sets and validated geometry.
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>

#include "function/function.hpp"
#include "geometry/box.hpp"
#include "logging/test_runner.hpp"
#include "set/affine_set.hpp"
#include "set/function_set.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;
using namespace std;

namespace Ariadne {

template<class X>
AffineConstraint<X> operator<=(const Dyadic& l, const Affine<X>& a) {
    return AffineConstraint<X>({l, a.value().precision()}, a, {+infty, a.value().precision()});
}
template<class X>
AffineConstraint<X> operator<=(const Affine<X>& a, const Dyadic& u) {
    return AffineConstraint<X>(-infty, a, u);
}
template<class X>
AffineConstraint<X> operator<=(const AffineConstraint<X>& ac, const Dyadic& u) {
    ARIADNE_ASSERT(decide(ac.upper_bound() == infty));
    return AffineConstraint<X>(ac.lower_bound(), ac.function(), {u, ac.lower_bound().precision()});
}
template<class X>
AffineConstraint<X> operator==(const Affine<X>& a, const Dyadic& u) {
    return AffineConstraint<X>(u, a, u);
}
} // namespace Ariadne

class TestAffineSet : public Ariadne::Testing::ContextFixture {
  private:
    ExactBoxType D;
    Vector<Affine<FloatDPBounds>> x;
    ValidatedAffineConstrainedImageSet set;

  public:
    explicit TestAffineSet(Ariadne::Testing::Context& context) : ContextFixture(context),
                                                                 set(ExactBoxType::unit_box(2), Matrix<FloatDP>(2, 2, dp), Vector<FloatDP>(2, dp)) {}

    Void test_basic_coverage() {
        RealBox real_box(2, RealInterval(-1, 1));
        ValidatedAffineConstrainedImageSet from_real(real_box);
        ARIADNE_TESTING_EXPECT_EQUAL(from_real.dimension(), 2u);
        ARIADNE_TESTING_PROVE(definitely(from_real.is_bounded()));
        ARIADNE_TESTING_EXPECT_EQUAL(from_real.number_of_constraints(), 0u);

        ExactBoxType exact_box({
            {-1.0_x, +1.0_x},
            {-2.0_x, +2.0_x}
        });
        ValidatedAffineConstrainedImageSet from_exact(exact_box);
        ARIADNE_TESTING_EXPECT_EQUAL(from_exact.domain(), exact_box);
        ARIADNE_TESTING_PROVE(definitely(from_exact.inside(ExactBoxType{
            {-2.0_x, +2.0_x},
            {-3.0_x, +3.0_x}
        })));

        Vector<ValidatedAffineModelDP> models = ValidatedAffineModelDP::scalings(ExactBoxType::unit_box(2), dp);
        ValidatedAffineConstrainedImageSet from_models(models);
        List<ValidatedAffineModelConstraintDP> constraints;
        ValidatedAffineConstrainedImageSet from_models_and_constraints(models, constraints);
        ARIADNE_TESTING_EXPECT_EQUAL(from_models.number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(from_models_and_constraints.number_of_constraints(), 0u);
        // The representation exposed to downstream algorithms retains its models.
        ARIADNE_TESTING_EXPECT_EQUAL(from_models.space_models().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(from_models_and_constraints.constraint_models().size(), 0u);
        std::unique_ptr<ValidatedAffineConstrainedImageSet> cloned(from_models.clone());
        ARIADNE_TESTING_EXPECT_EQUAL(cloned->dimension(), from_models.dimension());
        ARIADNE_TESTING_EXPECT_EQUAL(cloned->number_of_parameters(), from_models.number_of_parameters());

        auto y = ValidatedScalarMultivariateFunction::coordinates(2);
        ValidatedVectorMultivariateFunction h({y[0] + y[1], y[0] - y[1]});
        auto transformed = image(from_exact, h);
        ARIADNE_TESTING_EXPECT_EQUAL(transformed.dimension(), 2u);
    }

    Void test_pure_constraint() {
        D = ExactBoxType::unit_box(2);
        x = Affine<FloatDPBounds>::variables(2, dp);

        set = ValidatedAffineConstrainedImageSet(D, {x[0] + 2, x[0] + x[1]});

        set.new_parameter_constraint(x[0] + 0.5_x * x[1] + 0.75_x <= 0.0_x);
        set.new_parameter_constraint(-0.5_x * x[0] + 1.5_x * x[1] + 0.5_x <= 0.0_x);

        ARIADNE_TESTING_NOTE("set = " << (set));

        set = ValidatedAffineConstrainedImageSet(D, {x[0] + x[1] + 2, x[0] + x[1] + 2});
        set.new_parameter_constraint(-0.5_x * x[0] - 1.5_x * x[1] + 0.5_x <= 0.0_x);

        ARIADNE_TESTING_NOTE("set = " << (set));
    }

    Void test_nonunit_domain() {
        D = ExactBoxType({
            { 0.0_x, 2.0_x},
            {-3.0_x, 1.0_x}
        });
        x = Affine<FloatDPBounds>::variables(2, dp);

        set = ValidatedAffineConstrainedImageSet(D, {x[0], x[1]}, {x[0] + x[1] <= 1.0_x});
        ARIADNE_TESTING_EXPECT_EQUAL(set.domain(), D);
        ARIADNE_TESTING_EXPECT_EQUAL(set.number_of_constraints(), 1u);
        ARIADNE_TESTING_EXPECT_SAME(set.bounding_box(), UpperBoxType({
                                                            { 0.0_x, 2.0_x},
                                                            {-3.0_x, 1.0_x}
        }));

        ValidatedAffineConstrainedImageSet unconstrained(D, {x[0], x[1]});
        ARIADNE_TESTING_EXPECT_SAME(unconstrained.bounding_box(), UpperBoxType({
                                                                      { 0.0_x, 2.0_x},
                                                                      {-3.0_x, 1.0_x}
        }));
    }

    Void test_constrained_image() {
        D = ExactBoxType::unit_box(3);
        x = Affine<FloatDPBounds>::variables(3, dp);

        set = ValidatedAffineConstrainedImageSet(D, {2 * x[0] + 3 * x[1] + x[2], x[0] + x[1] - 2});

        set.new_parameter_constraint(1.0_x * x[0] + 0.5_x * x[1] + 0.25_x * x[2] <= 0.75_x);
        set.new_parameter_constraint(-2.0_x * x[0] - 3.0_x * x[1] - 1.0_x * x[2] <= 3.00_x);

        ARIADNE_TESTING_NOTE("set = " << (set));
    }

    Void test_separated() {
        D = ExactBoxType::unit_box(3);
        x = Affine<FloatDPBounds>::variables(3, dp);
        set = ValidatedAffineConstrainedImageSet(D, {-0.9375_x + 0.0625_x * x[0] + 0.5_x * x[1], 0.87890625_x - 0.1171875_x * x[0] + x[1] + 0.00390625_x * x[2]});
        set.new_parameter_constraint(-1.1875_x + 0.0625_x * x[0] + x[1] <= 0.0_x);
        ARIADNE_TESTING_NOTE("set = " << (set));
        ExactBoxType cell1({
            {  -1.0_x, -0.9375_x},
            {0.8125_x,   0.875_x}
        });
        ARIADNE_TESTING_NOTE("Constructed ExactBoxType cell1 = " << cell1); // subset
        ExactBoxType cell2({
            {-0.5625_x, -0.50_x},
            { 1.4375_x,   1.5_x}
        });
        ARIADNE_TESTING_NOTE("Constructed ExactBoxType cell2 = " << cell2); // overlaps
        ExactBoxType cell3({
            {-0.875_x,     -0.75_x},
            { 0.625_x, 0.7578125_x}
        });
        ARIADNE_TESTING_NOTE("Constructed ExactBoxType cell3 = " << cell3); // touches at (-0.875_x,0.7578125_x);
        ExactBoxType cell4({
            {-1.1875_x, -1.125_x},
            { 0.0625_x,  0.125_x}
        });
        ARIADNE_TESTING_NOTE("Constructed ExactBoxType cell4 = " << cell4); // almost touches
        ExactBoxType cell5({
            {-0.9375_x, -0.875_x},
            { 0.4375_x,    0.5_x}
        });
        ARIADNE_TESTING_NOTE("Constructed ExactBoxType cell5 = " << cell5); // disjoint
        ExactBoxType cell6({
            {-1.5_x, -1.375_x},
            { 0.5_x,  0.625_x}
        });
        ARIADNE_TESTING_NOTE("Constructed ExactBoxType cell6 = " << cell6); // disjoint; regression test

        ARIADNE_TESTING_PROVE(!definitely(set.separated(cell1)));
        ARIADNE_TESTING_PROVE(!definitely(set.separated(cell2)));
        ARIADNE_TESTING_PROVE(!definitely(set.separated(cell3)));
        // ARIADNE_TESTING_PROVE(!definitely(affine_set.overlaps(cell3)));
        ARIADNE_TESTING_PROVE(possibly(set.separated(cell4)));
        ARIADNE_TESTING_PROVE(definitely(set.separated(cell5)));
        ARIADNE_TESTING_PROVE(definitely(set.separated(cell6)));

        // The six separation outcomes above are the regression checks.
        // Rendering them here duplicated draw coverage already exercised by test_draw().
    }

    Void test_empty() {
        D = ExactBoxType::unit_box(2);
        x = Affine<FloatDPBounds>::variables(2, dp);

        ARIADNE_TESTING_NOTE("characteristics = " << get_characteristics(0.125_x * x[0] + 1.0_x));
        set = ValidatedAffineConstrainedImageSet(D, {0.125_x * x[0] + 1.0_x, 0.125_x * x[1] + 1.0_x});
        set.new_parameter_constraint(-1.0_x * x[1] <= -2.0_x);

        ARIADNE_TESTING_NOTE("set = " << (set));
        ARIADNE_TESTING_PROVE(set.is_empty());
    }

    Void test_logging_branches() {
        const unsigned int old_verbosity = Logger::instance().configuration().verbosity();
        Logger::instance().configuration().set_verbosity(10u);

        D = ExactBoxType::unit_box(2);
        x = Affine<FloatDPBounds>::variables(2, dp);
        set = ValidatedAffineConstrainedImageSet(D, {x[0], x[1]});

        ARIADNE_TESTING_PROVE(!definitely(set.separated(ExactBoxType{
            {-0.5_x, +0.5_x},
            {-0.5_x, +0.5_x}
        })));

        List<Pair<double, double>> vertices = set.boundary(0, 1);
        ARIADNE_TESTING_PROVE(!vertices.empty());

        // Exercise verbose tracing for the constrained simplex traversal,
        // whose per-pivot messages are not emitted for unconstrained zonotopes.
        set.new_parameter_constraint(x[0] <= 0.0_x);
        vertices = set.boundary(0, 1);
        ARIADNE_TESTING_PROVE(vertices.size() >= 4u);
        for (const auto& vertex: vertices) {
            ARIADNE_TESTING_PROVE(vertex.first <= 1e-6);
            ARIADNE_TESTING_PROVE(vertex.first >= -1.000001);
        }

        Logger::instance().configuration().set_verbosity(old_verbosity);
    }

    Void test_invalid() {
        {
            D = ExactBoxType::unit_box(1);
            x = Affine<FloatDPBounds>::variables(2, dp);
            ARIADNE_TESTING_EXPECT_EXCEPTION(ValidatedAffineConstrainedImageSet(D, {0.125_x * x[0] + 1.0_x, 0.125_x * x[1] + 1.0_x}));
        }
        {
            D = ExactBoxType::unit_box(2);
            x = Affine<FloatDPBounds>::variables(2, dp);
            auto a = Affine<FloatDPBounds>::variables(1, dp);
            ARIADNE_TESTING_EXPECT_EXCEPTION(ValidatedAffineConstrainedImageSet(D, {0.125_x * x[0] + 1.0_x, 0.125_x * x[1] + 1.0_x}, {a[0] <= -0.25_x}));

            set = ValidatedAffineConstrainedImageSet(D, {x[0], x[1]});
            ARIADNE_TESTING_EXPECT_EXCEPTION(set.new_parameter_constraint(a[0] <= -0.25_x));
            ARIADNE_TESTING_EXPECT_EXCEPTION(set.separated(ExactBoxType::unit_box(1)));
            ARIADNE_TESTING_EXPECT_EXCEPTION(set.inside(ExactBoxType::unit_box(1)));

            Vector<ValidatedAffineModelDP> bad_models = ValidatedAffineModelDP::scalings(ExactBoxType::unit_box(1), dp);
            List<ValidatedAffineModelConstraintDP> no_constraints;
            ARIADNE_TESTING_EXPECT_EXCEPTION(ValidatedAffineConstrainedImageSet(D, bad_models, no_constraints));

            ValidatedAffineModelConstraintDP bad_constraint(-1.0_x, bad_models[0], +1.0_x);
            List<ValidatedAffineModelConstraintDP> bad_constraints;
            bad_constraints.append(bad_constraint);
            Vector<ValidatedAffineModelDP> good_models = ValidatedAffineModelDP::scalings(D, dp);
            ARIADNE_TESTING_EXPECT_EXCEPTION(ValidatedAffineConstrainedImageSet(D, good_models, bad_constraints));
            ARIADNE_TESTING_EXPECT_EXCEPTION(set.new_constraint(bad_constraint));

            ValidatedAffineConstrainedImageSet empty_set(ExactBoxType(0u));
            ValidatedAffineModelDP zero_model(0u, dp);
            ValidatedAffineModelConstraintDP zero_model_constraint(-1.0_x, zero_model, +1.0_x);
            Affine<FloatDPBounds> zero_affine(0u, dp);
            ARIADNE_TESTING_EXPECT_EXCEPTION(empty_set.new_constraint(zero_model_constraint));
            ARIADNE_TESTING_EXPECT_EXCEPTION(empty_set.new_parameter_constraint(zero_affine <= 0.0_x));
            ARIADNE_TESTING_EXPECT_EXCEPTION(empty_set.number_of_parameters());

            Vector<Affine<FloatDPBounds>> empty_affine_models;
            Vector<ValidatedAffineModelDP> empty_models;
            ARIADNE_TESTING_EXPECT_EXCEPTION(ValidatedAffineConstrainedImageSet(D, empty_affine_models));
            ARIADNE_TESTING_EXPECT_EXCEPTION(ValidatedAffineConstrainedImageSet(D, empty_models, no_constraints));
            ARIADNE_TESTING_EXPECT_EXCEPTION(ValidatedAffineConstrainedImageSet{empty_models});

            Matrix<FloatDP> bad_matrix(1u, 1u, dp);
            Vector<FloatDP> bad_offset(2u, dp);
            ARIADNE_TESTING_EXPECT_EXCEPTION(ValidatedAffineConstrainedImageSet(D, bad_matrix, bad_offset));
            Matrix<FloatDP> empty_matrix(0u, 0u, dp);
            Vector<FloatDP> empty_offset(0u, dp);
            ARIADNE_TESTING_EXPECT_EXCEPTION(ValidatedAffineConstrainedImageSet(ExactBoxType(0u), empty_matrix, empty_offset));
        }
    }

    Void test_constrained_boundary() {
        // A square intersected with x<=0 is a rectangle [-1,0] x [-1,1].
        // In contrast to a zonotope, its boundary requires pivoting along
        // the affine inequality's slack-variable constraints.
        D = ExactBoxType::unit_box(2);
        x = Affine<FloatDPBounds>::variables(2, dp);
        set = ValidatedAffineConstrainedImageSet(D, {x[0], x[1]}, {x[0] <= 0.0_x});
        ARIADNE_TESTING_EXPECT_EQUAL(set.number_of_constraints(), 1u);

        List<Pair<double, double>> vertices = set.boundary(0, 1);
        ARIADNE_TESTING_PROVE(vertices.size() >= 4u);
        double min_x = +2.0, max_x = -2.0, min_y = +2.0, max_y = -2.0;
        for (const auto& vertex: vertices) {
            min_x = std::min(min_x, vertex.first);
            max_x = std::max(max_x, vertex.first);
            min_y = std::min(min_y, vertex.second);
            max_y = std::max(max_y, vertex.second);
        }
        ARIADNE_TESTING_PROVE(std::abs(min_x + 1.0) < 1e-6);
        ARIADNE_TESTING_PROVE(std::abs(max_x) < 1e-6);
        ARIADNE_TESTING_PROVE(std::abs(min_y + 1.0) < 1e-6);
        ARIADNE_TESTING_PROVE(std::abs(max_y - 1.0) < 1e-6);
    }

    Void test_outer_projection_enclosure() {
        // Certified outer projection is independent of the potentially
        // cycling exact-boundary simplex walk.
        D = ExactBoxType::unit_box(2);
        x = Affine<FloatDPBounds>::variables(2, dp);
        set = ValidatedAffineConstrainedImageSet(D, {x[0], x[1]}, {x[0] <= 0.0_x});
        const auto polygon = outer_project2d(set, 0u, 1u, 24u);
        ARIADNE_TESTING_PROVE(polygon.number_of_vertices() >= 3u);

        const auto contains = [&](const Rational& px, const Rational& py) {
            const auto& points = polygon.vertices();
            if (points.size() < 3u) return false;
            for (std::size_t i = 0u; i != points.size(); ++i) {
                const auto& a = points[i];
                const auto& b = points[(i + 1u) % points.size()];
                const Rational cross =
                    (b[0] - a[0]) * (py - a[1]) -
                    (b[1] - a[1]) * (px - a[0]);
                if (cross < Rational(0)) return false;
            }
            return true;
        };
        bool all_contained = true;
        for (int ix = -4; ix <= 0; ++ix) {
            for (int iy = -4; iy <= 4; ++iy) {
                const Rational px = Rational(ix) / Rational(4);
                const Rational py = Rational(iy) / Rational(4);
                if (!contains(px, py)) all_contained = false;
            }
        }
        ARIADNE_TESTING_PROVE(all_contained);
        ARIADNE_TESTING_EXPECT_EXCEPTION(outer_project2d(set, 0u, 0u));
        ARIADNE_TESTING_EXPECT_EXCEPTION(outer_project2d(set, 0u, 2u));
        ARIADNE_TESTING_EXPECT_EXCEPTION(outer_project2d(set, 0u, 1u, 0u));

        // No constraints: support bounds must still enclose a zonotope.
        set = ValidatedAffineConstrainedImageSet(D, {x[0], x[1]});
        const auto unconstrained = outer_project2d(set, 0u, 1u, 16u);
        ARIADNE_TESTING_PROVE(unconstrained.number_of_vertices() >= 3u);
        const auto& p = unconstrained.vertices();
        bool corners_inside = true;
        for (int sx : {-1, 1}) {
            for (int sy : {-1, 1}) {
                bool inside = true;
                for (std::size_t i = 0u; i != p.size(); ++i) {
                    const auto& a = p[i];
                    const auto& b = p[(i + 1u) % p.size()];
                    const Rational px(sx), py(sy);
                    if ((b[0] - a[0]) * (py - a[1]) -
                        (b[1] - a[1]) * (px - a[0]) < Rational(0)) {
                        inside = false;
                        break;
                    }
                }
                if (!inside) corners_inside = false;
            }
        }
        ARIADNE_TESTING_PROVE(corners_inside);
    }

    Void test_outer_projection_multiple_constraints() {
        // Synthetic three- and four-guard regressions; these do not replace
        // the original SPRE22 nc=3/nc=4 enclosure snapshots.
        D = ExactBoxType::unit_box(2);
        x = Affine<FloatDPBounds>::variables(2, dp);
        for (SizeType nc = 3u; nc != 5u; ++nc) {
            set = ValidatedAffineConstrainedImageSet(
                D, {x[0], x[1]},
                {x[0] + x[1] <= 1.0_x,
                 x[0] - x[1] <= 1.0_x,
                 -1.0_x * x[0] + x[1] <= 1.0_x});
            if (nc == 4u) {
                set.new_parameter_constraint(-1.0_x * x[0] - x[1] <= 1.0_x);
            }
            ARIADNE_TESTING_EXPECT_EQUAL(set.number_of_constraints(), nc);
            const auto polygon = outer_project2d(set, 0u, 1u, 16u);
            ARIADNE_TESTING_PROVE(polygon.number_of_vertices() >= 3u);
            const auto& points = polygon.vertices();
            bool all_contained = true;
            for (int ix = -4; ix <= 4; ++ix) {
                for (int iy = -4; iy <= 4; ++iy) {
                    if (!(ix + iy <= 4 && ix - iy <= 4 &&
                          -ix + iy <= 4 && (nc == 3u || -ix - iy <= 4))) continue;
                    const Rational px = Rational(ix) / Rational(4);
                    const Rational py = Rational(iy) / Rational(4);
                    for (std::size_t i = 0u; i != points.size(); ++i) {
                        const auto& a = points[i];
                        const auto& b = points[(i + 1u) % points.size()];
                        const Rational cross =
                            (b[0] - a[0]) * (py - a[1]) -
                            (b[1] - a[1]) * (px - a[0]);
                        if (cross < Rational(0)) all_contained = false;
                    }
                }
            }
            ARIADNE_TESTING_PROVE(all_contained);
        }
    }

    Void test_outer_projection_fallback_report() {
        D = ExactBoxType::unit_box(2);
        x = Affine<FloatDPBounds>::variables(2, dp);
        set = ValidatedAffineConstrainedImageSet(D, {x[0], x[1]}, {x[0] <= 0.0_x});

        // Avoid the simplex altogether. Every direction still has a sound
        // support certificate using the zero dual multiplier.
        const SizeType directions = 12u;
        const auto report = outer_project2d_with_report(
            set, 0u, 1u, directions, false);
        ARIADNE_TESTING_EXPECT_EQUAL(report.requested_directions, directions);
        ARIADNE_TESTING_EXPECT_EQUAL(report.certified_cuts, directions);
        ARIADNE_TESTING_EXPECT_EQUAL(report.simplex_cuts, 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(report.relaxed_cuts, directions);
        ARIADNE_TESTING_EXPECT_EQUAL(report.simplex_failures, 0u);
        ARIADNE_TESTING_PROVE(report.used_relaxation());
        // Twelve valid cuts can all be redundant: counting them is not
        // equivalent to proving that the box geometry became smaller.
        ARIADNE_TESTING_PROVE(report.box_only());
        ARIADNE_TESTING_EXPECT_EQUAL(report.polygon.number_of_vertices(), 4u);

        // An infeasible constraint does not license claiming mathematical
        // emptiness from an approximate LP solve or from a bounding box.
        set = ValidatedAffineConstrainedImageSet(D, {x[0], x[1]});
        set.new_parameter_constraint(-1.0_x * x[0] <= -2.0_x);
        ARIADNE_TESTING_PROVE(definitely(set.is_empty()));
        const auto empty_source = outer_project2d_with_report(
            set, 0u, 1u, directions, false);
        ARIADNE_TESTING_EXPECT_EQUAL(empty_source.relaxed_cuts, directions);
        ARIADNE_TESTING_PROVE(!empty_source.polygon.empty());
    }

    Void test_outer_projection_model_errors_and_guard_bounds() {
        // Both the image and the affine guard carry explicit nonzero errors.
        // The constructed points stay feasible for every permitted guard
        // error; the polygon must enclose the full image error interval.
        D = ExactBoxType::unit_box(2);
        Vector<ValidatedAffineModelDP> models =
            ValidatedAffineModelDP::scalings(D, dp);
        models[0].error() += FloatDPError(ExactDouble(0.125), dp);
        models[1].error() += FloatDPError(ExactDouble(0.25), dp);

        ValidatedAffineModelDP guard = models[0];
        guard.error() += FloatDPError(ExactDouble(0.125), dp);
        List<ValidatedAffineModelConstraintDP> guards;
        const FloatDPBounds low(ExactDouble(-2.0), ExactDouble(-1.875), dp);
        const FloatDPBounds high(ExactDouble(-0.25), ExactDouble(-0.125), dp);
        guards.append(ValidatedAffineModelConstraintDP(low, guard, high));
        set = ValidatedAffineConstrainedImageSet(D, models, guards);

        const auto result = outer_project2d_with_report(set, 0u, 1u, 24u, false);
        ARIADNE_TESTING_EXPECT_EQUAL(result.certified_cuts, 24u);
        const auto& polygon = result.polygon.vertices();
        ARIADNE_TESTING_PROVE(polygon.number_of_vertices() >= 3u);

        bool all_contained = true;
        for (int ix = -8; ix <= -7; ++ix) {
            for (int iy = -8; iy <= 8; ++iy) {
                // Deliberately check extremal output errors (not just model centres).
                for (int ex : {-1, 1}) {
                    for (int ey : {-1, 1}) {
                        const Rational px =
                            Rational(ix) / Rational(8) + Rational(ex) / Rational(8);
                        const Rational py =
                            Rational(iy) / Rational(8) + Rational(ey) / Rational(4);
                        for (std::size_t i = 0u; i != polygon.size(); ++i) {
                            const auto& a = polygon[i];
                            const auto& b = polygon[(i + 1u) % polygon.size()];
                            const Rational cross =
                                (b[0] - a[0]) * (py - a[1]) -
                                (b[1] - a[1]) * (px - a[0]);
                            if (cross < Rational(0)) all_contained = false;
                        }
                    }
                }
            }
        }
        ARIADNE_TESTING_PROVE(all_contained);
    }

    Void test_empty_constrained_boundary() {
        // The parameter domain is [-1,1]^2, but -x[1]<=-2 requires
        // x[1]>=2. The affine image is therefore empty and has no boundary.
        D = ExactBoxType::unit_box(2);
        x = Affine<FloatDPBounds>::variables(2, dp);
        set = ValidatedAffineConstrainedImageSet(D, {x[0], x[1]});
        set.new_parameter_constraint(-1.0_x * x[1] <= -2.0_x);

        ARIADNE_TESTING_PROVE(definitely(set.is_empty()));
        List<Pair<double, double>> vertices = set.boundary(0, 1);
        ARIADNE_TESTING_PROVE(vertices.empty());
    }

    Void test_large_zonotope_boundary() {
        // Keep n=501: the 501 distinct planar generators yield 2*n=1002
        // boundary vertices, deliberately exceeding the historical 1000-step
        // traversal cap.  Reducing n below 501 would weaken this regression.
        const SizeType n = 501u;
        D = ExactBoxType::unit_box(n);

        Affine<FloatDPBounds> fx(n, dp);
        Affine<FloatDPBounds> fy(n, dp);
        for (SizeType i = 0; i != n; ++i) {
            Nat j = static_cast<Nat>(i);
            fx[j] = FloatDPBounds(1, dp);
            fy[j] = FloatDPBounds(static_cast<Int>(i + 1u), dp);
        }

        set = ValidatedAffineConstrainedImageSet(D, {fx, fy});
        List<Pair<double, double>> vertices = set.boundary(0, 1);
        ARIADNE_TESTING_PROVE(vertices.size() > 1000u);
    }

    Void test_collinear_negative_zonotope_boundary() {
        // Two identical, negatively oriented generators must be traversed
        // from their upper parameter bounds.  Their equal angles also require
        // the deterministic generator-index tie-break in boundary sorting.
        D = ExactBoxType::unit_box(3);
        Affine<FloatDPBounds> fx(3, dp);
        Affine<FloatDPBounds> fy(3, dp);
        fx[0] = FloatDPBounds(-65536, dp);
        fx[1] = FloatDPBounds(-65536, dp);
        fx[2] = FloatDPBounds(+65536, dp);
        fy[0] = FloatDPBounds(+65536, dp);
        fy[1] = FloatDPBounds(+65536, dp);
        fy[2] = FloatDPBounds(+65536, dp);

        set = ValidatedAffineConstrainedImageSet(D, {fx, fy});
        List<Pair<double, double>> vertices = set.boundary(0, 1);
        ARIADNE_TESTING_EXPECT_EQUAL(vertices.size(), 10u);

        bool has_left = false, has_right = false, has_bottom = false, has_top = false;
        for (const auto& vertex: vertices) {
            has_left = has_left || vertex.first < -196607.0;
            has_right = has_right || vertex.first > 196607.0;
            has_bottom = has_bottom || vertex.second < -196607.0;
            has_top = has_top || vertex.second > 196607.0;
        }
        ARIADNE_TESTING_PROVE(has_left && has_right && has_bottom && has_top);
    }

    Void test_near_vertical_zonotope_boundary() {
        // The simplex minimiser accepts residuals below its progress threshold.
        // A near-vertical edge can therefore have a slightly negative outgoing
        // horizontal component at the reported minimiser.  Its atan2 angle must
        // be wrapped to [0,2*pi) rather than ordered before the other edges.
        D = ExactBoxType::unit_box(1);
        Affine<FloatDPBounds> fx(1, dp);
        Affine<FloatDPBounds> fy(1, dp);
        fx[0] = FloatDPBounds(-4e-13_pr, dp);
        fy[0] = FloatDPBounds(1, dp);

        set = ValidatedAffineConstrainedImageSet(D, {fx, fy});
        List<Pair<double, double>> vertices = set.boundary(0, 1);
        ARIADNE_TESTING_EXPECT_EQUAL(vertices.size(), 6u);

        bool narrow_x = true;
        bool has_low_y = false;
        bool has_high_y = false;
        for (const auto& vertex: vertices) {
            narrow_x = narrow_x && vertex.first > -1e-8 && vertex.first < +1e-8;
            has_low_y = has_low_y || vertex.second < -0.9;
            has_high_y = has_high_y || vertex.second > +0.9;
        }
        ARIADNE_TESTING_PROVE(narrow_x);
        ARIADNE_TESTING_PROVE(has_low_y);
        ARIADNE_TESTING_PROVE(has_high_y);
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_affine_sets");
    ARIADNE_TESTING_ADD_TEST(TestAffineSet, test_basic_coverage);
    ARIADNE_TESTING_ADD_TEST(TestAffineSet, test_empty);
    ARIADNE_TESTING_ADD_TEST(TestAffineSet, test_invalid);
    ARIADNE_TESTING_ADD_TEST(TestAffineSet, test_pure_constraint);
    ARIADNE_TESTING_ADD_TEST(TestAffineSet, test_constrained_image);
    ARIADNE_TESTING_ADD_TEST(TestAffineSet, test_nonunit_domain);
    ARIADNE_TESTING_ADD_TEST(TestAffineSet, test_separated);
    ARIADNE_TESTING_ADD_TEST(TestAffineSet, test_constrained_boundary);
    ARIADNE_TESTING_ADD_TEST(TestAffineSet, test_outer_projection_enclosure);
    ARIADNE_TESTING_ADD_TEST(TestAffineSet, test_outer_projection_multiple_constraints);
    ARIADNE_TESTING_ADD_TEST(TestAffineSet, test_outer_projection_fallback_report);
    ARIADNE_TESTING_ADD_TEST(TestAffineSet, test_outer_projection_model_errors_and_guard_bounds);
    ARIADNE_TESTING_ADD_TEST(TestAffineSet, test_empty_constrained_boundary);
    ARIADNE_TESTING_ADD_TEST(TestAffineSet, test_large_zonotope_boundary);
    ARIADNE_TESTING_ADD_TEST(TestAffineSet, test_collinear_negative_zonotope_boundary);
    ARIADNE_TESTING_ADD_TEST(TestAffineSet, test_near_vertical_zonotope_boundary);
    ARIADNE_TESTING_ADD_TEST(TestAffineSet, test_logging_branches);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
