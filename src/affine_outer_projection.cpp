/***************************************************************************
 *            src/affine_outer_projection.cpp
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

/** @file src/affine_outer_projection.cpp
 *  @brief Direct validated outer projection of constrained affine images.
 */

#include "set/affine_set.hpp"
#include "set/exact_support_polygon.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "algebra/matrix.hpp"
#include "algebra/vector.hpp"
#include "linear-programming/linear_simplex_algorithm.hpp"

namespace Ariadne {

OuterProjection2dResult outer_project2d_with_report(
    const ValidatedAffineConstrainedImageSet& set,
    SizeType xind, SizeType yind, SizeType directions, bool attempt_simplex) {
    if (xind >= set.dimension() || yind >= set.dimension() || xind == yind) {
        throw std::invalid_argument("Invalid or repeated affine projection coordinates");
    }
    if (directions == 0u) {
        throw std::invalid_argument("Outer projection requires at least one support direction");
    }

    const SizeType nx = set.number_of_parameters();
    const SizeType nc = set.number_of_constraints();
    const SizeType np = nx + nc + 2u;
    const auto& xa = set.space_models()[xind];
    const auto& ya = set.space_models()[yind];
    const auto& constraints = set.constraint_models();
    const auto box = set.bounding_box();
    const double xl = box[xind].lower_bound().raw().get_d();
    const double xu = box[xind].upper_bound().raw().get_d();
    const double yl = box[yind].lower_bound().raw().get_d();
    const double yu = box[yind].upper_bound().raw().get_d();
    if (!std::isfinite(xl) || !std::isfinite(xu) ||
        !std::isfinite(yl) || !std::isfinite(yu) ||
        xl > xu || yl > yu) {
        throw std::invalid_argument("Outer projection requires a finite validated bounding box");
    }

    const auto bd = [](double v) { return FloatDPBounds(ExactDouble(v), dp); };
    const auto from_float = [&](const FloatDP& v) { return bd(v.get_d()); };
    const FloatDPBounds zero = bd(0.0);
    const FloatDPBounds one = bd(1.0);
    const FloatDPBounds unit_interval(ExactDouble(-1.0), ExactDouble(1.0), dp);

    // An enclosing relaxation A*z=b, lower<=z<=upper. All coefficients,
    // bounds and affine-model errors are independently outward enclosed.
    std::vector<FloatDPBounds> lower(np, zero), upper(np, zero), rhs(nc, zero);
    std::vector<std::vector<FloatDPBounds>> coefficients(
        nc, std::vector<FloatDPBounds>(np, zero));
    for (SizeType j = 0u; j != nx; ++j) {
        lower[j] = -one;
        upper[j] = one;
    }
    for (SizeType i = 0u; i != nc; ++i) {
        const auto& constraint = constraints[i];
        const auto& fm = constraint.function();
        rhs[i] = from_float(fm.value());
        FloatDPBounds natural_range = from_float(fm.value());
        for (SizeType j = 0u; j != nx; ++j) {
            coefficients[i][j] = -from_float(fm.gradient(j));
            natural_range += from_float(fm.gradient(j)) * unit_interval;
        }
        coefficients[i][nx + i] = one;
        const FloatDPBounds error = bd(fm.error().raw().get_d());
        lower[nx + i] = bd(natural_range.lower_raw().get_d());
        upper[nx + i] = bd(natural_range.upper_raw().get_d());

        // The guard endpoints are validated outward bounds, not central values.
        const double cl = constraint.lower_bound().lower_raw().get_d();
        const double cu = constraint.upper_bound().upper_raw().get_d();
        if (std::isfinite(cl)) lower[nx + i] = max(lower[nx + i], bd(cl) - error);
        if (std::isfinite(cu)) upper[nx + i] = min(upper[nx + i], bd(cu) + error);
    }
    const FloatDPBounds err_x = bd(xa.error().raw().get_d());
    const FloatDPBounds err_y = bd(ya.error().raw().get_d());
    lower[nx + nc] = -err_x;
    upper[nx + nc] = err_x;
    lower[nx + nc + 1u] = -err_y;
    upper[nx + nc + 1u] = err_y;

    // The approximate simplex may suggest *any* dual multipliers: the
    // certificate below does not trust its feasibility or optimality.
    LinearSimplexSolver<FloatDPApproximation> solver;
    Matrix<FloatDPApproximation> A(nc, np, dp);
    Vector<FloatDPApproximation> b(nc, dp);
    Vector<FloatDPApproximation> l(np, dp), u(np, dp);
    Array<Slackness> vt;
    Array<SizeType> p(nc);
    Matrix<FloatDPApproximation> B(nc, nc, dp);
    Vector<FloatDPApproximation> x(np, dp), y(nc, dp);
    bool lp_ready = false;
    SizeType simplex_failures = 0u;

    if (nc != 0u && attempt_simplex) {
        for (SizeType j = 0u; j != nx; ++j) {
            l[j] = FloatDPApproximation(-1.0_x, dp);
            u[j] = FloatDPApproximation(+1.0_x, dp);
        }
        for (SizeType i = 0u; i != nc; ++i) {
            const auto& constraint = constraints[i];
            const auto& fm = constraint.function();
            for (SizeType j = 0u; j != nx; ++j) {
                A[i][j] = neg(numeric_cast<FloatDPApproximation>(fm.gradient(j)));
            }
            for (SizeType j = nx; j != np; ++j) A[i][j] = 0.0_x;
            A[i][nx + i] = 1.0_x;
            const FloatDPApproximation fe =
                numeric_cast<FloatDPApproximation>(fm.error().raw());
            b[i] = numeric_cast<FloatDPApproximation>(fm.value());
            l[nx + i] = numeric_cast<FloatDPApproximation>(constraint.lower_bound().value()) - fe;
            u[nx + i] = numeric_cast<FloatDPApproximation>(constraint.upper_bound().value()) + fe;
        }
        l[nx + nc] = static_cast<FloatDPApproximation>(-xa.error().raw());
        u[nx + nc] = static_cast<FloatDPApproximation>(+xa.error().raw());
        l[nx + nc + 1u] = static_cast<FloatDPApproximation>(-ya.error().raw());
        u[nx + nc + 1u] = static_cast<FloatDPApproximation>(+ya.error().raw());
        try {
            lp_ready = possibly(solver.hotstarted_feasible(
                l, u, A, b, vt, p, B, x, y));
        } catch (const std::exception&) {
            lp_ready = false;
            ++simplex_failures;
        }
    }

    std::vector<CertifiedSupportHalfplane2d> cuts;
    cuts.reserve(directions);
    SizeType simplex_cuts = 0u;
    SizeType relaxed_cuts = 0u;
    for (SizeType direction = 0u; direction != directions; ++direction) {
        // Each direction is the exact pair of binary64 values chosen here.
        const double angle = 2.0 * std::numbers::pi *
            static_cast<double>(direction) / static_cast<double>(directions);
        const double dx = std::cos(angle), dy = std::sin(angle);
        const FloatDPBounds bx = bd(dx), by = bd(dy);
        std::vector<FloatDPBounds> objective(np, zero);
        Vector<FloatDPApproximation> cost(np, FloatDPApproximation(0.0_x, dp));
        for (SizeType j = 0u; j != nx; ++j) {
            objective[j] = bx * from_float(xa.gradient(j)) +
                           by * from_float(ya.gradient(j));
            cost[j] = -(FloatDPApproximation(dx, dp) * FloatDPApproximation(xa.gradient(j)) +
                        FloatDPApproximation(dy, dp) * FloatDPApproximation(ya.gradient(j)));
        }
        objective[nx + nc] = bx;
        objective[nx + nc + 1u] = by;
        cost[nx + nc] = -FloatDPApproximation(dx, dp);
        cost[nx + nc + 1u] = -FloatDPApproximation(dy, dp);

        std::vector<FloatDPBounds> multipliers(nc, zero);
        bool used_simplex = false;
        if (lp_ready) {
            try {
                solver.hotstarted_minimise(cost, l, u, A, b, vt, p, B);
                for (SizeType i = 0u; i != nc; ++i) {
                    FloatDPApproximation estimate(0.0_x, dp);
                    for (SizeType k = 0u; k != nc; ++k) {
                        if (p[k] >= np) throw std::runtime_error("Invalid dual basis index");
                        estimate -= cost[p[k]] * B[k][i];
                    }
                    multipliers[i] = bd(estimate.raw().get_d());
                }
                used_simplex = true;
            } catch (const std::exception&) {
                lp_ready = false;
                ++simplex_failures;
                std::fill(multipliers.begin(), multipliers.end(), zero);
            }
        }

        // For any multiplier y: c*z <= y*b + sum(max((c-A^T*y)*l,
        //                                             (c-A^T*y)*u)).
        // Outward interval operations prove the bound independently of LP.
        FloatDPBounds support = bx * from_float(xa.value()) +
                                by * from_float(ya.value());
        for (SizeType i = 0u; i != nc; ++i) support += multipliers[i] * rhs[i];
        bool finite = true;
        for (SizeType j = 0u; j != np; ++j) {
            FloatDPBounds reduced = objective[j];
            for (SizeType i = 0u; i != nc; ++i) {
                reduced -= multipliers[i] * coefficients[i][j];
            }
            if (!std::isfinite(lower[j].lower_raw().get_d()) ||
                !std::isfinite(upper[j].upper_raw().get_d())) {
                finite = false;
                break;
            }
            support += max(reduced * lower[j], reduced * upper[j]);
        }
        if (!finite) continue;
        const double certified_upper = support.upper_raw().get_d();
        if (!std::isfinite(certified_upper)) continue;
        cuts.push_back({dx, dy, certified_upper});
        if (nc != 0u) {
            if (used_simplex) ++simplex_cuts;
            else ++relaxed_cuts;
        }
    }
    return {certified_support_polygon_2d({xl, xu, yl, yu}, cuts),
            directions, cuts.size(), simplex_cuts, relaxed_cuts, simplex_failures};
}

OuterPolygon2d outer_project2d(const ValidatedAffineConstrainedImageSet& set,
                               SizeType xind, SizeType yind, SizeType directions) {
    return outer_project2d_with_report(set, xind, yind, directions).polygon;
}

} // namespace Ariadne
