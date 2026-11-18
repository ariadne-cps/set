/***************************************************************************
 *            set/exact_support_polygon.hpp
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

/** @file set/exact_support_polygon.hpp
 *  @brief Exact intersection of certified planar support halfplanes.
 */

#ifndef ARIADNE_EXACT_SUPPORT_POLYGON_HPP
#define ARIADNE_EXACT_SUPPORT_POLYGON_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "numeric/builtin.hpp"
#include "numeric/rational.hpp"
#include "set/polytope.hpp"

namespace Ariadne {

/** @brief A halfplane nx*x+ny*y<=upper, interpreted as exact binary64 dyadics.
 *
 *  The calling algorithm must establish that upper is a conservative
 *  support bound. Exact clipping cannot certify an incorrectly supplied cut.
 */
struct CertifiedSupportHalfplane2d {
    double nx;
    double ny;
    double upper;
};

/** @brief Intersect a certified enclosing box and certified halfplanes exactly.
 *
 *  All binary64 coefficients are converted to exact rationals before any
 *  predicate or intersection. The result encloses the original projected
 *  set if and only if the caller supplied an enclosing box and sound support
 *  upper bounds. No heuristic clipping epsilon or display conversion occurs.
 *  Invalid or nonfinite inputs throw. Empty and degenerate intersections
 *  are returned explicitly.
 */
inline Polytope certified_support_polygon_2d(
    const std::array<double, 4>& bounds,
    const std::vector<CertifiedSupportHalfplane2d>& halfplanes)
{
    using Vertex = std::pair<Rational, Rational>;
    const auto exact = [](double value) -> Rational {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Support polygon requires finite binary64 coefficients");
        }
        return Rational(ExactDouble(value));
    };
    if (bounds[0] > bounds[1] || bounds[2] > bounds[3]) {
        throw std::invalid_argument("Support polygon initial box is reversed");
    }
    std::vector<Vertex> clipped = {
        {exact(bounds[0]), exact(bounds[2])},
        {exact(bounds[1]), exact(bounds[2])},
        {exact(bounds[1]), exact(bounds[3])},
        {exact(bounds[0]), exact(bounds[3])}
    };
    const Rational zero(0);
    for (const auto& plane : halfplanes) {
        const Rational a = exact(plane.nx);
        const Rational b = exact(plane.ny);
        const Rational c = exact(plane.upper);
        const auto distance = [&](const Vertex& point) {
            return a * point.first + b * point.second - c;
        };

        std::vector<Vertex> next;
        if (clipped.empty()) break;
        next.reserve(clipped.size() + 2u);
        Vertex previous = clipped.back();
        Rational previous_distance = distance(previous);
        for (const Vertex& current : clipped) {
            const Rational current_distance = distance(current);
            const bool was_inside = previous_distance <= zero;
            const bool is_inside = current_distance <= zero;
            if (was_inside != is_inside) {
                const Rational t = previous_distance / (previous_distance - current_distance);
                next.emplace_back(previous.first + t * (current.first - previous.first),
                                  previous.second + t * (current.second - previous.second));
            }
            if (is_inside) next.push_back(current);
            previous = current;
            previous_distance = current_distance;
        }
        clipped.swap(next);
    }

    // Use the generic exact V-polytope representation. Its reduction
    // removes duplicate/collinear generators and preserves degenerate sets.
    Polytope result(2u);
    for (const auto& vertex : clipped) {
        Vector<Rational> coordinates(2u, Rational(0));
        coordinates[0] = vertex.first;
        coordinates[1] = vertex.second;
        result.new_vertex(RationalPoint(coordinates));
    }
    result.reduce2d();
    return result;
}

} // namespace Ariadne

#endif // ARIADNE_EXACT_SUPPORT_POLYGON_HPP
