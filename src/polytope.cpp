/***************************************************************************
 *            src/polytope.cpp
 *
 *  Copyright  2007-20  Alberto Casagrande, Pieter Collins
 *  Copyright  2026     Luca Geretti
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

/** @file src/polytope.cpp
 *  @brief Exact vertex-hull operations and two-dimensional projection.
 */

#include "set/polytope.hpp"

#include <algorithm>
#include <cmath>
#include <utility>
#include <stdexcept>
#include <vector>

#include "algebra/vector.hpp"
#include "numeric/builtin.hpp"

namespace Ariadne {
namespace {

bool lexicographic_less_2d(const RationalPoint& lhs, const RationalPoint& rhs) {
    if (lhs[0] < rhs[0])
        return true;
    if (rhs[0] < lhs[0])
        return false;
    return lhs[1] < rhs[1];
}

bool same_point_2d(const RationalPoint& lhs, const RationalPoint& rhs) {
    return lhs[0] == rhs[0] && lhs[1] == rhs[1];
}

Rational oriented_area_twice(const RationalPoint& a,
                             const RationalPoint& b,
                             const RationalPoint& c) {
    return (b[0] - a[0]) * (c[1] - a[1]) -
           (b[1] - a[1]) * (c[0] - a[0]);
}

} // namespace

Polytope::Polytope(DimensionType dimension)
    : _dimension(dimension) {}

Polytope::Polytope(DimensionType dimension, const std::vector<VertexType>& points)
    : _dimension(dimension) {
    for (const auto& point: points) {
        this->new_vertex(point);
    }
}

void Polytope::new_vertex(const VertexType& point) {
    if (point.dimension() != _dimension) {
        throw std::invalid_argument("Polytope generator dimension does not match ambient dimension");
    }
    _vertices.push_back(point);
}

void Polytope::reduce2d() {
    if (_dimension != 2u) {
        throw std::invalid_argument("Polytope::reduce2d requires dimension two");
    }

    // Keep all original points immutable throughout the hull computation.
    // Sorting and maintaining indices avoids moving or assigning RationalPoint
    // instances while they may still be used by geometric predicates.
    std::vector<SizeType> order;
    order.reserve(_vertices.size());
    for (SizeType i = 0u; i != _vertices.size(); ++i) {
        if (_vertices[i].dimension() != 2u) {
            throw std::logic_error("Polytope contains a generator with invalid dimension");
        }
        order.push_back(i);
    }
    std::sort(order.begin(), order.end(),
              [this](SizeType lhs, SizeType rhs) {
                  return lexicographic_less_2d(_vertices[lhs], _vertices[rhs]);
              });
    order.erase(std::unique(order.begin(), order.end(),
                            [this](SizeType lhs, SizeType rhs) {
                                return same_point_2d(_vertices[lhs], _vertices[rhs]);
                            }),
                order.end());

    // The hull holds indices only. In particular, no point is modified by
    // either pass of the monotone-chain convex hull algorithm.
    std::vector<SizeType> hull;
    hull.reserve(2u * order.size());
    const Rational zero(0);

    for (const auto index: order) {
        // Keep the size guard separate: Rational comparisons return Ariadne
        // Boolean, whose overloaded conjunction does not short-circuit.
        while (hull.size() >= 2u) {
            if (oriented_area_twice(_vertices[hull[hull.size() - 2u]],
                                    _vertices[hull.back()], _vertices[index]) > zero)
                break;
            hull.pop_back();
        }
        hull.push_back(index);
    }
    const SizeType lower_size = hull.size();
    if (order.size() > 1u) {
        for (SizeType i = order.size() - 1u; i != 0u; --i) {
            const SizeType index = order[i - 1u];
            while (hull.size() > lower_size) {
                if (oriented_area_twice(_vertices[hull[hull.size() - 2u]],
                                        _vertices[hull.back()], _vertices[index]) > zero)
                    break;
                hull.pop_back();
            }
            hull.push_back(index);
        }
        hull.pop_back(); // Remove the repeated first vertex.
    }

    std::vector<VertexType> reduced;
    reduced.reserve(hull.size());
    for (const SizeType index: hull) {
        reduced.push_back(_vertices[index]);
    }
    _vertices.swap(reduced);
}

Polytope convex_hull(const Polytope& p1, const Polytope& p2) {
    if (p1.dimension() != p2.dimension()) {
        throw std::invalid_argument("Convex hull requires matching ambient dimensions");
    }
    Polytope result(p1.dimension());
    for (const auto& point: p1.vertices())
        result.new_vertex(point);
    for (const auto& point: p2.vertices())
        result.new_vertex(point);
    if (result.dimension() == 2u)
        result.reduce2d();
    return result;
}

Polytope project2d(const Polytope& p, SizeType xind, SizeType yind) {
    if (xind >= p.dimension() || yind >= p.dimension() || xind == yind) {
        throw std::invalid_argument("Invalid or repeated projection coordinates for Polytope");
    }
    Polytope result(2u);
    for (const auto& point: p.vertices()) {
        Vector<Rational> xy(2u, Rational(0));
        xy[0] = point[xind];
        xy[1] = point[yind];
        result.new_vertex(RationalPoint(xy));
    }
    result.reduce2d();
    return result;
}

RationalPoint mean_of_generators(const Polytope& p) {
    if (p.empty()) {
        throw std::invalid_argument("Cannot compute mean of an empty polytope's generators");
    }
    Vector<Rational> sum(p.dimension(), Rational(0));
    for (const auto& point: p.vertices()) {
        for (SizeType i = 0u; i != p.dimension(); ++i)
            sum[i] += point[i];
    }
    const Rational count(p.number_of_vertices());
    for (SizeType i = 0u; i != p.dimension(); ++i)
        sum[i] /= count;
    return RationalPoint(sum);
}

Polytope certified_support_polygon_2d(
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
