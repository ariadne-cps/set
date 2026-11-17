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
#include <stdexcept>
#include <vector>

#include "algebra/vector.hpp"

namespace Ariadne {
namespace {

bool lexicographic_less_2d(const RationalPoint& lhs, const RationalPoint& rhs) {
    if (lhs[0] < rhs[0]) return true;
    if (rhs[0] < lhs[0]) return false;
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

    // Sort references, not Point objects. Moving Point<Rational> instances
    // during sorting can leave a generator with an empty coordinate array.
    // The original generators remain unchanged until the new hull is ready.
    std::vector<const VertexType*> ordered;
    ordered.reserve(_vertices.size());
    for (const auto& point: _vertices) {
        if (point.dimension() != 2u) {
            throw std::logic_error("Polytope contains a generator with invalid dimension");
        }
        ordered.push_back(&point);
    }
    std::sort(ordered.begin(), ordered.end(),
              [](const VertexType* lhs, const VertexType* rhs) {
                  return lexicographic_less_2d(*lhs, *rhs);
              });
    ordered.erase(std::unique(ordered.begin(), ordered.end(),
                              [](const VertexType* lhs, const VertexType* rhs) {
                                  return same_point_2d(*lhs, *rhs);
                              }),
                  ordered.end());

    std::vector<VertexType> hull;
    hull.reserve(ordered.size());
    if (ordered.size() <= 2u) {
        for (const auto* point: ordered) hull.push_back(*point);
        _vertices.swap(hull);
        return;
    }

    const Rational zero(0);
    for (const auto* point: ordered) {
        while (hull.size() >= 2u &&
               oriented_area_twice(hull[hull.size()-2u], hull.back(), *point) <= zero) {
            hull.pop_back();
        }
        hull.push_back(*point);
    }
    const SizeType lower_size = hull.size();
    for (SizeType i = ordered.size() - 1u; i != 0u; --i) {
        const auto& point = *ordered[i-1u];
        while (hull.size() > lower_size &&
               oriented_area_twice(hull[hull.size()-2u], hull.back(), point) <= zero) {
            hull.pop_back();
        }
        hull.push_back(point);
    }
    hull.pop_back();
    _vertices.swap(hull);
}

Polytope convex_hull(const Polytope& p1, const Polytope& p2) {
    if (p1.dimension() != p2.dimension()) {
        throw std::invalid_argument("Convex hull requires matching ambient dimensions");
    }
    Polytope result(p1.dimension());
    for (const auto& point: p1.vertices()) result.new_vertex(point);
    for (const auto& point: p2.vertices()) result.new_vertex(point);
    if (result.dimension() == 2u) result.reduce2d();
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
        for (SizeType i = 0u; i != p.dimension(); ++i) sum[i] += point[i];
    }
    const Rational count(p.number_of_vertices());
    for (SizeType i = 0u; i != p.dimension(); ++i) sum[i] /= count;
    return RationalPoint(sum);
}

} // namespace Ariadne
