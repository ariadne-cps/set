/***************************************************************************
 *            set/polytope.hpp
 *
 *  Copyright  2005-20  Alberto Casagrande, Pieter Collins
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

/** @file set/polytope.hpp
 *  @brief Exact vertex-represented polytopes and two-dimensional projections.
 */

#ifndef ARIADNE_POLYTOPE_HPP
#define ARIADNE_POLYTOPE_HPP

#include <vector>

#include "geometry/point.hpp"
#include "numeric/rational.hpp"
#include "utility/typedefs.hpp"

namespace Ariadne {

/** @brief A bounded convex polytope represented by the convex hull of exact rational points.
 *
 *  The ambient dimension is retained even when no vertices are present. In
 *  dimensions greater than two, the stored points need not be extreme; no
 *  general vertex/facet conversion or arbitrary-dimensional hull reduction is
 *  implied by this representation.
 */
class Polytope {
  public:
    using VertexType = RationalPoint;
    using ConstIterator = std::vector<VertexType>::const_iterator;

    /** @brief Construct the empty zero-dimensional polytope. */
    Polytope() = default;
    /** @brief Construct an empty polytope with specified ambient dimension. */
    explicit Polytope(DimensionType dimension);
    /** @brief Construct the convex hull of the supplied points in a specified ambient dimension. */
    Polytope(DimensionType dimension, const std::vector<VertexType>& points);

    /** @brief Return the ambient dimension, including for empty polytopes. */
    DimensionType dimension() const { return _dimension; }
    /** @brief Return true if the generating set has no points. */
    bool empty() const { return _vertices.empty(); }
    /** @brief Return the number of stored generating points; some may be non-extreme. */
    SizeType number_of_vertices() const { return _vertices.size(); }
    /** @brief Return the selected generating point. */
    const VertexType& vertex(SizeType i) const { return _vertices.at(i); }
    /** @brief Return the stored generating points. */
    const std::vector<VertexType>& vertices() const { return _vertices; }
    /** @brief Return the first generating point. */
    ConstIterator vertices_begin() const { return _vertices.begin(); }
    /** @brief Return the one-past-the-last generating point. */
    ConstIterator vertices_end() const { return _vertices.end(); }

    /** @brief Adjoin a point to the generating set, preserving ambient dimension. */
    void new_vertex(const VertexType& point);
    /** @brief Remove non-extreme and duplicate generators in 2D using exact predicates.
     *
     *  Throws for dimensions other than two. Correctly handles the empty,
     *  singleton and collinear cases without constructing an artificial area.
     */
    void reduce2d();

  private:
    DimensionType _dimension = 0u;
    std::vector<VertexType> _vertices;
};

/** @brief Convex hull of two vertex-represented polytopes of the same ambient dimension. */
Polytope convex_hull(const Polytope& p1, const Polytope& p2);

/** @brief Exact convex hull of a coordinate projection, without an N-dimensional hull computation.
 *
 *  The two coordinates must be distinct and valid even if the polytope is
 *  empty. This operation is geometric, not a drawing operation.
 */
Polytope project2d(const Polytope& p, SizeType xind, SizeType yind);

/** @brief Arithmetic mean of the stored generators, not the polytope's volume centroid.
 *
 *  Throws if the polytope has no generators.
 */
RationalPoint mean_of_generators(const Polytope& p);

} // namespace Ariadne

#endif /* ARIADNE_POLYTOPE_HPP */
