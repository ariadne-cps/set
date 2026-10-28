/***************************************************************************
 *            tutorial_set.cpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <string>
#include "function/function.hpp"
#include "set/function_set.hpp"
#include "set/affine_set.hpp"
#include "set/curve.hpp"
#include "set/list_set.hpp"
#include "geometry/point.hpp"
#include "geometry/box.hpp"
using namespace Ariadne;
using std::cout;
namespace {
void section(const std::string& title) { cout << "\n== " << title << " ==\n"; }
template<class T> void show(const std::string& label,const T& value) { cout << label << ": " << value << "\n"; }
}
int main() {
    ExactBoxType workspace({{-2.0_x,2.0_x},{-1.5_x,1.5_x}});
    ExactBoxType observation({{-0.75_x,1.25_x},{-0.5_x,1.0_x}});
    section("Implicit regions: geometry defined by constraints");

    // Some sets are easier to describe by a condition than by a
    // parametrisation.  The following set is the unit disk, encoded as
    // 0 <= x^2+y^2 <= 1.  ConstraintSet represents the preimage itself.
    auto q = EffectiveScalarMultivariateFunction::coordinates(2);
    EffectiveVectorMultivariateFunction radius_squared({
        sqr(q[0]) + sqr(q[1])
    });
    ConstraintSet unit_disk(
        radius_squared,
        RealBox(1u, RealInterval(0, 1))
    );

    // Intersecting with a box turns an unbounded implicit description into
    // a bounded located set.  This makes bounding boxes and image operations
    // available while preserving the defining constraint.
    RealBox search_window({
        {-1.25_x, 1.25_x},
        {-1.25_x, 1.25_x}
    });
    BoundedConstraintSet bounded_disk = intersection(search_window, unit_disk);

    show("implicit unit disk", unit_disk);
    show("constraint count", unit_disk.number_of_constraints());
    show("bounded disk", bounded_disk);
    show("bounded disk domain", bounded_disk.domain());
    show("bounded disk bounding box", bounded_disk.bounding_box());
    show("disk overlaps observation", bounded_disk.overlaps(observation));

    // Images convert an implicitly defined feasible region into a
    // parametrised geometric set.  This is a common pattern after solving a
    // feasibility problem: preserve the original admissible region, then map
    // it into the variables that are useful downstream.
    EffectiveVectorMultivariateFunction shear({
        q[0] + q[1] / 2,
        q[1]
    });
    ConstrainedImageSet sheared_disk = image(bounded_disk, shear);
    show("sheared disk", sheared_disk);
    show("sheared disk bounding box", sheared_disk.bounding_box());

    auto sheared_parts = sheared_disk.split(0);
    show("first split image", sheared_parts.first.bounding_box());
    show("second split image", sheared_parts.second.bounding_box());

    section("Validated parametrised sets: dependence plus constraints");

    // Here the parameters u=(u0,u1) range over a box and are mapped through a
    // nonlinear function.  Unlike a plain bounding box, the representation
    // keeps the parameter dependence, so later algorithms can subdivide,
    // propagate constraints and build validated approximations.
    auto u = ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction parametrisation({
        u[0] + u[1],
        u[0] * u[1]
    });
    ExactBoxType parameter_domain = ExactBoxType::unit_box(2);
    ValidatedConstrainedImageSet image_set(parameter_domain, parametrisation);

    show("parameter domain", image_set.domain());
    show("parametrisation", image_set.function());
    show("image dimension", image_set.dimension());
    show("number of parameters", image_set.number_of_parameters());
    show("initial bounding box", image_set.bounding_box());

    // Constraints may live either in parameter space or in the image space.
    // Space constraints are automatically composed with the parametrisation.
    image_set.new_parameter_constraint(u[0] + u[1] <= ExactNumber(1.25_x));

    auto x = ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedExactConstraint state_constraint = x[0] - x[1] <= ExactNumber(1.5_x);
    image_set.new_space_constraint(state_constraint);

    show("constrained image set", image_set);
    show("constraint count", image_set.number_of_constraints());
    show("constraint bounds", image_set.constraint_bounds());
    show("bounding box before reduction", image_set.bounding_box());

    // reduce() contracts the parameter domain by propagating the constraints.
    // It does not change the represented set; it changes only the internal
    // domain used for subsequent computations.
    image_set.reduce();
    show("original parameter domain", image_set.domain());
    show("reduced parameter domain", image_set.reduced_domain());
    show("bounding box after reduction", image_set.bounding_box());
    show("is empty", image_set.is_empty());

    // Splitting and restriction are the two basic ways to localise a
    // nonlinear set.  Splitting is algorithmic; restriction is useful when
    // the user already knows which part of parameter space is relevant.
    auto image_parts = image_set.split(0);
    show("first validated split", image_parts.first.bounding_box());
    show("second validated split", image_parts.second.bounding_box());

    ExactBoxType left_parameter_domain({
        {-1.0_x, 0.0_x},
        {-1.0_x, 1.0_x}
    });
    ValidatedConstrainedImageSet restricted =
        image_set.restriction(left_parameter_domain);
    show("restricted parameter domain", restricted.domain());
    show("restricted image bounds", restricted.bounding_box());

    // Predicates act on the represented set rather than just on its bounding
    // box.  Results are validated logical values and may therefore be
    // indeterminate if the requested relation cannot be certified.
    ExactBoxType query_box({
        {-0.5_x, 0.5_x},
        {-0.5_x, 0.5_x}
    });
    show("overlaps query box", image_set.overlaps(query_box));
    show("separated from query box", image_set.separated(query_box));
    show("inside workspace", image_set.inside(workspace));

    section("Affine over-approximation: cheaper validated geometry");

    // An affine constrained image preserves correlations that a box loses,
    // while being cheaper to manipulate than a nonlinear constrained image.
    // It is therefore a useful hand-off representation for downstream
    // geometry algorithms.
    ValidatedAffineConstrainedImageSet affine_set(
        image_set.affine_over_approximation()
    );
    show("affine over-approximation", affine_set);
    show("affine dimension", affine_set.dimension());
    show("affine is bounded", affine_set.is_bounded());
    show("affine is empty", affine_set.is_empty());
    show("affine bounding box", affine_set.bounding_box());

    section("Interpolated curves: geometry from sampled trajectories");

    // InterpolatedCurve is useful when the data already comes as samples,
    // e.g. a simulation trace or a sequence of measured states.  The curve
    // retains the parameter values and linearly interpolates between points.
    FloatDPPoint p0({ExactDouble(0.0), ExactDouble(0.0)}, dp);
    FloatDPPoint p1({ExactDouble(0.5), ExactDouble(0.25)}, dp);
    FloatDPPoint p2({ExactDouble(1.0), ExactDouble(1.0)}, dp);

    InterpolatedCurve trajectory(FloatDP(0.0_x, dp), p0);
    trajectory.insert(FloatDP(0.5_x, dp), p1);
    trajectory.insert(FloatDP(1.0_x, dp), p2);

    show("trajectory", trajectory);
    show("stored samples", trajectory.size());
    show("trajectory dimension", trajectory.dimension());
    show("trajectory bounding box", trajectory.bounding_box());


    section("ListSet: finite collections of basic sets");
    ListSet<ExactBoxType> pieces;
    pieces.adjoin(ExactBoxType({{-1.0_x,0.0_x},{-1.0_x,0.0_x}}));
    pieces.adjoin(ExactBoxType({{0.0_x,1.0_x},{0.0_x,1.0_x}}));
    show("parts",pieces);
    show("count",pieces.size());
    show("bounding box",pieces.bounding_box());

    section("Smooth curves from functions");
    auto t=EffectiveScalarUnivariateFunction::coordinate();
    Curve smooth(EffectiveVectorUnivariateFunction({t,t*t}));
    FloatDP parameter(0.5_x,dp);
    show("curve dimension",smooth.dimension());
    show("curve value",smooth.value(parameter));
    show("curve tangent",smooth.tangent(parameter));

    cout << "\nSee tutorial_geometry for fundamentals and tutorial_paving for grids.\n";
    return 0;
}
