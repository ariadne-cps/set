#!/usr/bin/env python3
# tutorial_set.py
# Copyright  2026  Luca Geretti
#
# This file is part of Ariadne.
#
# Ariadne is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Ariadne is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.

"""Concrete sets and curves: companion to tutorial_set.cpp."""
from pyariadne import *

def section(title):
    print(f"\n== {title} ==")

def show(name, value):
    print(f"{name}: {value}")

def tutorial_set():
    section("Effective and validated constrained images")
    u0 = EffectiveScalarMultivariateFunction.coordinate(2, 0)
    u1 = EffectiveScalarMultivariateFunction.coordinate(2, 1)
    function = EffectiveVectorMultivariateFunction([u0 + u1, u0 * u1])
    domain = FloatDPExactBox([
        {exact(-1): exact(1)},
        {exact(-1): exact(1)},
    ])
    real_domain = RealBox([RealInterval(Real(-1), Real(1)),
                           RealInterval(Real(-1), Real(1))])
    effective = ConstrainedImageSet(real_domain, function)
    show("effective set", effective)
    show("effective domain", effective.domain())
    show("image dimension", effective.dimension())
    show("parameters", effective.number_of_parameters())
    show("effective bounding box", effective.bounding_box())
    effective_left, effective_right = effective.split(0)
    show("effective split 1", effective_left.bounding_box())
    show("effective split 2", effective_right.bounding_box())

    section("Validated dependence, reduction and restriction")
    image = ValidatedConstrainedImageSet(domain, function)
    show("validated set", image)
    show("validated domain", image.domain())
    show("validated bounding box", image.bounding_box())
    show("number of parameters", image.number_of_parameters())
    show("number of constraints", image.number_of_constraints())
    image.reduce()
    show("reduced domain", image.reduced_domain())
    half = FloatDPExactBox([
        {exact(-1): exact(0)}, {exact(-1): exact(1)}
    ])
    show("restricted image", image.restriction(half).bounding_box())
    first, second = image.split(0)
    show("split 1", first.bounding_box())
    show("split 2", second.bounding_box())

    section("Affine constrained approximation")
    affine = image.affine_over_approximation()
    show("affine image", affine)
    show("affine dimension", affine.dimension())
    show("affine bounds", affine.bounding_box())
    show("is bounded", affine.is_bounded())

    section("ListSet: finite unions of boxes")
    pieces = ListSetFloatDPExactBox()
    pieces.adjoin(domain)
    pieces.adjoin(half)
    show("number of box pieces", pieces.size())
    show("combined bounds", pieces.bounding_box())

    section("Smooth curves defined by functions")
    t = EffectiveScalarUnivariateFunction.coordinate()
    smooth = Curve(EffectiveVectorUnivariateFunction([t, t * t]))
    show("curve dimension", smooth.dimension())
    show("curve point", smooth.value(FloatDP(exact(0.5), DoublePrecision())))

    section("Interpolated sampled curve")
    dp = DoublePrecision()
    p0 = FloatDPPoint([
        FloatDP(exact(0), dp), FloatDP(exact(0), dp)
    ])
    p1 = FloatDPApproximationPoint([
        FloatDPApproximation(1, dp), FloatDPApproximation(1, dp)
    ])
    trajectory = InterpolatedCurve(FloatDP(exact(0), dp), p0)
    trajectory.insert(FloatDP(exact(1), dp), p1)
    show("trajectory", trajectory)
    show("samples", trajectory.size())
    show("dimension", trajectory.dimension())
    show("bounding box", trajectory.bounding_box())

    print("\nSee tutorial_geometry for primitives, tutorial_paving for grids.")

if __name__ == "__main__":
    tutorial_set()
