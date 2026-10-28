# test_set.py
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

from pyariadne import *

def test_validated_constrained_image_set_dimension():
    """The image dimension is independent of its parameter-domain dimension."""
    u0 = EffectiveScalarMultivariateFunction.coordinate(2, 0)
    u1 = EffectiveScalarMultivariateFunction.coordinate(2, 1)
    domain = FloatDPExactBox([
        {exact(-1): exact(1)},
        {exact(-1): exact(1)},
    ])
    mapping = EffectiveVectorMultivariateFunction([u0 + u1])
    image = ValidatedConstrainedImageSet(domain, mapping)
    assert image.number_of_parameters() == 2
    assert image.dimension() == 1



def test_new_set_tutorial_bindings():
    box = FloatDPExactBox([{exact(-1): exact(1)}])
    boxes = ListSetFloatDPExactBox()
    boxes.adjoin(box)
    assert boxes.size() == 1
    assert boxes.bounding_box().dimension() == 1
    t = EffectiveScalarUnivariateFunction.coordinate()
    smooth = Curve(EffectiveVectorUnivariateFunction([t, t * t]))
    assert smooth.dimension() == 2
    dp = DoublePrecision()
    assert len(str(smooth.value(FloatDP(exact(0.5), dp)))) > 0


def test_effective_constrained_image_set_real_box_domain():
    """A RealBox is constructed from intervals, not dicts keyed by Real."""
    real_domain = RealBox([
        RealInterval(Real(-1), Real(1)),
        RealInterval(Real(-1), Real(1)),
    ])
    u0 = EffectiveScalarMultivariateFunction.coordinate(2, 0)
    u1 = EffectiveScalarMultivariateFunction.coordinate(2, 1)
    mapping = EffectiveVectorMultivariateFunction([u0 + u1, u0 * u1])
    image = ConstrainedImageSet(real_domain, mapping)
    assert image.number_of_parameters() == 2
    assert image.dimension() == 2


def test_interpolated_curve_python_constructor_and_insert():
    """Python exposes a seeded constructor, not the C++ default constructor."""
    dp = DoublePrecision()
    p0 = FloatDPPoint([FloatDP(exact(0), dp), FloatDP(exact(0), dp)])
    p1 = FloatDPApproximationPoint([
        FloatDPApproximation(1, dp), FloatDPApproximation(1, dp)
    ])
    trajectory = InterpolatedCurve(FloatDP(exact(0), dp), p0)
    trajectory.insert(FloatDP(exact(1), dp), p1)
    assert trajectory.size() == 2
    assert trajectory.dimension() == 2
    assert trajectory.bounding_box().dimension() == 2
