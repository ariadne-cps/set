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
