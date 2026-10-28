// Coverage moved from function/test_multifunction.cpp without changing its
// construction or evaluation assertions. Concrete image types belong to set.
#include "algebra/sweeper.hpp"
#include "set/multifunction_image.hpp"
#include "function/taylor_function.hpp"
#include "macro/testing.hpp"

using namespace Ariadne;

Void test_multifunction_image()
{
    auto x0=EffectiveScalarMultivariateFunction::coordinate(3,0);
    auto x1=EffectiveScalarMultivariateFunction::coordinate(3,1);
    auto p=EffectiveScalarMultivariateFunction::coordinate(3,2);

    BoxDomainType dom({{-1,3},{-1,3}});
    BoxDomainType prms({{-1,1}});;

    EffectiveVectorMultivariateFunction f({1.5_dy-x0*x0-0.375_dy*x1+p,x0});

    Sweeper<FloatDP> swp;

    ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateMultifunctionPatch,mvfp,(f,prms));
    ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateMultifunction,mvf,(mvfp));

    ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateMultifunctionModel<DP>,mvfmdp,(dom,f,prms,swp));;
    ARIADNE_TEST_CONSTRUCT(ValidatedVectorMultivariateMultifunctionModel<DP>,mvfdp,(mvfmdp));;

    // The first two coordinates are inputs; the final coordinate is uncertainty.
    ARIADNE_TEST_EQUALS(mvfp.argument_size(),2u);
    ARIADNE_TEST_EQUALS(mvfp.error_domain().dimension(),1u);
    ARIADNE_TEST_EQUALS(mvfmdp.argument_size(),2u);
    ARIADNE_TEST_EQUALS(mvfmdp.domain().dimension(),2u);
    ARIADNE_TEST_EQUALS(mvfmdp.error_domain().dimension(),1u);

    Vector<ValidatedNumber> v(2,1);
    ARIADNE_TEST_ASSIGN_CONSTRUCT(ValidatedLocatedSet<RealVector>,res,mvfp(v));
    ARIADNE_TEST_ASSIGN(res,mvf(v));
    ARIADNE_TEST_ASSIGN(res,mvfmdp(v));
    ARIADNE_TEST_ASSIGN(res,mvfdp(v));
}

Int main() {
    ARIADNE_TEST_CALL(test_multifunction_image());
    return ARIADNE_TEST_FAILURES;
}
