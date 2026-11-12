/***************************************************************************
 *            tests/test_multifunction_image.cpp
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

/** @file tests/test_multifunction_image.cpp
 *  @brief Test set-valued multifunction patches and models.
 */

// Coverage moved from function/test_multifunction.cpp without changing its
// construction or evaluation assertions. Concrete image types belong to set.
#include "algebra/sweeper.hpp"
#include "function/taylor_function.hpp"
#include "set/multifunction_image.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

class TestMultifunctionImage : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    Void test_multifunction_image() {
        auto x0 = EffectiveScalarMultivariateFunction::coordinate(3, 0);
        auto x1 = EffectiveScalarMultivariateFunction::coordinate(3, 1);
        auto p = EffectiveScalarMultivariateFunction::coordinate(3, 2);

        BoxDomainType dom({
            {-1, 3},
            {-1, 3}
        });
        BoxDomainType prms({
            {-1, 1}
        });
        ;

        EffectiveVectorMultivariateFunction f({1.5_dy - x0 * x0 - 0.375_dy * x1 + p, x0});

        Sweeper<FloatDP> swp;

        ValidatedVectorMultivariateMultifunctionPatch mvfp(f, prms);
        ARIADNE_TESTING_NOTE("Constructed ValidatedVectorMultivariateMultifunctionPatch mvfp = " << mvfp);
        ValidatedVectorMultivariateMultifunction mvf(mvfp);
        ARIADNE_TESTING_NOTE("Constructed ValidatedVectorMultivariateMultifunction mvf = " << mvf);

        ValidatedVectorMultivariateMultifunctionModel<DP> mvfmdp(dom, f, prms, swp);
        ARIADNE_TESTING_NOTE("Constructed ValidatedVectorMultivariateMultifunctionModel<DP> mvfmdp = " << mvfmdp);
        ;
        ValidatedVectorMultivariateMultifunctionModel<DP> mvfdp(mvfmdp);
        ARIADNE_TESTING_NOTE("Constructed ValidatedVectorMultivariateMultifunctionModel<DP> mvfdp = " << mvfdp);
        ;

        // The first two coordinates are inputs; the final coordinate is uncertainty.
        ARIADNE_TESTING_EXPECT_EQUAL(mvfp.argument_size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(mvfp.error_domain().dimension(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(mvfmdp.argument_size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(mvfmdp.domain().dimension(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(mvfmdp.error_domain().dimension(), 1u);

        Vector<ValidatedNumber> v(2, 1);
        ValidatedLocatedSet<RealVector> res = mvfp(v);
        ARIADNE_TESTING_NOTE("Assigned ValidatedLocatedSet<RealVector> res = " << res);
        res = (mvf(v));
        ARIADNE_TESTING_NOTE("Assigned res = " << res);
        res = (mvfmdp(v));
        ARIADNE_TESTING_NOTE("Assigned res = " << res);
        res = (mvfdp(v));
        ARIADNE_TESTING_NOTE("Assigned res = " << res);
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_multifunction_image");
    ARIADNE_TESTING_ADD_TEST(TestMultifunctionImage, test_multifunction_image);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
