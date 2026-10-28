/***************************************************************************
 *            test_function_set_sensitivity.cpp
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

#include "set/function_set.hpp"
#include "function/taylor_function.hpp"
#include "algebra/sweeper.hpp"
#include "algebra/matrix.hpp"
#include "macro/testing.hpp"

using namespace Ariadne;

// Analysis helpers consumed by downstream pavers; currently defined in set.
namespace Ariadne {
Matrix<FloatDPError> nonlinearities_zeroth_order(
    const ValidatedVectorMultivariateTaylorFunctionModelDP&, const ExactBoxType&);
Pair<SizeType,FloatDPError> nonlinearity_index_and_error(
    const ValidatedVectorMultivariateTaylorFunctionModelDP&, const ExactBoxType);
Pair<SizeType,FloatDPError> nonlinearity_index_and_error(
    const ValidatedVectorMultivariateFunction&, const ExactBoxType&);
Pair<SizeType,FloatDPError> lipschitz_index_and_error(
    const ValidatedVectorMultivariateFunction&, const ExactBoxType&);
}

Void test_nonlinearity_dominant_coordinate() {
    ExactBoxType domain{{-1.0_x,1.0_x},{-2.0_x,2.0_x}};
    auto x=ValidatedScalarMultivariateFunction::coordinate(2u,0u);
    auto y=ValidatedScalarMultivariateFunction::coordinate(2u,1u);
    ValidatedVectorMultivariateFunction polynomial({x*x, x+3*y*y});

    ValidatedVectorMultivariateTaylorFunctionModelDP model(
        domain,polynomial,TrivialSweeper<FloatDP>(dp));
    Matrix<FloatDPError> nonlinearities=nonlinearities_zeroth_order(model,domain);
    ARIADNE_TEST_EQUALS(nonlinearities.row_size(),2u);
    ARIADNE_TEST_EQUALS(nonlinearities.column_size(),2u);
    ARIADNE_TEST_ASSERT(nonlinearities[0][0].raw().get_d()>0.0);
    ARIADNE_TEST_ASSERT(nonlinearities[1][1].raw().get_d()>
                        nonlinearities[0][0].raw().get_d());

    auto direct=nonlinearity_index_and_error(model,domain);
    ARIADNE_TEST_EQUALS(direct.first,1u);
    ARIADNE_TEST_ASSERT(direct.second.raw().get_d()>0.0);

    // The non-model overload must construct a Taylor model when needed.
    auto converted=nonlinearity_index_and_error(polynomial,domain);
    ARIADNE_TEST_EQUALS(converted.first,1u);
    ARIADNE_TEST_ASSERT(converted.second.raw().get_d()>0.0);

    // A function wrapping an existing model must take the fast dispatch path.
    ValidatedConstrainedImageSet modeled_set(
        domain,ValidatedVectorMultivariateFunctionModelDP(model));
    auto modeled=nonlinearity_index_and_error(modeled_set.function(),domain);
    ARIADNE_TEST_EQUALS(modeled.first,1u);
}

Void test_nonlinearity_ignores_smaller_rows() {
    // The second output has no higher-order terms: the earlier nonzero
    // row remains the maximal one, even though another output is visited.
    ExactBoxType domain{{-1.0_x,1.0_x},{-1.0_x,1.0_x}};
    auto x=ValidatedScalarMultivariateFunction::coordinate(2u,0u);
    auto y=ValidatedScalarMultivariateFunction::coordinate(2u,1u);
    ValidatedVectorMultivariateFunction f({x*x,y});
    ValidatedVectorMultivariateTaylorFunctionModelDP model(
        domain,f,TrivialSweeper<FloatDP>(dp));
    Matrix<FloatDPError> contributions=nonlinearities_zeroth_order(model,domain);
    ARIADNE_TEST_ASSERT(contributions[0][0].raw().get_d()>0.0);
    ARIADNE_TEST_EQUALS(contributions[1][0].raw().get_d(),0.0);
    ARIADNE_TEST_EQUALS(contributions[1][1].raw().get_d(),0.0);
    auto nonlinearity=nonlinearity_index_and_error(model,domain);
    ARIADNE_TEST_EQUALS(nonlinearity.first,0u);
    ARIADNE_TEST_ASSERT(nonlinearity.second.raw().get_d()>0.0);
}

Void test_lipschitz_dominant_coordinate() {
    ExactBoxType domain{{-1.0_x,1.0_x},{-2.0_x,2.0_x}};
    auto x=ValidatedScalarMultivariateFunction::coordinate(2u,0u);
    auto y=ValidatedScalarMultivariateFunction::coordinate(2u,1u);
    ValidatedVectorMultivariateFunction f({x+2*y,3*y});

    auto sensitivity=lipschitz_index_and_error(f,domain);
    ARIADNE_TEST_EQUALS(sensitivity.first,1u);
    ARIADNE_TEST_ASSERT(sensitivity.second.raw().get_d()>0.0);

    // No variable has positive sensitivity for a constant mapping.
    auto constant=ValidatedVectorMultivariateFunction::zeros(2u,2u);
    auto zero=lipschitz_index_and_error(constant,domain);
    ARIADNE_TEST_EQUALS(zero.first,2u);
}

Int main() {
    ARIADNE_TEST_CALL(test_nonlinearity_dominant_coordinate());
    ARIADNE_TEST_CALL(test_nonlinearity_ignores_smaller_rows());
    ARIADNE_TEST_CALL(test_lipschitz_dominant_coordinate());
    return ARIADNE_TEST_FAILURES;
}
