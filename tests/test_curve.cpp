/***************************************************************************
 *            test_curve.cpp
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

#include <cmath>
#include <memory>
#include <sstream>

#include "set/curve.hpp"

#include "geometry/point.hpp"
#include "geometry/box.hpp"
#include "ariadne-testing.hpp"

using namespace Ariadne;

Void test_function_curve() {
    auto x=EffectiveScalarUnivariateFunction::coordinate();
    EffectiveVectorUnivariateFunction f({x,x*x});
    Curve curve(f);
    ARIADNE_TESTING_EXPECT_EQUAL(curve.dimension(),2u);
    ARIADNE_TESTING_PROVE(curve.smoothness()==SMOOTH);

    Curve::ParameterType s(2.0_x,dp);
    auto pt=curve.value(s);
    ARIADNE_TESTING_EXPECT_EQUAL(pt.size(),2u);
    ARIADNE_TESTING_PROVE(std::abs(pt[0].raw().get_d()-2.0)<1e-8);
    ARIADNE_TESTING_PROVE(std::abs(pt[1].raw().get_d()-4.0)<1e-8);

    auto tangent=curve.tangent(s);
    ARIADNE_TESTING_EXPECT_EQUAL(tangent.size(),2u);
    ARIADNE_TESTING_PROVE(std::abs(tangent[0].raw().get_d()-1.0)<1e-8);
    ARIADNE_TESTING_PROVE(std::abs(tangent[1].raw().get_d()-4.0)<1e-8);

    Curve copy(curve);
    ARIADNE_TESTING_EXPECT_EQUAL(copy.dimension(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(copy.value(s).size(),2u);
    std::unique_ptr<Curve> cloned(curve.clone());
    ARIADNE_TESTING_EXPECT_EQUAL(cloned->dimension(),2u);
    ARIADNE_TESTING_PROVE(std::abs(cloned->value(s)[1].raw().get_d()-4.0)<1e-8);

    std::ostringstream repr;
    repr<<curve;
    ARIADNE_TESTING_PROVE(repr.str().find("Curve( function=")!=std::string::npos);
}

Void test_interpolated_curve() {
    using P=InterpolatedCurve::PointType;
    P a({FloatDPApproximation(1,dp),FloatDPApproximation(2,dp)});
    P b({FloatDPApproximation(3,dp),FloatDPApproximation(4,dp)});
    P c({FloatDPApproximation(-1,dp),FloatDPApproximation(5,dp)});

    InterpolatedCurve curve(a);
    ARIADNE_TESTING_EXPECT_EQUAL(curve.size(),1u);
    ARIADNE_TESTING_EXPECT_EQUAL(curve.dimension(),2u);
    curve.insert(InterpolatedCurve::ParameterType(1.0_x,dp),b);
    curve.insert(Dyadic(2),c);
    ARIADNE_TESTING_EXPECT_EQUAL(curve.size(),3u);
    // The curve's end() intentionally returns the last sample, not past-the-end.
    ARIADNE_TESTING_EXPECT_EQUAL(curve.end()->first,InterpolatedCurve::ParameterType(2.0_x,dp));
    ARIADNE_TESTING_PROVE(decide(curve.end()->second==c));

    auto box=curve.bounding_box();
    ARIADNE_TESTING_EXPECT_EQUAL(box.dimension(),2u);
    ARIADNE_TESTING_NOTE("box = " << box);
    ARIADNE_TESTING_PROVE(definitely(subset(box,ExactBoxType({{-2.0_x,4.0_x},{1.0_x,6.0_x}}))));
    // All three points must contribute to the range, including the final one.
    ARIADNE_TESTING_PROVE(std::abs(box[0].lower_bound().raw().get_d()+1.0)<1e-8);
    ARIADNE_TESTING_PROVE(std::abs(box[0].upper_bound().raw().get_d()-3.0)<1e-8);
    ARIADNE_TESTING_PROVE(std::abs(box[1].lower_bound().raw().get_d()-2.0)<1e-8);
    ARIADNE_TESTING_PROVE(std::abs(box[1].upper_bound().raw().get_d()-5.0)<1e-8);

    std::unique_ptr<InterpolatedCurve> cloned(curve.clone());
    ARIADNE_TESTING_EXPECT_EQUAL(cloned->size(),3u);
    ARIADNE_TESTING_EXPECT_EQUAL(cloned->dimension(),2u);

    std::ostringstream repr;
    repr<<curve;
    ARIADNE_TESTING_PROVE(repr.str().find("InterpolatedCurve( size=3")!=std::string::npos);
}

;

class TestCurve : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    Void test_curve_overloads_and_drawing() {
    using P=InterpolatedCurve::PointType;
    P a({FloatDPApproximation(1,dp),FloatDPApproximation(2,dp)});
    P b({FloatDPApproximation(3,dp),FloatDPApproximation(4,dp)});
    InterpolatedCurve empty;
    ARIADNE_TESTING_EXPECT_EQUAL(empty.size(),0u);
    empty.insert(Dyadic(0),a);
    ARIADNE_TESTING_EXPECT_EQUAL(empty.size(),1u);
    InterpolatedCurve parameter(FloatDP(2,dp),a);
    ARIADNE_TESTING_EXPECT_EQUAL(parameter.begin()->first,FloatDP(2,dp));
    InterpolatedCurve dyadic(Dyadic(3),a);
    ARIADNE_TESTING_EXPECT_EQUAL(dyadic.begin()->first,FloatDP(3,dp));
    Vector<FloatDP> raw({1.0_x,2.0_x},dp);
    InterpolatedCurve from_raw(FloatDP(0,dp),raw);
    from_raw.insert(FloatDP(1,dp),raw);
    ARIADNE_TESTING_EXPECT_EQUAL(from_raw.size(),2u);
    P wrong_dimension({FloatDPApproximation(1,dp)});
    ARIADNE_TESTING_EXPECT_EXCEPTION(from_raw.insert(FloatDP(2,dp),wrong_dimension));
    ARIADNE_TESTING_EXPECT_EXCEPTION(from_raw.insert(FloatDP(2,dp),Vector<FloatDP>({1.0_x},dp)));
    ARIADNE_TESTING_PROVE(decide(from_raw.begin()->second==a));
    InterpolatedCurve open(a,b);
    open.insert(Dyadic(2),a);
    std::ostringstream repr;
    open._write(repr);
    ARIADNE_TESTING_PROVE(repr.str().find("InterpolatedCurve")!=String::npos);
    auto x=EffectiveScalarUnivariateFunction::coordinate();
    Curve function_curve(EffectiveVectorUnivariateFunction({x,x*x}));
    const CurveInterface& interface=function_curve;
    repr<<interface;
    ARIADNE_TESTING_PROVE(repr.str().find("Curve( function=")!=String::npos);
    Vector<FloatMPApproximation> mp({1.5_x},MultiplePrecision(128));
    ARIADNE_TESTING_EXPECT_EQUAL(mp[0].get_d(),1.5);
}

};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_curve");
    ARIADNE_TESTING_ADD_TEST(TestCurve, test_curve_overloads_and_drawing);
    ARIADNE_TESTING_ADD_TEST(TestCurve, test_function_curve);
    ARIADNE_TESTING_ADD_TEST(TestCurve, test_interpolated_curve);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
