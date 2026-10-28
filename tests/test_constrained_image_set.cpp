/***************************************************************************
 *            test_constrained_image_set.cpp
 *
 *  Copyright  2009-20  Pieter Collins
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
#include <limits>
#include "function/function.hpp"
#include "function/taylor_model.hpp"
#include "algebra/algebra.hpp"
#include "geometry/box.hpp"
#include "set/affine_set.hpp"
#include "set/function_set.hpp"
#include "function/formula.hpp"

#include "macro/testing.hpp"

using namespace Ariadne;
using namespace std;

class TestConstrainedImageSet
{
  private:
  public:
    TestConstrainedImageSet() = default;

    Void test() {
        ExactIntervalType ivl(-4.0_x,4.0_x);
        ivl=ExactIntervalType{-4.0_x,4.0_x};
        ARIADNE_TEST_CALL(test_constructor());
        ARIADNE_TEST_CALL(test_domain());
        ARIADNE_TEST_CALL(test_geometry());
        ARIADNE_TEST_CALL(test_separated());
        ARIADNE_TEST_CALL(test_split());
        ARIADNE_TEST_CALL(test_affine_approximation());
    }

    Void test_constructor() {
        List<EffectiveScalarMultivariateFunction> s=EffectiveScalarMultivariateFunction::coordinates(3);
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);

        RealBox d(3,RealInterval(-1,+2));
        EffectiveConstrainedImageSet set1(d,{s[0],s[0]*s[0]/4+s[1]+s[2]/2});
        set1.new_parameter_constraint(0<=s[0]+s[1]<=1);
        set1.new_space_constraint(x[0]+x[1]<=2);
        set1.apply({x[0]+x[1],x[0]-x[1]*x[1]});
        ARIADNE_TEST_PRINT(set1);

        EffectiveConstrainedImageSet set2(d,{s[0],s[0]*s[0]/4+s[1]+s[2]/2},{0<=s[0]+s[1]<=1,x[0]+x[1]<=2});
        ARIADNE_TEST_PRINT(set2);
    }

    Void test_domain() {
        List<EffectiveScalarMultivariateFunction> s=EffectiveScalarMultivariateFunction::coordinates(3);
        RealBox d(3,RealInterval(Decimal(-1.1),Decimal(+2.1)));
        EffectiveConstrainedImageSet set(d,{s[0],s[0]*s[0]/4+s[1]+s[2]/2});
        ARIADNE_TEST_EQUALS(set.number_of_parameters(),3);

        set.new_parameter_constraint(0<=s[0]+s[1]<=1);

        ExactBoxType bx(cast_exact_box(set.bounding_box()));
        LowerKleenean overlaps = set.overlaps(bx);
        ValidatedLowerKleenean check_overlaps=overlaps.check(Effort::get_default());
        ARIADNE_TEST_ASSERT(definitely(check_overlaps));
    }

    Void test_geometry() {
        List<EffectiveScalarMultivariateFunction> p=EffectiveScalarMultivariateFunction::coordinates(1);
        List<EffectiveScalarMultivariateFunction> s=EffectiveScalarMultivariateFunction::coordinates(3);
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        ExactBoxType box1(2);
        ExactBoxType box2(2);
        ExactBoxType box3(2);

        // Test the polytope
        EffectiveConstrainedImageSet polytope({{-2,+2},{-2,+2}},{x[0],x[1]});
        polytope.new_parameter_constraint(x[0]+Real(1.5_x)*+x[1]<=1);
        box1=ExactBoxType( {ExactIntervalType(1.0_x,2.0_x),ExactIntervalType(0.5_x,1.0_x)} );
        box2=ExactBoxType( {ExactIntervalType(0.0_x,1.0_x),ExactIntervalType(0.5_x,1.0_x)} );
        ARIADNE_TEST_ASSERT(definitely(polytope.separated(box1),Effort::get_default()));
        //ARIADNE_TEST_ASSERT(polytope.overlaps(box2));

        // Test the unit disc
        EffectiveConstrainedImageSet disc({RealInterval(-2,+2),RealInterval(-2,+2)},{x[0],x[1]});
        disc.new_parameter_constraint(x[0]*x[0]+x[1]*x[1]<=1);
        box1=ExactBoxType( {ExactIntervalType(-0.5_x,0.5_x),ExactIntervalType(0.25_x,0.5_x)} );
        box2=ExactBoxType( {ExactIntervalType(1,2),ExactIntervalType(0.5_x,1)} );
        box3=ExactBoxType( {ExactIntervalType(0.75_x,2),ExactIntervalType(-1.0_x,-0.5_x)} );
        //ARIADNE_TEST_ASSERT(disc.overlaps(box1));
        ARIADNE_TEST_ASSERT(definitely(disc.separated(box2),Effort::get_default()));
        //ARIADNE_TEST_ASSERT(disc.overlaps(box3));

        // Test a one-dimensional parabolic set
        EffectiveConstrainedImageSet parabola({RealInterval(-1,+1)},{p[0],p[0]*p[0]});
        box1=ExactBoxType( {ExactIntervalType(0,0.5_x),ExactIntervalType(0.5_x,1)} );
        box2=ExactBoxType( {ExactIntervalType(0.75_x,2),ExactIntervalType(0.5_x,1)} );
        ARIADNE_TEST_PRINT(parabola);
        ARIADNE_TEST_ASSERT(parabola.separated(box1));
        //ARIADNE_TEST_ASSERT(parabola.overlaps(box2));

        // Test whether the second iterate of the Henon map intersects a box
        Real half(0.5_x);
        RealBox d(2,RealInterval(-half,+half));
        Real a(1.5_x); Real b(0.375_x);
        EffectiveVectorMultivariateFunction h={a-x[0]*x[0]-b*x[1],x[0]};
        EffectiveVectorMultivariateFunction f=compose(h,h);
        EffectiveConstrainedImageSet set(d,f);
        set.new_parameter_constraint(0<=x[0]+x[1]<=1);

        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_ASSERT(set.separated(ExactBoxType{{1.0_x,1.25_x},{1.0_x,1.25_x}}));
        //ARIADNE_TEST_ASSERT(set.overlaps(ExactBoxType(2, -1.0_x,-0.875_x, 1.375_x,1.625_x)));

        ARIADNE_TEST_PRINT((ExactPointType({0.375_x,-0.375_x},dp)));
        ARIADNE_TEST_PRINT(f(ExactPointType({0.375_x,-0.375_x},dp)));
    }

    Void test_separated() {
        List<EffectiveScalarMultivariateFunction> s=EffectiveScalarMultivariateFunction::coordinates(3);
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);

        RealBox d(3,RealInterval(Decimal(-1.1),Decimal(+2.1)));
//        RealBox d(3,RealInterval(Decimal(-0.1015625_x),Decimal(+2.1015625_x)));
        EffectiveConstrainedImageSet set(d,{s[0],s[0]*s[0]/4+s[1]+s[2]/2});
        set.new_parameter_constraint(0<=s[0]+s[1]<=1);
        ExactBoxType b(cast_exact_box(set.bounding_box()));
        List<ExactBoxType> stack(1u,b);
        while(!stack.empty()) {
            ExactBoxType bx=stack.back();
            stack.pop_back();
            if(possibly(bx.radius()>=0.125_x)) {
                Pair<ExactBoxType,ExactBoxType> sbx=bx.split();
                stack.append(sbx.first); stack.append(sbx.second);
            } else {
                ValidatedLowerKleenean overlaps = check(set.overlaps(bx),Effort::get_default());
                LowerKleenean separated = set.separated(bx);
            if(definitely(overlaps)) {
                } else if(definitely(separated,Effort::get_default())) {
                } else {
                }
            }
        }


        // Regression test
        EffectiveConstrainedImageSet quadratic_set(RealBox({{-1.0_x,1.0_x},{-1.0_x,1.0_x}}));
        quadratic_set.apply( {2*x[0]+x[1]+x[0]*x[0]/4,x[0]+x[1]} );
        ExactBoxType box{{0.750000_x,1.00000_x},{0.00000_x,0.250000_x}};
        ARIADNE_TEST_ASSERT(not definitely(set.separated(box).check(Effort::get_default())));
    }




    Void test_split() {
        EffectiveScalarMultivariateFunction o=EffectiveScalarMultivariateFunction::constant(3,1);
        EffectiveScalarMultivariateFunction s0=EffectiveScalarMultivariateFunction::coordinate(3,0);
        EffectiveScalarMultivariateFunction s1=EffectiveScalarMultivariateFunction::coordinate(3,1);
        EffectiveScalarMultivariateFunction s2=EffectiveScalarMultivariateFunction::coordinate(3,2);
        EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction x1=EffectiveScalarMultivariateFunction::coordinate(2,1);
        EffectiveVectorMultivariateFunction translation;
        RealBox d(3,RealInterval(-1,+1));
        EffectiveConstrainedImageSet set(d,{s0,s1+s2*s2/2});
        set.new_parameter_constraint(s0+Real(0.75_x)*s1+s2<=Real(0));

        EffectiveConstrainedImageSet subset1,subset2;
        make_lpair(subset1,subset2)=set.split(0);
        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_PRINT(subset1);
        ARIADNE_TEST_PRINT(subset2);
        ARIADNE_TEST_PRINT(set.split(1));

        EffectiveConstrainedImageSet subset11,subset12,subset21,subset22;
        make_lpair(subset11,subset12)=subset1.split(0);
        make_lpair(subset21,subset22)=subset2.split(0);
        ARIADNE_TEST_PRINT(subset11);
        translation={x0-Real(2.5_x),x1};
        subset11.apply(translation);
        ARIADNE_TEST_PRINT(subset11);

        translation={x0-Real(2.5_x),x1};
        set.apply(translation);
    }

    Void test_affine_approximation() {
        // Test conversionn is exact for the affine set -2<x<1; 0<y<2 3x+y<1
        List<EffectiveScalarMultivariateFunction> s=EffectiveScalarMultivariateFunction::coordinates(2);
        RealBox d={RealInterval(-2,1),RealInterval(0,2)};
        EffectiveConstrainedImageSet set(d,{s[0],s[1]});
        set.new_parameter_constraint(3*s[0]+s[1]<=1);
        ValidatedAffineConstrainedImageSet affine_set=set.affine_approximation();
        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_PRINT(affine_set);
        //ARIADNE_TEST_PRINT(set.affine_approximation());
    }




};

class TestValidatedConstrainedImageSet {
  public:
    Void test() const {
        ARIADNE_TEST_CALL(test_construct());
        ARIADNE_TEST_CALL(test_geometry());
    }

    Void test_construct() const {
        List<ValidatedScalarMultivariateFunction> s=ValidatedScalarMultivariateFunction::coordinates(3);
        List<ValidatedScalarMultivariateFunction> x=ValidatedScalarMultivariateFunction::coordinates(2);

        ExactBoxType d({{-1.0_x,2.0_x},{1.0_x,3.0_x},{-2.0_x,-1.0_x}});

        ValidatedConstrainedImageSet set1(d);
        ARIADNE_TEST_EQUALS(set1.dimension(),3);

        ValidatedConstrainedImageSet set2(d,{s[0],s[0]*s[0]/4+s[1]+s[2]/2});
        set2.new_parameter_constraint(s[0]*s[1]>=0);
        set2.new_space_constraint(x[0]+x[1]>=0);
        ARIADNE_TEST_EQUALS(set2.dimension(),2);
        set2.apply({x[0]+x[1],x[0]-x[1]*x[1]});
    }

    Void test_geometry() const {

        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        ExactBoxType box1(2);
        ExactBoxType box2(2);
        ExactBoxType box3(2);

        ValidatedConstrainedImageSet idisc(ExactBoxType({{-2.0_x,+2.0_x},{-2.0_x,+2.0_x}}),{x[0],x[1]});
        idisc.new_parameter_constraint(x[0]*x[0]+x[1]*x[1]<=1);
        box1=ExactBoxType( {ExactIntervalType(-0.5_x,0.5_x),ExactIntervalType(0.25_x,0.75_x)} );
        box1=ExactBoxType( {ExactIntervalType(-0.5_x,0.75_x),ExactIntervalType(0.25_x,0.75_x)} );
        box2=ExactBoxType( {ExactIntervalType(1,2),ExactIntervalType(0.5_x,1.0_x)} );
        box3=ExactBoxType( {ExactIntervalType(0.75_x,2.0_x),ExactIntervalType(-1.0_x,-0.5_x)} );
        ARIADNE_TEST_ASSERT(idisc.overlaps(box1));
        ARIADNE_TEST_ASSERT(idisc.separated(box2));
        ARIADNE_TEST_ASSERT(idisc.overlaps(box3));
    }
};


Int main()
{
    TestConstrainedImageSet().test();
    TestValidatedConstrainedImageSet().test();
    return ARIADNE_TEST_FAILURES;
}
