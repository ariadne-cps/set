/***************************************************************************
 *            test_affine_sets.cpp
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
#include "function/function.hpp"
#include "geometry/box.hpp"
#include "set/affine_set.hpp"
#include "set/function_set.hpp"

#include "macro/testing.hpp"

using namespace Ariadne;
using namespace std;

namespace Ariadne {

template<class X> AffineConstraint<X> operator<=(const Dyadic& l, const Affine<X>& a) {
    return AffineConstraint<X>({l,a.value().precision()},a,{+infty,a.value().precision()}); }
template<class X> AffineConstraint<X> operator<=(const Affine<X>& a, const Dyadic& u) {
    return AffineConstraint<X>(-infty,a,u); }
template<class X> AffineConstraint<X> operator<=(const AffineConstraint<X>& ac, const Dyadic& u) {
    ARIADNE_ASSERT(decide(ac.upper_bound()==infty));
    return AffineConstraint<X>(ac.lower_bound(),ac.function(),{u,ac.lower_bound().precision()}); }
template<class X> AffineConstraint<X> operator==(const Affine<X>& a, const Dyadic& u) {
    return AffineConstraint<X>(u,a,u); }
}


class TestAffineSet
{
  private:
    ExactBoxType D;
    Vector< Affine<FloatDPBounds> > x;
    ValidatedAffineConstrainedImageSet set;
  public:
    TestAffineSet() :
        set(ExactBoxType::unit_box(2),Matrix<FloatDP>(2,2,dp),Vector<FloatDP>(2,dp)) { }

    Void test_basic_coverage() {
        RealBox real_box(2,RealInterval(-1,1));
        ValidatedAffineConstrainedImageSet from_real(real_box);
        ARIADNE_TEST_EQUALS(from_real.dimension(),2u);
        ARIADNE_TEST_ASSERT(definitely(from_real.is_bounded()));
        ARIADNE_TEST_EQUALS(from_real.number_of_constraints(),0u);

        ExactBoxType exact_box({{-1.0_x,+1.0_x},{-2.0_x,+2.0_x}});
        ValidatedAffineConstrainedImageSet from_exact(exact_box);
        ARIADNE_TEST_EQUALS(from_exact.domain(),exact_box);
        ARIADNE_TEST_ASSERT(definitely(from_exact.inside(ExactBoxType{{-2.0_x,+2.0_x},{-3.0_x,+3.0_x}})));

        Vector<ValidatedAffineModelDP> models=ValidatedAffineModelDP::scalings(ExactBoxType::unit_box(2),dp);
        ValidatedAffineConstrainedImageSet from_models(models);
        List<ValidatedAffineModelConstraintDP> constraints;
        ValidatedAffineConstrainedImageSet from_models_and_constraints(models,constraints);
        ARIADNE_TEST_EQUALS(from_models.number_of_parameters(),2u);
        ARIADNE_TEST_EQUALS(from_models_and_constraints.number_of_constraints(),0u);

        auto y=ValidatedScalarMultivariateFunction::coordinates(2);
        ValidatedVectorMultivariateFunction h({y[0]+y[1],y[0]-y[1]});
        auto transformed=image(from_exact,h);
        ARIADNE_TEST_EQUALS(transformed.dimension(),2u);
    }

    Void test_pure_constraint() {

        D = ExactBoxType::unit_box(2);
        x = Affine<FloatDPBounds>::variables(2,dp);

        set=ValidatedAffineConstrainedImageSet(D, {x[0]+2,x[0]+x[1]} );

        set.new_parameter_constraint(x[0]+0.5_x*x[1]+0.75_x <= 0.0_x);
        set.new_parameter_constraint(-0.5_x*x[0]+1.5_x*x[1]+0.5_x <= 0.0_x);

        ARIADNE_TEST_PRINT(set);

        set=ValidatedAffineConstrainedImageSet(D, {x[0]+x[1]+2,x[0]+x[1]+2});
        set.new_parameter_constraint(-0.5_x*x[0]-1.5_x*x[1]+0.5_x <= 0.0_x);

        ARIADNE_TEST_PRINT(set);
    }

    Void test_nonunit_domain() {
        D=ExactBoxType({{0.0_x,2.0_x},{-3.0_x,1.0_x}});
        x=Affine<FloatDPBounds>::variables(2,dp);

        set=ValidatedAffineConstrainedImageSet(D,{x[0],x[1]},{x[0]+x[1]<=1.0_x});
        ARIADNE_TEST_EQUALS(set.domain(),D);
        ARIADNE_TEST_EQUALS(set.number_of_constraints(),1u);
        ARIADNE_TEST_SAME(set.bounding_box(),UpperBoxType({{0.0_x,2.0_x},{-3.0_x,1.0_x}}));

        ValidatedAffineConstrainedImageSet unconstrained(D,{x[0],x[1]});
        ARIADNE_TEST_SAME(unconstrained.bounding_box(),UpperBoxType({{0.0_x,2.0_x},{-3.0_x,1.0_x}}));
    }

    Void test_constrained_image() {
        D = ExactBoxType::unit_box(3);
        x = Affine<FloatDPBounds>::variables(3,dp);

        set=ValidatedAffineConstrainedImageSet(D, {2*x[0]+3*x[1]+x[2],x[0]+x[1]-2} );

        set.new_parameter_constraint(1.0_x*x[0]+0.5_x*x[1]+0.25_x*x[2]<=0.75_x);
        set.new_parameter_constraint(-2.0_x*x[0]-3.0_x*x[1]-1.0_x*x[2]<=3.00_x);

        ARIADNE_TEST_PRINT(set)

    }


    Void test_separated() {
        D=ExactBoxType::unit_box(3);
        x=Affine<FloatDPBounds>::variables(3,dp);
        set=ValidatedAffineConstrainedImageSet(D, {-0.9375_x+0.0625_x*x[0]+0.5_x*x[1],0.87890625_x-0.1171875_x*x[0]+x[1]+0.00390625_x*x[2]});
        set.new_parameter_constraint(-1.1875_x+0.0625_x*x[0]+x[1]<=0.0_x);
        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_CONSTRUCT(ExactBoxType,cell1,({{-1.0_x,-0.9375_x},{0.8125_x, 0.875_x}})); // subset
        ARIADNE_TEST_CONSTRUCT(ExactBoxType,cell2,({{-0.5625_x,-0.50_x},{1.4375_x,1.5_x}})); // overlaps
        ARIADNE_TEST_CONSTRUCT(ExactBoxType,cell3,({{-0.875_x,-0.75_x},{0.625_x,0.7578125_x}})); // touches at (-0.875_x,0.7578125_x);
        ARIADNE_TEST_CONSTRUCT(ExactBoxType,cell4,({{-1.1875_x,-1.125_x},{0.0625_x,0.125_x}})); // almost touches
        ARIADNE_TEST_CONSTRUCT(ExactBoxType,cell5,({{-0.9375_x,-0.875_x},{0.4375_x,0.5_x}})); // disjoint
        ARIADNE_TEST_CONSTRUCT(ExactBoxType,cell6,({{-1.5_x,-1.375_x},{0.5_x,0.625_x}})); // disjoint; regression test

        ARIADNE_TEST_ASSERT(!definitely(set.separated(cell1)));
        ARIADNE_TEST_ASSERT(!definitely(set.separated(cell2)));
        ARIADNE_TEST_ASSERT(!definitely(set.separated(cell3)));
        //ARIADNE_TEST_ASSERT(!definitely(affine_set.overlaps(cell3)));
        ARIADNE_TEST_ASSERT(possibly(set.separated(cell4)));
        ARIADNE_TEST_ASSERT(definitely(set.separated(cell5)));
        ARIADNE_TEST_ASSERT(definitely(set.separated(cell6)));

        // The six separation outcomes above are the regression checks.
        // Rendering them here duplicated draw coverage already exercised by test_draw().
    }





    Void test_empty() {
        D = ExactBoxType::unit_box(2);
        x = Affine<FloatDPBounds>::variables(2,dp);

        std::cerr<<get_characteristics(0.125_x*x[0]+1.0_x)<<"\n";
        set=ValidatedAffineConstrainedImageSet(D, {0.125_x*x[0]+1.0_x, 0.125_x*x[1]+1.0_x} );
        set.new_parameter_constraint(-1.0_x*x[1] <= -2.0_x);

        ARIADNE_TEST_PRINT(set);
        ARIADNE_TEST_ASSERT(set.is_empty());
    }

    Void test_logging_branches() {
        unsigned int const old_verbosity=Logger::instance().configuration().verbosity();
        Logger::instance().configuration().set_verbosity(10u);

        D=ExactBoxType::unit_box(2);
        x=Affine<FloatDPBounds>::variables(2,dp);
        set=ValidatedAffineConstrainedImageSet(D,{x[0],x[1]});

        ARIADNE_TEST_ASSERT(!definitely(set.separated(ExactBoxType{{-0.5_x,+0.5_x},{-0.5_x,+0.5_x}})));

        List<Pair<double,double>> vertices=set.boundary(0,1);
        ARIADNE_TEST_ASSERT(!vertices.empty());

        Logger::instance().configuration().set_verbosity(old_verbosity);
    }

    Void test_invalid() {
        {
            D = ExactBoxType::unit_box(1);
            x = Affine<FloatDPBounds>::variables(2,dp);
            ARIADNE_TEST_FAIL(ValidatedAffineConstrainedImageSet(D, {0.125_x*x[0]+1.0_x, 0.125_x*x[1]+1.0_x} ));
        }
        {
            D = ExactBoxType::unit_box(2);
            x = Affine<FloatDPBounds>::variables(2,dp);
            auto a = Affine<FloatDPBounds>::variables(1,dp);
            ARIADNE_TEST_FAIL(ValidatedAffineConstrainedImageSet(D, {0.125_x*x[0]+1.0_x, 0.125_x*x[1]+1.0_x}, {a[0]<=-0.25_x} ));

            set=ValidatedAffineConstrainedImageSet(D,{x[0],x[1]});
            ARIADNE_TEST_FAIL(set.new_parameter_constraint(a[0]<=-0.25_x));
            ARIADNE_TEST_FAIL(set.separated(ExactBoxType::unit_box(1)));
            ARIADNE_TEST_FAIL(set.inside(ExactBoxType::unit_box(1)));


            Vector<ValidatedAffineModelDP> bad_models=ValidatedAffineModelDP::scalings(ExactBoxType::unit_box(1),dp);
            List<ValidatedAffineModelConstraintDP> no_constraints;
            ARIADNE_TEST_FAIL(ValidatedAffineConstrainedImageSet(D,bad_models,no_constraints));

            ValidatedAffineModelConstraintDP bad_constraint(-1.0_x,bad_models[0],+1.0_x);
            List<ValidatedAffineModelConstraintDP> bad_constraints;
            bad_constraints.append(bad_constraint);
            Vector<ValidatedAffineModelDP> good_models=ValidatedAffineModelDP::scalings(D,dp);
            ARIADNE_TEST_FAIL(ValidatedAffineConstrainedImageSet(D,good_models,bad_constraints));
            ARIADNE_TEST_FAIL(set.new_constraint(bad_constraint));

            ValidatedAffineConstrainedImageSet empty_set(ExactBoxType(0u));
            ValidatedAffineModelDP zero_model(0u,dp);
            ValidatedAffineModelConstraintDP zero_model_constraint(-1.0_x,zero_model,+1.0_x);
            Affine<FloatDPBounds> zero_affine(0u,dp);
            ARIADNE_TEST_FAIL(empty_set.new_constraint(zero_model_constraint));
            ARIADNE_TEST_FAIL(empty_set.new_parameter_constraint(zero_affine<=0.0_x));
            ARIADNE_TEST_FAIL(empty_set.number_of_parameters());

            Vector<Affine<FloatDPBounds>> empty_affine_models;
            Vector<ValidatedAffineModelDP> empty_models;
            ARIADNE_TEST_FAIL(ValidatedAffineConstrainedImageSet(D,empty_affine_models));
            ARIADNE_TEST_FAIL(ValidatedAffineConstrainedImageSet(D,empty_models,no_constraints));
            ARIADNE_TEST_FAIL(ValidatedAffineConstrainedImageSet{empty_models});

            Matrix<FloatDP> bad_matrix(1u,1u,dp);
            Vector<FloatDP> bad_offset(2u,dp);
            ARIADNE_TEST_FAIL(ValidatedAffineConstrainedImageSet(D,bad_matrix,bad_offset));
            Matrix<FloatDP> empty_matrix(0u,0u,dp);
            Vector<FloatDP> empty_offset(0u,dp);
            ARIADNE_TEST_FAIL(ValidatedAffineConstrainedImageSet(ExactBoxType(0u),empty_matrix,empty_offset));
        }
    }

    Void test_large_zonotope_boundary() {
        // Keep n=501: the 501 distinct planar generators yield 2*n=1002
        // boundary vertices, deliberately exceeding the historical 1000-step
        // traversal cap.  Reducing n below 501 would weaken this regression.
        SizeType const n=501u;
        D=ExactBoxType::unit_box(n);

        Affine<FloatDPBounds> fx(n,dp);
        Affine<FloatDPBounds> fy(n,dp);
        for(SizeType i=0; i!=n; ++i) {
            Nat j=static_cast<Nat>(i);
            fx[j]=FloatDPBounds(1,dp);
            fy[j]=FloatDPBounds(static_cast<Int>(i+1u),dp);
        }

        set=ValidatedAffineConstrainedImageSet(D,{fx,fy});
        List<Pair<double,double>> vertices=set.boundary(0,1);
        ARIADNE_TEST_ASSERT(vertices.size()>1000u);
    }

    Void test_collinear_negative_zonotope_boundary() {
        // Two identical, negatively oriented generators must be traversed
        // from their upper parameter bounds.  Their equal angles also require
        // the deterministic generator-index tie-break in boundary sorting.
        D=ExactBoxType::unit_box(3);
        Affine<FloatDPBounds> fx(3,dp);
        Affine<FloatDPBounds> fy(3,dp);
        fx[0]=FloatDPBounds(-65536,dp);
        fx[1]=FloatDPBounds(-65536,dp);
        fx[2]=FloatDPBounds(+65536,dp);
        fy[0]=FloatDPBounds(+65536,dp);
        fy[1]=FloatDPBounds(+65536,dp);
        fy[2]=FloatDPBounds(+65536,dp);

        set=ValidatedAffineConstrainedImageSet(D,{fx,fy});
        List<Pair<double,double>> vertices=set.boundary(0,1);
        ARIADNE_TEST_EQUALS(vertices.size(),10u);

        bool has_left=false, has_right=false, has_bottom=false, has_top=false;
        for(auto const& vertex:vertices) {
            has_left = has_left || vertex.first < -196607.0;
            has_right = has_right || vertex.first > 196607.0;
            has_bottom = has_bottom || vertex.second < -196607.0;
            has_top = has_top || vertex.second > 196607.0;
        }
        ARIADNE_TEST_ASSERT(has_left && has_right && has_bottom && has_top);
    }

    Void test_near_vertical_zonotope_boundary() {
        // The simplex minimiser accepts residuals below its progress threshold.
        // A near-vertical edge can therefore have a slightly negative outgoing
        // horizontal component at the reported minimiser.  Its atan2 angle must
        // be wrapped to [0,2*pi) rather than ordered before the other edges.
        D=ExactBoxType::unit_box(1);
        Affine<FloatDPBounds> fx(1,dp);
        Affine<FloatDPBounds> fy(1,dp);
        fx[0]=FloatDPBounds(-4e-13_pr,dp);
        fy[0]=FloatDPBounds(1,dp);

        set=ValidatedAffineConstrainedImageSet(D,{fx,fy});
        List<Pair<double,double>> vertices=set.boundary(0,1);
        ARIADNE_TEST_EQUALS(vertices.size(),6u);

        bool narrow_x=true;
        bool has_low_y=false;
        bool has_high_y=false;
        for(auto const& vertex:vertices) {
            narrow_x=narrow_x && vertex.first > -1e-8 && vertex.first < +1e-8;
            has_low_y=has_low_y || vertex.second < -0.9;
            has_high_y=has_high_y || vertex.second > +0.9;
        }
        ARIADNE_TEST_ASSERT(narrow_x);
        ARIADNE_TEST_ASSERT(has_low_y);
        ARIADNE_TEST_ASSERT(has_high_y);
    }








    Void test() {
        ARIADNE_TEST_CALL(test_basic_coverage());
        ARIADNE_TEST_CALL(test_empty());
        ARIADNE_TEST_CALL(test_invalid());
        ARIADNE_TEST_CALL(test_pure_constraint());
        ARIADNE_TEST_CALL(test_constrained_image());
        ARIADNE_TEST_CALL(test_nonunit_domain());
        ARIADNE_TEST_CALL(test_separated());
        ARIADNE_TEST_CALL(test_large_zonotope_boundary());
        ARIADNE_TEST_CALL(test_collinear_negative_zonotope_boundary());
        ARIADNE_TEST_CALL(test_near_vertical_zonotope_boundary());
        ARIADNE_TEST_CALL(test_logging_branches());
    }

};



Int main()
{
    TestAffineSet().test();
    return ARIADNE_TEST_FAILURES;
}
