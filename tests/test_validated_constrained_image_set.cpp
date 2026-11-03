/***************************************************************************
 *            test_validated_constrained_image_set.cpp
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
#include "set/affine_set.hpp"
#include "function/taylor_function.hpp"
#include "ariadne-testing.hpp"
#include <memory>
#include <sstream>
using namespace Ariadne;

class TestValidatedConstrainedImageSetSuite : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    Void test_default_image_sets() {
    ConstrainedImageSet effective;
    ValidatedConstrainedImageSet validated;
    ARIADNE_TESTING_EXPECT_EQUAL(effective.dimension(),0u);
    ARIADNE_TESTING_EXPECT_EQUAL(validated.dimension(),0u);
    ARIADNE_TESTING_EXPECT_EQUAL(effective.number_of_parameters(),0u);
    ARIADNE_TESTING_EXPECT_EQUAL(validated.number_of_parameters(),0u);
    ARIADNE_TESTING_EXPECT_EQUAL(effective.number_of_constraints(),0u);
    ARIADNE_TESTING_EXPECT_EQUAL(validated.number_of_constraints(),0u);
    ARIADNE_TESTING_EXPECT_EQUAL(effective.bounding_box().dimension(),0u);
    ARIADNE_TESTING_EXPECT_EQUAL(validated.bounding_box().dimension(),0u);
}

Void test_validated_image_operations() {
    auto x=ValidatedScalarMultivariateFunction::coordinate(2u,0u);
    auto y=ValidatedScalarMultivariateFunction::coordinate(2u,1u);
    ValidatedVectorMultivariateFunction mapping({x+y,x-y});
    ExactBoxType domain{{-1.0_x,1.0_x},{-2.0_x,2.0_x}};
    ValidatedConstrainedImageSet set(domain,mapping);
    set.new_parameter_constraint(x>=0);
    set.new_space_constraint(x<=ValidatedNumber(2));
    ARIADNE_TESTING_EXPECT_EQUAL(set.constraints().size(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(set.constraint(0u).function().argument_size(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(set.constraint_function().result_size(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(set.constraint_bounds().dimension(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(set.number_of_parameters(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(set.dimension(),2u);
    std::unique_ptr<ValidatedConstrainedImageSet> copy(set.clone());
    ARIADNE_TESTING_EXPECT_EQUAL(copy->number_of_constraints(),2u);
    auto pieces=set.split();
    ARIADNE_TESTING_EXPECT_EQUAL(pieces.first.number_of_constraints(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(pieces.second.number_of_constraints(),2u);
    ARIADNE_TESTING_PROVE(same(pieces.first.domain(),ExactBoxType{{-1.0_x,1.0_x},{-2.0_x,0.0_x}}));
    auto explicit_pieces=set.split(0u);
    ARIADNE_TESTING_PROVE(same(explicit_pieces.second.domain(),ExactBoxType{{0.0_x,1.0_x},{-2.0_x,2.0_x}}));
    auto restricted=set.restriction(ExactBoxType{{0.0_x,0.5_x},{-1.0_x,1.0_x}});
    ARIADNE_TESTING_EXPECT_EQUAL(restricted.number_of_constraints(),2u);
    auto affine=set.affine_approximation();
    auto enclosing_affine=set.affine_over_approximation();
    ARIADNE_TESTING_EXPECT_EQUAL(affine.dimension(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(enclosing_affine.dimension(),2u);
    set.reduce();
    ARIADNE_TESTING_PROVE(set.reduced_domain()[0].lower_bound()>=FloatDP(0,dp));
    ARIADNE_TESTING_PROVE(!definitely(set.is_empty()));
    ValidatedConstrainedImageSet impossible(domain,mapping,List<ValidatedConstraint>{x>=3});
    ARIADNE_TESTING_PROVE(definitely(impossible.is_empty()));
    ValidatedConstrainedImageSet box_set(domain);
    ARIADNE_TESTING_PROVE(definitely(box_set.inside(ExactBoxType{{-3.0_x,3.0_x},{-3.0_x,3.0_x}})));
    ARIADNE_TESTING_PROVE(definitely(box_set.separated(ExactBoxType{{3.0_x,4.0_x},{3.0_x,4.0_x}})));
    ARIADNE_TESTING_PROVE(definitely(box_set.overlaps(ExactBoxType{{-0.5_x,0.5_x},{-0.5_x,0.5_x}})));
    auto transformed=image(box_set,mapping);
    ARIADNE_TESTING_EXPECT_EQUAL(transformed.dimension(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(transformed.bounding_box()[0].lower_bound().raw(),FloatDP(-3,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(transformed.bounding_box()[0].upper_bound().raw(),FloatDP(3,dp));
    transformed.apply(ValidatedVectorMultivariateFunction::identity(2u));
    ARIADNE_TESTING_EXPECT_EQUAL(transformed.dimension(),2u);
    auto wrong=ValidatedScalarMultivariateFunction::coordinate(1u,0u);
    ARIADNE_TESTING_EXPECT_EXCEPTION(ValidatedConstrainedImageSet(ExactBoxType{{0.0_x,1.0_x}},mapping));
    ARIADNE_TESTING_EXPECT_EXCEPTION(ValidatedConstrainedImageSet(ExactBoxType{{0.0_x,1.0_x}},mapping,List<ValidatedConstraint>{}));
    ARIADNE_TESTING_EXPECT_EXCEPTION(set.new_parameter_constraint(wrong<=ValidatedNumber(1)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(set.new_space_constraint(wrong<=ValidatedNumber(1)));
    auto taylor=ValidatedVectorMultivariateTaylorFunctionModelDP(domain,mapping,ThresholdSweeper<FloatDP>(dp,1e-12));
    ValidatedConstrainedImageSet model_set(domain,ValidatedVectorMultivariateFunctionModelDP(taylor));
    ARIADNE_TESTING_EXPECT_EQUAL(model_set.dimension(),2u);
    ARIADNE_TESTING_EXPECT_EXCEPTION(ValidatedConstrainedImageSet(ExactBoxType{{0.0_x,1.0_x}},ValidatedVectorMultivariateFunctionModelDP(taylor)));
    std::ostringstream repr;
    repr<<set;
    ARIADNE_TESTING_PROVE(repr.str().find("ValidatedConstrainedImageSet")!=String::npos);
}

Void test_one_sided_satisfaction() {
    auto x=EffectiveScalarMultivariateFunction::coordinate(1u,0u);
    EffectiveVectorMultivariateFunction identity({x});
    ConstrainedImageSet effective(RealBox{{-2,2}},identity,List<EffectiveConstraint>{-0.5_x<=x<=0.5_x});
    ValidatedConstrainedImageSet validated(ExactBoxType{{-2.0_x,2.0_x}},identity,List<ValidatedConstraint>{-0.5_x<=x<=0.5_x});
    ARIADNE_TESTING_PROVE(definitely(effective.satisfies(x>=-1,Effort(1))));
    ARIADNE_TESTING_PROVE(definitely(validated.satisfies(x>=-1)));
    ARIADNE_TESTING_PROVE(definitely(effective.satisfies(x<=1,Effort(1))));
    ARIADNE_TESTING_PROVE(definitely(validated.satisfies(x<=1)));
    ARIADNE_TESTING_PROVE(definitely(effective.satisfies(-1<=x<=1,Effort(1))));
    ARIADNE_TESTING_PROVE(definitely(validated.satisfies(-1<=x<=1)));
    ARIADNE_TESTING_PROVE(definitely(!effective.satisfies(x>=1,Effort(1))));
    ARIADNE_TESTING_PROVE(definitely(!validated.satisfies(x>=1)));
    ARIADNE_TESTING_PROVE(definitely(!effective.satisfies(x<=-1,Effort(1))));
    ARIADNE_TESTING_PROVE(definitely(!validated.satisfies(x<=-1)));
    ARIADNE_TESTING_PROVE(definitely(effective.satisfies(x<=3,Effort(1))));
    ARIADNE_TESTING_PROVE(definitely(validated.satisfies(x<=3)));
    ConstrainedImageSet box(RealBox{{-2,2}});
    auto two_variables=EffectiveVectorMultivariateFunction::identity(2u);
    ARIADNE_TESTING_EXPECT_EXCEPTION(ConstrainedImageSet(RealBox{{0,1}},two_variables));
    ARIADNE_TESTING_EXPECT_EXCEPTION(ConstrainedImageSet(RealBox{{0,1}},two_variables,List<EffectiveConstraint>{}));
    ARIADNE_TESTING_EXPECT_EXCEPTION(box.new_parameter_constraint(EffectiveScalarMultivariateFunction::coordinate(2u,0u)<=1));
    ARIADNE_TESTING_EXPECT_EXCEPTION(box.new_space_constraint(EffectiveScalarMultivariateFunction::coordinate(2u,0u)<=1));
    std::unique_ptr<ConstrainedImageSet> clone(box.clone());
    ARIADNE_TESTING_EXPECT_EQUAL(clone->dimension(),1u);
    auto split_box=box.split();
    ARIADNE_TESTING_EXPECT_EQUAL(split_box.first.dimension(),1u);
    ARIADNE_TESTING_PROVE(definitely(box.inside(ExactBoxType{{-3.0_x,3.0_x}},Effort(1))));
    ARIADNE_TESTING_PROVE(definitely(box.inside(ExactBoxType{{-3.0_x,3.0_x}}),Effort(1)));
    ARIADNE_TESTING_EXPECT_EQUAL(image(box,identity).dimension(),1u);
}
Void test_join_image_sets() {
    auto x=ValidatedScalarMultivariateFunction::coordinate(1u,0u);
    ValidatedVectorMultivariateFunction f1({x,x*x});
    ValidatedVectorMultivariateFunction f2({x+ValidatedNumber(1),x*x+ValidatedNumber(FloatDPBounds(0.125_x,0.25_x,dp))});
    ExactBoxType left{{-1.0_x,0.0_x}},right{{0.0_x,1.0_x}};
    ValidatedConstrainedImageSet first(left,f1),second(right,f2);
    ValidatedVectorMultivariateTaylorFunctionModelDP tm1(left,f1,TrivialSweeper<FloatDP>(dp));
    ValidatedVectorMultivariateTaylorFunctionModelDP tm2(right,f2,TrivialSweeper<FloatDP>(dp));
    auto modeled_join=join(ValidatedConstrainedImageSet(left,ValidatedVectorMultivariateFunctionModelDP(tm1)),ValidatedConstrainedImageSet(right,ValidatedVectorMultivariateFunctionModelDP(tm2)));
    ARIADNE_TESTING_PROVE(same(modeled_join.domain(),ExactBoxType{{-1.0_x,1.0_x}}));
    auto combined=join(first,second);
    ARIADNE_TESTING_PROVE(same(combined.domain(),ExactBoxType{{-1.0_x,1.0_x}}));
    for(Int k=-4;k!=5;++k) {
        ExactDouble t=ExactDouble(static_cast<double>(k)/4.0);
        ExactBoxType point{{t,t}};
        const auto& original=(k<=0 ? f1 : f2);
        ARIADNE_TESTING_PROVE(definitely(subset(cast_exact_box(apply(original,point)),cast_exact_box(apply(combined.function(),point)))));
    }
    ARIADNE_TESTING_PROVE(same(join(second,first).domain(),combined.domain()));
    ARIADNE_TESTING_PROVE(same(join(first,first).domain(),left));
    ARIADNE_TESTING_EXPECT_EXCEPTION(join(first,ValidatedConstrainedImageSet(ExactBoxType{{1.0_x,2.0_x}},f1)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(join(ValidatedConstrainedImageSet(ExactBoxType{{-2.0_x,-1.0_x}},f1),second));
    ValidatedConstrainedImageSet constrained(left,f1,List<ValidatedConstraint>{x<=ValidatedNumber(0)});
    ARIADNE_TESTING_EXPECT_EXCEPTION(join(constrained,second));
    ARIADNE_TESTING_EXPECT_EXCEPTION(join(first,constrained));
    ARIADNE_TESTING_EXPECT_EXCEPTION(join(first,ValidatedConstrainedImageSet(left,ValidatedVectorMultivariateFunction({x}))));
    auto xy=ValidatedVectorMultivariateFunction::identity(2u);
    ARIADNE_TESTING_EXPECT_EXCEPTION(join(first,ValidatedConstrainedImageSet(ExactBoxType{{-1.0_x,1.0_x},{-1.0_x,1.0_x}},xy)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(join(ValidatedConstrainedImageSet(ExactBoxType{{-1.0_x,0.0_x},{-1.0_x,0.0_x}},xy),ValidatedConstrainedImageSet(ExactBoxType{{0.0_x,1.0_x},{0.0_x,1.0_x}},xy)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(first.restriction(ExactBoxType{{-2.0_x,0.0_x}}));
    ARIADNE_TESTING_EXPECT_EXCEPTION(first.split(1u));
    ARIADNE_TESTING_EXPECT_EXCEPTION(ValidatedConstrainedImageSet().split());
    ConstrainedImageSet effective;
    ARIADNE_TESTING_EXPECT_EXCEPTION(effective.split());
    ARIADNE_TESTING_EXPECT_EXCEPTION(effective.split(0u));
    auto e=EffectiveScalarMultivariateFunction::coordinate(1u,0u);
    BoundedConstraintSet bounded(RealBox{{-1,1}},List<EffectiveConstraint>{e<=0});
    ARIADNE_TESTING_EXPECT_EXCEPTION(image(bounded,EffectiveVectorMultivariateFunction::identity(2u)));
    ConstrainedImageSet converted(bounded);
    ARIADNE_TESTING_EXPECT_EQUAL(converted.number_of_constraints(),1u);
    auto predicate=converted.inside(ExactBoxType{{-2.0_x,2.0_x}});
    std::unique_ptr<Detail::LogicalInterface> predicate_copy(predicate.repr().reference()._copy());
    ARIADNE_TESTING_PROVE(Detail::definitely(predicate_copy->_check(Effort(1))));
    std::ostringstream text; predicate_copy->_write(text);
    ARIADNE_TESTING_PROVE(!text.str().empty());
}

    Void test_split_chooses_widest_parameter() {
    // With three coordinates, the largest domain radius belongs to coordinate 1.
    ConstrainedImageSet effective(RealBox{{-1,1},{-3,3},{-2,2}});
    auto split_effective=effective.split();
    ARIADNE_TESTING_EXPECT_EQUAL(split_effective.first.number_of_parameters(),3u);
    ARIADNE_TESTING_EXPECT_EQUAL(split_effective.second.number_of_parameters(),3u);
    ARIADNE_TESTING_PROVE(decide(split_effective.first.domain()[1].upper_bound()<Real(1)));
    ARIADNE_TESTING_PROVE(decide(split_effective.second.domain()[1].lower_bound()>Real(-1)));

    // With the first coordinate widest, neither later coordinate should be chosen.
    ConstrainedImageSet first_widest(RealBox{{-4,4},{-1,1}});
    auto split_first=first_widest.split();
    ARIADNE_TESTING_PROVE(decide(split_first.first.domain()[0].upper_bound()<Real(1)));
    ARIADNE_TESTING_PROVE(decide(split_first.second.domain()[0].lower_bound()>Real(-1)));

    // A tie in a validated domain should be resolved by the first coordinate.
    ValidatedConstrainedImageSet validated(ExactBoxType::unit_box(2u));
    auto split_validated=validated.split();
    ARIADNE_TESTING_PROVE(same(split_validated.first.domain(),
        ExactBoxType{{-1.0_x,0.0_x},{-1.0_x,1.0_x}}));
    ARIADNE_TESTING_PROVE(same(split_validated.second.domain(),
        ExactBoxType{{0.0_x,1.0_x},{-1.0_x,1.0_x}}));
}

};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_validated_constrained_image_set");
    ARIADNE_TESTING_ADD_TEST(TestValidatedConstrainedImageSetSuite, test_split_chooses_widest_parameter);
    ARIADNE_TESTING_ADD_TEST(TestValidatedConstrainedImageSetSuite, test_join_image_sets);
    ARIADNE_TESTING_ADD_TEST(TestValidatedConstrainedImageSetSuite, test_default_image_sets);
    ARIADNE_TESTING_ADD_TEST(TestValidatedConstrainedImageSetSuite, test_validated_image_operations);
    ARIADNE_TESTING_ADD_TEST(TestValidatedConstrainedImageSetSuite, test_one_sided_satisfaction);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
