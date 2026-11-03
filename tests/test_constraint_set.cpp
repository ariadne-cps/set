/***************************************************************************
 *            test_constraint_set.cpp
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
#include "ariadne-testing.hpp"
#include <memory>
#include <sstream>
using namespace Ariadne;

class TestConstraintSet : public Ariadne::Testing::ContextFixture {
    Void check_predicate_clone(const LowerKleenean& predicate) {
    ARIADNE_TESTING_PROVE(definitely(predicate,Effort(1)));
    std::unique_ptr<Detail::LogicalInterface> copied(predicate.repr().reference()._copy());
    ARIADNE_TESTING_PROVE(Detail::definitely(copied->_check(Effort(1))));
    std::ostringstream repr;
    predicate.repr().reference()._write(repr);
    copied->_write(repr);
    ARIADNE_TESTING_PROVE(!repr.str().empty());
}

  public:
    using ContextFixture::ContextFixture;

    Void test_constraint_sets() {
    auto x=EffectiveScalarMultivariateFunction::coordinate(2u,0u);
    auto y=EffectiveScalarMultivariateFunction::coordinate(2u,1u);
    EffectiveVectorMultivariateFunction identity({x,y});
    ConstraintSet constraint(identity,RealBox{{-1,1},{-1,1}});
    ConstraintSet list_constraint(List<EffectiveConstraint>{-1<=x<=1});
    ARIADNE_TESTING_EXPECT_EXCEPTION(ConstraintSet(identity,RealBox{{-1,1}}));
    ARIADNE_TESTING_EXPECT_EXCEPTION(ConstraintSet(List<EffectiveConstraint>{x<=1,EffectiveScalarMultivariateFunction::coordinate(1u,0u)<=1}));
    ARIADNE_TESTING_EXPECT_EXCEPTION(BoundedConstraintSet(RealBox{{0,1}},identity,RealBox{{-1,1},{-1,1}}));
    ConstraintSet empty_constraint(List<EffectiveConstraint>{});
    ARIADNE_TESTING_EXPECT_EQUAL(empty_constraint.dimension(),0u);
    ARIADNE_TESTING_EXPECT_EQUAL(constraint.dimension(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(constraint.number_of_constraints(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(constraint.constraints().size(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(constraint.constraint(0u).function().argument_size(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(constraint.codomain().dimension(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(constraint.constraint_function().result_size(),2u);
    std::unique_ptr<ConstraintSet> clone(constraint.clone());
    ARIADNE_TESTING_EXPECT_EQUAL(clone->dimension(),2u);
    ExactBoxType centre{{-0.25_x,0.25_x},{-0.25_x,0.25_x}};
    ExactBoxType far{{3.0_x,4.0_x},{3.0_x,4.0_x}};
    check_predicate_clone(constraint.covers(centre));
    check_predicate_clone(constraint.overlaps(centre));
    check_predicate_clone(constraint.separated(far));
    ARIADNE_TESTING_PROVE(definitely(constraint.covers(centre,Effort(1))));
    ARIADNE_TESTING_PROVE(definitely(constraint.overlaps(centre,Effort(1))));
    ARIADNE_TESTING_PROVE(definitely(constraint.separated(far,Effort(1))));
    ARIADNE_TESTING_EXPECT_EQUAL(intersection(constraint,list_constraint).number_of_constraints(),3u);
    RealBox domain{{-2,2},{-3,3}};
    BoundedConstraintSet box_only(domain);
    BoundedConstraintSet bounded(domain,identity,RealBox{{-1,1},{-1,1}});
    BoundedConstraintSet bounded_list(domain,List<EffectiveConstraint>{-1<=x<=1});
    ARIADNE_TESTING_EXPECT_EQUAL(bounded.number_of_constraints(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(bounded.constraints().size(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(bounded.constraint(1u).function().argument_size(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(bounded.codomain().dimension(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(bounded.constraint_function().result_size(),2u);
    std::unique_ptr<BoundedConstraintSet> bounded_clone(bounded.clone());
    ARIADNE_TESTING_EXPECT_EQUAL(bounded_clone->dimension(),2u);
    check_predicate_clone(bounded.covers(centre));
    check_predicate_clone(bounded.overlaps(centre));
    check_predicate_clone(bounded.separated(far));
    ExactBoxType large{{-4.0_x,4.0_x},{-4.0_x,4.0_x}};
    check_predicate_clone(bounded.inside(large));
    ARIADNE_TESTING_PROVE(!definitely(bounded.covers(far,Effort(1))));
    ARIADNE_TESTING_PROVE(!definitely(bounded.overlaps(far,Effort(1))));
    ARIADNE_TESTING_PROVE(definitely(bounded.separated(ExactBoxType{{1.5_x,1.75_x},{1.5_x,1.75_x}},Effort(1))));
    ARIADNE_TESTING_PROVE(definitely(box_only.overlaps(centre,Effort(1))));
    ARIADNE_TESTING_PROVE(definitely(box_only.inside(large,Effort(1))));
    ARIADNE_TESTING_PROVE(definitely(subset(ExactBoxType{{-2.0_x,2.0_x},{-3.0_x,3.0_x}},cast_exact_box(bounded.bounding_box()))));
    RealBox clipped_domain{{-0.5_x,0.5_x},{-2,2}};
    auto clipped=intersection(bounded_list,clipped_domain);
    auto reverse_clipped=intersection(clipped_domain,bounded_list);
    ARIADNE_TESTING_EXPECT_EQUAL(clipped.dimension(),2u);
    ARIADNE_TESTING_PROVE(decide(clipped.domain()==clipped_domain));
    ARIADNE_TESTING_PROVE(decide(reverse_clipped.domain()==clipped_domain));
    auto combined=intersection(bounded_list,BoundedConstraintSet(clipped_domain));
    ARIADNE_TESTING_PROVE(decide(combined.domain()==clipped_domain));
    ARIADNE_TESTING_EXPECT_EQUAL(combined.number_of_constraints(),1u);
    auto combined_constraints=intersection(bounded_list,constraint);
    auto reverse_combined_constraints=intersection(constraint,bounded_list);
    ARIADNE_TESTING_PROVE(decide(combined_constraints.domain()==domain));
    ARIADNE_TESTING_PROVE(decide(reverse_combined_constraints.domain()==domain));
    ARIADNE_TESTING_EXPECT_EQUAL(combined_constraints.number_of_constraints(),3u);
    ARIADNE_TESTING_EXPECT_EQUAL(intersection(constraint,domain).dimension(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(intersection(domain,constraint).dimension(),2u);
    auto mapped=image(bounded_list,identity);
    ARIADNE_TESTING_EXPECT_EQUAL(mapped.number_of_constraints(),1u);
    ARIADNE_TESTING_EXPECT_EQUAL(mapped.dimension(),2u);
    ConstrainedImageSet converted(bounded_list);
    ARIADNE_TESTING_EXPECT_EQUAL(converted.number_of_constraints(),1u);
    ARIADNE_TESTING_EXPECT_EQUAL(class_name<ConstraintSet>(),String("EffectiveConstraintSet"));
    ARIADNE_TESTING_EXPECT_EQUAL(class_name<BoundedConstraintSet>(),String("EffectiveBoundedConstraintSet"));
    ARIADNE_TESTING_EXPECT_EQUAL(class_name<ConstrainedImageSet>(),String("EffectiveConstrainedImageSet"));
    ARIADNE_TESTING_EXPECT_EQUAL(class_name<ValidatedConstrainedImageSet>(),String("ValidatedConstrainedImageSet"));
    ARIADNE_TESTING_EXPECT_EQUAL(class_name<ValidatedAffineConstrainedImageSet>(),String("ValidatedAffineConstrainedImageSet"));
    std::ostringstream repr;
    repr<<constraint<<bounded<<bounded_list;
    ARIADNE_TESTING_PROVE(repr.str().find("BoundedConstraintSet")!=String::npos);
}

};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_constraint_set");
    ARIADNE_TESTING_ADD_TEST(TestConstraintSet, test_constraint_sets);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
