/***************************************************************************
 *            test_list_set.cpp
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

// ListSet belongs to the advanced set module, not geometry.
#include "set/list_set.hpp"
#include "ariadne-testing.hpp"
#include <memory>
#include <sstream>
#include <type_traits>
using namespace Ariadne;

class TestListSet : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    Void test_list_set() {
    ExactBoxType a{{0.0_x,1.0_x},{0.0_x,1.0_x}}, b{{2.0_x,3.0_x},{2.0_x,3.0_x}};
    ListSet<ExactBoxType> sets(2u);
    ARIADNE_TESTING_PROVE(sets.empty());
    ARIADNE_TESTING_EXPECT_EQUAL(sets.dimension(),0u);
    ARIADNE_TESTING_EXPECT_EQUAL(sets.bounding_box().dimension(),0u);
    ARIADNE_TESTING_PROVE(sets.bounding_boxes().empty());
    std::ostringstream empty_text; empty_text<<sets;
    ARIADNE_TESTING_EXPECT_EQUAL(empty_text.str(),String("ListSet( )"));
    sets.adjoin(a); sets.adjoin(b);
    ARIADNE_TESTING_EXPECT_EQUAL(sets.size(),2u);
    ARIADNE_TESTING_PROVE(sets.capacity()>=sets.size());
    ARIADNE_TESTING_EXPECT_EQUAL(sets.dimension(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(sets[1],b);
    ARIADNE_TESTING_PROVE(same(cast_exact_box(sets.bounding_box()),ExactBoxType{{0.0_x,3.0_x},{0.0_x,3.0_x}}));
    ARIADNE_TESTING_EXPECT_EQUAL(sets.bounding_boxes().size(),2u);
    std::unique_ptr<ListSet<ExactBoxType>> copy(sets.clone());
    ARIADNE_TESTING_EXPECT_EQUAL(*copy,sets);
    auto summary=sets.summary();
    ARIADNE_TESTING_EXPECT_EQUAL(summary.size,2u);
    ARIADNE_TESTING_EXPECT_EQUAL(summary.dimension,2u);
    std::ostringstream text; sets._write(text); text<<summary;
    ARIADNE_TESTING_PROVE(text.str().find("s=2, d=2")!=String::npos);
    const auto& const_sets=sets;
    ListSet<ExactBoxType> from_iterators(const_sets.begin(),const_sets.end());
    ARIADNE_TESTING_EXPECT_EQUAL(from_iterators,sets);
    ListSet<ExactBoxType> from_vector(std::vector<ExactBoxType>{a,b});
    ARIADNE_TESTING_EXPECT_EQUAL(from_vector,sets);
    ListSet<UpperBoxType> converted(sets);
    ARIADNE_TESTING_EXPECT_EQUAL(converted.size(),2u);
    copy->erase(copy->begin());
    ARIADNE_TESTING_EXPECT_EQUAL(copy->pop(),b);
    ARIADNE_TESTING_PROVE(copy->empty());
    copy->push_back(a); copy->clear();
    copy->adjoin(sets);
    copy->adjoin(ListSet<ExactBoxType>());
    copy->adjoin(ListSet<ExactBoxType>(a));
    ARIADNE_TESTING_EXPECT_EQUAL(copy->size(),3u);
    ARIADNE_TESTING_EXPECT_EXCEPTION(sets.adjoin(ExactBoxType{{0.0_x,1.0_x}}));
    ARIADNE_TESTING_EXPECT_EXCEPTION(sets.adjoin(ListSet<ExactBoxType>(ExactBoxType{{0.0_x,1.0_x}})));
}

};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_list_set");
    // Use the official fixture variant lifecycle to avoid GCC's false
    // -Wmaybe-uninitialized at the pointer-to-member registration boundary.
    auto variant_1 = Ariadne::Testing::fixture_variant<TestListSet>(ariadne_test_suite, "TestListSet");
    ARIADNE_TESTING_ADD_VARIANT_TEST(variant_1, test_list_set);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
