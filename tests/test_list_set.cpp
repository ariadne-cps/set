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
#include "macro/testing.hpp"
#include <memory>
#include <sstream>
#include <type_traits>
using namespace Ariadne;

Void test_list_set() {
    ExactBoxType a{{0.0_x,1.0_x},{0.0_x,1.0_x}}, b{{2.0_x,3.0_x},{2.0_x,3.0_x}};
    ListSet<ExactBoxType> sets(2u);
    ARIADNE_TEST_ASSERT(sets.empty());
    ARIADNE_TEST_EQUALS(sets.dimension(),0u);
    ARIADNE_TEST_EQUALS(sets.bounding_box().dimension(),0u);
    ARIADNE_TEST_ASSERT(sets.bounding_boxes().empty());
    std::ostringstream empty_text; empty_text<<sets;
    ARIADNE_TEST_EQUALS(empty_text.str(),String("ListSet( )"));
    sets.adjoin(a); sets.adjoin(b);
    ARIADNE_TEST_EQUALS(sets.size(),2u);
    ARIADNE_TEST_ASSERT(sets.capacity()>=sets.size());
    ARIADNE_TEST_EQUALS(sets.dimension(),2u);
    ARIADNE_TEST_EQUALS(sets[1],b);
    ARIADNE_TEST_ASSERT(same(cast_exact_box(sets.bounding_box()),ExactBoxType{{0.0_x,3.0_x},{0.0_x,3.0_x}}));
    ARIADNE_TEST_EQUALS(sets.bounding_boxes().size(),2u);
    std::unique_ptr<ListSet<ExactBoxType>> copy(sets.clone());
    ARIADNE_TEST_EQUALS(*copy,sets);
    auto summary=sets.summary();
    ARIADNE_TEST_EQUALS(summary.size,2u);
    ARIADNE_TEST_EQUALS(summary.dimension,2u);
    std::ostringstream text; sets._write(text); text<<summary;
    ARIADNE_TEST_ASSERT(text.str().find("s=2, d=2")!=String::npos);
    const auto& const_sets=sets;
    ListSet<ExactBoxType> from_iterators(const_sets.begin(),const_sets.end());
    ARIADNE_TEST_EQUALS(from_iterators,sets);
    ListSet<ExactBoxType> from_vector(std::vector<ExactBoxType>{a,b});
    ARIADNE_TEST_EQUALS(from_vector,sets);
    ListSet<UpperBoxType> converted(sets);
    ARIADNE_TEST_EQUALS(converted.size(),2u);
    copy->erase(copy->begin());
    ARIADNE_TEST_EQUALS(copy->pop(),b);
    ARIADNE_TEST_ASSERT(copy->empty());
    copy->push_back(a); copy->clear();
    copy->adjoin(sets);
    copy->adjoin(ListSet<ExactBoxType>());
    copy->adjoin(ListSet<ExactBoxType>(a));
    ARIADNE_TEST_EQUALS(copy->size(),3u);
    ARIADNE_TEST_FAIL(sets.adjoin(ExactBoxType{{0.0_x,1.0_x}}));
    ARIADNE_TEST_FAIL(sets.adjoin(ListSet<ExactBoxType>(ExactBoxType{{0.0_x,1.0_x}})));
}


Int main() {
    ARIADNE_TEST_CALL(test_list_set());
    return ARIADNE_TEST_FAILURES;
}
