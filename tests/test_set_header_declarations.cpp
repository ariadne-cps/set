// Set headers must be self-contained without paving declarations or includes.
#include "set/function_set.hpp"
#include "set/affine_set.hpp"
#include <type_traits>
#include <utility>

using namespace Ariadne;
static_assert(std::is_same_v<decltype(std::declval<const ValidatedConstrainedImageSet&>().dimension()),DimensionType>);
static_assert(std::is_same_v<decltype(std::declval<const ValidatedAffineConstrainedImageSet&>().bounding_box()),UpperBoxType>);

int main() { return 0; }
