/** @file tests/fixtures/affine_projection_reader.hpp
 *  @brief Test-only reader for a bit-exact high-dimensional constrained-image record.
 *
 *  This reads a fixed binary record with no runtime external dependency.
 */
#ifndef ARIADNE_SET_TEST_AFFINE_PROJECTION_READER_HPP
#define ARIADNE_SET_TEST_AFFINE_PROJECTION_READER_HPP

#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "function/constraint.hpp"
#include "function/taylor_function.hpp"
#include "set/function_set.hpp"

namespace Ariadne::AffineProjectionFixture {

// A high-bit tag is a run of zero bytes; other tags precede literal bytes.
// Encoding does not alter any of the original binary64 coefficients.
inline std::vector<std::uint8_t> decode_rle_base64(std::string_view encoded) {
    std::vector<std::uint8_t> packed;
    unsigned accumulator = 0u;
    unsigned bits = 0u;
    for (char ch: encoded) {
        unsigned digit;
        if (ch >= 'A' && ch <= 'Z')
            digit = static_cast<unsigned>(ch - 'A');
        else if (ch >= 'a' && ch <= 'z')
            digit = static_cast<unsigned>(ch - 'a') + 26u;
        else if (ch >= '0' && ch <= '9')
            digit = static_cast<unsigned>(ch - '0') + 52u;
        else if (ch == '+')
            digit = 62u;
        else if (ch == '/')
            digit = 63u;
        else if (ch == '=')
            break;
        else
            throw std::runtime_error("Invalid captured projection fixture encoding");
        accumulator = (accumulator << 6u) | digit;
        bits += 6u;
        if (bits >= 8u) {
            bits -= 8u;
            packed.push_back(static_cast<std::uint8_t>((accumulator >> bits) & 255u));
        }
    }

    std::vector<std::uint8_t> bytes;
    for (std::size_t pos = 0u; pos < packed.size();) {
        const unsigned marker = packed[pos++];
        const std::size_t count = (marker & 127u) + 1u;
        if ((marker & 128u) != 0u) {
            bytes.insert(bytes.end(), count, std::uint8_t{0});
        } else {
            if (count > packed.size() - pos)
                throw std::runtime_error("Truncated captured projection RLE");
            for (std::size_t k = 0u; k != count; ++k)
                bytes.push_back(packed[pos++]);
        }
    }
    return bytes;
}

class BinaryReader {
  public:
    explicit BinaryReader(const std::vector<std::uint8_t>& bytes) : _bytes(bytes) {}

    std::uint64_t read_uint(unsigned width) {
        if (width > 8u || width > _bytes.size() - _offset)
            throw std::runtime_error("Truncated captured projection record");
        std::uint64_t value = 0u;
        for (unsigned k = 0u; k < width; ++k)
            value |= std::uint64_t(_bytes[_offset++]) << (8u * k);
        return value;
    }

    FloatDP read_float(bool allow_infinity = false) {
        const double value = std::bit_cast<double>(read_uint(8u));
        if (std::isnan(value) || (!allow_infinity && !std::isfinite(value)))
            throw std::runtime_error("Nonfinite captured projection coefficient");
        if (std::isinf(value))
            return FloatDP::inf(std::signbit(value) ? Sign::NEGATIVE : Sign::POSITIVE, dp);
        return FloatDP(ExactDouble(value), dp);
    }

    FloatDPBounds read_bounds() {
        const FloatDP lower = read_float(true);
        const FloatDP upper = read_float(true);
        return FloatDPBounds(lower, upper);
    }

    bool at_end() const { return _offset == _bytes.size(); }

  private:
    const std::vector<std::uint8_t>& _bytes;
    std::size_t _offset = 0u;
};

inline ValidatedTaylorModelDP read_taylor_model(BinaryReader& reader, SizeType dimension) {
    const SizeType arguments = reader.read_uint(4u);
    if (arguments != dimension || arguments > 4096u)
        throw std::runtime_error("captured projection model dimension mismatch");
    // The captured critical record uses ThresholdSweeper, tag 1.
    if (reader.read_uint(1u) != 1u)
        throw std::runtime_error("Unsupported captured projection sweeper");
    const FloatDP threshold = reader.read_float();
    const FloatDP error = reader.read_float();
    const SizeType terms = reader.read_uint(4u);
    if (terms > 10000000u)
        throw std::runtime_error("Excessive captured projection term count");
    Expansion<MultiIndex, FloatDP> expansion(arguments, dp);
    expansion.reserve(terms);
    for (SizeType i = 0u; i != terms; ++i) {
        MultiIndex index(arguments);
        for (SizeType j = 0u; j != arguments; ++j) {
            const std::uint64_t exponent = reader.read_uint(4u);
            if (exponent > static_cast<std::uint64_t>(
                               std::numeric_limits<MultiIndexData::IndexType>::max()))
                throw std::runtime_error("captured projection exponent exceeds MultiIndex index range");
            index[j] = static_cast<MultiIndexData::IndexType>(exponent);
        }
        expansion.append(index, reader.read_float());
    }
    return ValidatedTaylorModelDP(
        expansion, FloatDPError(error),
        Sweeper<FloatDP>(new ThresholdSweeper<FloatDP>(dp, threshold)));
}

struct CapturedSet {
    ValidatedConstrainedImageSet set;
    SizeType x;
    SizeType y;
};

// Decodes one original record including the 8-byte length prefix.
inline CapturedSet decode_high_constraint_record(std::string_view encoded) {
    const auto data = decode_rle_base64(encoded);
    BinaryReader reader(data);
    const std::uint64_t length = reader.read_uint(8u);
    if (length + 8u != data.size() || reader.read_uint(8u) != 11591u)
        throw std::runtime_error("captured projection critical record index/length mismatch");
    const SizeType projection_n = reader.read_uint(4u);
    const SizeType projection_x = reader.read_uint(4u);
    const SizeType projection_y = reader.read_uint(4u);
    // This is a legacy v1 capture: the old writer added one extra
    // unprojected time component. The original projection is unchanged.
    if (projection_n != 8u || projection_x != 1u || projection_y != 2u)
        throw std::runtime_error("Unexpected captured projection");
    const SizeType dimension = reader.read_uint(4u);
    if (dimension != 26u)
        throw std::runtime_error("Unexpected captured projection parameter count");
    ExactBoxType domain(dimension);
    for (SizeType i = 0u; i != dimension; ++i) {
        const FloatDP lower = reader.read_float();
        const FloatDP upper = reader.read_float();
        domain[i] = ExactIntervalType(lower, upper);
    }
    const SizeType count = reader.read_uint(4u);
    if (count != projection_n + 1u)
        throw std::runtime_error("Expected legacy captured projection extra time component");
    Vector<ValidatedTaylorModelDP> models(
        count, ValidatedTaylorModelDP(dimension, TrivialSweeper<FloatDP>(dp)));
    for (SizeType i = 0u; i != count; ++i)
        models[i] = read_taylor_model(reader, dimension);
    const ValidatedVectorMultivariateTaylorFunctionModelDP mapping(domain, models);

    const SizeType nc = reader.read_uint(4u);
    if (nc != 4u)
        throw std::runtime_error("Unexpected captured projection constraint count");
    List<ValidatedConstraint> constraints;
    for (SizeType i = 0u; i != nc; ++i) {
        const FloatDPBounds lower = reader.read_bounds();
        const auto model = read_taylor_model(reader, dimension);
        const FloatDPBounds upper = reader.read_bounds();
        constraints.append(ValidatedConstraint(
            lower, ValidatedScalarMultivariateTaylorFunctionModelDP(domain, model), upper));
    }
    if (!reader.at_end())
        throw std::runtime_error("Trailing captured projection record bytes");
    return {ValidatedConstrainedImageSet(
                domain, ValidatedVectorMultivariateFunctionModelDP(mapping), constraints),
            projection_x, projection_y};
}

} // namespace Ariadne::AffineProjectionFixture

#endif // ARIADNE_SET_TEST_AFFINE_PROJECTION_READER_HPP
