/***************************************************************************
 *            tests/fixtures/spre22_affine_subset.hpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Set.
 *
 *  Ariadne Set is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Set is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Set. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file tests/fixtures/spre22_affine_subset.hpp
 *  @brief Bit-exact SPRE22 enclosure captures 1 and 17 for Set regression.
 *  @details Extracted unchanged from Paving fixture at commit
 *  f7ad74742c07785b180916d45c664d6415760ccb, generated from the
 *  real SPRE22 evolver (horizon 1/16, max jumps 1, 48 enclosures).
 *  These are nonlinear Taylor sets with one constraint each; they are
 *  NOT the original affine nc=3/nc=4 cycling examples.
 *  Fixture is test-only and does not introduce a Paving dependency.
 */

#ifndef ARIADNE_SET_SPRE22_AFFINE_SUBSET_HPP
#define ARIADNE_SET_SPRE22_AFFINE_SUBSET_HPP

#include <bit>
#include <cstdint>
#include <vector>

#include "function/constraint.hpp"
#include "function/taylor_function.hpp"
#include "set/function_set.hpp"

namespace Ariadne::SPRE22Fixture {
inline std::vector<ValidatedConstrainedImageSet> affine_projection_samples() {
    std::vector<ValidatedConstrainedImageSet> sets;
    sets.reserve(2u);
    {
        ExactBoxType domain(6u);
        domain[0u] = ExactIntervalType(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x0ULL})), dp), FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x0ULL})), dp));
        domain[1u] = ExactIntervalType(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xc08ce80000000000ULL})), dp), FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xc08c840000000000ULL})), dp));
        domain[2u] = ExactIntervalType(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xc07a900000000000ULL})), dp), FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xc079c80000000000ULL})), dp));
        domain[3u] = ExactIntervalType(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x0ULL})), dp), FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x4014000000000000ULL})), dp));
        domain[4u] = ExactIntervalType(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x0ULL})), dp), FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x4014000000000000ULL})), dp));
        domain[5u] = ExactIntervalType(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x0ULL})), dp), FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f957ce2511a7bffULL})), dp));
        Vector<ValidatedTaylorModelDP> state_models(5u, ValidatedTaylorModelDP(6u, Sweeper<FloatDP>(new TrivialSweeper<FloatDP>(dp))));
        {
            ValidatedTaylorModelDP tm(6u, Sweeper<FloatDP>(new ThresholdSweeper<FloatDP>(dp, FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e80000000000000ULL})), dp))));
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f857ce2511a7bffULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f857ce2511a7bffULL})), dp);
            tm.set_error(FloatDPError(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3c60076a6f34f740ULL})), dp)));
            state_models[0u] = std::move(tm);
        }
        {
            ValidatedTaylorModelDP tm(6u, Sweeper<FloatDP>(new ThresholdSweeper<FloatDP>(dp, FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e80000000000000ULL})), dp))));
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 4u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e85a3637821728aULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3ed0dbd84af8afcdULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e8b0011a4efc363ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbefa0b551f9be422ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3eb2ff2903349777ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf395c9e89e83a2cULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbef42a642f9b0ae9ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f63f8628ba973dfULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3ec39c246066c0c8ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f9a0e0cb883ca8aULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e828fb8f55629d3ULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf047b646489da34ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f9fede74eea1d24ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3eb3d0782a2223e3ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f9a748cf0301afeULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x4018fffada66e270ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xc08cb5c52f655309ULL})), dp);
            tm.set_error(FloatDPError(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3ebd542b73150cbbULL})), dp)));
            state_models[1u] = std::move(tm);
        }
        {
            ValidatedTaylorModelDP tm(6u, Sweeper<FloatDP>(new ThresholdSweeper<FloatDP>(dp, FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e80000000000000ULL})), dp))));
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3ed0e06b9e9352d0ULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e8f2e0358cd5d4fULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbee7ecfbebffdbc0ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf396373472314f8ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbeb2f4c21280b64eULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbef743086286ce47ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f52210918519521ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f9a0dd536fa46b6ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbec3917deee0259cULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf07a0926c913660ULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbe801171e4c57a21ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f9d29481e580a0aULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f9a74710ad0bc3fULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbeb3c5bc8daaa00cULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x4018fffa100fe404ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xc07a2b8fef2b821eULL})), dp);
            tm.set_error(FloatDPError(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3ebba350676c5b86ULL})), dp)));
            state_models[2u] = std::move(tm);
        }
        {
            ValidatedTaylorModelDP tm(6u, Sweeper<FloatDP>(new ThresholdSweeper<FloatDP>(dp, FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e80000000000000ULL})), dp))));
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 4u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbe9f1dc49a04cebdULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbee8468ddc4fdb30ULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbea370b78401824bULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f12b9efd86def85ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbedbe57bdb78cb5bULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f5242c33e27aa70ULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbe99af7262447336ULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f0d3ee285064fbeULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf7cf487d2e5f0d0ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f2c54e681c68de7ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbfb2e4babb7abd3bULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3ee9ff97ccfc1308ULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf6e0b8e7ce1ec63ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3fddbdfe00766a0dULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f2d379ba1a47914ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x4003668bc018ee0cULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3eead055d62d3b91ULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf6e81c1126e45bbULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x4007c65fb6137a4aULL})), dp);
            tm.set_error(FloatDPError(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3ed3cc5a4f12522dULL})), dp)));
            state_models[3u] = std::move(tm);
        }
        {
            ValidatedTaylorModelDP tm(6u, Sweeper<FloatDP>(new ThresholdSweeper<FloatDP>(dp, FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e80000000000000ULL})), dp))));
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 4u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbe8bdc447e446484ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbee84d18044e7017ULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbea673202c5d2d56ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f0157678958b420ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f5247b845d26c01ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3edbda235ad20449ULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f10e2e70150e742ULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e963d194b03e899ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf6a920fbd949e91ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbfb2e9d1784b2968ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf2c456502fcab61ULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf71544469382ffdULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbee68184cd7ba03fULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3fcb001fa16a4a42ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x4003666269eeeb52ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf2d27be1a344f44ULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf7198839e3ed683ULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbee736401efab9c0ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x4005b6b7f149d7a3ULL})), dp);
            tm.set_error(FloatDPError(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3ed1f89c2cb51459ULL})), dp)));
            state_models[4u] = std::move(tm);
        }
        ValidatedVectorMultivariateTaylorFunctionModelDP state_map(domain, state_models);
        List<ValidatedConstraint> constraints;
        {
            ValidatedTaylorModelDP tm(6u, Sweeper<FloatDP>(new ThresholdSweeper<FloatDP>(dp, FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e80000000000000ULL})), dp))));
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f857ce2511a7bffULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf857ce2511a7bffULL})), dp);
            tm.set_error(FloatDPError(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3c64076a6f34f740ULL})), dp)));
            const ValidatedScalarMultivariateTaylorFunctionModelDP fn(domain, tm);
            constraints.append(ValidatedConstraint(ValidatedNumber(FloatDPBounds(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xfff0000000000000ULL})), dp), FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xfff0000000000000ULL})), dp))), ValidatedScalarMultivariateFunction(fn), ValidatedNumber(FloatDPBounds(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x0ULL})), dp), FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x0ULL})), dp)))));
        }
        sets.emplace_back(domain, ValidatedVectorMultivariateFunctionModelDP(state_map), constraints);
    }
    {
        ExactBoxType domain(6u);
        domain[0u] = ExactIntervalType(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x0ULL})), dp), FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x0ULL})), dp));
        domain[1u] = ExactIntervalType(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xc08ce80000000000ULL})), dp), FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xc08c840000000000ULL})), dp));
        domain[2u] = ExactIntervalType(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xc07a900000000000ULL})), dp), FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xc079c80000000000ULL})), dp));
        domain[3u] = ExactIntervalType(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x0ULL})), dp), FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x4014000000000000ULL})), dp));
        domain[4u] = ExactIntervalType(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x0ULL})), dp), FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x4014000000000000ULL})), dp));
        domain[5u] = ExactIntervalType(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x0ULL})), dp), FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f957ce2511a7c30ULL})), dp));
        Vector<ValidatedTaylorModelDP> state_models(5u, ValidatedTaylorModelDP(6u, Sweeper<FloatDP>(new TrivialSweeper<FloatDP>(dp))));
        {
            ValidatedTaylorModelDP tm(6u, Sweeper<FloatDP>(new ThresholdSweeper<FloatDP>(dp, FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e80000000000000ULL})), dp))));
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f857ce2511a7c30ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3fa01da9bcd3dd0cULL})), dp);
            tm.set_error(FloatDPError(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3c8557ce2511a7c3ULL})), dp)));
            state_models[0u] = std::move(tm);
        }
        {
            ValidatedTaylorModelDP tm(6u, Sweeper<FloatDP>(new ThresholdSweeper<FloatDP>(dp, FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e80000000000000ULL})), dp))));
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 4u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e845d2ce892100fULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3ecfb9fa0914a8faULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e8afffba32d5b19ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbef884e6aec1719aULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3eb1dfae8a75f866ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf37dfd9459f6d16ULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbef2f97ac0fda986ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f62ca7e44d6bbf9ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3edbaf2e20bb2f7cULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f9884173a8d72abULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e9a3827ebf4aea7ULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf1dd08512708398ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3fa4cf060bd3f1a4ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3ee5749c9f1cae12ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3fb340b25dea828fULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e9a7ad202fd5a50ULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x4018ffd29a2aa473ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xc08cb531e54c4be1ULL})), dp);
            tm.set_error(FloatDPError(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3ece49f5ef4fc16eULL})), dp)));
            state_models[1u] = std::move(tm);
        }
        {
            ValidatedTaylorModelDP tm(6u, Sweeper<FloatDP>(new ThresholdSweeper<FloatDP>(dp, FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e80000000000000ULL})), dp))));
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3ecfc20ceb31b45cULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e8f2de607f99e1aULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbee68375a46a9ee0ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf37e6469fabe32aULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbeb1d632331f9324ULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbef5e2ec35027c9fULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f510c12dbcff49aULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f988375a38270cdULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbedba0120ed1efa5ULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf2132250eed83b4ULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbe96b262e35fa0f4ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3fa0c75621547887ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3fb34074da942d8fULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbee568ea8e7ac1deULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x4018ffcba1b97845ULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbe96ec26b2945c58ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xc07a2a94e166bcd9ULL})), dp);
            tm.set_error(FloatDPError(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3ecc937098494c8bULL})), dp)));
            state_models[2u] = std::move(tm);
        }
        {
            ValidatedTaylorModelDP tm(6u, Sweeper<FloatDP>(new ThresholdSweeper<FloatDP>(dp, FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e80000000000000ULL})), dp))));
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 4u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbe9d48aaeb59fa4cULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbee6d7938342b782ULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbea370a7abad4b05ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f11a13a79aa6eceULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbed97054fe25a174ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f51309de7d56724ULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbe99af5a3beaacd9ULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f0b87cdcd865552ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf7b4265b5a6d903ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f28fbd34026dc9dULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbfb1c91240d4156aULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3ee6ed9b9cd93c16ULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf6c454373f3c231ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3fdbfc611440e858ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f449dfbade78bd4ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x40024133d4b6cd9cULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f02eba10aa77639ULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf863384a26b92b9ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x400efd1b3a4635dfULL})), dp);
            tm.set_error(FloatDPError(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3ee40a3c0aad51a6ULL})), dp)));
            state_models[3u] = std::move(tm);
        }
        {
            ValidatedTaylorModelDP tm(6u, Sweeper<FloatDP>(new ThresholdSweeper<FloatDP>(dp, FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e80000000000000ULL})), dp))));
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 4u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbe8a374e512a7d01ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbee6dd57dc6593c7ULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbea6730b10d1a044ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 3u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f0051681802a37cULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f5135490c85eeb0ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3ed96652cc137b83ULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f0fcac0b0841f18ULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e963d072aeb6810ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 2u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf6900924177f0a9ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbfb1cddbf8d7b0bbULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf28ee26cd0463c0ULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf704df5c10bb4baULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbee3d9043c68da2bULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3fc963a189347c3fULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 1u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x400240bb80f5d2d8ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 1u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf4492bb349c20fcULL})), dp);
            tm[MultiIndex{0u, 0u, 1u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf899c122bd137a7ULL})), dp);
            tm[MultiIndex{0u, 1u, 0u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf0060f1f9bbfbc8ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x4008fcb14bbf680eULL})), dp);
            tm.set_error(FloatDPError(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3ee23af00cc2d0d1ULL})), dp)));
            state_models[4u] = std::move(tm);
        }
        ValidatedVectorMultivariateTaylorFunctionModelDP state_map(domain, state_models);
        List<ValidatedConstraint> constraints;
        {
            ValidatedTaylorModelDP tm(6u, Sweeper<FloatDP>(new ThresholdSweeper<FloatDP>(dp, FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3e80000000000000ULL})), dp))));
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 1u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3f857ce2511a7c30ULL})), dp);
            tm[MultiIndex{0u, 0u, 0u, 0u, 0u, 0u}] = FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xbf857ce2511a7c30ULL})), dp);
            tm.set_error(FloatDPError(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x3c655f3894469f0cULL})), dp)));
            const ValidatedScalarMultivariateTaylorFunctionModelDP fn(domain, tm);
            constraints.append(ValidatedConstraint(ValidatedNumber(FloatDPBounds(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xfff0000000000000ULL})), dp), FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0xfff0000000000000ULL})), dp))), ValidatedScalarMultivariateFunction(fn), ValidatedNumber(FloatDPBounds(FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x0ULL})), dp), FloatDP(cast_exact(std::bit_cast<double>(std::uint64_t{0x0ULL})), dp)))));
        }
        sets.emplace_back(domain, ValidatedVectorMultivariateFunctionModelDP(state_map), constraints);
    }
    return sets;
}
} // namespace Ariadne::SPRE22Fixture

#endif // ARIADNE_SET_SPRE22_AFFINE_SUBSET_HPP
