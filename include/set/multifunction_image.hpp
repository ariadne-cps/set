/***************************************************************************
 *            set/multifunction_image.hpp
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

/** @file set/multifunction_image.hpp
 *  @brief Multifunction patches and models producing constrained-image sets.
 */

// Multifunctions with constrained-image results belong to set, not function.
#ifndef ARIADNE_SET_MULTIFUNCTION_IMAGE_HPP
#define ARIADNE_SET_MULTIFUNCTION_IMAGE_HPP

#include "function/function_model.hpp"
#include "function/multifunction.hpp"
#include "set/function_set.hpp"

namespace Ariadne {

using ValidatedImageSet = ValidatedConstrainedImageSet;

/** @brief Forward-declare a multifunction patch parameterized by function signature and resulting set. */
template<class P, class SIG, template<class, class> class SET = LocatedSet>
class MultifunctionPatch;
template<class P>
using VectorMultivariateMultifunctionPatch = MultifunctionPatch<P, RealVector(RealVector)>;
using ValidatedVectorMultivariateMultifunctionPatch = VectorMultivariateMultifunctionPatch<ValidatedTag>;

/** @brief Represent a vector multifunction using a function with bounded error parameters. */
template<class P>
class MultifunctionPatch<P, RealVector(RealVector)> {
    using ARG = RealVector;
    using RES = RealVector;
    using SIG = RES(ARG);
    Function<P, SIG> _f;
    BoxDomainType _pdom;

  public:
    /** @brief Construct a multifunction patch from a function and its parameter domain. */
    MultifunctionPatch(VectorMultivariateFunction<P> f, BoxDomainType pdom) : _f(f), _pdom(pdom) {}

    /** @brief Return the number of free arguments, excluding error parameters. */
    SizeType argument_size() const { return this->_f.argument_size() - this->_pdom.dimension(); }
    /** @brief Return the box of uncertainty parameters. */
    BoxDomainType error_domain() const { return this->_pdom; }

    /** @brief Evaluate the multifunction at an argument to obtain a validated image set. */
    ValidatedImageSet operator()(const Vector<ValidatedNumber>& x) const;

    /** @brief Write a textual representation of the multifunction. */
    friend OutputStream& operator<<(OutputStream& os, const MultifunctionPatch<P, SIG>& fm) {
        return os << "MultifunctionPatch(function=" << fm._f << ", parameter_domain=" << fm._pdom << ")";
    }
};

template<class P>
auto MultifunctionPatch<P, RealVector(RealVector)>::operator()(const Vector<ValidatedNumber>& x) const -> ValidatedImageSet {
    auto fx = ValidatedFunction<SIG>::constant(error_domain().dimension(), x);
    auto fid = ValidatedFunction<SIG>::identity(error_domain().dimension());
    Function<P, SIG> fim = compose(_f, join(fx, fid));
    return ValidatedImageSet(this->error_domain(), fim);
}

static_assert(AMultifunction<MultifunctionPatch<ValidatedTag, RealVector(RealVector)>, ValidatedTag, RealVector(RealVector), LocatedSet>);

template<class PR>
class Sweeper;

/** @brief Forward-declare a multifunction model with a numeric precision type. */
template<class P, class SIG, class PR>
class MultifunctionModel;
template<class P, class PR>
using VectorMultivariateMultifunctionModel = MultifunctionModel<P, RealVector(RealVector), PR>;
template<class PR>
using ValidatedVectorMultivariateMultifunctionModel = VectorMultivariateMultifunctionModel<ValidatedTag, PR>;

/** @brief Represent a vector multifunction with a validated function model. */
template<class P, class PR>
class MultifunctionModel<P, RealVector(RealVector), PR> {
    using ARG = RealVector;
    using RES = RealVector;
    using SIG = RES(ARG);
    using FLT = RawFloatType<PR>;
    FunctionModel<P, SIG, PR> _f;
    SizeType _as;

  public:
    /** @brief Construct a multifunction model from its argument domain, function, parameter box and sweeper. */
    MultifunctionModel(BoxDomainType dom, VectorMultivariateFunction<P> f, BoxDomainType params, Sweeper<FLT> swp);

    /** @brief Return the number of free arguments, excluding error parameters. */
    SizeType argument_size() const { return this->_as; }
    /** @brief Return the admissible domain of the free arguments. */
    BoxDomainType domain() const { return project(this->_f.domain(), Range(0u, this->argument_size())); }
    /** @brief Return the box of uncertainty parameters. */
    BoxDomainType error_domain() const { return project(this->_f.domain(), Range(this->argument_size(), this->_f.argument_size())); }

    /** @brief Evaluate the multifunction at an argument to obtain a validated image set. */
    ValidatedImageSet operator()(const Vector<ValidatedNumber>& x) const;

    /** @brief Write a textual representation of the multifunction. */
    friend OutputStream& operator<<(OutputStream& os, const MultifunctionModel<P, SIG, PR>& fm) {
        return os << "MultifunctionModel(model=" << fm._f << ", argument_size=" << fm._as << ")";
    }
};

template<class P, class PR>
auto MultifunctionModel<P, RealVector(RealVector), PR>::operator()(const Vector<ValidatedNumber>& x) const -> ValidatedImageSet {
    auto fx = factory(this->_f).create_constants(error_domain(), x);
    auto fid = factory(this->_f).create_identity(error_domain());
    FunctionModel<P, SIG, PR> fim = compose(_f, join(fx, fid));
    return ValidatedImageSet(this->error_domain(), fim);
}

static_assert(AMultifunction<MultifunctionModel<ValidatedTag, RealVector(RealVector), DoublePrecision>, ValidatedTag, RealVector(RealVector), LocatedSet>);

} // namespace Ariadne

#endif // ARIADNE_SET_MULTIFUNCTION_IMAGE_HPP
