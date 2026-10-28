/***************************************************************************
 *            multifunction_image.hpp
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

// Multifunctions with constrained-image results belong to set, not function.
#ifndef ARIADNE_SET_MULTIFUNCTION_IMAGE_HPP
#define ARIADNE_SET_MULTIFUNCTION_IMAGE_HPP

#include "function/multifunction.hpp"
#include "function/function_model.hpp"
#include "set/function_set.hpp"

namespace Ariadne {

using ValidatedImageSet = ValidatedConstrainedImageSet;

template<class P, class SIG, template<class,class>class SET=LocatedSet> class MultifunctionPatch;
template<class P> using VectorMultivariateMultifunctionPatch = MultifunctionPatch<P,RealVector(RealVector)>;
using ValidatedVectorMultivariateMultifunctionPatch = VectorMultivariateMultifunctionPatch<ValidatedTag>;

template<class P> class MultifunctionPatch<P,RealVector(RealVector)> {
    using ARG=RealVector; using RES=RealVector; using SIG=RES(ARG);
    Function<P,SIG> _f;
    BoxDomainType _pdom;
  public:
    MultifunctionPatch(VectorMultivariateFunction<P> f, BoxDomainType pdom) : _f(f), _pdom(pdom) { }

    SizeType argument_size() const { return this->_f.argument_size()-this->_pdom.dimension(); }
    BoxDomainType error_domain() const { return this->_pdom; }

    ValidatedImageSet operator() (Vector<ValidatedNumber> const& x) const;

    friend OutputStream& operator<<(OutputStream& os, MultifunctionPatch<P,SIG> const& fm) {
        return os << "MultifunctionPatch(function="<<fm._f<<", parameter_domain="<<fm._pdom<<")"; }
};

template<class P> auto
MultifunctionPatch<P,RealVector(RealVector)>::operator() (Vector<ValidatedNumber> const& x) const -> ValidatedImageSet
{
    auto fx=ValidatedFunction<SIG>::constant(error_domain().dimension(),x);
    auto fid=ValidatedFunction<SIG>::identity(error_domain().dimension());
    Function<P,SIG> fim=compose(_f,join(fx,fid));
    return ValidatedImageSet(this->error_domain(),fim);
}

static_assert(AMultifunction<MultifunctionPatch<ValidatedTag,RealVector(RealVector)>,ValidatedTag,RealVector(RealVector),LocatedSet>);



template<class PR> class Sweeper;

template<class P, class SIG, class PR> class MultifunctionModel;
template<class P, class PR> using VectorMultivariateMultifunctionModel = MultifunctionModel<P,RealVector(RealVector),PR>;
template<class PR> using ValidatedVectorMultivariateMultifunctionModel = VectorMultivariateMultifunctionModel<ValidatedTag,PR>;

template<class P, class PR> class MultifunctionModel<P,RealVector(RealVector),PR> {
    using ARG=RealVector; using RES=RealVector; using SIG=RES(ARG);
    using FLT=RawFloatType<PR>;
    FunctionModel<P,SIG,PR> _f;
    SizeType _as;
  public:
    MultifunctionModel(BoxDomainType dom, VectorMultivariateFunction<P> f, BoxDomainType params, Sweeper<FLT> swp);

    SizeType argument_size() const { return this->_as; }
    BoxDomainType domain() const { return project(this->_f.domain(),Range(0u,this->argument_size())); }
    BoxDomainType error_domain() const { return project(this->_f.domain(),Range(this->argument_size(),this->_f.argument_size())); }

    ValidatedImageSet operator() (Vector<ValidatedNumber> const& x) const;

    friend OutputStream& operator<<(OutputStream& os, MultifunctionModel<P,SIG,PR> const& fm) {
        return os << "MultifunctionModel(model="<<fm._f<<", argument_size="<<fm._as<<")"; }
};


template<class P, class PR> auto
MultifunctionModel<P,RealVector(RealVector),PR>::operator() (Vector<ValidatedNumber> const& x) const -> ValidatedImageSet
{
    auto fx=factory(this->_f).create_constants(error_domain(),x);
    auto fid=factory(this->_f).create_identity(error_domain());
    FunctionModel<P,SIG,PR> fim=compose(_f,join(fx,fid));
    return ValidatedImageSet(this->error_domain(),fim);
}

static_assert(AMultifunction<MultifunctionModel<ValidatedTag,RealVector(RealVector),DoublePrecision>,ValidatedTag,RealVector(RealVector),LocatedSet>);


} // namespace Ariadne

#endif // ARIADNE_SET_MULTIFUNCTION_IMAGE_HPP
