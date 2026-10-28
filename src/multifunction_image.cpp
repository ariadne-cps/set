/***************************************************************************
 *            multifunction_image.cpp
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

// Implementation and explicit instantiations for set-valued multifunctions.
#include "set/multifunction_image.hpp"
#include "function/taylor_model.hpp"
#include "algebra/algebra.hpp"
#include "function/taylor_function.hpp"

namespace Ariadne {

template<class P, class PR> MultifunctionModel<P,RealVector(RealVector),PR>::
    MultifunctionModel(BoxDomainType dom, VectorMultivariateFunction<P> f, BoxDomainType params, Sweeper<FLT> swp)
        : _f(VectorMultivariateTaylorFunctionModel<P,FLT>(product(dom,params),f,swp)), _as(dom.dimension()) { }

template class MultifunctionPatch<ValidatedTag,RealVector(RealVector)>;
template class MultifunctionModel<ValidatedTag,RealVector(RealVector),DoublePrecision>;

} // namespace Ariadne
