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
