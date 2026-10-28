#ifndef ARIADNE_SET_DETAIL_LINEAR_PROGRAM_STATE_HPP
#define ARIADNE_SET_DETAIL_LINEAR_PROGRAM_STATE_HPP

#include "solving/linear_simplex_algorithm.hpp"
#include "utility/array.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"

namespace Ariadne {

template<class X, class XX>
struct LinearProgram {
    LinearProgram() : LinearProgram(0u,0u,X(dp)) { }
    LinearProgram(SizeType m, SizeType n, X const& z);
    Matrix<X> A;
    Vector<X> b;
    Vector<X> c;
    Vector<X> l;
    Vector<X> u;
    Array<Slackness> vt;
    Array<SizeType> p;
    Matrix<XX> B;
    Vector<XX> x;
    Vector<XX> y;
    Vector<XX> z;
};

template<class X, class XX> LinearProgram<X,XX>::LinearProgram(SizeType m, SizeType n, X const& zero)
    : A(m,n,zero), b(m,zero), c(n,zero), l(n,zero), u(n,zero)
    , vt(n), p(n), B(m,m,zero)
    , x(n,zero), y(m,zero),z(n,zero) { }

} // namespace Ariadne

#endif
