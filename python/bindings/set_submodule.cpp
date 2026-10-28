/***************************************************************************
 *            set_submodule.cpp
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

#include "pybind11.hpp"
#include "function-utilities.hpp"
#include "numeric_submodule.hpp"
#include "interval-utilities.hpp"
#include "geometry/geometry.hpp"
#include "set/curve.hpp"
#include "set/function_set.hpp"
#include "set/affine_set.hpp"
#include "set/list_set.hpp"

namespace Ariadne {



template<class F> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Bounds<F>>& x);

template<> struct PythonTemplateName<Point> { static std::string get() { return "Point"; } };
template<> struct PythonTemplateName<Box> { static std::string get() { return "Box"; } };

template<class X> struct PythonClassName<Point<X>> {
    std::string get() const { return python_template_class_name<X>("Point"); } };
template<class X> struct PythonClassName<Box<X>> {
    std::string get() const { return python_template_class_name<X>("Box"); } };

template<class UB> struct PythonClassName<Box<Interval<UB>>> {
    std::string get() const { return python_class_name<UB>()+"Box"; } };
template<> struct PythonClassName<Box<Interval<FloatDP>>> {
    std::string get() const { return "FloatDPExactBox"; } };
template<> struct PythonClassName<Box<Interval<FloatDPUpperBound>>> {
    std::string get() const { return "FloatDPUpperBox"; } };
template<> struct PythonClassName<Box<Interval<FloatDPLowerBound>>> {
    std::string get() const { return "FloatDPLowerBox"; } };
template<> struct PythonClassName<Box<Interval<FloatDPApproximation>>> {
    std::string get() const { return "FloatDPApproximateBox"; } };


template<class UB> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Box<UB>>& repr) {
    Box<UB> const& bx=repr.reference();
    os << python_class_name<Box<UB>>() << "([";
    for(SizeType i=0; i!=bx.dimension(); ++i) {
        if(i!=0) { os << ','; }
        os << python_literal(bx[i]);
    }
    os << "])";
    return os;
}

template<class P, class T> class OpenSetWrapper;
template<class P, class T> class ClosedSetWrapper;
template<class P, class T> class OvertSetWrapper;
template<class P, class T> class BoundedSetWrapper;
template<class P, class T> class CompactSetWrapper;
template<class P, class T> class RegularSetWrapper;
template<class P, class T> class LocatedSetWrapper;

template<class T> class OpenSetWrapper<EffectiveTag,T>
  : public pybind11::wrapper<OpenSetInterface<EffectiveTag,T>>
{
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    OpenSetInterface<EffectiveTag,T>* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    LowerKleenean covers(const BasicSetType& r) const { return this->get_override("covers")(r); }
    LowerKleenean overlaps(const BasicSetType& r) const { return this->get_override("overlaps")(r); }
    OutputStream& _write(OutputStream& os) const { return this->get_override("_write")(os); }
};

template<class T> class ClosedSetWrapper<EffectiveTag,T>
  : public pybind11::wrapper<ClosedSetInterface<EffectiveTag,T>>
{
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    ClosedSetInterface<EffectiveTag,T>* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    LowerKleenean separated(const BasicSetType& r) const { return this->get_override("separated")(r); }
    OutputStream& _write(OutputStream& os) const { return this->get_override("_write")(os); }
};


template<class T> class OvertSetWrapper<EffectiveTag,T>
  : public pybind11::wrapper<OvertSetInterface<EffectiveTag,T>>
{
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    OvertSetInterface<EffectiveTag,T>* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    LowerKleenean overlaps(const BasicSetType& r) const { return this->get_override("overlaps")(r); }
    OutputStream& _write(OutputStream& os) const { return this->get_override("_write")(os); }
};


template<class T> class BoundedSetWrapper<EffectiveTag,T>
  : public pybind11::wrapper<BoundedSetInterface<EffectiveTag,T>>
{
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::BoundingSetType BoundingSetType;
    BoundedSetInterface<EffectiveTag,T>* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    LowerKleenean inside(const BasicSetType& r) const { return this->get_override("inside")(r); }
    BoundingSetType bounding_box() const { return this->get_override("bounding_box")(); }
    OutputStream& _write(OutputStream& os) const { return this->get_override("_write")(os); }
};

template<class T> class CompactSetWrapper<EffectiveTag,T>
  : public pybind11::wrapper<CompactSetInterface<EffectiveTag,T>>
{
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::BoundingSetType BoundingSetType;
    CompactSetInterface<EffectiveTag,T>* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    LowerKleenean separated(const BasicSetType& r) const { return this->get_override("separated")(r); }
    LowerKleenean inside(const BasicSetType& r) const { return this->get_override("inside")(r); }
    LowerKleenean is_bounded() const { return this->get_override("is_bounded")(); }
    BoundingSetType bounding_box() const { return this->get_override("bounding_box")(); }
    OutputStream& _write(OutputStream& os) const { return this->get_override("_write")(os); }
};

template<class T> class RegularSetWrapper<EffectiveTag,T>
  : public pybind11::wrapper<RegularSetInterface<EffectiveTag,T>>
{
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    RegularSetInterface<EffectiveTag,T>* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    LowerKleenean overlaps(const BasicSetType& r) const { return this->get_override("overlaps")(r); }
    LowerKleenean covers(const BasicSetType& r) const { return this->get_override("covers")(r); }
    LowerKleenean separated(const BasicSetType& r) const { return this->get_override("separated")(r); }
    OutputStream& _write(OutputStream& os) const { return this->get_override("_write")(os); }
};

template<class T> class LocatedSetWrapper<EffectiveTag,T>
  : public pybind11::wrapper<LocatedSetInterface<EffectiveTag,T>>
{
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::BoundingSetType BoundingSetType;
    LocatedSetInterface<EffectiveTag,T>* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    LowerKleenean overlaps(const BasicSetType& r) const { return this->get_override("overlaps")(r); }
    LowerKleenean separated(const BasicSetType& r) const { return this->get_override("separated")(r); }
    LowerKleenean inside(const BasicSetType& r) const { return this->get_override("inside")(r); }
    LowerKleenean is_bounded() const { return this->get_override("is_bounded")(); }
    BoundingSetType bounding_box() const { return this->get_override("bounding_box")(); }
    OutputStream& _write(OutputStream& os) const { return this->get_override("_write")(os); }
};




template<class T> class OpenSetWrapper<ValidatedTag,T>
  : public pybind11::wrapper<OpenSetInterface<ValidatedTag,T>>
{
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    OpenSetInterface<ValidatedTag,T>* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    ValidatedLowerKleenean covers(const BasicSetType& r) const { return this->get_override("covers")(r); }
    ValidatedLowerKleenean overlaps(const BasicSetType& r) const { return this->get_override("overlaps")(r); }
    OutputStream& _write(OutputStream& os) const { return this->get_override("_write")(os); }
};

template<class T> class ClosedSetWrapper<ValidatedTag,T>
  : public pybind11::wrapper<ClosedSetInterface<ValidatedTag,T>>
{
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    ClosedSetInterface<ValidatedTag,T>* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    ValidatedLowerKleenean separated(const BasicSetType& r) const { return this->get_override("separated")(r); }
    OutputStream& _write(OutputStream& os) const { return this->get_override("_write")(os); }
};


template<class T> class OvertSetWrapper<ValidatedTag,T>
  : public pybind11::wrapper<OvertSetInterface<ValidatedTag,T>>
{
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    OvertSetInterface<ValidatedTag,T>* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    ValidatedLowerKleenean overlaps(const BasicSetType& r) const { return this->get_override("overlaps")(r); }
    OutputStream& _write(OutputStream& os) const { return this->get_override("_write")(os); }
};


template<class T> class BoundedSetWrapper<ValidatedTag,T>
  : public pybind11::wrapper<BoundedSetInterface<ValidatedTag,T>>
{
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::BoundingSetType BoundingSetType;
    BoundedSetInterface<ValidatedTag,T>* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    ValidatedLowerKleenean inside(const BasicSetType& r) const { return this->get_override("inside")(r); }
    BoundingSetType bounding_box() const { return this->get_override("bounding_box")(); }
    OutputStream& _write(OutputStream& os) const { return this->get_override("_write")(os); }
};


template<class T> class CompactSetWrapper<ValidatedTag,T>
  : public pybind11::wrapper<CompactSetInterface<ValidatedTag,T>>
{
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::BoundingSetType BoundingSetType;
    CompactSetInterface<ValidatedTag,T>* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    ValidatedLowerKleenean separated(const BasicSetType& r) const { return this->get_override("separated")(r); }
    ValidatedLowerKleenean inside(const BasicSetType& r) const { return this->get_override("inside")(r); }
    ValidatedLowerKleenean is_bounded() const { return this->get_override("is_bounded")(); }
    BoundingSetType bounding_box() const { return this->get_override("bounding_box")(); }
    OutputStream& _write(OutputStream& os) const { return this->get_override("_write")(os); }
};

template<class T> class RegularSetWrapper<ValidatedTag,T>
  : public pybind11::wrapper<RegularSetInterface<ValidatedTag,T>>
{
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    RegularSetInterface<ValidatedTag,T>* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    ValidatedLowerKleenean overlaps(const BasicSetType& r) const { return this->get_override("overlaps")(r); }
    ValidatedLowerKleenean covers(const BasicSetType& r) const { return this->get_override("covers")(r); }
    ValidatedLowerKleenean separated(const BasicSetType& r) const { return this->get_override("separated")(r); }
    OutputStream& _write(OutputStream& os) const { return this->get_override("_write")(os); }
};

template<class T> class LocatedSetWrapper<ValidatedTag,T>
  : public pybind11::wrapper<LocatedSetInterface<ValidatedTag,T>>
{
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::BoundingSetType BoundingSetType;
    LocatedSetInterface<ValidatedTag,T>* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    ValidatedLowerKleenean overlaps(const BasicSetType& r) const { return this->get_override("overlaps")(r); }
    ValidatedLowerKleenean separated(const BasicSetType& r) const { return this->get_override("separated")(r); }
    ValidatedLowerKleenean inside(const BasicSetType& r) const { return this->get_override("inside")(r); }
    ValidatedLowerKleenean is_bounded() const { return this->get_override("is_bounded")(); }
    BoundingSetType bounding_box() const { return this->get_override("bounding_box")(); }
    OutputStream& _write(OutputStream& os) const { return this->get_override("_write")(os); }
};


} // namespace Ariadne

using namespace Ariadne;

Void export_set_interface(pybind11::module& module) {
    using P=EffectiveTag; using T=RealVector;
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetInterfaceBase<T>::BoundingSetType BoundingSetType;

    pybind11::class_<OvertSetInterface<P,T>, OvertSetWrapper<P,T>> overt_set_interface_class(module,"OvertSet");
    overt_set_interface_class.def("overlaps",(LowerKleenean(OvertSetInterface<P,T>::*)(const BasicSetType& bx)const) &OvertSetInterface<P,T>::overlaps);

    pybind11::class_<OpenSetInterface<P,T>, OpenSetWrapper<P,T>> open_set_interface_class(module,"OpenSet",overt_set_interface_class);
    open_set_interface_class.def("covers",(LowerKleenean(OpenSetInterface<P,T>::*)(const BasicSetType& bx)const) &OpenSetInterface<P,T>::covers);

    pybind11::class_<ClosedSetInterface<P,T>, ClosedSetWrapper<P,T>> closed_set_interface_class(module,"ClosedSet");
    closed_set_interface_class.def("separated",(LowerKleenean(ClosedSetInterface<P,T>::*)(const BasicSetType& bx)const) &ClosedSetInterface<P,T>::separated);

    pybind11::class_<BoundedSetInterface<P,T>, BoundedSetWrapper<P,T>> bounded_set_interface_class(module,"BoundedSet");
    bounded_set_interface_class.def("inside",(LowerKleenean(BoundedSetInterface<P,T>::*)(const BasicSetType& bx)const) &BoundedSetInterface<P,T>::inside);
    bounded_set_interface_class.def("bounding_box", (BoundingSetType(BoundedSetInterface<P,T>::*)()const)&BoundedSetInterface<P,T>::bounding_box);

    pybind11::class_<CompactSetInterface<P,T>, CompactSetWrapper<P,T>, ClosedSetInterface<P,T>, BoundedSetInterface<P,T>> compact_set_interface_class(module,"CompactSet", pybind11::multiple_inheritance());

    pybind11::class_<RegularSetInterface<P,T>, RegularSetWrapper<P,T>, OpenSetInterface<P,T>,ClosedSetInterface<P,T>> regular_set_interface_class(module,"RegularSet", pybind11::multiple_inheritance());
    pybind11::class_<LocatedSetInterface<P,T>, LocatedSetWrapper<P,T>, OvertSetInterface<P,T>,CompactSetInterface<P,T>> located_set_interface_class(module,"LocatedSet", pybind11::multiple_inheritance());


    pybind11::class_<OvertSetInterface<ValidatedTag,T>, OvertSetWrapper<ValidatedTag,T>> validated_overt_set_interface_class(module,"ValidatedOvertSet");
    validated_overt_set_interface_class.def("overlaps",(ValidatedLowerKleenean(OvertSetInterface<ValidatedTag,T>::*)(const BasicSetType& bx)const) &OvertSetInterface<ValidatedTag,T>::overlaps);

    pybind11::class_<OpenSetInterface<ValidatedTag,T>, OpenSetWrapper<ValidatedTag,T>> validated_open_set_interface_class(module,"ValidatedOpenSet",validated_overt_set_interface_class);
    validated_open_set_interface_class.def("covers",(ValidatedLowerKleenean(OpenSetInterface<ValidatedTag,T>::*)(const BasicSetType& bx)const) &OpenSetInterface<ValidatedTag,T>::covers);

    pybind11::class_<ClosedSetInterface<ValidatedTag,T>, ClosedSetWrapper<ValidatedTag,T>> validated_closed_set_interface_class(module,"ValidatedClosedSet");
    validated_closed_set_interface_class.def("separated",(ValidatedLowerKleenean(ClosedSetInterface<ValidatedTag,T>::*)(const BasicSetType& bx)const) &ClosedSetInterface<ValidatedTag,T>::separated);

    pybind11::class_<BoundedSetInterface<ValidatedTag,T>, BoundedSetWrapper<ValidatedTag,T>> validated_bounded_set_interface_class(module,"ValidatedBoundedSet");
    validated_bounded_set_interface_class.def("inside",(ValidatedLowerKleenean(BoundedSetInterface<ValidatedTag,T>::*)(const BasicSetType& bx)const) &BoundedSetInterface<ValidatedTag,T>::inside);
    validated_bounded_set_interface_class.def("bounding_box", (BoundingSetType(BoundedSetInterface<ValidatedTag,T>::*)()const)&BoundedSetInterface<ValidatedTag,T>::bounding_box);

    pybind11::class_<CompactSetInterface<ValidatedTag,T>, CompactSetWrapper<ValidatedTag,T>, ClosedSetInterface<ValidatedTag,T>, BoundedSetInterface<ValidatedTag,T>> validated_compact_set_interface_class(module,"ValidatedCompactSet", pybind11::multiple_inheritance());

    pybind11::class_<RegularSetInterface<ValidatedTag,T>, RegularSetWrapper<ValidatedTag,T>, OpenSetInterface<ValidatedTag,T>, ClosedSetInterface<ValidatedTag,T>> validated_regular_set_interface_class(module,"ValidatedRegularSet", pybind11::multiple_inheritance());
    pybind11::class_<LocatedSetInterface<ValidatedTag,T>, LocatedSetWrapper<ValidatedTag,T>, OvertSetInterface<ValidatedTag,T>,CompactSetInterface<ValidatedTag,T>> validated_located_set_interface_class(module,"ValidatedLocatedSet", pybind11::multiple_inheritance());
}
/*

Void export_zonotope(pybind11::module& module)
{
    pybind11::class_<Zonotope,pybind11::bases<CompactSetInterface<T>,OpenSetInterface<T>>> zonotope_class(module,"Zonotope");
    zonotope_class.def(pybind11::init<Zonotope>());
    zonotope_class.def(pybind11::init<Vector<FloatDP>,Matrix<FloatDP>,Vector<FloatDPError>>());
    zonotope_class.def(pybind11::init<Vector<FloatDP>,Matrix<FloatDP>>());
    zonotope_class.def(pybind11::init<BasicSetType>());
    zonotope_class.def("centre",&Zonotope::centre);
    zonotope_class.def("generators",&Zonotope::generators);
    zonotope_class.def("error",&Zonotope::error);
    zonotope_class.def("contains",&Zonotope::contains);
    zonotope_class.def("split", (ListSet<Zonotope>(*)(const Zonotope&)) &split);
    zonotope_class.def("__str__",&__cstr__<Zonotope>);

    module.def("contains", (ValidatedKleenean(*)(const Zonotope&,const ExactPoint&)) &contains);
    module.def("separated", (ValidatedKleenean(*)(const Zonotope&,const BasicSetType&)) &separated);
    module.def("overlaps", (ValidatedKleenean(*)(const Zonotope&,const BasicSetType&)) &overlaps);
    module.def("separated", (ValidatedKleenean(*)(const Zonotope&,const Zonotope&)) &separated);

    module.def("polytope", (Polytope(*)(const Zonotope&)) &polytope);
    module.def("orthogonal_approximation", (Zonotope(*)(const Zonotope&)) &orthogonal_approximation);
    module.def("orthogonal_over_approximation", (Zonotope(*)(const Zonotope&)) &orthogonal_over_approximation);
    module.def("error_free_over_approximation", (Zonotope(*)(const Zonotope&)) &error_free_over_approximation);

//    module.def("image", (Zonotope(*)(const Zonotope&, const ValidatedVectorMultivariateFunction&)) &image);
}

Void export_polytope(pybind11::module& module)
{
    pybind11::class_<Polytope,pybind11::bases<LocatedSetInterface<T>>> polytope_class(module,"Polytope");
    polytope_class.def(pybind11::init<Polytope>());
    polytope_class.def(pybind11::init<Int>());
    polytope_class.def("new_vertex",&Polytope::new_vertex);
    polytope_class.def("__iter__",boost::python::range(&Polytope::vertices_begin,&Polytope::vertices_end));
    polytope_class.def(self_ns::str(self));
}

*/

Void export_curve(pybind11::module& module)
{
    pybind11::class_<Curve> smooth_curve_class(module,"Curve");
    smooth_curve_class.def(pybind11::init<EffectiveVectorUnivariateFunction>());
    smooth_curve_class.def("dimension", &Curve::dimension);
    smooth_curve_class.def("value", &Curve::value);
    smooth_curve_class.def("tangent", &Curve::tangent);
    smooth_curve_class.def("__str__", &__cstr__<Curve>);
    pybind11::class_<InterpolatedCurve> interpolated_curve_class(module,"InterpolatedCurve");
    interpolated_curve_class.def(pybind11::init<InterpolatedCurve>());
    interpolated_curve_class.def(pybind11::init<FloatDP,FloatDPPoint>());
    interpolated_curve_class.def("insert", (Void(InterpolatedCurve::*)(const FloatDP&, const Point<FloatDPApproximation>&)) &InterpolatedCurve::insert);
    interpolated_curve_class.def("__iter__", [](InterpolatedCurve const& c){return pybind11::make_iterator(c.begin(),c.end());});
    interpolated_curve_class.def("__str__", &__cstr__<InterpolatedCurve>);
    interpolated_curve_class.def("size", &InterpolatedCurve::size);
    interpolated_curve_class.def("dimension", &InterpolatedCurve::dimension);
    interpolated_curve_class.def("bounding_box", &InterpolatedCurve::bounding_box);


}



Void export_affine_set(pybind11::module& module)
{
    pybind11::class_<ValidatedAffineConstrainedImageSet,pybind11::bases<ValidatedEuclideanCompactSetInterface>>
        affine_set_class(module,"ValidatedAffineConstrainedImageSet", pybind11::multiple_inheritance());
    affine_set_class.def(pybind11::init<ValidatedAffineConstrainedImageSet>());
    affine_set_class.def(pybind11::init<RealBox>());
    affine_set_class.def(pybind11::init<ExactBoxType>());
    affine_set_class.def(pybind11::init<Vector<ExactIntervalType>, Matrix<FloatDP>, Vector<FloatDP> >());
    affine_set_class.def("new_parameter_constraint", (Void(ValidatedAffineConstrainedImageSet::*)(const Constraint<Affine<FloatDPBounds>,FloatDPBounds>&)) &ValidatedAffineConstrainedImageSet::new_parameter_constraint);
    affine_set_class.def("new_constraint", (Void(ValidatedAffineConstrainedImageSet::*)(const Constraint<AffineModel<ValidatedTag,FloatDP>,FloatDPBounds>&)) &ValidatedAffineConstrainedImageSet::new_constraint);
    affine_set_class.def("dimension", &ValidatedAffineConstrainedImageSet::dimension);
    affine_set_class.def("is_bounded", &ValidatedAffineConstrainedImageSet::is_bounded);
    affine_set_class.def("is_empty", &ValidatedAffineConstrainedImageSet::is_empty);
    affine_set_class.def("bounding_box", &ValidatedAffineConstrainedImageSet::bounding_box);
    affine_set_class.def("separated", &ValidatedAffineConstrainedImageSet::separated);
    affine_set_class.def("boundary", &ValidatedAffineConstrainedImageSet::boundary);
    affine_set_class.def("__str__",&__cstr__<ValidatedAffineConstrainedImageSet>);

    module.def("image", (ValidatedAffineConstrainedImageSet(*)(ValidatedAffineConstrainedImageSet const&,ValidatedVectorMultivariateFunction const&)) &_image_);
}

Void export_constraint_set(pybind11::module& module)
{
//    from_python< List<EffectiveConstraint> >();

    pybind11::class_<ConstraintSet,pybind11::bases<EffectiveEuclideanRegularSetInterface,EffectiveEuclideanOpenSetInterface> >
        constraint_set_class(module,"ConstraintSet", pybind11::multiple_inheritance());
    constraint_set_class.def(pybind11::init<ConstraintSet>());
    constraint_set_class.def(pybind11::init< List<EffectiveConstraint> >());
    constraint_set_class.def("dimension", &ConstraintSet::dimension);
    constraint_set_class.def("number_of_constraints", &ConstraintSet::number_of_constraints);
    constraint_set_class.def("__str__", &__cstr__<ConstraintSet>);

//    pybind11::class_<BoundedConstraintSet,pybind11::bases<DrawableWrapper> >
    pybind11::class_<BoundedConstraintSet,pybind11::bases<EffectiveEuclideanRegularSetInterface,EffectiveEuclideanLocatedSetInterface> >
        bounded_constraint_set_class(module,"BoundedConstraintSet", pybind11::multiple_inheritance());
    bounded_constraint_set_class.def(pybind11::init<BoundedConstraintSet>());
    bounded_constraint_set_class.def(pybind11::init< RealBox, List<EffectiveConstraint> >());
    bounded_constraint_set_class.def("dimension", &BoundedConstraintSet::dimension);
    bounded_constraint_set_class.def("domain", &BoundedConstraintSet::domain);
    bounded_constraint_set_class.def("bounding_box", &BoundedConstraintSet::bounding_box);
    bounded_constraint_set_class.def("number_of_constraints", &BoundedConstraintSet::number_of_constraints);
    bounded_constraint_set_class.def("__str__", &__cstr__<BoundedConstraintSet>);

    module.def("intersection", (ConstraintSet(*)(ConstraintSet const&,ConstraintSet const&)) &_intersection_);
    module.def("intersection", (BoundedConstraintSet(*)(ConstraintSet const&, RealBox const&)) &_intersection_);
    module.def("intersection", (BoundedConstraintSet(*)(RealBox const&, ConstraintSet const&)) &_intersection_);

    module.def("intersection", (BoundedConstraintSet(*)(BoundedConstraintSet const&, BoundedConstraintSet const&)) &_intersection_);
    module.def("intersection", (BoundedConstraintSet(*)(BoundedConstraintSet const&, RealBox const&)) &_intersection_);
    module.def("intersection", (BoundedConstraintSet(*)(RealBox const&, BoundedConstraintSet const&)) &_intersection_);
    module.def("intersection", (BoundedConstraintSet(*)(ConstraintSet const&, BoundedConstraintSet const&)) &_intersection_);
    module.def("intersection", (BoundedConstraintSet(*)(BoundedConstraintSet const&, ConstraintSet const&)) &_intersection_);

    module.def("image", (ConstrainedImageSet(*)(BoundedConstraintSet const&, EffectiveVectorMultivariateFunction const&)) &_image_);

}


Void export_constrained_image_set(pybind11::module& module)
{
//    from_python< List<ValidatedConstraint> >();

    pybind11::class_<ConstrainedImageSet,pybind11::bases<EffectiveEuclideanLocatedSetInterface> >
        constrained_image_set_class(module,"ConstrainedImageSet");
    constrained_image_set_class.def(pybind11::init<ConstrainedImageSet>());
    constrained_image_set_class.def(pybind11::init<BoundedConstraintSet>());
    constrained_image_set_class.def(pybind11::init<RealBox,EffectiveVectorMultivariateFunction>());
    constrained_image_set_class.def(pybind11::init<RealBox,EffectiveVectorMultivariateFunction,List<EffectiveConstraint> >());
    constrained_image_set_class.def("dimension", &ConstrainedImageSet::dimension);
    constrained_image_set_class.def("domain", &ConstrainedImageSet::domain);
    constrained_image_set_class.def("bounding_box", &ConstrainedImageSet::bounding_box);
    constrained_image_set_class.def("number_of_parameters", &ConstrainedImageSet::number_of_parameters);
    constrained_image_set_class.def("split", (Pair<ConstrainedImageSet,ConstrainedImageSet>(ConstrainedImageSet::*)(SizeType)const) &ConstrainedImageSet::split);
    constrained_image_set_class.def("split", (Pair<ConstrainedImageSet,ConstrainedImageSet>(ConstrainedImageSet::*)()const) &ConstrainedImageSet::split);
//    	constrained_image_set_class.def("affine_over_approximation", &ValidatedConstrainedImageSet::affine_over_approximation);
    constrained_image_set_class.def("__str__",&__cstr__<ConstrainedImageSet>);

//    pybind11::class_<ValidatedConstrainedImageSet,pybind11::bases<CompactSetInterface<T>> >
    pybind11::class_<ValidatedConstrainedImageSet,pybind11::bases<ValidatedEuclideanLocatedSetInterface> >
        validated_constrained_image_set_class(module,"ValidatedConstrainedImageSet", pybind11::multiple_inheritance());
    validated_constrained_image_set_class.def(pybind11::init<ValidatedConstrainedImageSet>());
    validated_constrained_image_set_class.def(pybind11::init<ExactBoxType>());
    validated_constrained_image_set_class.def(pybind11::init<ExactBoxType,EffectiveVectorMultivariateFunction>());
    validated_constrained_image_set_class.def(pybind11::init<ExactBoxType,ValidatedVectorMultivariateFunction>());
    validated_constrained_image_set_class.def(pybind11::init<ExactBoxType,ValidatedVectorMultivariateFunction,List<ValidatedConstraint> >());
    validated_constrained_image_set_class.def(pybind11::init<ExactBoxType,ValidatedVectorMultivariateFunctionModelDP>());
    validated_constrained_image_set_class.def("dimension", &ValidatedConstrainedImageSet::dimension);
    validated_constrained_image_set_class.def("domain", &ValidatedConstrainedImageSet::domain);
    validated_constrained_image_set_class.def("reduce", &ValidatedConstrainedImageSet::reduce);
    validated_constrained_image_set_class.def("reduced_domain", &ValidatedConstrainedImageSet::reduced_domain);
    validated_constrained_image_set_class.def("restriction", &ValidatedConstrainedImageSet::restriction);
    validated_constrained_image_set_class.def("function", &ValidatedConstrainedImageSet::function);
    validated_constrained_image_set_class.def("constraint", &ValidatedConstrainedImageSet::constraint);
    validated_constrained_image_set_class.def("number_of_parameters", &ValidatedConstrainedImageSet::number_of_parameters);
    validated_constrained_image_set_class.def("number_of_constraints", &ValidatedConstrainedImageSet::number_of_constraints);
    validated_constrained_image_set_class.def("apply", &ValidatedConstrainedImageSet::apply);
    validated_constrained_image_set_class.def("new_space_constraint", (Void(ValidatedConstrainedImageSet::*)(const ValidatedConstraint&))&ValidatedConstrainedImageSet::new_space_constraint);
    validated_constrained_image_set_class.def("new_parameter_constraint", (Void(ValidatedConstrainedImageSet::*)(const ValidatedConstraint&))&ValidatedConstrainedImageSet::new_parameter_constraint);
    validated_constrained_image_set_class.def("affine_approximation", &ValidatedConstrainedImageSet::affine_approximation);
    validated_constrained_image_set_class.def("affine_over_approximation", &ValidatedConstrainedImageSet::affine_over_approximation);
    validated_constrained_image_set_class.def("bounding_box", &ValidatedConstrainedImageSet::bounding_box);
    validated_constrained_image_set_class.def("inside", &ValidatedConstrainedImageSet::inside);
    validated_constrained_image_set_class.def("separated", &ValidatedConstrainedImageSet::separated);
    validated_constrained_image_set_class.def("overlaps", &ValidatedConstrainedImageSet::overlaps);
    validated_constrained_image_set_class.def("split", (Pair<ValidatedConstrainedImageSet,ValidatedConstrainedImageSet>(ValidatedConstrainedImageSet::*)()const) &ValidatedConstrainedImageSet::split);
    validated_constrained_image_set_class.def("split", (Pair<ValidatedConstrainedImageSet,ValidatedConstrainedImageSet>(ValidatedConstrainedImageSet::*)(SizeType)const) &ValidatedConstrainedImageSet::split);
    validated_constrained_image_set_class.def("__str__", &__cstr__<ValidatedConstrainedImageSet>);
    validated_constrained_image_set_class.def("__repr__", &__cstr__<ValidatedConstrainedImageSet>);

    //module.def("product", (ValidatedConstrainedImageSet(*)(const ValidatedConstrainedImageSet&,const BasicSetType&)) &product);
}
void export_list_set(pybind11::module& module) {
    using SetOfBoxes = ListSet<ExactBoxType>;
    pybind11::class_<SetOfBoxes> box_list_class(module, "ListSetFloatDPExactBox");
    box_list_class.def(pybind11::init<>());
    box_list_class.def("adjoin", (Void(SetOfBoxes::*)(const ExactBoxType&)) &SetOfBoxes::adjoin);
    box_list_class.def("size", &SetOfBoxes::size);
    box_list_class.def("dimension", &SetOfBoxes::dimension);
    box_list_class.def("bounding_box", &SetOfBoxes::bounding_box);
    box_list_class.def("__len__", &SetOfBoxes::size);
    box_list_class.def("__str__", &__cstr__<SetOfBoxes>);
}
void set_submodule(pybind11::module& module) {
    export_list_set(module);
    export_set_interface(module);
    export_curve(module);
    export_affine_set(module);
    export_constraint_set(module);
    export_constrained_image_set(module);
}
