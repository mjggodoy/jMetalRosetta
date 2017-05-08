// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   rosetta_binders.hpp
/// @brief  Support for custom binders for some Rosetta template classes
/// @author Sergey Lyskov


#ifndef _INCLUDED_rosetta_binders_hpp_
#define _INCLUDED_rosetta_binders_hpp_

#include <pybind11/operators.h>
#include <pybind11/stl_bind.h>

#include <utility/vectorL.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>

#include <type_traits>
#include <sstream>

#include <set>
#include <map>

namespace rosetta_binders {


// template <typename T, class Allocator>
// void vector_binder(pybind11::module &m, char const *name, char const * /*allocator name*/) {
// 	pybind11::vector_binder<T, Allocator, std::shared_ptr< std::vector<T, Allocator> > >(m, name);
// }


template<typename T>
//constexpr auto has_equal_operator(int) -> decltype(std::declval<T>() == std::declval<T>(), std::declval<T>() != std::declval<T>(), bool()) { return true; }
constexpr auto has_equal_operator(int) -> decltype(std::is_same<std::true_type, decltype(std::declval<T>() == std::declval<T>())>::value  and
												   std::is_same<std::true_type, decltype(std::declval<T>() != std::declval<T>())>::value) { return true; }
template<typename T>
constexpr bool has_equal_operator(...) { return false; }

template<typename T, typename SFINAE = void>
struct has_equal_operator_s {
	static const bool value = false;
};
template<typename T>
struct has_equal_operator_s<T>
{
	static const bool value = has_equal_operator<T>(0);
};
template <typename A>
struct has_equal_operator_s< std::vector<A> >
{
	static const bool value = has_equal_operator_s<A>::value;
};
// template <typename A>
// struct has_equal_operator_s< std::deque<A> >
// {
// 	static const bool value = has_equal_operator_s<A>::value;
// };
template <typename A>
struct has_equal_operator_s< std::set<A> >
{
	static const bool value = has_equal_operator_s<A>::value;
};
template <typename A, typename B>
struct has_equal_operator_s< std::pair<A,B> >
{
	static const bool value = has_equal_operator_s<A>::value and has_equal_operator_s<B>::value;
};
template <typename A, typename B>
struct has_equal_operator_s< std::map<A,B> >
{
	static const bool value = has_equal_operator_s<A>::value and has_equal_operator_s<B>::value;
};


// Rosetta specific types
template <typename A>
struct has_equal_operator_s< utility::vector0<A> >
{
	static const bool value = has_equal_operator_s<A>::value;
};
template <typename A>
struct has_equal_operator_s< utility::vector1<A> >
{
	static const bool value = has_equal_operator_s<A>::value;
};



namespace has_insertion_operator_implementation {
enum class False {};
struct any_type {
    template<typename T> any_type(T const&);
};
False operator<<(std::ostream const&, any_type const&);
}
template<typename T>
constexpr bool has_insertion_operator() {
	using namespace has_insertion_operator_implementation;
	return std::is_same< decltype(std::declval<std::ostream&>() << std::declval<T>()), std::ostream & >::value;
}
template<typename T>
struct has_insertion_operator_s {
	static const bool value = has_insertion_operator<T>();
};


// template <typename T, typename Iterator, typename... Extra> pybind11::iterator make_iterator(Iterator first, Iterator last, Extra&&... extra) {
//     typedef pybind11::detail::iterator_state<Iterator> state;
//     if (!pybind11::detail::get_type_info(typeid(state))) {
//         pybind11::class_<state>(pybind11::handle(), "")
//             .def("__iter__", [](state &s) -> state& { return s; })
//             .def("__next__", [](state &s) -> T {
//                 if (s.it == s.end)
//                     throw pybind11::stop_iteration();
//                 return *s.it++;
//             }, pybind11::return_value_policy::reference_internal, std::forward<Extra>(extra)...);
//     }
//     return (pybind11::iterator) pybind11::cast(state { first, last });
// }
// template <typename T, typename Type, typename... Extra> pybind11::iterator make_iterator(Type &value, Extra&&... extra) {
//     return make_iterator<T>(std::begin(value), std::end(value), extra...);
// }



template<typename Vector, platform::SSize L, typename T, typename Allocator>
class utility_vector_binder
{
	using SizeType = typename Vector::size_type;
	using SSizeType = typename Vector::ssize_type;
    using ItType   = typename Vector::iterator;
	using Class_ = pybind11::class_<Vector, std::shared_ptr< Vector > >;


	// template<typename U = T, typename std::enable_if< std::is_constructible<U>{} >::type * = nullptr>
	// void maybe_constructible(Class_ &cl) {
	// 	cl.def(pybind11::init<>());
	// }
	// template<typename U = T, typename std::enable_if< !std::is_constructible<U>{} >::type * = nullptr>
	// void maybe_constructible(Class_ &cl) {}

	template<typename U = T, typename std::enable_if< std::is_default_constructible<U>::value >::type * = nullptr>
	void maybe_default_constructible(Class_ &cl) {
		cl.def(pybind11::init<SizeType>());
		cl.def("resize", (void (Vector::*)(SizeType count)) &Vector::resize, "changes the number of elements stored");
	}
	template<typename U = T, typename std::enable_if< !std::is_default_constructible<U>::value >::type * = nullptr>
	void maybe_default_constructible(Class_ &) {}


	template<typename U = T, typename std::enable_if< std::is_copy_constructible<U>::value >::type * = nullptr>
	void maybe_copy_constructible(Class_ &cl) {
		cl.def(pybind11::init< Vector const &>());
	}
	template<typename U = T, typename std::enable_if< !std::is_copy_constructible<U>::value >::type * = nullptr>
	void maybe_copy_constructible(Class_ &) {}


	template<typename U = T, typename std::enable_if< has_equal_operator_s<U>::value >::type * = nullptr>
	void maybe_has_equal_operator(Class_ &cl) {
	    cl.def(pybind11::self == pybind11::self);
	    cl.def(pybind11::self != pybind11::self);

		cl.def("count", [](Vector const &v, T const & value) { return std::count(v.begin(), v.end(), value); }, "counts the elements that are equal to value");

		cl.def("remove", [](Vector &v, T const &t) {
				auto p = std::find(v.begin(), v.end(), t);
				if(p != v.end()) v.erase(p);
				else throw pybind11::value_error();
			}, "Remove the first item from the list whose value is x. It is an error if there is no such item.");

		cl.def("__contains__", [](Vector const &v, T const &t) { return std::find(v.begin(), v.end(), t) != v.end(); }, "return true if item in the container");
	}
	template<typename U = T, typename std::enable_if< !has_equal_operator_s<U>::value >::type * = nullptr>
	void maybe_has_equal_operator(Class_ &) {}


	template<typename U = T, typename std::enable_if< has_insertion_operator_s<U>::value >::type * = nullptr>
	void maybe_has_insertion_operator(std::string const &name, Class_ &cl) {
		cl.def("__repr__", [name](Vector &v) {
				std::ostringstream s;
				s << "vector" << L << '_' << name << '[';
				if( v.size() ) {
					for(SizeType i=v.l(); i<=v.u(); ++i) {
						s << v[i];
						if(i != v.u()) s << ", ";
					}
				}
				s << ']';
 				return s.str();
			});

	}
	template<typename U = T, typename std::enable_if< !has_insertion_operator_s<U>::value >::type * = nullptr>
	void maybe_has_insertion_operator(std::string const &, Class_ &) {}


public:
	utility_vector_binder(pybind11::module &m, std::string const &lower_index, std::string const &name, std::string const & /*allocator name*/) {
		Class_ cl(m, std::string("vector"+ lower_index +'_' + name).c_str() );

		cl.def(pybind11::init<>());

		//maybe_constructible(cl);
		maybe_default_constructible(cl);
		maybe_copy_constructible(cl);

		// Element access
		cl.def("front", [](Vector &v) -> T {
				if(v.size()) return v.front();
				else throw pybind11::index_error();
			}, "access the first element");
		cl.def("back", [](Vector &v) -> T {
				if(v.size()) return v.back();
				else throw pybind11::index_error();
			}, "access the last element ");
		// Not needed, the operator[] is already providing bounds checking cl.def("at", (T& (Vector::*)(SizeType i)) &Vector::at, "access specified element with bounds checking");

		// Capacity, C++ style
		cl.def("max_size",      &Vector::max_size,      "returns the maximum possible number of elements");
		cl.def("reserve",       &Vector::reserve,       "reserves storage");
		cl.def("capacity",      &Vector::capacity,      "returns the number of elements that can be held in currently allocated storage");
		cl.def("shrink_to_fit", &Vector::shrink_to_fit, "reduces memory usage by freeing unused memory");

		// Modifiers, C++ style
		cl.def("clear", &Vector::clear, "clears the contents");

		// Modifiers, Python style
		cl.def("append", [](Vector &v, const T &value) { v.push_back(value); }, "adds an element to the end");
		//cl.def("insert", [](Vector &v, SizeType i, const T&t) {v.insert(v.begin()+i, t);}, "insert an item at a given position");
		cl.def("extend", [](Vector &v, Vector &src) { v.reserve( v.size() + src.size() ); v.insert(v.end(), src.begin(), src.end()); }, "extend the list by appending all the items in the given vector");
		cl.def("pop", [](Vector &v) {
				if(v.size()) {
					T t = v.back();
					v.pop_back();
					return t;
				}
				else throw pybind11::index_error();
			}, "remove and return last item");

		cl.def("pop", [](Vector &v, SizeType i) {
				if( v.empty()  or i > v.u()  or  i < v.l() ) throw pybind11::index_error();
				T t = v[i];
				v.erase(v.begin() + i);
				return t;
			}, "remove and return item at index");

		cl.def("erase", [](Vector &v, SizeType i) {
				if( v.empty()  or i > v.u()  or  i < v.l() ) throw pybind11::index_error();
				v.erase(v.begin() + i);
			}, "erases element at index");


		// Python friendly bindings
		// #ifdef PYTHON_ABI_VERSION // Python 3+
		// 	cl.def("__bool__",    [](Vector &v) -> bool { return v.size() != 0; }); // checks whether the container has any elements in it
		// #else
		// 	cl.def("__nonzero__", [](Vector &v) -> bool { return v.size() != 0; }); // checks whether the container has any elements in it
		// #endif

		cl.def("__bool__", [](const Vector &v) -> bool {
				return !v.empty();
			},
			"Check whether the list is nonempty");


		cl.def("__getitem__", [](Vector const &v, SSizeType i) -> T {
				if( v.empty()  or i > v.u()  or  i < v.l() ) throw pybind11::index_error();
				return v[i];
			});

		cl.def("__setitem__", [](Vector &v, SSizeType i, T const & t) {
				if( v.empty()  or  i > v.u()  or  i < v.l() ) throw pybind11::index_error();
				v[i] = t;
			});

		cl.def("__len__", [](Vector const &v) { return v.size(); } ); // workaround for ld: warning: direct access in ... means the weak symbol cannot be overridden at runtime. This was likely caused by different translation units being compiled with different visibility settings.

		cl.def("__iter__", [](Vector &v) {
				return pybind11::make_iterator<pybind11::return_value_policy::reference_internal, ItType, ItType, T>(v.begin(), v.end());
			},
			pybind11::keep_alive<0, 1>() /* Essential: keep list alive while iterator exists */
			);

		/// Slicing protocol
		// cl.def("__setitem__", [](Vector &v, pybind11::slice slice,  Vector const &value) {
		// 		pybind11::ssize_t start, stop, step, slicelength;
		// 		if(!slice.compute(v.size(), &start, &stop, &step, &slicelength))
		// 			throw pybind11::error_already_set();
		// 		if((size_t) slicelength != value.size())
		// 			throw std::runtime_error("Left and right hand size of slice assignment have different sizes!");
		// 		for(int i=0; i<slicelength; ++i) {
		// 			v[start] = value[i]; start += step;
		// 		}
		// 	});

		cl.def("l", (SizeType (Vector::*)() const) &Vector::l, "lower index");
		cl.def("u", (SizeType (Vector::*)() const) &Vector::u, "upper index");

		// fixme later: add separate binding for to/from std::vector  cl.def("swap",   &Vector::swap, "swaps the contents");


		// Comparisons
		maybe_has_equal_operator(cl);

		// Printing
		maybe_has_insertion_operator(name, cl);
	}

	utility_vector_binder(pybind11::module &m, std::string const &name, std::string const & allocator) : utility_vector_binder(m, std::to_string(L), name, allocator) {}
};


template<platform::SSize L, typename T, typename A> using vectorL_binder = utility_vector_binder<utility::vectorL<L, T, A>, L, T, A>;
template< typename T, typename A > using vector0_binder = utility_vector_binder<utility::vector0<T, A>, 0, T, A>;
template< typename T, typename A > using vector1_binder = utility_vector_binder<utility::vector1<T, A>, 1, T, A>;


template< typename T >
void xyzVector_add_on_binder(pybind11::class_<numeric::xyzVector<T>, std::shared_ptr< numeric::xyzVector<T> > > &cl) {
	using Vector = numeric::xyzVector<T>;
	cl.def_property("x", (T const & (Vector::*)() const) &Vector::x, (void (numeric::xyzVector<T>::*)(const T &)) &Vector::x);
	cl.def_property("y", (T const & (Vector::*)() const) &Vector::y, (void (numeric::xyzVector<T>::*)(const T &)) &Vector::y);
	cl.def_property("z", (T const & (Vector::*)() const) &Vector::z, (void (numeric::xyzVector<T>::*)(const T &)) &Vector::z);

	cl.def(pybind11::self +  pybind11::self);
	cl.def(pybind11::self -  pybind11::self);
	cl.def(pybind11::self += pybind11::self);
	cl.def(pybind11::self -= pybind11::self);

	cl.def(pybind11::self +  T());
	cl.def(pybind11::self -  T());
    cl.def(pybind11::self *= T());
    cl.def(pybind11::self /= T());
}






//using Vector = utility::vectorL<L, T, Allocator>;
//template< typename T, typename Allocator > using vector0_binder = vectorL_binder<0, T, Allocator>;
//template< typename T, typename Allocator > using vector1_binder = vectorL_binder<1, T, Allocator>;


// template< typename T, typename Allocator >
// class vector0_binder : public vectorL_binder<0, T, Allocator>
// {
// public:
// 	using vectorL_binder<0, T, Allocator>::vectorL_binder;
// 	vector0_binder(pybind11::module &m, std::string const &name, std::string const & allocator) : vector0_binder(m, "0", name, allocator) {}
// };


// template< typename T, typename Allocator >
// class vector1_binder : public vectorL_binder<1, T, Allocator>
// {
// public:
// 	using vectorL_binder<1, T, Allocator>::vectorL_binder;
// 	vector1_binder(pybind11::module &m, std::string const &name, std::string const & allocator) : vector1_binder(m, "1", name, allocator) {}
// };



		// cl.def("vector", (const class std::vector<class core::scoring::rna::data::RNA_Reactivity, class std::allocator<class core::scoring::rna::data::RNA_Reactivity> > & (utility::vectorL<1,core::scoring::rna::data::RNA_Reactivity,std::allocator<core::scoring::rna::data::RNA_Reactivity>>::*)() const) &utility::vectorL<1, core::scoring::rna::data::RNA_Reactivity, std::allocator<core::scoring::rna::data::RNA_Reactivity> >::vector, "doc");
		// cl.def("vector", (class std::vector<class core::scoring::rna::data::RNA_Reactivity, class std::allocator<class core::scoring::rna::data::RNA_Reactivity> > & (utility::vectorL<1,core::scoring::rna::data::RNA_Reactivity,std::allocator<core::scoring::rna::data::RNA_Reactivity>>::*)()) &utility::vectorL<1, core::scoring::rna::data::RNA_Reactivity, std::allocator<core::scoring::rna::data::RNA_Reactivity> >::vector, "doc");
		// cl.def("append", (class utility::vectorL<1, class core::scoring::rna::data::RNA_Reactivity, class std::allocator<class core::scoring::rna::data::RNA_Reactivity> > & (utility::vectorL<1,core::scoring::rna::data::RNA_Reactivity,std::allocator<core::scoring::rna::data::RNA_Reactivity>>::*)(const class utility::vectorL<1, class core::scoring::rna::data::RNA_Reactivity, class std::allocator<class core::scoring::rna::data::RNA_Reactivity> > &)) &utility::vectorL<1, core::scoring::rna::data::RNA_Reactivity, std::allocator<core::scoring::rna::data::RNA_Reactivity> >::append, "doc", pybind11::arg("v"));
		// cl.def("add_back", (class utility::vectorL<1, class core::scoring::rna::data::RNA_Reactivity, class std::allocator<class core::scoring::rna::data::RNA_Reactivity> > & (utility::vectorL<1,core::scoring::rna::data::RNA_Reactivity,std::allocator<core::scoring::rna::data::RNA_Reactivity>>::*)(const class core::scoring::rna::data::RNA_Reactivity &)) &utility::vectorL<1, core::scoring::rna::data::RNA_Reactivity, std::allocator<core::scoring::rna::data::RNA_Reactivity> >::add_back, "doc", pybind11::arg("t"));
		// cl.def("remove_back", (class utility::vectorL<1, class core::scoring::rna::data::RNA_Reactivity, class std::allocator<class core::scoring::rna::data::RNA_Reactivity> > & (utility::vectorL<1,core::scoring::rna::data::RNA_Reactivity,std::allocator<core::scoring::rna::data::RNA_Reactivity>>::*)()) &utility::vectorL<1, core::scoring::rna::data::RNA_Reactivity, std::allocator<core::scoring::rna::data::RNA_Reactivity> >::remove_back, "doc");
		// cl.def("shrink", (void (utility::vectorL<1,core::scoring::rna::data::RNA_Reactivity,std::allocator<core::scoring::rna::data::RNA_Reactivity>>::*)()) &utility::vectorL<1, core::scoring::rna::data::RNA_Reactivity, std::allocator<core::scoring::rna::data::RNA_Reactivity> >::shrink, "doc");
		// cl.def("has", (bool (utility::vectorL<1,core::scoring::rna::data::RNA_Reactivity,std::allocator<core::scoring::rna::data::RNA_Reactivity>>::*)(const unsigned long) const) &utility::vectorL<1, core::scoring::rna::data::RNA_Reactivity, std::allocator<core::scoring::rna::data::RNA_Reactivity> >::has, "doc", pybind11::arg("i"));
		// cl.def("at", (const class core::scoring::rna::data::RNA_Reactivity & (utility::vectorL<1,core::scoring::rna::data::RNA_Reactivity,std::allocator<core::scoring::rna::data::RNA_Reactivity>>::*)(const unsigned long) const) &utility::vectorL<1, core::scoring::rna::data::RNA_Reactivity, std::allocator<core::scoring::rna::data::RNA_Reactivity> >::at, "doc", pybind11::arg("i"));
		// cl.def("at", (class core::scoring::rna::data::RNA_Reactivity & (utility::vectorL<1,core::scoring::rna::data::RNA_Reactivity,std::allocator<core::scoring::rna::data::RNA_Reactivity>>::*)(const unsigned long)) &utility::vectorL<1, core::scoring::rna::data::RNA_Reactivity, std::allocator<core::scoring::rna::data::RNA_Reactivity> >::at, "doc", pybind11::arg("i"));
		// cl.def("swap", (void (utility::vectorL<1,core::scoring::rna::data::RNA_Reactivity,std::allocator<core::scoring::rna::data::RNA_Reactivity>>::*)(class utility::vectorL<1, class core::scoring::rna::data::RNA_Reactivity, class std::allocator<class core::scoring::rna::data::RNA_Reactivity> > &)) &utility::vectorL<1, core::scoring::rna::data::RNA_Reactivity, std::allocator<core::scoring::rna::data::RNA_Reactivity> >::swap, "doc", pybind11::arg("v"));



} // namespace rosetta_binders

#endif // _INCLUDED_rosetta_binders_hpp_
