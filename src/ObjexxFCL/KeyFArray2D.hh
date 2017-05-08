#ifndef INCLUDED_ObjexxFCL_KeyFArray2D_hh
#define INCLUDED_ObjexxFCL_KeyFArray2D_hh


// KeyFArray2D: Key-Access Fortran-Compatible 2D Array
//
// Project: Objexx Fortran Compatibility Library (ObjexxFCL)
//
// Version: 3.0.0
//
// Language: C++
//
// Copyright (c) 2000-2009 Objexx Engineering, Inc. All Rights Reserved.
// Use of this source code or any derivative of it is restricted by license.
// Licensing is available from Objexx Engineering, Inc.:  http://objexx.com  Objexx@objexx.com


// ObjexxFCL Headers
#include <ObjexxFCL/KeyFArray2D.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArrayInitializer.hh>


namespace ObjexxFCL {


/// @brief KeyFArray2D: Key-Access Fortran-Compatible 2D Array
template< typename T >
class KeyFArray2D :
	public FArray2< T >,
	public ObserverMulti
{


private: // Types


	typedef  FArray2< T >  Super;
	typedef  internal::InitializerSentinel  InitializerSentinel;


private: // Friend


	template< typename > friend class KeyFArray2D;
	friend class FArray2D< T >;
	friend class FArray2P< T >;
	friend class FArray2A< T >;


public: // Types


	typedef  typename Super::Base  Base;
	typedef  typename Base::Section  Section;
	typedef  typename Super::IR  SIR;
	typedef  DynamicIndexRange  IR;

	// STL Style
	typedef  typename Base::value_type  value_type;
	typedef  typename Base::reference  reference;
	typedef  typename Base::const_reference  const_reference;
	typedef  typename Base::pointer  pointer;
	typedef  typename Base::const_pointer  const_pointer;
	typedef  typename Base::size_type  size_type;
	typedef  typename Base::difference_type  difference_type;

	// C++ Style
	typedef  typename Base::Value  Value;
	typedef  typename Base::Reference  Reference;
	typedef  typename Base::ConstReference  ConstReference;
	typedef  typename Base::Pointer  Pointer;
	typedef  typename Base::ConstPointer  ConstPointer;
	typedef  typename Base::Size  Size;
	typedef  typename Base::Difference  Difference;

	typedef  FArrayInitializer< T, ObjexxFCL::KeyFArray2D >  Initializer;
	typedef  typename Initializer::function_type  InitializerFunction;

	using Super::array_;
	using Super::array_size_;
	using Super::sarray_;
	using Super::shift_;
	using Super::shift_set;
	using Super::size_;
	using Super::size_of;
	using Super::s1_;


public: // Creation


	/// @brief Default Constructor
	inline
	KeyFArray2D()
	{
		insert_as_observer();
	}


	/// @brief Copy Constructor
	inline
	KeyFArray2D( KeyFArray2D const & a ) :
		Super( a ),
		ObserverMulti(),
		I1_( a.I1_ ),
		I2_( a.I2_ )
	{
		insert_as_observer();
	}


	/// @brief Copy Constructor Template
	template< typename U >
	inline
	explicit
	KeyFArray2D( KeyFArray2D< U > const & a ) :
		Super( a ),
		I1_( a.I1_ ),
		I2_( a.I2_ )
	{
		insert_as_observer();
	}


	/// @brief Super Constructor Template
	template< typename U >
	inline
	explicit
	KeyFArray2D( FArray2< U > const & a ) :
		Super( a ),
		I1_( a.I1() ),
		I2_( a.I2() )
	{
		insert_as_observer();
	}


	/// @brief IndexRange Constructor
	inline
	KeyFArray2D( IR const & I1_a, IR const & I2_a ) :
		Super( size_of( I1_a, I2_a ) ),
		I1_( I1_a ),
		I2_( I2_a )
	{
		setup_real();
		insert_as_observer();
	}


	/// @brief IndexRange + Initializer Value Constructor
	inline
	KeyFArray2D( IR const & I1_a, IR const & I2_a, T const & t ) :
		Super( size_of( I1_a, I2_a ), InitializerSentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		initializer_( t )
	{
		setup_real();
		initialize();
		insert_as_observer();
	}


	/// @brief IndexRange + Initializer Function Constructor
	inline
	KeyFArray2D( IR const & I1_a, IR const & I2_a, InitializerFunction const & function_a ) :
		Super( size_of( I1_a, I2_a ), InitializerSentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		initializer_( function_a )
	{
		setup_real();
		initialize();
		insert_as_observer();
	}


	/// @brief Super + IndexRange Constructor Template
	template< typename U >
	inline
	KeyFArray2D( FArray2< U > const & a, IR const & I1_a, IR const & I2_a ) :
		Super( size_of( I1_a, I2_a ) ),
		I1_( I1_a ),
		I2_( I2_a )
	{
		setup_real();
		if ( dimensions_initialized() ) {
			if ( a.dimensions_initialized() ) { // Copy array data where overlap
				int const b1( std::max( I1_.l(), a.l1() ) ), e1( std::min( I1_.u(), a.u1() ) );
				int const b2( std::max( I2_.l(), a.l2() ) ), e2( std::min( I2_.u(), a.u2() ) );
				for ( int i2 = b2; i2 <= e2; ++i2 ) {
					for ( int i1 = b1; i1 <= e1; ++i1 ) {
						operator ()( i1, i2 ) = T( a( i1, i2 ) );
					}
				}
			}
		}
		insert_as_observer();
	}


	/// @brief Super + IndexRange + Fill Value Constructor Template
	template< typename U >
	inline
	KeyFArray2D( FArray2< U > const & a, IR const & I1_a, IR const & I2_a, T const & t ) :
		Super( size_of( I1_a, I2_a ) ),
		I1_( I1_a ),
		I2_( I2_a )
	{
		setup_real();
		if ( dimensions_initialized() ) {
			(*this) = t; // Initialize array with fill value
			if ( a.dimensions_initialized() ) { // Copy array data where overlap
				int const b1( std::max( I1_.l(), a.l1() ) ), e1( std::min( I1_.u(), a.u1() ) );
				int const b2( std::max( I2_.l(), a.l2() ) ), e2( std::min( I2_.u(), a.u2() ) );
				for ( int i2 = b2; i2 <= e2; ++i2 ) {
					for ( int i1 = b1; i1 <= e1; ++i1 ) {
						operator ()( i1, i2 ) = T( a( i1, i2 ) );
					}
				}
			}
		}
		insert_as_observer();
	}


	/// @brief Diagonal Matrix Named Constructor
	inline
	static
	KeyFArray2D
	diag( IR const & I_a, T const & d )
	{
		KeyFArray2D D( I_a, I_a );
		D.to_diag( d );
		return D;
	}


	/// @brief Identity Matrix Named Constructor
	inline
	static
	KeyFArray2D
	identity( IR const & I_a )
	{
		KeyFArray2D D( I_a, I_a );
		D.to_diag( T( 1 ) );
		return D;
	}


	/// @brief Destructor
	inline
	virtual
	~KeyFArray2D()
	{}


public: // Assignment


	/// @brief Copy Assignment
	inline
	KeyFArray2D &
	operator =( KeyFArray2D const & a )
	{
		if ( this != &a ) {
			if ( ! equal_dimension( a ) ) dimension( a );
			Base::operator =( a );
		}
		return *this;
	}


	/// @brief Super Assignment
	inline
	KeyFArray2D &
	operator =( Super const & a )
	{
		if ( this != &a ) {
			if ( ! equal_dimension( a ) ) dimension( a );
			Base::operator =( a );
		}
		return *this;
	}


	/// @brief Super Assignment Template
	template< typename U >
	inline
	KeyFArray2D &
	operator =( FArray2< U > const & a )
	{
		if ( ! equal_dimension( a ) ) dimension( a );
		Base::operator =( a );
		return *this;
	}


	/// @brief += Array Template
	template< typename U >
	inline
	KeyFArray2D &
	operator +=( FArray2< U > const & a )
	{
		Super::operator +=( a );
		return *this;
	}


	/// @brief *= Array Template
	template< typename U >
	inline
	KeyFArray2D &
	operator *=( FArray2< U > const & a )
	{
		Super::operator *=( a );
		return *this;
	}


	/// @brief -= Array Template
	template< typename U >
	inline
	KeyFArray2D &
	operator -=( FArray2< U > const & a )
	{
		Super::operator -=( a );
		return *this;
	}


	/// @brief = Value
	inline
	KeyFArray2D &
	operator =( T const & t )
	{
		Super::operator =( t );
		return *this;
	}


	/// @brief += Value
	inline
	KeyFArray2D &
	operator +=( T const & t )
	{
		Super::operator +=( t );
		return *this;
	}


	/// @brief -= Value
	inline
	KeyFArray2D &
	operator -=( T const & t )
	{
		Super::operator -=( t );
		return *this;
	}


	/// @brief *= Value
	inline
	KeyFArray2D &
	operator *=( T const & t )
	{
		Super::operator *=( t );
		return *this;
	}


	/// @brief /= Value
	inline
	KeyFArray2D &
	operator /=( T const & t )
	{
		Super::operator /=( t );
		return *this;
	}


public: // Subscript


	/// @brief array( i1, i2 ) const
	template< typename K1, typename K2 >
	inline
	T const &
	operator ()( K1 const & i1, K2 const & i2 ) const
	{
		assert( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) );
		return sarray_[ ( i2 * s1_ ) + i1 ];
	}


	/// @brief array( i1, i2 )
	template< typename K1, typename K2 >
	inline
	T &
	operator ()( K1 const & i1, K2 const & i2 )
	{
		assert( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) );
		return sarray_[ ( i2 * s1_ ) + i1 ];
	}


	/// @brief Section Starting at array( i1, i2 )
	template< typename K1, typename K2 >
	inline
	Section const
	a( K1 const & i1, K2 const & i2 ) const
	{
		assert( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) );
		size_type const offset( ( ( i2 * s1_ ) + i1 ) - shift_ );
		return Section( array_size_ - offset, array_ + offset );
	}


	/// @brief Linear Index
	template< typename K1, typename K2 >
	inline
	size_type
	index( K1 const & i1, K2 const & i2 ) const
	{
		assert( ( I1_.initialized() ) && ( I2_.initialized() ) );
		return ( ( ( i2 * s1_ ) + i1 ) - shift_ );
	}


	/// @brief array[ i ] const: Linear Subscript
	inline
	T const &
	operator []( size_type const i ) const
	{
		assert( i < size_ );
		return array_[ i ];
	}


	/// @brief array[ i ]: Linear Subscript
	inline
	T &
	operator []( size_type const i )
	{
		assert( i < size_ );
		return array_[ i ];
	}


public: // Predicate


	/// @brief Dimensions Initialized?
	inline
	bool
	dimensions_initialized() const
	{
		return ( ( I1_.initialized() ) && ( I2_.initialized() ) );
	}


	/// @brief Contains Indexed Element?
	template< typename K1, typename K2 >
	inline
	bool
	contains( K1 const & i1, K2 const & i2 ) const
	{
		return ( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) );
	}


	/// @brief Initializer Active?
	inline
	bool
	initializer_active() const
	{
		return initializer_.is_active();
	}


public: // Inspector


	/// @brief IndexRange of Dimension 1
	inline
	IR const &
	I1() const
	{
		return I1_;
	}


	/// @brief Lower Index of Dimension 1
	inline
	int
	l1() const
	{
		return I1_.l();
	}


	/// @brief Upper Index of Dimension 1
	inline
	int
	u1() const
	{
		return I1_.u();
	}


	/// @brief IndexRange of Dimension 2
	inline
	IR const &
	I2() const
	{
		return I2_;
	}


	/// @brief Lower Index of Dimension 2
	inline
	int
	l2() const
	{
		return I2_.l();
	}


	/// @brief Upper Index of Dimension 2
	inline
	int
	u2() const
	{
		return I2_.u();
	}


	/// @brief Size of Dimension 2
	inline
	size_type
	size2() const
	{
		return I2_.size();
	}


public: // Modifier


	/// @brief Clear
	inline
	KeyFArray2D &
	clear()
	{
		Super::clear();
		I1_.clear_no_notify();
		I2_.clear_no_notify();
		initializer_.clear();
		notify();
		return *this;
	}


	/// @brief Dimension by IndexRange
	inline
	KeyFArray2D &
	dimension( IR const & I1_a, IR const & I2_a )
	{
		initializer_.clear();
		I1_.assign_no_notify( I1_a );
		I2_.assign_no_notify( I2_a );
		dimension_real();
		notify();
		return *this;
	}


	/// @brief Dimension by IndexRange + Initializer Value
	inline
	KeyFArray2D &
	dimension( IR const & I1_a, IR const & I2_a, T const & t )
	{
		initializer_ = t;
		I1_.assign_no_notify( I1_a );
		I2_.assign_no_notify( I2_a );
		dimension_real();
		initialize();
		notify();
		return *this;
	}


	/// @brief Dimension by IndexRange + Initializer Function
	inline
	KeyFArray2D &
	dimension( IR const & I1_a, IR const & I2_a, InitializerFunction const & function_a )
	{
		initializer_ = function_a;
		I1_.assign_no_notify( I1_a );
		I2_.assign_no_notify( I2_a );
		dimension_real();
		initialize();
		notify();
		return *this;
	}


	/// @brief Dimension by Array Template
	template< typename U >
	inline
	KeyFArray2D &
	dimension( FArray2< U > const & a )
	{
		initializer_.clear();
		I1_.assign_no_notify( a.I1() );
		I2_.assign_no_notify( a.I2() );
		dimension_real();
		notify();
		return *this;
	}


	/// @brief Dimension by Array + Initializer Value Template
	template< typename U >
	inline
	KeyFArray2D &
	dimension( FArray2< U > const & a, T const & t )
	{
		initializer_ = t;
		I1_.assign_no_notify( a.I1() );
		I2_.assign_no_notify( a.I2() );
		dimension_real();
		initialize();
		notify();
		return *this;
	}


	/// @brief Dimension by Array + Initializer Function Template
	template< typename U >
	inline
	KeyFArray2D &
	dimension( FArray2< U > const & a, InitializerFunction const & function_a )
	{
		initializer_ = function_a;
		I1_.assign_no_notify( a.I1() );
		I2_.assign_no_notify( a.I2() );
		dimension_real();
		initialize();
		notify();
		return *this;
	}


	/// @brief Data-Preserving Redimension by IndexRange
	inline
	KeyFArray2D &
	redimension( IR const & I1_a, IR const & I2_a )
	{
		KeyFArray2D( *this, I1_a, I2_a ).swap( *this );
		return *this;
	}


	/// @brief Data-Preserving Redimension by IndexRange + Fill Value
	inline
	KeyFArray2D &
	redimension( IR const & I1_a, IR const & I2_a, T const & t )
	{
		KeyFArray2D( *this, I1_a, I2_a, t ).swap( *this );
		return *this;
	}


	/// @brief Data-Preserving Redimension by Array Template
	template< typename U >
	inline
	KeyFArray2D &
	redimension( FArray2< U > const & a )
	{
		KeyFArray2D( *this, a.I1(), a.I2() ).swap( *this );
		return *this;
	}


	/// @brief Data-Preserving Redimension by Array + Fill Value Template
	template< typename U >
	inline
	KeyFArray2D &
	redimension( FArray2< U > const & a, T const & t )
	{
		KeyFArray2D( *this, a.I1(), a.I2(), t ).swap( *this );
		return *this;
	}


	/// @brief Set Initializer Value
	inline
	KeyFArray2D &
	initializer( T const & t )
	{
		initializer_ = t;
		return *this;
	}


	/// @brief Set Initializer Function
	inline
	KeyFArray2D &
	initializer( InitializerFunction const & function_a )
	{
		initializer_ = function_a;
		return *this;
	}


	/// @brief Clear Initializer
	inline
	KeyFArray2D &
	initializer_clear()
	{
		initializer_.clear();
		return *this;
	}


	/// @brief Initialize
	inline
	KeyFArray2D &
	initialize()
	{
		if ( ( initializer_.is_active() ) && ( dimensions_initialized() ) ) {
			if ( initializer_.is_value() ) {
				(*this) = initializer_.value();
			} else if ( initializer_.is_function() ) {
				initializer_.function()( *this );
			}
		}
		return *this;
	}


	/// @brief Swap
	inline
	KeyFArray2D &
	swap( KeyFArray2D & v )
	{
		swap2DB( v );
		I1_.swap_no_notify( v.I1_ );
		I2_.swap_no_notify( v.I2_ );
		std::swap( initializer_, v.initializer_ );
		notify(); // So proxy FArrays can reattach
		v.notify(); // So proxy FArrays can reattach
		return *this;
	}


public: // Observer Modifier


	/// @brief Update
	inline
	void
	update()
	{
		dimension_real();
		initialize();
	}


	/// @brief Update for Destruction of a Subject
	inline
	void
	destructed( Subject const & )
	{}


public: // Friend


	/// @brief Swap
	friend
	inline
	void
	swap( KeyFArray2D & a, KeyFArray2D & b )
	{
		a.swap( b );
	}


public: // Generator


	/// @brief -Array
	friend
	inline
	KeyFArray2D
	operator -( KeyFArray2D const & a )
	{
		KeyFArray2D r( a );
		r *= T( -1 );
		return r;
	}


	/// @brief Array + Array
	friend
	inline
	KeyFArray2D
	operator +( KeyFArray2D const & a, KeyFArray2D const & b )
	{
		KeyFArray2D r( a );
		r += b;
		return r;
	}


	/// @brief Array - Array
	friend
	inline
	KeyFArray2D
	operator -( KeyFArray2D const & a, KeyFArray2D const & b )
	{
		KeyFArray2D r( a );
		r -= b;
		return r;
	}


	/// @brief Array + Value
	friend
	inline
	KeyFArray2D
	operator +( KeyFArray2D const & a, T const & t )
	{
		KeyFArray2D r( a );
		r += t;
		return r;
	}


	/// @brief Value + Array
	friend
	inline
	KeyFArray2D
	operator +( T const & t, KeyFArray2D const & a )
	{
		KeyFArray2D r( a );
		r += t;
		return r;
	}


	/// @brief Array - Value
	friend
	inline
	KeyFArray2D
	operator -( KeyFArray2D const & a, T const & t )
	{
		KeyFArray2D r( a );
		r -= t;
		return r;
	}


	/// @brief Value - Array
	friend
	inline
	KeyFArray2D
	operator -( T const & t, KeyFArray2D const & a )
	{
		KeyFArray2D r( a );
		r *= T( -1 );
		r += t;
		return r;
	}


	/// @brief Array * Value
	friend
	inline
	KeyFArray2D
	operator *( KeyFArray2D const & a, T const & t )
	{
		KeyFArray2D r( a );
		r *= t;
		return r;
	}


	/// @brief Value * Array
	friend
	inline
	KeyFArray2D
	operator *( T const & t, KeyFArray2D const & a )
	{
		KeyFArray2D r( a );
		r *= t;
		return r;
	}


	/// @brief Array / Value
	friend
	inline
	KeyFArray2D
	operator /( KeyFArray2D const & a, T const & t )
	{
		KeyFArray2D r( a );
		r /= t;
		return r;
	}


	/// @brief Transposed
	friend
	inline
	KeyFArray2D
	transposed( KeyFArray2D const & a )
	{
		assert( a.square() );
		KeyFArray2D aT( a.I1(), a.I2() );
		for ( int i = a.l1(), ie = a.u1(); i <= ie; ++i ) {
			for ( int j = a.l2(), je = a.u2(); j <= je; ++j ) {
				aT( i, j ) = a( j, i );
			}
		}
		return aT;
	}


protected: // Functions


	/// @brief Dimension by IndexRanges
	inline
	void
	dimension_assign( SIR const & I1_a, SIR const & I2_a )
	{
		initializer_.clear();
		I1_.assign_no_notify( I1_a );
		I2_.assign_no_notify( I2_a );
		dimension_real();
		notify();
	}


private: // Functions


	/// @brief Setup for IndexRange Constructor
	inline
	void
	setup_real()
	{
		s1_ = I1_.size();
		if ( dimensions_initialized() ) {
			shift_set( ( I2_.lz() * s1_ ) + I1_.lz() );
		} else {
			shift_set( 0 );
		}
	}


	/// @brief Dimension by Current IndexRanges
	inline
	void
	dimension_real()
	{
		s1_ = I1_.size();
		if ( dimensions_initialized() ) {
			resize( size_of( s1_, I2_.size() ) );
			shift_set( ( I2_.lz() * s1_ ) + I1_.lz() );
#ifdef OBJEXXFCL_FARRAY_SIZE_REPORT
			size_report();
#endif // OBJEXXFCL_FARRAY_SIZE_REPORT
		} else {
			Base::clear();
		}
	}


	/// @brief Insert as Observer of the IndexRanges
	inline
	void
	insert_as_observer()
	{
		I1_.insert_observer( *this );
		I2_.insert_observer( *this );
#ifdef OBJEXXFCL_FARRAY_SIZE_REPORT
		size_report();
#endif // OBJEXXFCL_FARRAY_SIZE_REPORT
	}

	/* // unused private
	/// @brief Remove as Observer of the IndexRanges
	inline
	void
	remove_as_observer()
	{
		I1_.remove_observer( *this );
		I2_.remove_observer( *this );
	}*/


#ifdef OBJEXXFCL_FARRAY_SIZE_REPORT
	/// @brief Report size if at least value defined for OBJEXXFCL_FARRAY_SIZE_REPORT
	/// @note  Size is based on sizeof( T ) so T-controlled heap memory is not counted
	inline
	void
	size_report() const
	{
		if ( size_ * sizeof( T ) >= OBJEXXFCL_FARRAY_SIZE_REPORT ) {
			std::cout << "  Index ranges: " << I1_ << ' ' << I2_ << std::endl;
		}
	}
#endif // OBJEXXFCL_FARRAY_SIZE_REPORT


private: // Data


	/// @brief Dimension 1 index range
	IR I1_;

	/// @brief Dimension 2 index range
	IR I2_;

	/// @brief Array initializer
	Initializer initializer_;


}; // KeyFArray2D


/// @brief Swap
template< typename T >
void
swap( KeyFArray2D< T > & a, KeyFArray2D< T > & b );


/// @brief -Array
template< typename T >
KeyFArray2D< T >
operator -( KeyFArray2D< T > const & a );


/// @brief Array + Array
template< typename T >
KeyFArray2D< T >
operator +( KeyFArray2D< T > const & a, KeyFArray2D< T > const & b );


/// @brief Array - Array
template< typename T >
KeyFArray2D< T >
operator -( KeyFArray2D< T > const & a, KeyFArray2D< T > const & b );


/// @brief Array + Value
template< typename T >
KeyFArray2D< T >
operator +( KeyFArray2D< T > const & a, T const & t );


/// @brief Value + Array
template< typename T >
KeyFArray2D< T >
operator +( T const & t, KeyFArray2D< T > const & a );


/// @brief Array - Value
template< typename T >
KeyFArray2D< T >
operator -( KeyFArray2D< T > const & a, T const & t );


/// @brief Value - Array
template< typename T >
KeyFArray2D< T >
operator -( T const & t, KeyFArray2D< T > const & a );


/// @brief Array * Value
template< typename T >
KeyFArray2D< T >
operator *( KeyFArray2D< T > const & a, T const & t );


/// @brief Value * Array
template< typename T >
KeyFArray2D< T >
operator *( T const & t, KeyFArray2D< T > const & a );


/// @brief Array / Value
template< typename T >
KeyFArray2D< T >
operator /( KeyFArray2D< T > const & a, T const & t );


/// @brief Transposed
template< typename T >
KeyFArray2D< T >
transposed( KeyFArray2D< T > const & a );


} // namespace ObjexxFCL


#ifndef NO_STD_SWAP_OVERLOADS


// std::swap Overloads for Efficiency
//
// Technically you cannot add template functions overloads to namespace std
// but this works with most compilers and makes it much faster if someone uses
// std::swap instead of swap or ObjexxFCL::swap.  The legal alternative would be
// to add specializations of swap for each anticipated instantiation.


namespace std {


/// @brief std::swap( KeyFArray2D, KeyFArray2D )
template< typename T >
inline
void
swap( ObjexxFCL::KeyFArray2D< T > & a, ObjexxFCL::KeyFArray2D< T > & b )
{
	a.swap( b );
}


} // namespace std


#endif // NO_STD_SWAP_OVERLOADS


#endif // INCLUDED_ObjexxFCL_KeyFArray2D_HH
