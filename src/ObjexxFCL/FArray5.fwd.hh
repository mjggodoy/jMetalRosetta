#ifndef INCLUDED_ObjexxFCL_FArray5_fwd_hh
#define INCLUDED_ObjexxFCL_FArray5_fwd_hh


// FArray5 Forward Declarations
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


// C++ Headers
#include <cstddef>
#include <string>


namespace ObjexxFCL {


// Forward Declarations
template< typename > class FArray5;
class byte;
typedef  byte  sbyte;
class ubyte;
class Fstring;


// Types
typedef  FArray5< bool >                FArray5_bool;
typedef  FArray5< byte >                FArray5_byte;
typedef  FArray5< sbyte >               FArray5_sbyte;
typedef  FArray5< ubyte >               FArray5_ubyte;
typedef  FArray5< short int >           FArray5_short;
typedef  FArray5< int >                 FArray5_int;
typedef  FArray5< long int >            FArray5_long;
typedef  FArray5< unsigned short int >  FArray5_ushort;
typedef  FArray5< unsigned int >        FArray5_uint;
typedef  FArray5< unsigned long int >   FArray5_ulong;
typedef  FArray5< std::size_t >         FArray5_size_t;
typedef  FArray5< std::size_t >         FArray5_size;
typedef  FArray5< float >               FArray5_float;
typedef  FArray5< double >              FArray5_double;
typedef  FArray5< long double >         FArray5_longdouble;
typedef  FArray5< char >                FArray5_char;
typedef  FArray5< unsigned char >       FArray5_uchar;
typedef  FArray5< signed char >         FArray5_schar;
typedef  FArray5< std::string >         FArray5_string;
typedef  FArray5< Fstring >             FArray5_Fstring;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray5_fwd_HH
