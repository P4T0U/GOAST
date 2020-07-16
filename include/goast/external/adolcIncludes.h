// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Includes to use the ADOL-C autodiff library
 * \author Sassen
 * \warning This is still work in progress!
 */

#ifndef _ADOLCINCLUDES_H
#define _ADOLCINCLUDES_H

#ifdef GOAST_WITH_ADOLC

#include <goast/Core/EigenIncludes.h>

#include <adolc/adolc.h>
#ifdef _OPENMP
#include <adolc/adolc_openmp.h>
#endif

// Make the adouble type Eigen compatible
inline const adouble &conj( const adouble &x ) { return x; }
inline const adouble &real( const adouble &x ) { return x; }
inline adouble imag( const adouble & ) { return 0.; }
inline adouble abs( const adouble &x ) { return fabs( x ); }
inline adouble abs2( const adouble &x ) { return x * x; }

namespace Eigen {
  template<>
  struct NumTraits<adouble>
          : NumTraits<double> {
    typedef adouble Real;
    typedef adouble NonInteger;
    typedef adouble Nested;
    enum {
      IsComplex = 0,
      IsInteger = 0,
      IsSigned = 1,
      RequireInitialization = 1,
      ReadCost = 1,
      AddCost = 3,
      MulCost = 3
    };
  };
}

#ifdef ADOLC_TRACELESS
#include <adolc/adtl.h>

// Make the adtl::adouble type Eigen-compatible
namespace adtl {
  inline const adouble &conj( const adouble &x ) { return x; }
  inline const adouble &real( const adouble &x ) { return x; }
  inline adouble imag( const adouble & ) { return 0.; }
  inline adouble abs( const adouble &x ) { return fabs( x ); }
  inline adouble abs2( const adouble &x ) { return x * x; }
}

namespace Eigen {
  template<>
  struct NumTraits<adtl::adouble>
          : NumTraits<double> {
    typedef adtl::adouble Real;
    typedef adtl::adouble NonInteger;
    typedef adtl::adouble Nested;
    enum {
      IsComplex = 0,
      IsInteger = 0,
      IsSigned = 1,
      RequireInitialization = 1,
      ReadCost = 1,
      AddCost = 10,
      MulCost = 10
    };
  };
}

#endif // ADOLC_TRACELESS

#endif // GOAST_WITH_ADOLC

#endif //_ADOLCINCLUDES_H
