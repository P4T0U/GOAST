// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPTIMIZATION_MOSEKINTERFACE_H
#define OPTIMIZATION_MOSEKINTERFACE_H

#ifdef GOAST_WITH_MOSEK

#include <exception>

#include "../optInterface.h"

#include "mosek.h"

/**
 * \brief Custom exception for error occuring in MOSEK
 * \author Sassen
 */
class MOSEKException : public std::exception {
  const MSKrescodee _r;
  std::string _msg;
public:
  explicit MOSEKException( const MSKrescodee r, MSKenv_t &env, MSKtask_t &task ) : _r( r ) {
    MSK_deletetask( &task );
    MSK_deleteenv( &env );

    char symname[MSK_MAX_STR_LEN];
    char desc[MSK_MAX_STR_LEN];
    MSK_getcodedesc( _r, symname, desc );

    std::stringstream s;
    s << "MOSEK ERROR (" << _r << "): " << symname << " - '" << desc << "'" << std::endl;
    _msg = s.str();
  }

  const char *what() const noexcept override {
    return _msg.c_str();
  }
};

/***
 * \brief MOSEK solver for quadratic programming
 * \author Sassen
 * \tparam ConfiguratorType Container with used datatypes
 *
 * This solves a convex quadratic optimization problem with linear and constant bounds, i.e.
 *
 * minimize 1/2 * x^T * Q * x + c^T * x + cf
 * subject to lc ≤ Ax ≤ uc
 *            lx ≤ x ≤ ux
 */
template<typename ConfiguratorType>
class MOSEKQuadraticSolver : public OptimizationBase<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;

  // Problem description
  const MatrixType &_Q;
  const VectorType &_c;
  const RealType &_cf;

  const MatrixType &_A;
  const VectorType &_linearLowerBounds;
  const VectorType &_linearUpperBounds;

  const VectorType &_constantLowerBounds;
  const VectorType &_constantUpperBounds;

  // MOSEK Parameters
  std::map<MSKiparame, int> integerParameters;
  std::map<MSKdparame, double> doubleParameters;

  bool _quiet;

public:
  MOSEKQuadraticSolver( const MatrixType &Q,
                        const VectorType &c,
                        const RealType &cf,
                        const MatrixType &A,
                        const VectorType &linearLowerBounds, const VectorType &linearUpperBounds,
                        const VectorType &constantLowerBounds, const VectorType &constantUpperBounds,
                        bool quiet = false ) : _Q( Q ), _c( c ), _cf( cf ), _A( A ),
                                               _linearLowerBounds( linearLowerBounds ),
                                               _linearUpperBounds( linearUpperBounds ),
                                               _constantLowerBounds( constantLowerBounds ),
                                               _constantUpperBounds( constantUpperBounds ), _quiet( quiet ) {
    assert( Q.rows() == Q.cols());
    // \todo Check if Q is symmetric

  }

  void solve( const VectorType & /*start*/, VectorType &s ) const override {
    solve( s );
  }

  void solve( VectorType &z_k ) const {
    MSKenv_t env = nullptr;
    MSKtask_t task = nullptr;
    MSKrescodee r;

    // Size of problem
    const int numVariables = _Q.rows();
    const int numConstraints = _A.rows();

    if ( z_k.size() != numVariables )
      z_k.resize( numVariables );

    // Convert Eigen matrices to MOSEK compatible datatypes

    // Q
    std::vector<int> Qi( _Q.nonZeros(), -1 );
    std::vector<int> Qj( _Q.nonZeros(), -1 );
    std::vector<RealType> Qv( _Q.nonZeros(), 0. );

    int index = 0;
    for ( int k = 0; k < _Q.outerSize(); ++k ) {
      for ( typename MatrixType::InnerIterator it( _Q, k ); it; ++it ) {
        if ( it.row() >= it.col()) {
          Qi[index] = it.row();
          Qj[index] = it.col();
          Qv[index] = it.value();
          index++;
        }
      }
    }
    Qi.resize( index );
    Qj.resize( index );
    Qv.resize( index );

    // Prepare MOSEK

    // Create the MOSEK environment.
    r = MSK_makeenv( &env, nullptr );
    if ( r != MSK_RES_OK ) throw MOSEKException( r, env, task );

    r = MSK_maketask( env, numConstraints, numVariables, &task );
    if ( r != MSK_RES_OK ) throw MOSEKException( r, env, task );

    if ( !_quiet ) {
      r = MSK_linkfunctotaskstream( task, MSK_STREAM_LOG, nullptr, printstr );
      if ( r != MSK_RES_OK ) throw MOSEKException( r, env, task );
    }

    r = MSK_appendcons( task, numConstraints );
    if ( r != MSK_RES_OK ) throw MOSEKException( r, env, task );

    r = MSK_appendvars( task, numVariables );
    if ( r != MSK_RES_OK ) throw MOSEKException( r, env, task );

    r = MSK_putcfix( task, _cf );
    if ( r != MSK_RES_OK ) throw MOSEKException( r, env, task );

    for ( int j = 0; j < numVariables; ++j ) {
      /* Set the linear term c_j in the objective.*/
      r = MSK_putcj( task, j, _c[j] );
      if ( r != MSK_RES_OK ) throw MOSEKException( r, env, task );

      /* Set the bounds on variable j. */
      MSKboundkeye boundType;
      if ( _constantLowerBounds[j] == _constantUpperBounds[j] )
        boundType = MSK_BK_FX;
      else if ( _constantLowerBounds[j] <= -MSK_INFINITY && _constantUpperBounds[j] >= MSK_INFINITY )
        boundType = MSK_BK_FR;
      else if ( _constantLowerBounds[j] <= -MSK_INFINITY )
        boundType = MSK_BK_UP;
      else if ( _constantUpperBounds[j] >= MSK_INFINITY )
        boundType = MSK_BK_LO;
      else
        boundType = MSK_BK_RA;

      r = MSK_putvarbound( task, j, boundType, _constantLowerBounds[j], _constantUpperBounds[j] );
      if ( r != MSK_RES_OK ) throw MOSEKException( r, env, task );
    }

    // Set A, how depends on storage order:
    if ( _A.Options == 0 ) { // column-major
      for ( int j = 0; j < numVariables; ++j ) {
        r = MSK_putacol( task, j,
                         _A.outerIndexPtr()[j + 1] - _A.outerIndexPtr()[j],
                         &_A.innerIndexPtr()[_A.outerIndexPtr()[j]],
                         &_A.valuePtr()[_A.outerIndexPtr()[j]] );
        if ( r != MSK_RES_OK ) throw MOSEKException( r, env, task );
      }
    }
    else { // row-major
      for ( int i = 0; i < numConstraints; ++i ) {
        r = MSK_putarow( task, i,
                         _A.outerIndexPtr()[i + 1] - _A.outerIndexPtr()[i],
                         &_A.innerIndexPtr()[_A.outerIndexPtr()[i]],
                         &_A.valuePtr()[_A.outerIndexPtr()[i]] );
        if ( r != MSK_RES_OK ) throw MOSEKException( r, env, task );
      }
    }


    /* Set the bounds on constraints.
       for i=1, ...,NUMCON : blc[i] <= constraint i <= buc[i] */
    for ( int i = 0; i < numConstraints; ++i ) {
      MSKboundkeye boundType;
      if ( _linearLowerBounds[i] == _linearUpperBounds[i] )
        boundType = MSK_BK_FX;
      else if ( _linearLowerBounds[i] <= -MSK_INFINITY && _linearUpperBounds[i] >= MSK_INFINITY )
        boundType = MSK_BK_FR;
      else if ( _linearLowerBounds[i] <= -MSK_INFINITY )
        boundType = MSK_BK_UP;
      else if ( _linearUpperBounds[i] >= MSK_INFINITY )
        boundType = MSK_BK_LO;
      else
        boundType = MSK_BK_RA;

      r = MSK_putconbound( task, i, boundType, _linearLowerBounds[i], _linearUpperBounds[i] );
      if ( r != MSK_RES_OK ) throw MOSEKException( r, env, task );
    }

    r = MSK_putmaxnumqnz( task, _Q.nonZeros());
    if ( r != MSK_RES_OK ) throw MOSEKException( r, env, task );

    r = MSK_putqobj( task, Qv.size(), &Qi[0], &Qj[0], &Qv[0] );
    if ( r != MSK_RES_OK ) throw MOSEKException( r, env, task );

    MSKrescodee trmcode;

    /* Run optimizer */
    r = MSK_optimizetrm( task, &trmcode );

    if ( !_quiet )
      MSK_solutionsummary( task, MSK_STREAM_MSG );

    if ( r != MSK_RES_OK ) throw MOSEKException( r, env, task );

    MSKsolstae solsta;
    int j;

    MSK_getsolsta( task, MSK_SOL_ITR, &solsta );

    switch ( solsta ) {
      case MSK_SOL_STA_OPTIMAL:
        MSK_getxx( task, MSK_SOL_ITR, z_k.data());
        break;

      case MSK_SOL_STA_DUAL_INFEAS_CER:
      case MSK_SOL_STA_PRIM_INFEAS_CER:
        std::cerr << "MOSEKInterface: Primal or dual infeasibility certificate found." << std::endl;
        break;

      case MSK_SOL_STA_UNKNOWN:
        std::cerr << "The status of the solution could not be determined. Termination code: "
                  << std::to_string( trmcode ) << std::endl;
        break;

      default:
        std::cerr << "MOSEKInterface: Untreated solution status." << std::endl;
        break;
    }

    MSK_deletetask( &task );
    MSK_deleteenv( &env );
  }

  void setParameter( const std::string &name, std::string value ) override {
    throw std::runtime_error( "MOSEKQuadraticSolver::setParameter(): This class has no string parameters." );
  }

  void setParameter( const std::string &name, RealType value ) override {
    throw std::runtime_error( "MOSEKQuadraticSolver::setParameter(): Unknown parameter '" + name + "'." );
  }

  void setParameter( const std::string &name, int value ) override {
    if ( name == "print_level" )
      _quiet = (value < 5);
    else
      throw std::runtime_error( "MOSEKQuadraticSolver::setParameter(): Unknown parameter '" + name + "'." );
  }

protected:
  static void MSKAPI printstr( void *handle, const char str[] ) {
    printf( "%s", str );
  }
};

#endif // GOAST_WITH_MOSEK
#endif // OPTIMIZATION_MOSEKINTERFACE_H
