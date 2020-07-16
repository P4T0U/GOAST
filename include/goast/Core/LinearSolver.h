// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include <Eigen/UmfPackSupport>
#include <Eigen/CholmodSupport>
#include <Eigen/SparseLU>

enum LINEAR_SOLVER_TYPE {
  UMFPACK_LU_FACT = 1,
  CHOLMOD_LLT = 2,
  EIGEN_CG = 3,
  EIGEN_BICGSTAB = 4,
  EIGEN_LU = 5
};

/**
 * \brief Runtime configurable interface to linear solver from Eigen
 * \tparam ConfiguratorType Container with datatypes
 */
template<typename ConfiguratorType>
class LinearSolver {

protected:
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;

  int _solverType;
  bool _isPrepared;

  Eigen::UmfPackLU<MatrixType> _LUSolver;
  Eigen::SparseLU<MatrixType> _eLUSolver;
  Eigen::CholmodSupernodalLLT<MatrixType> _CholmodSolver;
//  Eigen::CholmodDecomposition<MatrixType> _CholmodSolver;
  Eigen::ConjugateGradient<MatrixType, Eigen::Lower | Eigen::Upper> _CGSolver;
  Eigen::BiCGSTAB<MatrixType> _BiCGSolver;

public:
  LinearSolver( LINEAR_SOLVER_TYPE linearSolver = UMFPACK_LU_FACT ) : _solverType( linearSolver ), _isPrepared(false) {}

  //
  void prepareSolver( const MatrixType &systemMatrix ) {
    switch ( _solverType ) {
      case UMFPACK_LU_FACT: {
        _LUSolver.compute( systemMatrix );
        if ( _LUSolver.info() != Eigen::Success )
          throw BasicException( "LinearSolver::prepareSolver: UmfPackLU solver failed!" );
        break;
      }

      case CHOLMOD_LLT: {
        _CholmodSolver.compute( systemMatrix );
        if ( _CholmodSolver.info() != Eigen::Success )
          throw BasicException( "LinearSolver::prepareSolver: CHOLMOD_LLT solver failed!" );
        break;
      }

      case EIGEN_LU: {
        _eLUSolver.compute( systemMatrix );
        if ( _eLUSolver.info() != Eigen::Success )
          throw BasicException( "LinearSolver::prepareSolver: Eigen LU solver failed!" );
        break;
      }

      case EIGEN_CG: {
        _CGSolver.compute( systemMatrix );
        break;
      }

      case EIGEN_BICGSTAB : {
        _BiCGSolver.compute( systemMatrix );
      }

      default:
        throw BasicException( "LinearSolver::prepareSolver: unknown solver type!" );
        break;
    }

    _isPrepared = true;
  }

  //
  void backSubstitute( const VectorType &rhs, VectorType &solution ) const {

    if( !_isPrepared )
      throw BasicException("LinearSolver::backSubstitute(): solver has not been prepared yet!");

    switch ( _solverType ) {

      case UMFPACK_LU_FACT: {
        solution = _LUSolver.solve( rhs );
        break;
      }

      case CHOLMOD_LLT: {
        solution = _CholmodSolver.solve( rhs );
        break;
      }

      case EIGEN_LU: {
        solution = _eLUSolver.solve( rhs );
        break;
      }

      case EIGEN_CG: {
        solution = _CGSolver.solve( rhs );
        break;
      }

      case EIGEN_BICGSTAB : {
        solution = _BiCGSolver.solve( rhs );
      }

      default:
        throw BasicException( "LinearSolver::backSubstitute: unknown solver type!" );
        break;
    }
  }

  // solve Ax = b for x
  void solve( const MatrixType &A, const VectorType &b, VectorType &x ) {
    prepareSolver( A );
    backSubstitute( b, x );
  }

};

#endif //LINEARSOLVER_H
