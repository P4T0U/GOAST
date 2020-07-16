// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPTIMIZATION_TRLIBINTERFACE_H
#define OPTIMIZATION_TRLIBINTERFACE_H

#ifdef GOAST_WITH_TRLIB

#include <cstdio>
extern "C"
{
#include <trlib.h>
}

#include "optInterface.h"
#include "optUtils.h"

/**
 *
 * \tparam ConfiguratorType
 * \todo Proper interface, especially expose parameters and tolerances
 * \todo Add warm start capability
 */
template<typename ConfiguratorType>
class trlibSolver : public OptimizationBase<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;

  const MatrixType &_H;
  const VectorType &_c;

  const RealType _radius;
  RealType _tolerance;
  int _maxIterations;

  bool _quiet;

public:
  trlibSolver( const MatrixType &H, const VectorType &c, const RealType radius,
               const RealType tolerance, const int maxIterations = 1000, bool quiet = false )
          : _H( H ), _c( c ), _tolerance( tolerance ), _radius( radius ),
            _maxIterations( maxIterations ), _quiet( quiet ) {}

  void solve( const VectorType & /*start*/, VectorType &s ) const override {
    solve(s);
  }

  void solve( VectorType &s ) const {
    s.resize( _c.size());
    s.setZero();

    // Setup memory
    trlib_int_t iwork_size, fwork_size, h_pointer;
    trlib_krylov_memory_size( _maxIterations, &iwork_size, &fwork_size, &h_pointer );
    auto *iwork = new trlib_int_t[iwork_size];
    auto *fwork = new trlib_flt_t[fwork_size];

    // Parameters
    trlib_int_t equality = 0; // set to 1 if trust region constraint should be enforced as equality
    trlib_int_t itmax_lanczos = 100; // maximum number of Lanczos type iterations.
    trlib_int_t ctl_invariant = TRLIB_CLC_NO_EXP_INV;
    trlib_int_t refine = 1; // set to 1 if iterative refinement should be used on solving linear systems, otherwise to 0

    trlib_int_t verbose = 1; // determines the verbosity level of output that is written to fout
    if (_quiet)
      verbose = 0;
    trlib_int_t unicode = 0;  //  set to 1 if fout can handle unicode, otherwise to 0
    trlib_flt_t tol_rel_i = -2.0; // relative stopping tolerance for interior solution
    trlib_flt_t tol_abs_i = 0.0; // absolute stopping tolerance for interior solution
    trlib_flt_t tol_rel_b = -3.0; // relative stopping tolerance for boundary solution
    trlib_flt_t tol_abs_b = 0.0; // absolute stopping tolerance for boundary solution
    trlib_flt_t obj_lo = -1e20; //  lower bound on objective, returns if a point is found with function value <= obj_lo
    trlib_int_t convexify = 1; // set to 1 if you like to monitor if the tridiagonal solution and the backtransformed solution match and if not resolve with a convexified problem
    trlib_int_t earlyterm = 1; // set to 1 if you like to terminate in the boundary case if it unlikely that much progress will be made fast but no convergence is reached, else set to 0

    trlib_int_t ret = 0; // return value

    // DOF
    trlib_int_t n = _c.size();

    trlib_int_t init = TRLIB_CLS_INIT, inc = 1, itp1 = 0;
    trlib_int_t iter;

    trlib_krylov_prepare_memory( _maxIterations, fwork );

    RealType g_dot_g = 0.0, v_dot_g = 0.0, p_dot_Hp = 0.0, flt1, flt2, flt3;
    trlib_int_t action, ityp;

    VectorType g( _c );
    VectorType &v = g; //! \todo Include preconditioner
    VectorType p( -g );
    VectorType gm( n );
    gm.setZero();

    VectorType Hp;
    VectorType temp( n );

    FullMatrixType Q( n, _maxIterations + 1 );
    Q.setZero();
    char prefix = 0;

    while ( true ) {
      ret = trlib_krylov_min( init, _radius, equality, _maxIterations, itmax_lanczos,
                              tol_rel_i, tol_abs_i, tol_rel_b, tol_abs_b,
                              TRLIB_EPS * TRLIB_EPS, obj_lo, ctl_invariant, convexify, earlyterm,
                              g_dot_g, v_dot_g, p_dot_Hp, iwork, fwork,
                              refine, verbose, unicode, &prefix, stdout, NULL,
                              &action, &iter, &ityp, &flt1, &flt2, &flt3 );
      this->status.Iteration = iter;
      init = 0;
      switch ( action ) {
        case TRLIB_CLA_INIT:
          // s = 0
          s.setZero();
          // g = g_0
          g = _c;
          // g_ = 0
          gm.setZero();
          // v = M^{-1} * g
          //! \todo Include preconditioner
          // p = -v
          p = -v;
          // Hp = H*p
          Hp = _H * p;
          // g_dot_g = <g,g>
          g_dot_g = g.squaredNorm();
          // v_dot_g = <v,g>
          v_dot_g = v.dot( g );
          // p_dot_Hp = <p, Hp>
          p_dot_Hp = p.dot( Hp );
          // q_0 = 1/sqrt(v_dot_g) * v
          Q.col( 0 ) = v / std::sqrt( v_dot_g );
          break;
        case TRLIB_CLA_RETRANSF: {
          Eigen::Map<VectorType> h_i( fwork + h_pointer, iter + 1 );
          // s = Q_i * h_i
          s = Q.block( 0, 0, n, iter + 1 ) * h_i;

          break;
        }
        case TRLIB_CLA_UPDATE_STATIO:
          if ( ityp == TRLIB_CLT_CG) {
            s += flt1 * p;
          }
          break;
        case TRLIB_CLA_UPDATE_GRAD:
          if ( ityp == TRLIB_CLT_CG) {
            // q_i = flt2 * v
            Q.col( iter ) = flt2 * v;
            // g_ = g
            gm = g;
            // g = g+flt1* Hp
            g += flt1 * Hp;
          }
          if ( ityp == TRLIB_CLT_L) {
            // s = Hp + flt1 *g + flt2*g_
            s = Hp + flt1 * g + flt2 * gm;
            // g_ = flt3 * g
            gm = flt3 * g;
            // g = s
            g = s;
          }
          // v = M^{-1} * g
          //! \todo Include preconditioner
          // g_dot_g = <g,g>
          g_dot_g = g.squaredNorm();
          // v_dot_g = <v,g>
          v_dot_g = v.dot( g );
          break;
        case TRLIB_CLA_UPDATE_DIR:
          if ( ityp == TRLIB_CLT_CG) {
            // p = - v + flt2 * p with flt1 = -1
            p = -v + flt2 * p;
          }
          if ( ityp == TRLIB_CLT_L) {
            // p = flt1 * v
            p = flt1 * v;
            // (optionally) M-orthonormalize against Q_{i-1}
          }

          // Hp = H*p
          Hp = _H * p;
          // p_dot_Hp = <p, Hp>
          p_dot_Hp = p.dot( Hp );

          if ( ityp == TRLIB_CLT_L) {
            // q_i = p
            Q.col( iter ) = p;
          }
          break;
        case TRLIB_CLA_CONV_HARD:
          // v_dot_g = <H*s + g_0 +flt1 * M*s, M^{-1}(H*s+g_0)+flt1 * s
          //! \todo Include preconditioner
          temp = _H * s;
          temp += _c + flt1 * s;
          v_dot_g = temp.squaredNorm();
          break;
        case TRLIB_CLA_NEW_KRYLOV:
          std::cerr << "TRLIB: Hit invariant Krylov subspace. Please implement proper reorthogonalization!"
                    << std::endl;
          break;
        case TRLIB_CLA_OBJVAL:
          // g_dot_g = 1/2 * <s,Hs> + <g,s>
          temp = 0.5 * _H * s;
          temp += _c;

          g_dot_g = s.dot( temp );
          break;
        default:
          break;

      }

      if ( ret < 10 ) { break; }
    }

    this->status.reasonOfTermination = ret;
    this->status.Residual = v_dot_g;

    delete[] iwork;
    delete[] fwork;

  }

  void setParameter( const std::string &name, std::string value ) override {
    throw std::runtime_error( "trlibSolver::setParameter(): This class has no string parameters." );
  }

  void setParameter( const std::string &name, RealType value ) override {
    if ( name == "tolerance" )
      _tolerance = value;
    else
      throw std::runtime_error( "trlibSolver::setParameter(): Unknown parameter '" + name + "'." );
  }

  void setParameter( const std::string &name, int value ) override {
    if ( name == "maximum_iterations" )
      _maxIterations = value;
    else if ( name == "print_level" )
      _quiet = (value < 5);
    else
      throw std::runtime_error( "trlibSolver::setParameter(): Unknown parameter '" + name + "'." );
  }
};

#endif // GOAST_WITH_TRLIB

#endif //OPTIMIZATION_TRLIBINTERFACE_H
