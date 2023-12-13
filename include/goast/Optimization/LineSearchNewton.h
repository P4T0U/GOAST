// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Header for instances of the unconstrained trust-region Newton method
 * \author Sassen
 */

#ifndef OPTIMIZATION_LINESEARCHNEWTON_H
#define OPTIMIZATION_LINESEARCHNEWTON_H

#include "optInterface.h"
#include "optParameters.h"
#include "optUtils.h"

#include "interfaces/CholmodInterface.h"

template<typename ConfiguratorType>
class LineSearchNewton : public OptimizationBase<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::TensorType TensorType;

  // Epsilon for safe-guarding computation of reduction against numerically instabilities
  const RealType eps_red = 100 * std::numeric_limits<RealType>::epsilon();

  // Functionals
  const BaseOp<VectorType, RealType> &_F;
  const BaseOp<VectorType, VectorType> &_DF;
  const BaseOp<VectorType, MatrixType> &_D2F;


  // Linear solver
  mutable Eigen::CholmodSupernodalLLT<MatrixType> _cholmodSolver;
  RealType _dir_beta = 1e-3;
  RealType _tau_factor = 2.;

  // Stopping criteria
  int _maxIterations;
  RealType _stopEpsilon;

  // Line search
  TIMESTEP_CONTROLLER _linesearchMethod = ARMIJO;
  RealType _sigma = 0.1;
  RealType _beta = 0.9;
  RealType _start_tau = 1.0;
  RealType _tau_min = 1e-12;
  RealType _tau_max = 4.;

  bool _reducedDirection = false;

  const std::vector<int> *_bdryMask;

  QUIET_MODE _quietMode;
  
public:
  LineSearchNewton( const BaseOp<VectorType, RealType> &F,
                    const BaseOp<VectorType, VectorType> &DF,
                    const BaseOp<VectorType, MatrixType> &D2F,
                    const OptimizationParameters<ConfiguratorType> &optPars ) 
  : _F( F ), _DF( DF ), _D2F( D2F ), 
    _maxIterations( optPars.getNewtonIterations() ),
    _stopEpsilon( optPars.getStopEpsilon() ), 
    _bdryMask( nullptr ), 
    _quietMode( optPars.getQuietMode() ) { }
                                           
  LineSearchNewton( const BaseOp<VectorType, RealType> &F,
                    const BaseOp<VectorType, VectorType> &DF,
                    const BaseOp<VectorType, MatrixType> &D2F,
                    const RealType stopEpsilon = 1e-8,
                    const int maxIterations = 100,
                    QUIET_MODE quietMode = SUPERQUIET ) : _F( F ), _DF( DF ), _D2F( D2F ), _maxIterations( maxIterations ),
                                           _stopEpsilon( stopEpsilon ), _bdryMask( nullptr ), _quietMode( quietMode ) { }

  LineSearchNewton( const BaseOp<VectorType, RealType> &F,
                    const BaseOp<VectorType, VectorType> &DF,
                    const BaseOp<VectorType, MatrixType> &D2F,
                    const RealType stopEpsilon = 1e-8,
                    const int maxIterations = 100,
                    bool quiet = false ) : _F( F ), _DF( DF ), _D2F( D2F ), _maxIterations( maxIterations ),
                                           _stopEpsilon( stopEpsilon ), _bdryMask( nullptr ),
                                           _quietMode( quiet ? SUPERQUIET  : SHOW_ALL ) { }

  void setBoundaryMask( const std::vector<int> &Mask ) {
    _bdryMask = &Mask;
  }


  void solve( const VectorType &x_0, VectorType &x_k ) const override {
    this->status.Iteration = 0;
    this->status.totalTime = 0.;
    this->status.additionalTimings["Evaluation"] = 0.;
    this->status.additionalTimings["Direction"] = 0.;
    this->status.additionalTimings["LineSearch"] = 0.;

    auto t_start_setup = std::chrono::high_resolution_clock::now();

    int numDofs = x_0.size();
    int numDofs_red = x_0.size();

    x_k = x_0;

    VectorType p_k( x_0.size());
    p_k.setZero();

    RealType alpha_k  = _start_tau;
    RealType tau_k = 0;    

    VectorType tmp_x_k( x_0.size());
    RealType F_k, tmp_F_k;

    VectorType grad_F_k;
    MatrixType Hess_F_k;

    // initial evaluation of F, DF and D^2F
    auto t_start_eval = std::chrono::high_resolution_clock::now();
    _F.apply( x_k, F_k );
    // compute gradient and check norm
    _DF.apply( x_k, grad_F_k );
    if ( _bdryMask )
      applyMaskToVector<VectorType>( *_bdryMask, grad_F_k );

    RealType gradNorm = grad_F_k.norm();
    if ( gradNorm < _stopEpsilon ) {
      auto t_end_eval = std::chrono::high_resolution_clock::now();
      this->status.additionalTimings["Evaluation"] += std::chrono::duration<RealType, std::milli>(
              t_end_eval - t_start_eval ).count();

      if ( _quietMode == SHOW_ALL || _quietMode == SHOW_TERMINATION_INFO  )
        std::cout << " -- LSN -- Initial gradient norm below epsilon." << std::endl;

      return;
    }
    // compute Hessian
    _D2F.apply( x_k, Hess_F_k );
    if ( _bdryMask )
      applyMaskToSymmetricMatrix<MatrixType>( *_bdryMask, Hess_F_k );

    auto t_end_eval = std::chrono::high_resolution_clock::now();
    this->status.additionalTimings["Evaluation"] += std::chrono::duration<RealType, std::milli>(t_end_eval - t_start_eval ).count();

    // Sparse identity matrix, because modifying the diagonal of a sparse matrix is not as elegant otherwise
//    MatrixType I( n, n );
//    for ( int i = 0; i < n; i++ ) {
//      if ( _bdryMask )
//        if ( std::find( _bdryMask->begin(), _bdryMask->end(), i ) != _bdryMask->end())
//          continue;
//      I.insert( i, i ) = 1.;
//    }
//    I.makeCompressed();

    // permute vertex indices when using reduced directions (here the boundary indices are moved to the end)
    VectorType diagonal( x_0.size());
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, typename MatrixType::StorageIndex> bndPermutation(numDofs);
    bndPermutation.setIdentity();
    if (_reducedDirection &&  _bdryMask ) {
      numDofs_red -= _bdryMask->size();

      int nonBndIdx = 0;
      int bndIdx = 0;
      for ( int i = 0; i < numDofs; i++ ) {
        if ( std::find( _bdryMask->begin(), _bdryMask->end(), i ) != _bdryMask->end()) {
          bndPermutation.indices()[numDofs_red + bndIdx] = i;
          bndIdx++;
        }
        else {
          bndPermutation.indices()[nonBndIdx] = i;
          nonBndIdx++;
        }
      }
    }

    VectorType grad_F_k_red(numDofs_red);
    MatrixType Hess_F_k_red(numDofs_red, numDofs_red);
    VectorType p_k_red(numDofs_red);

    if ( _reducedDirection ) {
      Hess_F_k_red = (bndPermutation.inverse() * Hess_F_k * bndPermutation).block( 0, 0, numDofs_red, numDofs_red );
//      grad_F_k_red = (bndPermutation.inverse() * grad_F_k).segment( 0, numDofs_red );
      diagonal.resize(numDofs_red);
    }

    // Modified Hessian
//    MatrixType H_lambda( Hess_F_k + I );

    StepsizeControl<ConfiguratorType> stepsizeControl( _F, _DF, _linesearchMethod, _sigma, _beta, _start_tau, _tau_min, _tau_max );

    auto t_end_setup = std::chrono::high_resolution_clock::now();
    this->status.additionalTimings["Setup"] = std::chrono::duration<RealType, std::milli>( t_end_setup - t_start_setup ).count();

    // initialize cholmod
    auto t_start = std::chrono::high_resolution_clock::now();

    if (_reducedDirection ) {
      _cholmodSolver.analyzePattern( Hess_F_k_red );
    }
    else {
      _cholmodSolver.analyzePattern( Hess_F_k );
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    this->status.additionalTimings["Direction"] += std::chrono::duration<RealType, std::milli>( t_end - t_start ).count();

    _cholmodSolver.cholmod().print = 0;

    // print header for console output
    if ( _quietMode == SHOW_ALL )
        printHeader();
    if ( _quietMode == SHOW_ALL )
      printConsoleOutput<std::chrono::duration<RealType, std::milli>>(0, F_k, grad_F_k.norm(), 0, 0, 0, t_end_setup - t_start_setup, t_end - t_start, std::chrono::duration<RealType, std::milli>(0 ), t_end_eval - t_start_eval);

    // start Newton iteration
    for ( int k = 1; k <= _maxIterations; k++ ) {
      this->status.Iteration = k;

      auto t_start = std::chrono::high_resolution_clock::now();

      // Step 1: Compute descent direction
      auto t_start_dir = std::chrono::high_resolution_clock::now();

      int j = 0;
      if (_reducedDirection ) {
        Hess_F_k_red = (bndPermutation.inverse() * Hess_F_k * bndPermutation).block(0,0,numDofs_red,numDofs_red);
        grad_F_k_red = (bndPermutation.inverse() * grad_F_k).segment(0, numDofs_red);


        diagonal = Hess_F_k_red.diagonal();
        RealType min_diag = diagonal.minCoeff();

        if ( min_diag > -eps_red )
          tau_k = 0;
        else
          tau_k = -min_diag + _dir_beta;

        while ( true ) {
          _cholmodSolver.setShift( tau_k );
          _cholmodSolver.factorize( Hess_F_k_red );

          if ( _cholmodSolver.info() == Eigen::Success )
            break;

          tau_k = std::max( _tau_factor * tau_k, _dir_beta );
          j++;
        }

        // Actually compute reduced descent direction
        p_k_red = -_cholmodSolver.solve( grad_F_k_red );

        p_k.segment( 0, numDofs_red ) = p_k_red;
        p_k.segment( numDofs_red, numDofs - numDofs_red ).setZero();

        p_k = bndPermutation * p_k;
      }
      else {
        diagonal = Hess_F_k.diagonal();
        RealType min_diag = diagonal.minCoeff();

        if ( min_diag > -eps_red )
          tau_k = 0;
        else
          tau_k = -min_diag + _dir_beta; // std::max( -min_diag + _dir_beta, tau_k / _tau_factor ); //

        while ( true ) {
//        H_lambda = Hess_F_k + tau_k * I;
//        _cholmodSolver.factorize( H_lambda );

          _cholmodSolver.setShift( tau_k );
          _cholmodSolver.factorize( Hess_F_k );

          if ( _cholmodSolver.info() == Eigen::Success )
            break;

          tau_k = std::max( _tau_factor * tau_k, _dir_beta );
          j++;
        }

        // Actually compute descent direction
        p_k = -_cholmodSolver.solve( grad_F_k );
      }

//      std::cout << "pk Norm: " << std::scientific << p_k.norm() << std::endl;
//      std::cout << "pk_red Norm: " << (bndPermutation.inverse() * p_k).segment(0, numDofs_red).norm() << std::endl;
//      std::cout << "pk_nonred Norm: " << (bndPermutation.inverse() * p_k).segment(numDofs_red, numDofs - numDofs_red).norm() << std::endl;

      auto t_end_dir = std::chrono::high_resolution_clock::now();
      this->status.additionalTimings["Direction"] += std::chrono::duration<RealType, std::milli>(t_end_dir - t_start_dir ).count();

      // Step 2: Line search
      auto t_start_ls = std::chrono::high_resolution_clock::now();
      alpha_k = stepsizeControl.getStepsize( x_k, grad_F_k, p_k, alpha_k, F_k );
      auto t_end_ls = std::chrono::high_resolution_clock::now();
      this->status.additionalTimings["LineSearch"] += std::chrono::duration<RealType, std::milli>(t_end_ls - t_start_ls ).count();

      // step size zero -> TERMINATE!
      if ( alpha_k == 0. ) {
        auto t_end = std::chrono::high_resolution_clock::now();
        this->status.totalTime += std::chrono::duration<RealType, std::milli>( t_end - t_start ).count();

        if ( _quietMode == SHOW_ALL ) {
            printConsoleOutput(k, F_k, grad_F_k.norm(), alpha_k, tau_k, j, t_end - t_start, t_end_dir - t_start_dir, t_end_ls - t_start_ls, t_end_eval - t_start_eval);
            std::cout << " -- LSN -- Iter " << std::setw(3) << k << ": " << "Step size too small." << std::endl;
        }
        break;
      }

      // actual update of position
      x_k += alpha_k * p_k;

      // start new evaluation
      t_start_eval = std::chrono::high_resolution_clock::now();
      _F.apply( x_k, F_k );
      _DF.apply( x_k, grad_F_k );
      _D2F.apply( x_k, Hess_F_k );
      if ( _bdryMask )
        applyMaskToSymmetricMatrixAndVector<MatrixType, VectorType>( *_bdryMask, Hess_F_k, grad_F_k );
      gradNorm = grad_F_k.norm();
      t_end_eval = std::chrono::high_resolution_clock::now();
      this->status.additionalTimings["Evaluation"] += std::chrono::duration<RealType, std::milli>(t_end_eval - t_start_eval ).count();

      // stop total timing
      auto t_end = std::chrono::high_resolution_clock::now();
      this->status.totalTime += std::chrono::duration<RealType, std::milli>( t_end - t_start ).count();

      // gradient norm small enough -> TERMINATE!
      if ( gradNorm < _stopEpsilon ) {
        if ( _quietMode == SHOW_ALL ) {
            printConsoleOutput( k, F_k, gradNorm, alpha_k, tau_k, j, t_end - t_start, t_end_dir - t_start_dir, t_end_ls - t_start_ls, t_end_eval - t_start_eval);
            std::cout << " -- LSN -- Iter " << std::setw(3) << k << ": " << "Gradient norm below epsilon." << std::endl;
        }
        break;
      }

      // console output for that iteration
      if ( _quietMode == SHOW_ALL )
          printConsoleOutput( k, F_k, gradNorm, alpha_k, tau_k, j, t_end - t_start, t_end_dir - t_start_dir, t_end_ls - t_start_ls, t_end_eval - t_start_eval);
                  
    } // end of Newton iteration

    // final console output
    if ( (_quietMode == SHOW_ALL) || (_quietMode == SHOW_TERMINATION_INFO) || ( (_quietMode == SHOW_ONLY_IF_FAILED) && (gradNorm > _stopEpsilon) ) ){
      std::cout << " -- LSN -- Final   : " << std::scientific << std::setprecision( 6 )
                << F_k
                << " || " << gradNorm
                << " ||===|| " << std::fixed << std::setprecision( 2 ) << this->status.totalTime
                << " || " << this->status.additionalTimings["Direction"]
                << " || " << this->status.additionalTimings["LineSearch"]
                << " || " << this->status.additionalTimings["Evaluation"]
                << std::endl;
    }
  }

  void setParameter( const std::string &name, std::string value ) override {
    throw std::runtime_error( "LineSearchNewton::setParameter(): Unknown string parameter '" + name + "'." );
  }

  void setParameter( const std::string &name, RealType value ) override {
    if ( name == "sigma" )
      _sigma = value;
    else if ( name == "beta" )
      _beta = value;
    else if ( name == "initial_stepsize" )
      _start_tau = value;
    else if ( name == "minimal_stepsize" )
      _tau_min = value;
    else if ( name == "maximal_stepsize" )
      _tau_max = value;
    else if ( name == "direction_beta" )
      _dir_beta = value;
    else if ( name == "tau_increase" )
      _tau_factor = value;
    else if ( name == "tolerance" )
      _stopEpsilon = value;
    else
      throw std::runtime_error( "LineSearchNewton::setParameter(): Unknown real parameter '" + name + "'." );
  }

  void setParameter( const std::string &name, int value ) override {
    if ( name == "maximum_iterations" )
      _maxIterations = value;
    else if ( name == "stepsize_control" )
      _linesearchMethod = static_cast<TIMESTEP_CONTROLLER> (value);
    else if ( name == "print_level" )
      _quietMode =  static_cast<QUIET_MODE> (value);
    else if ( name == "reduced_direction" )
      _reducedDirection = static_cast<bool> (value);
    else if ( name == "gradient_iterations" )
      return;
    else if ( name == "BFGS_iterations" )
      return;
    else
      throw std::runtime_error( "LineSearchNewton::setParameter(): Unknown integer parameter '" + name + "'." );
  }

protected:
    void printHeader() const {
            std::cout << " -- LSN --         : " << std::scientific << std::setprecision( 6 )
                      << "  F[x_k]  "  << " || " << "  |DF[x_k]|  " << " ||===|| "
                      << " stepsize " << " || " << "  shift  " << " || " << "num LS" << " ||===|| "
                      << "total time" << " || " << "time dir" << " || " << "time ls" << " || " << "time eval"
                      << std::endl;
    }

    template<typename IntType>
    void printConsoleOutput( int iteration,
                             RealType F,
                             RealType NormDF,
                             RealType stepsize,
                             RealType shift,
                             int numLSEvals,
                             IntType DeltaTotalTime,
                             IntType DeltaDirTime,
                             IntType DeltaLSTime,
                             IntType DeltaEvalTime  ) const {
            std::cout << " -- LSN -- Iter " << std::setw( 3 ) << iteration << ": " << std::scientific << std::setprecision( 6 )
                      << F
                      << " || " << NormDF
                      << " ||===|| " << std::setw( 13 ) << stepsize
                      << " || " << shift
                      << " || " << numLSEvals
                      << std::setprecision( 2 ) << std::fixed
                      << " ||===|| " << std::setw( 6 ) << std::chrono::duration<RealType, std::milli>( DeltaTotalTime ).count()
                      << " || " << std::setw( 6 ) << std::chrono::duration<RealType, std::milli>( DeltaDirTime ).count()
                      << " || " << std::setw( 6 ) << std::chrono::duration<RealType, std::milli>(  DeltaLSTime ).count()
                      << " || " << std::setw( 6 ) << std::chrono::duration<RealType, std::milli>( DeltaEvalTime ).count()
                      << std::endl;
  }

};

#endif //OPTIMIZATION_LINESEARCHNEWTON_H
