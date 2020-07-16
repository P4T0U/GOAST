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

#ifndef OPTIMIZATION_TRUSTREGIONNEWTON_H
#define OPTIMIZATION_TRUSTREGIONNEWTON_H

#include "optInterface.h"
#include "optUtils.h"
#include "SteihaugCG.h"
#include "MoreSorensen.h"
#include "interfaces/trlibInterface.h"

template<typename ConfiguratorType>
class TrustRegionNewton : public OptimizationBase<ConfiguratorType> {
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

  RealType _maxRadius;
  RealType _minRadius = 1e-8;
  RealType _initRadius;
  RealType _eta;
  RealType _minStepsize = 1e-9;
  int _maxIterations;
  int _cgIterations;

  RealType _stopEpsilon;

  bool _diagonalPreconditioning = true;

  bool _quiet;

  std::string _subproblemSolver = "SteihaugCG";

  const std::vector<int> *_bdryMask;

  // Trust-region subproblem solver parameters
  std::map<std::string, int> _trsolverIntParameters;
  std::map<std::string, RealType> _trsolverRealParameters;
  std::map<std::string, std::string> _trsolverStringParameters;


public:
  TrustRegionNewton( const BaseOp<VectorType, RealType> &F,
                     const BaseOp<VectorType, VectorType> &DF,
                     const BaseOp<VectorType, MatrixType> &D2F,
                     const RealType initRadius,
                     const RealType maxRadius,
                     const RealType stopEpsilon = 1e-8,
                     const int maxIterations = 100,
                     const int cgIterations = 100,
                     const RealType eta = 0.25,
                     bool quiet = false ) : _F( F ), _DF( DF ), _D2F( D2F ), _maxRadius( maxRadius ),
                                            _initRadius( initRadius ), _eta( eta ),
                                            _maxIterations( maxIterations ), _cgIterations( cgIterations ),
                                            _stopEpsilon( stopEpsilon ), _bdryMask( nullptr ), _quiet( quiet ) {}

  void setBoundaryMask( const std::vector<int> &Mask ) {
    _bdryMask = &Mask;
  }


  void solve( const VectorType &x_0, VectorType &x_k ) const override {
    this->status.Iteration = 0;
    this->status.totalTime = 0.;
    this->status.additionalTimings["Preconditioner"] = 0.;
    this->status.additionalTimings["Subproblem"] = 0.;
    this->status.additionalTimings["Evaluation"] = 0.;

    RealType trRadius = _initRadius;
    RealType eta_k, eps_k, rho_k;

    int n = x_0.size();

    x_k = x_0;

    VectorType p_k( x_0.size());
    p_k.setZero();

    VectorType tmp_x_k( x_0.size());
    RealType F_k, tmp_F_k;
    RealType m_red, f_red; // predicted and actual reduction

    VectorType grad_F_k;
    MatrixType Hess_F_k;
    _F.apply( x_k, F_k );
    _DF.apply( x_k, grad_F_k );
    _D2F.apply( x_k, Hess_F_k );
    if ( _bdryMask )
      applyMaskToSymmetricMatrixAndVector<MatrixType, VectorType>( *_bdryMask, Hess_F_k, grad_F_k );

    RealType initGradNorm = grad_F_k.norm();

    VectorType diagonal( x_0.size());
    MatrixType Dinv( x_0.size(), x_0.size());
    for ( int i = 0; i < n; i++ ) {
      Dinv.coeffRef( i, i ) = 1;
    }

    VectorType c( n );
    MatrixType H( n, n );

    auto t_start_eval = std::chrono::high_resolution_clock::now();
    auto t_end_eval = std::chrono::high_resolution_clock::now();


    for ( int k = 0; k < _maxIterations; k++ ) {
      // Step 1: Solve trust-region subproblem
      this->status.Iteration = k;

      auto t_start = std::chrono::high_resolution_clock::now();

      // Diagonal Preconditioning
      auto t_start_pre = std::chrono::high_resolution_clock::now();

      if ( _diagonalPreconditioning ) {
        diagonal = Hess_F_k.diagonal();
        for ( int i = 0; i < n; i++ ) {
          if ( std::abs( diagonal[i] ) > eps_red )
            Dinv.coeffRef( i, i ) = 1 / diagonal[i];
          else
            Dinv.coeffRef( i, i ) = 0;
        }
      }

      c.noalias() = Dinv * grad_F_k;
      H = Dinv * Hess_F_k * Dinv;

      auto t_end_pre = std::chrono::high_resolution_clock::now();
      this->status.additionalTimings["Preconditioner"] += std::chrono::duration<RealType, std::milli>(
              t_end_pre - t_start_pre ).count();

      // Compute forcing sequence / epsilon
      eta_k = std::min( 0.5, std::sqrt( c.norm()));
      eps_k = eta_k * c.norm();

      // Apply trust-region subproblem solver
      auto t_start_inner = std::chrono::high_resolution_clock::now();
      SolverStatus<ConfiguratorType> trsolverStatus = solveTrustRegionSubproblem( H, c, trRadius, eps_k, p_k );
      auto t_end_inner = std::chrono::high_resolution_clock::now();
      this->status.additionalTimings["Subproblem"] += std::chrono::duration<RealType, std::milli>(
              t_end_inner - t_start_inner ).count();

      RealType pkn = p_k.norm();

      p_k = Dinv * p_k;

      // Step 2: Determine reduction ration
      tmp_x_k = x_k + p_k; // temporary new iterate
      _F.apply( tmp_x_k, tmp_F_k ); // temporary new function value

      m_red = -grad_F_k.dot( p_k ) - 0.5 * p_k.dot( Hess_F_k * p_k ); // predicted reduction i.e. in the quadratic model
      f_red = (F_k - tmp_F_k);

      if ((std::abs( f_red ) < eps_red && std::abs( m_red ) < eps_red) || std::abs( f_red - m_red ) < eps_red ) {
        if ( !_quiet )
          std::cout << " -- TRN -- Iter " << std::setw( 3 ) << k << ": " << "Cutoff active in rho-computation"
                    << std::endl;
        rho_k = 1.;
      }
      else
        rho_k = f_red / m_red; // actual over predicted reduction

      // Step 3: Update trust region radius
      if ( rho_k < 0.25 )
        trRadius = trRadius / 4.;
      else if ( rho_k > 0.75 && std::abs( pkn - trRadius ) <= eps_red )
        trRadius = std::min( 2 * trRadius, _maxRadius );

      // Step 4: Accept or decline new iterate
      if ( rho_k > _eta ) {
        x_k = tmp_x_k;
        F_k = tmp_F_k;

        t_start_eval = std::chrono::high_resolution_clock::now();
        _DF.apply( x_k, grad_F_k );
        _D2F.apply( x_k, Hess_F_k );
        if ( _bdryMask )
          applyMaskToSymmetricMatrixAndVector<MatrixType, VectorType>( *_bdryMask, Hess_F_k, grad_F_k );
        t_end_eval = std::chrono::high_resolution_clock::now();
        this->status.additionalTimings["Evaluation"] += std::chrono::duration<RealType, std::milli>(
                t_end_eval - t_start_eval ).count();

        if ( grad_F_k.norm() < _stopEpsilon ) {
          if ( !_quiet )
            std::cout << " -- TRN -- Iter " << std::setw( 3 ) << k << ": " << "Gradient norm below epsilon."
                    << std::endl;
          auto t_end = std::chrono::high_resolution_clock::now();
          this->status.totalTime += std::chrono::duration<RealType, std::milli>( t_end - t_start ).count();
          break;
        }
      }

      if ( trRadius < _minRadius ) {
        if ( !_quiet )
          std::cout << " -- TRN -- Iter " << std::setw( 3 ) << k << ": " << "Trust region too small." << std::endl;
        auto t_end = std::chrono::high_resolution_clock::now();
        this->status.totalTime += std::chrono::duration<RealType, std::milli>( t_end - t_start ).count();
        break;
      }

      if ( p_k.template lpNorm<Eigen::Infinity>() < _minStepsize ) {
        if ( !_quiet )
          std::cout << " -- TRN -- Iter " << std::setw( 3 ) << k << ": " << "Step size too small." << std::endl;
        auto t_end = std::chrono::high_resolution_clock::now();
        this->status.totalTime += std::chrono::duration<RealType, std::milli>( t_end - t_start ).count();
        break;
      }


      auto t_end = std::chrono::high_resolution_clock::now();
      this->status.totalTime += std::chrono::duration<RealType, std::milli>( t_end - t_start ).count();

      if ( !_quiet )
        std::cout << " -- TRN -- Iter " << std::setw( 3 ) << k << ": " << std::scientific << std::setprecision(6)
                  << F_k
                  << " || " << grad_F_k.norm()
                  << " ||===|| " << std::setw( 13 ) <<  rho_k
                  << " || " << p_k.template lpNorm<Eigen::Infinity>()
                  << " || " << trRadius
                  << " ||===|| " << trsolverStatus.Residual
                  << " || " << std::fixed << std::setw( 2 ) << trsolverStatus.reasonOfTermination
                  << " || " << std::setw( 4 ) << trsolverStatus.Iteration
                  << std::setprecision(2)
                  << " ||===|| " << std::setw( 6 ) << std::chrono::duration<RealType, std::milli>( t_end - t_start ).count()
                  << " || " << std::setw( 6 ) << std::chrono::duration<RealType, std::milli>( t_end_inner - t_start_inner ).count()
                  << " || " << std::setw( 6 ) << std::chrono::duration<RealType, std::milli>( t_end_pre - t_start_pre ).count()
                  << " || " << std::setw( 6 ) << std::chrono::duration<RealType, std::milli>( t_end_eval - t_start_eval ).count()
                  << std::endl;
    }

    if ( !_quiet )
      std::cout << " -- TRN -- Final   : " << std::scientific << std::setprecision(6)
                << F_k
                << " || " << grad_F_k.norm()
                << " ||===|| " << std::setw( 13 ) << rho_k
                << " || " << p_k.template lpNorm<Eigen::Infinity>()
                << " || " << trRadius
                << " ||===|| " << std::fixed << std::setprecision(2) << this->status.totalTime
                << " || " << this->status.additionalTimings["Subproblem"]
                << std::endl;
  }

  SolverStatus<ConfiguratorType> solveTrustRegionSubproblem( const MatrixType &H, const VectorType &c,
                                                             const RealType &trRadius, const RealType &eps_k,
                                                             VectorType &p ) const {
    if ( _subproblemSolver == "SteihaugCG" ) {
      SteihaugCGMethod<ConfiguratorType> trSolver( H, c, trRadius, eps_k, _cgIterations, true );
      trSolver.setParameters( _trsolverIntParameters );
      trSolver.setParameters( _trsolverRealParameters );
      trSolver.setParameters( _trsolverStringParameters );
      trSolver.setParameter( "tolerance", eps_k ); // Just to be safe

      trSolver.solve( p );

      return trSolver.status;
    }
    else if ( _subproblemSolver == "MoreSorensen" ) {
      MoreSorensenMethod<ConfiguratorType> trSolver( H, c, trRadius, eps_k, 100, 0.01, 0.1, 0.2, true );
      trSolver.setParameters( _trsolverIntParameters );
      trSolver.setParameters( _trsolverRealParameters );
      trSolver.setParameters( _trsolverStringParameters );
      trSolver.setParameter( "tolerance", eps_k ); // Just to be safe

      trSolver.setBoundaryMask( *_bdryMask );

      trSolver.solve( p );

      return trSolver.status;
    }
#ifdef GOAST_WITH_TRLIB
    else if ( _subproblemSolver == "trlib" ) {
      trlibSolver<ConfiguratorType> trSolver( H, c, trRadius, eps_k, _cgIterations, true );
      trSolver.setParameters( _trsolverIntParameters );
      trSolver.setParameters( _trsolverRealParameters );
      trSolver.setParameters( _trsolverStringParameters );
      trSolver.setParameter( "tolerance", eps_k ); // Just to be safe

      trSolver.solve( p );

      return trSolver.status;
    }
#endif
    else
      throw std::runtime_error( "TrustRegionNewton::solveTrustRegionSubproblem(): Unknown subproblem method!" );
  }

  void setParameter( const std::string &name, std::string value ) override {
    if ( name == "subproblem_solver" )
      _subproblemSolver = value;
    else if ( name.rfind( "trsolver__", 0 ) == 0 )
      _trsolverStringParameters[name.substr( 10, std::string::npos )] = value;
    else
      throw std::runtime_error( "TrustRegionNewton::setParameter(): Unknown parameter '" + name + "'." );
  }

  void setParameter( const std::string &name, RealType value ) override {
    if ( name == "maximal_radius" )
      _maxRadius = value;
    else if ( name == "initial_radius" )
      _initRadius = value;
    else if ( name == "minimal_radius" )
      _minRadius = value;
    else if ( name == "minimal_stepsize" )
      _minStepsize = value;
    else if ( name == "accept_reduction_ratio" )
      _eta = value;
    else if ( name == "tolerance" )
      _stopEpsilon = value;
    else if ( name.rfind( "trsolver__", 0 ) == 0 )
      _trsolverRealParameters[name.substr( 10, std::string::npos )] = value;
    else
      throw std::runtime_error( "TrustRegionNewton::setParameter(): Unknown parameter '" + name + "'." );
  }

  void setParameter( const std::string &name, int value ) override {
    if ( name == "maximum_iterations" )
      _maxIterations = value;
    else if ( name == "cg_iterations" )
      _cgIterations = value;
    else if ( name == "print_level" )
      _quiet = (value < 5);
    else if ( name == "diagonal_preconditioning" )
      _diagonalPreconditioning = static_cast<bool>(value);
    else if ( name.rfind( "trsolver__", 0 ) == 0 )
      _trsolverIntParameters[name.substr( 10, std::string::npos )] = value;
    else
      throw std::runtime_error( "TrustRegionNewton::setParameter(): Unknown parameter '" + name + "'." );
  }

};

#endif //OPTIMIZATION_TRUSTREGIONNEWTON_H
