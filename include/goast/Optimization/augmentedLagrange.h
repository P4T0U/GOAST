// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by sassen on 02.05.19.
//

#ifndef OPTIMIZATION_AUGMENTEDLAGRANGE_H
#define OPTIMIZATION_AUGMENTEDLAGRANGE_H

#include <iostream>
#include <chrono>
#include <ctime>
#include <string>

#include <goast/Core/Auxiliary.h>
#include <goast/Core/BaseOpInterface.h>

#include <goast/external/ipoptBoxConstraintSolver.h>

#include "TrustRegionNewton.h"
#include "LineSearchNewton.h"


template<typename ConfiguratorType>
class AugmentedLagrangeFunctional
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {

  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;

  const RealType &_alpha;
  const RealType &_mu;
  const VectorType &_lambda;
  const BaseOp<VectorType, RealType> &_F;
  const BaseOp<VectorType, VectorType> &_G;


public:
  AugmentedLagrangeFunctional( const BaseOp<VectorType, RealType> &F,
                               const BaseOp<VectorType, VectorType> &G,
                               const VectorType &lambda,
                               const RealType &mu,
                               const RealType &alpha )
          : _F( F ), _G( G ), _lambda( lambda ), _mu( mu ), _alpha( alpha ) {
    if ( _lambda.size() != _G.getTargetDimension())
      throw std::length_error( "AugmentedLagrangeFunctional: Wrong number of lagrange multipliers!" );
  }


  void apply( const VectorType &Arg, RealType &Dest ) const override {
    _F.apply( Arg, Dest );

    Dest *= _alpha;

    VectorType linearPart;
    _G.apply( Arg, linearPart );

    Dest -= _lambda.transpose() * linearPart;
    Dest += _mu / 2. * linearPart.squaredNorm();

    if ( std::isnan( Dest ))
      Dest = std::numeric_limits<RealType>::infinity();
  }

};

template<typename ConfiguratorType>
class AugmentedLagrangeDerivative : public BaseOp<typename ConfiguratorType::VectorType> {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;

  const BaseOp<VectorType, RealType> &_F;
  const BaseOp<VectorType, VectorType> &_DF;
  const BaseOp<VectorType, VectorType> &_G;
  const BaseOp<VectorType, MatrixType> &_DG;
  const RealType &_alpha;
  const RealType &_mu;
  const VectorType &_lambda;

public:
  AugmentedLagrangeDerivative( const BaseOp<VectorType, RealType> &F,
                               const BaseOp<VectorType, VectorType> &DF,
                               const BaseOp<VectorType, VectorType> &G,
                               const BaseOp<VectorType, MatrixType> &DG,
                               const VectorType &lambda,
                               const RealType &mu,
                               const RealType &alpha )
          : _F( F ), _DF( DF ), _G( G ), _DG( DG ), _lambda( lambda ), _mu( mu ), _alpha( alpha ) {
    if ( _lambda.size() != _G.getTargetDimension())
      throw std::length_error( "AugmentedLagrangeDerivative: Wrong number of lagrange multipliers!" );
  }

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    _DF.apply( Arg, Dest );

    Dest.array() *= _alpha;

    VectorType linearPart;
    _G.apply( Arg, linearPart );

    MatrixType JacobianOfLinearPart;
    _DG.apply( Arg, JacobianOfLinearPart );

    if ( _G.getTargetDimension() != JacobianOfLinearPart.rows())
      throw std::length_error( "AugmentedLagrangeDerivative: Wrong number of rows!" );

    for ( int i = 0; i < JacobianOfLinearPart.rows(); i++ )
      Dest -= _lambda[i] * JacobianOfLinearPart.row( i );

    Dest += _mu * JacobianOfLinearPart.transpose() * linearPart;
  }

};


template<typename ConfiguratorType>
class AugmentedLagrangeHessian
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;
  typedef typename ConfiguratorType::TensorType TensorType;

  const BaseOp<VectorType, RealType> &_F;
  const BaseOp<VectorType, VectorType> &_DF;
  const BaseOp<VectorType, MatrixType> &_D2F;

  const BaseOp<VectorType, VectorType> &_G;
  const BaseOp<VectorType, MatrixType> &_DG;
  const BaseOp<VectorType, TensorType> &_D2G;

  const RealType &_alpha;
  const RealType &_mu;
  const VectorType &_lambda;

public:
  mutable std::map<std::string, RealType> timings;

  AugmentedLagrangeHessian( const BaseOp<VectorType, RealType> &F,
                            const BaseOp<VectorType, VectorType> &DF,
                            const BaseOp<VectorType, MatrixType> &D2F,
                            const BaseOp<VectorType, VectorType> &G,
                            const BaseOp<VectorType, MatrixType> &DG,
                            const BaseOp<VectorType, TensorType> &D2G,
                            const VectorType &lambda,
                            const RealType &mu,
                            const RealType &alpha )
          : _F( F ), _DF( DF ), _D2F( D2F ), _G( G ), _DG( DG ), _D2G( D2G ), _lambda( lambda ), _mu( mu ),
            _alpha( alpha ) {
    if ( _lambda.size() != _G.getTargetDimension())
      throw std::length_error( "AugmentedLagrangeHessian: Wrong number of lagrange multipliers!" );
  }

  void apply( const VectorType &Arg, MatrixType &Dest ) const override {
    timings.clear();

    // Overall triplet list
    TripletListType tripletList, energyTripletList;

    // Energy contribution
    auto t_start = std::chrono::high_resolution_clock::now();

    _D2F.pushTriplets( Arg, energyTripletList );

    tripletList.reserve( energyTripletList.size());
    for ( const auto &trip : energyTripletList )
      tripletList.emplace_back( trip.row(), trip.col(), _alpha * trip.value());

    auto t_end = std::chrono::high_resolution_clock::now();

    timings["EnergyTriplets"] = std::chrono::duration<RealType, std::milli>( t_end - t_start ).count();

    // Lagrangian contribution
//    _D2G.pushTriplets( Arg, tripletList, -_lambda );

    // Penalty contribution
    t_start = std::chrono::high_resolution_clock::now();
    VectorType linearPart;
    _G.apply( Arg, linearPart );
    t_end = std::chrono::high_resolution_clock::now();

    timings["ConstraintEvaluation"] = std::chrono::duration<RealType, std::milli>( t_end - t_start ).count();

    t_start = std::chrono::high_resolution_clock::now();
    VectorType factor( -_lambda + _mu * linearPart );
    _D2G.pushTriplets( Arg, tripletList, factor );
    t_end = std::chrono::high_resolution_clock::now();

    timings["ConstraintHessian"] = std::chrono::duration<RealType, std::milli>( t_end - t_start ).count();

    // Assemble Dest
    t_start = std::chrono::high_resolution_clock::now();
    Dest.resize( Arg.size(), Arg.size());
    Dest.setZero();
    Dest.setFromTriplets( tripletList.begin(), tripletList.end());
    t_end = std::chrono::high_resolution_clock::now();

    timings["PartialAssemble"] = std::chrono::duration<RealType, std::milli>( t_end - t_start ).count();

    // DG(x) * DG(x)
    t_start = std::chrono::high_resolution_clock::now();
    MatrixType JacobianOfLinearPart;
    _DG.apply( Arg, JacobianOfLinearPart );
    t_end = std::chrono::high_resolution_clock::now();

    timings["ConstraintGradient"] = std::chrono::duration<RealType, std::milli>( t_end - t_start ).count();

    t_start = std::chrono::high_resolution_clock::now();
    MatrixType JacobianSquared( Arg.size(), Arg.size());
    JacobianSquared = _mu * (JacobianOfLinearPart.transpose() * JacobianOfLinearPart);
    t_end = std::chrono::high_resolution_clock::now();

    timings["SquaringJacobian"] = std::chrono::duration<RealType, std::milli>( t_end - t_start ).count();

    t_start = std::chrono::high_resolution_clock::now();
    Dest += JacobianSquared;
    t_end = std::chrono::high_resolution_clock::now();

    timings["Adding"] = std::chrono::duration<RealType, std::milli>( t_end - t_start ).count();

  }

  virtual void pushTriplets( const VectorType &Arg, TripletListType &Dest ) const override {
    TripletListType energyTripletList;

    // Energy contribution
    _D2F.pushTriplets( Arg, energyTripletList );

    Dest.reserve( Dest.size() + energyTripletList.size());
    for ( const auto &trip : energyTripletList )
      Dest.emplace_back( trip.row(), trip.col(), _alpha * trip.value());

    VectorType linearPart;
    _G.apply( Arg, linearPart );

    VectorType factor = -_lambda + _mu * linearPart;
    _D2G.pushTriplets( Arg, Dest, factor );

    // DG(x) * DG(x)
    MatrixType JacobianOfLinearPart;
    _DG.apply( Arg, JacobianOfLinearPart );

    MatrixType JacobianSquared( Arg.size(), Arg.size());
    JacobianSquared = _mu * (JacobianOfLinearPart.transpose() * JacobianOfLinearPart);

    Dest.reserve( Dest.size() + JacobianSquared.nonZeros());

    for ( int k = 0; k < JacobianSquared.outerSize(); ++k )
      for ( typename MatrixType::InnerIterator it( JacobianSquared, k ); it; ++it )
        Dest.emplace_back( it.row(), it.col(), it.value());

  }

};

template<typename ConfiguratorType>
class AugmentedLagrangeMethod : public OptimizationBase<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::TensorType TensorType;

  const BaseOp<VectorType, RealType> &_F;
  const BaseOp<VectorType, VectorType> &_DF;
  const BaseOp<VectorType, MatrixType> &_D2F;

  const BaseOp<VectorType, VectorType> &_G;
  const BaseOp<VectorType, MatrixType> &_DG;
  const BaseOp<VectorType, TensorType> &_D2G;

  int _maxIterations;
  bool _quiet;
  const std::vector<int> *_bdryMask;
  const VectorType *_lowerBounds;
  const VectorType *_upperBounds;

  // Augmented Lagrange parameters
  RealType _mu_0 = 10;

  RealType _eta_0 = 1 / std::pow( _mu_0, 0.1 );
  RealType _eta_star;

  RealType _tau_0 = 1 / _mu_0;
  RealType _tau_star;

  RealType _maxPenalty = 1e10;

  RealType _etaDecreaseExponent = 0.3;
  RealType _penaltyIncrease = 5.;

  bool _relativeStoppingCriterion = true;

  RealType _alpha = 1.;

  std::string _innerMethod = "LSN";

  // Inner solver parameters
  std::map<std::string, int> _innerIntParameters;
  std::map<std::string, RealType> _innerRealParameters;
  std::map<std::string, std::string> _innerStringParameters;

  // Lagrange multiplier
  mutable VectorType _Lambda;

public:
  AugmentedLagrangeMethod( const BaseOp<VectorType, RealType> &_F, const BaseOp<VectorType, VectorType> &_DF,
                           const BaseOp<VectorType, MatrixType> &_D2F, const BaseOp<VectorType, VectorType> &_G,
                           const BaseOp<VectorType, MatrixType> &_DG, const BaseOp<VectorType, TensorType> &_D2G,
                           int _maxIterations, RealType stopEpsilon, bool _quiet )
          : _F( _F ), _DF( _DF ), _D2F( _D2F ), _G( _G ), _DG( _DG ), _D2G( _D2G ), _maxIterations( _maxIterations ),
            _eta_star( stopEpsilon ), _tau_star( stopEpsilon ), _quiet( _quiet ),
            _bdryMask( nullptr ), _lowerBounds( nullptr ), _upperBounds( nullptr ),
            _Lambda( VectorType::Constant( _G.getTargetDimension(), 0. )) {

  }

  // set boundary mask
  void setBoundaryMask( const std::vector<int> &Mask ) {
    _bdryMask = &Mask;
  }

  void setVariableBounds( const VectorType &lowerBounds, const VectorType &upperBounds ) override {
    _lowerBounds = &lowerBounds;
    _upperBounds = &upperBounds;
  }

  VectorType getLambda() const {
    return _Lambda;
  }

  void setLambda( const VectorType &Lambda ) {
    if ( Lambda.size() != _G.getTargetDimension())
      throw std::length_error( "AugmentedLagrangeMethod::setLambda: Lagrange multipliers have the wrong size!" );
    _Lambda = Lambda;
  }

  void solve( const VectorType &x_0, VectorType &x_k ) const {
    this->status.Iteration = 0;
    this->status.totalTime = 0.;
    this->status.additionalTimings["InnerOptimization"] = 0.;
    this->status.additionalTimings["InnerPre"] = 0.;
    this->status.additionalTimings["TrustRegion"] = 0.;
    this->status.additionalIterations["Inner"] = 0;


    // Initial value
    x_k = x_0;

    // Parameters
    RealType eta_k = _eta_0;
    RealType tau_k = _tau_0;
    RealType mu_k = _mu_0;

    // Lagrange Multipliers
    VectorType &Lambda = _Lambda;
//    Lambda = VectorType::Constant( _G.getTargetDimension(), 0. );

    // Functionals for inner solver
    AugmentedLagrangeFunctional<ConfiguratorType> aLF( _F, _G, Lambda, mu_k, _alpha );
    AugmentedLagrangeDerivative<ConfiguratorType> aLG( _F, _DF, _G, _DG, Lambda, mu_k, _alpha );
    AugmentedLagrangeHessian<ConfiguratorType> aLH( _F, _DF, _D2F, _G, _DG, _D2G, Lambda, mu_k, _alpha );

    RealType zero = 0;
    AugmentedLagrangeDerivative<ConfiguratorType> LG( _F, _DF, _G, _DG, Lambda, zero, _alpha );

    // Variables for values at iterates
    VectorType constraintViolation( _G.getTargetDimension());
    RealType functionValue;
    VectorType functionGradient;
    RealType ALValue;
    VectorType ALGrad;
    MatrixType ALHess;

    aLF.apply( x_k, ALValue );
    aLG.apply( x_k, ALGrad );
    if ( _bdryMask )
      applyMaskToVector( *_bdryMask, ALGrad );
    RealType initALGradNorm = _relativeStoppingCriterion ? ALGrad.template lpNorm<Eigen::Infinity>() : 1.;

    aLH.apply( x_k, ALHess );
    if ( _bdryMask )
      applyMaskToSymmetricMatrix( *_bdryMask, ALHess );
    _G.apply( x_k, constraintViolation );
    _F.apply( x_k, functionValue );
    _DF.apply( x_k, functionGradient );

//    FullMatrixType mat(ALHess);
//    Eigen::SelfAdjointEigenSolver<FullMatrixType> es;
//    es.compute(mat);
//    std::cerr << " -- Eigenvalues (min/max): " << es.eigenvalues()[0] << " / " << es.eigenvalues()[x_k.size() - 1] << std::endl;

//    Lambda = -mu_k * constraintViolation;


    std::vector<RealType> linearLowerBounds( x_k.size(), -2.e+19 );
    std::vector<RealType> linearUpperBounds( x_k.size(), 2.e+19 );

    if ( _lowerBounds && _upperBounds ) {
      for ( int i = 0; i < x_k.size(); i++ ) {
        linearLowerBounds[i] = (*_lowerBounds)[i];
        linearUpperBounds[i] = (*_upperBounds)[i];
      }
    }

    if ( _bdryMask ) {
      for ( int idx : *_bdryMask ) {
        linearLowerBounds[idx] = x_0[idx];
        linearUpperBounds[idx] = x_0[idx];
      }
    }

    if ( !_quiet ) {
      std::cout << "=======================================================================" << std::endl;
      std::cout << "Start Augmented Lagrange method with " << _maxIterations << " iterations and eps = " << _eta_star
                << "."
                << std::endl;
      std::cout << "=======================================================================" << std::endl;
      std::cout << "#iter -- F -- ||DF|| -- ||G||_inf -- ALF -- ||DL||_inf -- mu -- eta -- tau" << std::endl;
    }


    int iterations = 0;
    if ( !_quiet )
      std::cout << "AL -- Iter " << std::setw( 3 ) << iterations << ": "
                << std::scientific << functionValue
                << " || " << constraintViolation.template lpNorm<Eigen::Infinity>()
                << " || " << ALGrad.template lpNorm<Eigen::Infinity>() << (_relativeStoppingCriterion ? "*" : "")
                << " ||===|| " << mu_k
                << " || " << eta_k
                << " || " << tau_k
                << std::endl;

    while ( iterations < _maxIterations &&
            (constraintViolation.template lpNorm<Eigen::Infinity>() > _eta_star ||
             ALGrad.template lpNorm<Eigen::Infinity>() / initALGradNorm > _tau_star)) {
      iterations++;
      this->status.Iteration++;

      auto t_start = std::chrono::high_resolution_clock::now();

      // Solve Subproblem
      auto t_start_inner = std::chrono::high_resolution_clock::now();
      SolverStatus<ConfiguratorType> innerStatus = solveInnerProblem( x_k, aLF, aLG, aLH, tau_k,
                                                                      linearLowerBounds, linearUpperBounds );
      auto t_end_inner = std::chrono::high_resolution_clock::now();

      this->status.additionalIterations["Inner"] += innerStatus.Iteration;
      for ( const auto &timing : innerStatus.additionalTimings )
        this->status.additionalTimings["Inner__" + timing.first] += timing.second;


      this->status.additionalTimings["InnerOptimization"] += std::chrono::duration<RealType, std::milli>(
              t_end_inner - t_start_inner ).count();

      // Evaluate functionals at new iterate
      aLF.apply( x_k, ALValue );

      _G.apply( x_k, constraintViolation );
      _F.apply( x_k, functionValue );
      _DF.apply( x_k, functionGradient );

      aLG.apply( x_k, ALGrad );
      if ( _bdryMask )
        applyMaskToVector( *_bdryMask, ALGrad );

      // Check how good constraints are fulfilled
      if ( constraintViolation.template lpNorm<Eigen::Infinity>() <= eta_k ) {
        // Update multipliers, tighten tolerances
        for ( int i = 0; i < _G.getTargetDimension(); i++ )
          Lambda[i] = Lambda[i] - mu_k * constraintViolation[i];
        eta_k = std::max( eta_k / std::pow( mu_k, _etaDecreaseExponent ), _eta_star );
        tau_k = std::max( tau_k / mu_k, _tau_star );
      }
      else if ( mu_k != _maxPenalty ) {
        // Increase penalty parameter, tighten tolerances
        mu_k = std::min( _maxPenalty, _penaltyIncrease * mu_k );
        eta_k = std::max( 1 / std::pow( mu_k, 0.1 ), _eta_star );
        tau_k = std::max( 1 / mu_k, _tau_star );
      }

      auto t_end = std::chrono::high_resolution_clock::now();
      this->status.totalTime += std::chrono::duration<RealType, std::milli>( t_end - t_start ).count();

      if ( !_quiet )
        std::cout << "AL -- Iter " << std::setw( 3 ) << iterations << ": "
                  << std::scientific << std::setprecision(6) << functionValue
                  << " || " << constraintViolation.template lpNorm<Eigen::Infinity>()
                  << " || " << ALGrad.template lpNorm<Eigen::Infinity>() / initALGradNorm
                  << " ||===|| " << mu_k
                  << " || " << eta_k
                  << " || " << tau_k
                  << std::fixed << std::setprecision( 2 )
                  << " ||===|| " << std::chrono::duration<double, std::milli>( t_end - t_start ).count()
                  << " || " << std::chrono::duration<double, std::milli>( t_end_inner - t_start_inner ).count()
                  << std::endl;
    }


    if ( !_quiet ) {
      std::cout << "=======================================================================" << std::endl;
      std::cout << "Finished Augmented Lagrange method after " << iterations << " steps." << std::endl;
      std::cout << "AL -- Iter " << std::setw( 3 ) << iterations << ": "
                << std::scientific << functionValue
                << " || " << constraintViolation.template lpNorm<Eigen::Infinity>()
                << " || " << ALGrad.template lpNorm<Eigen::Infinity>() / initALGradNorm
                << " ||===|| " << mu_k
                << " || " << eta_k
                << " || " << tau_k
                << std::fixed
                << " ||===|| " << this->status.totalTime
                << " || " << this->status.additionalTimings["InnerOptimization"]
                << std::endl;
      std::cout << "=======================================================================" << std::endl << std::endl;
    }
  }

  SolverStatus<ConfiguratorType> solveInnerProblem( VectorType &x_k,
                                                    const AugmentedLagrangeFunctional<ConfiguratorType> &aLF,
                                                    const AugmentedLagrangeDerivative<ConfiguratorType> &aLG,
                                                    const AugmentedLagrangeHessian<ConfiguratorType> &aLH,
                                                    const RealType &tau_k,
                                                    const std::vector<RealType> &linearLowerBounds,
                                                    const std::vector<RealType> &linearUpperBounds ) const {
    if ( _innerMethod == "TRN" ) {
      TrustRegionNewton<ConfiguratorType> innerOpt( aLF, aLG, aLH, 10., std::numeric_limits<RealType>::infinity(),
                                                    tau_k, 300, x_k.size());
      innerOpt.setParameters( _innerIntParameters );
      innerOpt.setParameters( _innerStringParameters );
      innerOpt.setParameters( _innerRealParameters );
      innerOpt.setParameter( "tolerance", tau_k ); // Just to be safe

      innerOpt.setBoundaryMask( *_bdryMask );

      innerOpt.solve( x_k, x_k );

      return innerOpt.status;
    }
    else if ( _innerMethod == "LSN" ) {
      int LSN_gditer = _innerIntParameters.count( "gradient_iterations" ) ?
                       _innerIntParameters.at( "gradient_iterations" ) : 0;
      int LSN_bfgsiter = _innerIntParameters.count( "BFGS_iterations" ) ?
                         _innerIntParameters.at( "BFGS_iterations" ) : 0;
      TIMESTEP_CONTROLLER LSN_timestep = static_cast<TIMESTEP_CONTROLLER > (
              _innerIntParameters.count( "stepsize_control" ) ? _innerIntParameters.at( "stepsize_control" ) : 1);
      bool LSN_quiet = _innerIntParameters.count( "print_level" ) ?
                       (_innerIntParameters.at( "print_level" ) < 5) : _quiet;

      if ( LSN_gditer > 0 ) {
        GradientDescent<ConfiguratorType> innerGD( aLF, aLG, LSN_gditer, tau_k, LSN_timestep, LSN_quiet );
        innerGD.setBoundaryMask( *_bdryMask );
        innerGD.solve( x_k, x_k );
      }

      if ( LSN_bfgsiter > 0 ) {
        QuasiNewtonBFGS<ConfiguratorType> innerBFGS( aLF, aLG, LSN_bfgsiter, tau_k, LSN_timestep, 50, LSN_quiet );
        innerBFGS.setBoundaryMask( *_bdryMask );
        innerBFGS.solve( x_k, x_k );
      }


      LineSearchNewton<ConfiguratorType> innerOpt( aLF, aLG, aLH, tau_k, 300, _quiet );

      innerOpt.setParameters( _innerIntParameters );
      innerOpt.setParameters( _innerStringParameters );
      innerOpt.setParameters( _innerRealParameters );
      innerOpt.setParameter( "tolerance", tau_k ); // Just to be safe

      innerOpt.setBoundaryMask( *_bdryMask );


      innerOpt.solve( x_k, x_k );

      return innerOpt.status;
    }

    else if ( _innerMethod == "CN" ) {
      int LSN_gditer = _innerIntParameters.count( "gradient_iterations" ) ?
                       _innerIntParameters.at( "gradient_iterations" ) : 500;
      int LSN_bfgsiter = _innerIntParameters.count( "BFGS_iterations" ) ?
                         _innerIntParameters.at( "BFGS_iterations" ) : 500;
      int LSN_niter = _innerIntParameters.count( "newton_iterations" ) ?
                      _innerIntParameters.at( "newton_iterations" ) : 500;
      int LSN_solver = _innerIntParameters.count( "linear_solver" ) ?
                       _innerIntParameters.at( "linear_solver" ) : 1;
      QUIET_MODE LSN_quiet = _innerIntParameters.count( "print_level" ) ?
                             (_innerIntParameters.at( "print_level" ) < 5 ? SUPERQUIET : SHOW_ALL) :
                             (_quiet ? SUPERQUIET : SHOW_ALL);
      TIMESTEP_CONTROLLER LSN_timestep = static_cast<TIMESTEP_CONTROLLER > (
              _innerIntParameters.count( "stepsize_control" ) ? _innerIntParameters.at( "stepsize_control" ) : 1);

      GradientDescent<ConfiguratorType> innerGD( aLF, aLG, LSN_gditer, tau_k, LSN_timestep, LSN_quiet );
      innerGD.setBoundaryMask( *_bdryMask );
      innerGD.solve( x_k, x_k );

      QuasiNewtonBFGS<ConfiguratorType> innerBFGS( aLF, aLG, LSN_bfgsiter, tau_k, LSN_timestep, 50, LSN_quiet );
      innerBFGS.setBoundaryMask( *_bdryMask );
      innerBFGS.solve( x_k, x_k );

      NewtonOptimizationMethod<ConfiguratorType> innerOpt( aLF, aLG, aLH, LSN_niter, tau_k, LSN_timestep, LSN_quiet );

      innerOpt.setBoundaryMask( *_bdryMask );
      innerOpt.setSolver( static_cast<LINEAR_SOLVER_TYPE >(LSN_solver));

      innerOpt.solve( x_k, x_k );

      return innerOpt.status;
    }
#ifdef GOAST_WITH_IPOPT
    else if ( _innerMethod == "Ipopt" ) {
      int ipopt_iter = _innerIntParameters.count( "maximum_iterations" ) ?
                       _innerIntParameters.at( "maximum_iterations" ) : 500;
      int ipopt_solver = _innerIntParameters.count( "linear_solver" ) ?
                         _innerIntParameters.at( "linear_solver" ) : 0;
      int ipopt_printlevel = _innerIntParameters.count( "print_level" ) ?
                             _innerIntParameters.at( "print_level" ) : (_quiet ? 0 : 5);

      IpoptBoxConstraintSecondOrderSolver<ConfiguratorType, true> innerOpt( aLF, aLG, aLH,
                                                                            ipopt_iter,
                                                                            tau_k * 1e2,
                                                                            linearLowerBounds,
                                                                            linearUpperBounds,
                                                                            ipopt_solver,
                                                                            ipopt_printlevel );
      innerOpt.solve( x_k, x_k );

      return innerOpt.status;
    }
#endif
    else {
      throw std::runtime_error( "AugmentedLagrangeMethod::solveInnerProblem(): Unknown inner method!" );
    }
  }

  void setParameter( const std::string &name, std::string value ) override {
    if ( name == "inner_method" )
      _innerMethod = value; //! \todo Add safeguard against undefined methods
    else if ( name.rfind( "inner__", 0 ) == 0 )
      _innerStringParameters[name.substr( 7, std::string::npos )] = value;
    else
      throw std::logic_error( "AugmentedLagrangeMethod::setParameter(): This class has no string parameters." );
  }

  void setParameter( const std::string &name, RealType value ) override {
    if ( name == "maximal_penalty" )
      _maxPenalty = value;
    else if ( name == "initial_penalty" ) {
      _mu_0 = value;
      _eta_0 = 1 / std::pow( _mu_0, 0.1 );
      _tau_0 = 1 / _mu_0;
    }
    else if ( name == "penalty_factor" )
      _penaltyIncrease = value;
    else if ( name == "constraint_decrease_exponent" )
      _etaDecreaseExponent = value;
    else if ( name == "constraint_tolerance" )
      _eta_star = value;
    else if ( name == "optimality_tolerance" )
      _tau_star = value;
    else if ( name == "functional_scaling" )
      _alpha = value;
    else if ( name.rfind( "inner__", 0 ) == 0 )
      _innerRealParameters[name.substr( 7, std::string::npos )] = value;
    else
      throw std::runtime_error( "AugmentedLagrangeMethod::setParameter(): Unknown parameter '" + name + "'." );
  }

  void setParameter( const std::string &name, int value ) override {
    if ( name == "maximum_iterations" )
      _maxIterations = value;
    else if ( name == "print_level" )
      _quiet = (value < 5);
    else if ( name == "relative_stopping" )
      _relativeStoppingCriterion = static_cast<bool>(value);
    else if ( name.rfind( "inner__", 0 ) == 0 )
      _innerIntParameters[name.substr( 7, std::string::npos )] = value;
    else
      throw std::runtime_error( "AugmentedLagrangeMethod::setParameter(): Unknown parameter '" + name + "'." );
  }
};

#endif //OPTIMIZATION_AUGMENTEDLAGRANGE_H
