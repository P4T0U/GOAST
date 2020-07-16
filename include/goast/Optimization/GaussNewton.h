// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Header for Gauss-Newton method
 * \author Heeren
 *
 * \todo Convert to new interface
 * \todo Documentation!
 */

#ifndef REDUCEDBASIS_GAUSSNEWTON_H
#define REDUCEDBASIS_GAUSSNEWTON_H

#include "optInterface.h"
#include "stepsizeControl.h"

//! \brief  Stepsize control for Gauss-Newton
//! \author Sassen
template<typename ConfiguratorType>
class StepsizeControlForGaussNewton : public StepsizeControlInterface<ConfiguratorType> {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;

  const BaseOp<VectorType, VectorType> &_E;
  const BaseOp<VectorType, MatrixType> &_DE;

public:
  StepsizeControlForGaussNewton( const BaseOp<VectorType, VectorType> &E,
                                 const BaseOp<VectorType, MatrixType> &DE,
                                 TIMESTEP_CONTROLLER timestepController,
                                 RealType sigma = 0.1,
                                 RealType beta = 0.9,
                                 RealType startTau = 1.0,
                                 RealType tauMin = 1e-12,
                                 RealType tauMax = 4. ) : StepsizeControlInterface<ConfiguratorType>( sigma, beta,
                                                                                                      timestepController,
                                                                                                      startTau, tauMin,
                                                                                                      tauMax ),
                                                          _E( E ), _DE( DE ) {}

  //Returns the scalar objective function evaluated at CurrentPosition.
  RealType evaluateLinesearchFnc( const VectorType &CurrentPosition ) const {
    VectorType f;
    _E.apply( CurrentPosition, f );
    return 0.5 * f.squaredNorm();
  }

  //Returns the dot product of the energy derivative and the descent direction.
  RealType evaluateLinesearchGrad( const VectorType &Position, const VectorType &DescentDir ) const {
    VectorType linearPart;
    _E.apply( Position, linearPart );

    MatrixType JacobianOfLinearPartTransp;
    _DE.applyTransposed( Position, JacobianOfLinearPartTransp );

    VectorType tmp = JacobianOfLinearPartTransp * linearPart;

    if ( this->_bdryMask )
      applyMaskToVector( *this->_bdryMask, tmp );

    return tmp.dot( DescentDir );
  }

};


/**
 * \brief Implements the Gauss-Newton algorithm to optimize a scalar-valued function \f$E[x] = 1/2 f(x)^T f(x)\f$ for \f$x \in R^n\f$
 * \author Heeren
 *
 * Here \f$ f: R^n \to R^m\f$, \f$J = DE[x] = Df(x)^T f(x)\f$ and \f$ m \f$ is the dimension of the range of \f$ f \f$.
 *
 * The Hessian of E is approximated by \f$ D^2E[x] \approx Df(x)^T Df(x)\f$, assuming either \f$D^2f\f$ or \f$f\f$ being small
 */
template<typename ConfiguratorType>
class GaussNewtonAlgorithm {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::VectorType VectorType;

  //! stores the evaluation of DF and is deleted in the destructor.
  mutable MatrixType _MatDfTransp;

  const int _dimRangeF;
  const BaseOp<VectorType> &_F;
  const BaseOp<VectorType, MatrixType> &_DF;
  int _maxIterations;
  RealType _stopCriterion, _sigma, _tauMin, _tauMax;
  bool _verbose;
  std::vector<int> _bdryMask;
  StepsizeControlForGaussNewton<ConfiguratorType> _stepsizeControl;

public:
  mutable SolverStatus<ConfiguratorType> status;

  GaussNewtonAlgorithm( int DimRangeF,
                        const BaseOp<VectorType> &F,
                        const BaseOp<VectorType, MatrixType> &DF,
                        const int MaxIterations = 50,
                        const RealType stopCriterion = 1e-8,
                        const TIMESTEP_CONTROLLER TimestepController = SIMPLE_TIMESTEP_CONTROL,
                        const bool quiet = false,
                        const RealType sigma = 0.1,
                        const RealType tauMin = 1.e-6,
                        const RealType tauMax = 4. )
          : _dimRangeF( DimRangeF ),
            _F( F ),
            _DF( DF ),
            _maxIterations( MaxIterations ),
            _stopCriterion( stopCriterion ),
            _sigma( sigma ),
            _tauMin( tauMin ),
            _tauMax( tauMax ),
            _verbose( !quiet ),
            _stepsizeControl( _F, _DF, TimestepController, sigma, 0.9, 1., tauMin, tauMax ) {}

  void setQuiet() { _verbose = false; }

  void setBoundaryMask( const std::vector<int> &mask ) {
    _bdryMask.resize( mask.size());
    _bdryMask = mask;
  }

  void solve( const VectorType &Arg, VectorType &Dest ) const {
    status.Iteration = 0;
    status.totalTime = 0.;
    status.additionalTimings["Gradient"] = 0.;
    status.additionalTimings["Hessian"] = 0.;
    status.additionalTimings["LinearSolver"] = 0.;
    status.additionalTimings["Stepsize"] = 0.;
    status.additionalTimings["Update"] = 0.;

    auto total_start = std::chrono::high_resolution_clock::now();

    Dest.resize( Arg.size());
    VectorType f( _dimRangeF );
    VectorType direction( Dest.size());
    _MatDfTransp.resize( Arg.size(), _dimRangeF );

    int iteration = 0;
    RealType energyOld = std::numeric_limits<RealType>::infinity();
    RealType energy = std::numeric_limits<RealType>::infinity();

    Dest = Arg;

    _F.apply( Dest, f );
    energy = 0.5 * f.squaredNorm();
    status.Residual = energy;

    if ( _verbose ) std::cerr << "Initial energy " << energy << std::endl;

    if ( energy < 1e-15 ) return;

    RealType tau = 1.;
    RealType res = 1.;

    while ((iteration < _maxIterations) && (res > _stopCriterion)) {

      auto t_start = std::chrono::high_resolution_clock::now();

      // J^T = (Df)^T
      auto t_start_grad = std::chrono::high_resolution_clock::now();
      _DF.applyTransposed( Dest, _MatDfTransp );

      VectorType currGradient = _MatDfTransp * f;
      if ( _bdryMask.size() > 0 )
        applyMaskToVector( _bdryMask, currGradient );

      auto t_end_grad = std::chrono::high_resolution_clock::now();
      status.additionalTimings["Gradient"] += std::chrono::duration<double, std::milli>( t_end_grad - t_start_grad ).count();

      // solve J^T J d = J^T f for direction d, with J = Df
      // Caution: direction is negative direction due to "wrong" sign of rhs!

      computeNegativeDirection( _MatDfTransp, currGradient, direction );


      if ( checkForNANsAndINFs( direction ) )
        throw BasicException("GaussNewtonAlgorithmBase::solve(): NANs in direction!");

      // compute descent direction as negative gradient and stepsize
      direction *= -1;
      auto t_start_step = std::chrono::high_resolution_clock::now();
      tau = _stepsizeControl.getStepsize( Dest, currGradient, direction, std::min( 2 * tau, 1. ), energy );
      auto t_end_step = std::chrono::high_resolution_clock::now();
      status.additionalTimings["Stepsize"] += std::chrono::duration<double, std::milli>( t_end_step - t_start_step ).count();

      // update
      auto t_start_update = std::chrono::high_resolution_clock::now();
      Dest += tau * direction;
      _F.apply( Dest, f );
      energy = 0.5 * f.squaredNorm();
      status.Residual = energy;
      auto t_end_update = std::chrono::high_resolution_clock::now();
      status.additionalTimings["Update"] += std::chrono::duration<double, std::milli>( t_end_update - t_start_update ).count();

      // terminate if stepsize is zero or energy haven't decreased
      if ((energyOld < energy) || tau < 1e-15 )
        break;

      // recompute residual and update of variables
      res = computeResidual( energyOld, energy, tau, direction );
      energyOld = energy;
      ++iteration;
      status.Iteration++;

      auto t_end = std::chrono::high_resolution_clock::now();
      if ( _verbose )
        std::cout << "Step " << iteration << ", tau " << tau << ", energy " << energy << ", res = " << res;
      if ( _verbose )
        std::cout << ", t = " << std::chrono::duration<double, std::ratio<1> >( t_end - t_start ).count() << "s"
                                                                                                          << std::endl;

      status.totalTime += std::chrono::duration<double, std::milli>( t_end - t_start ).count();

    }

    auto total_end = std::chrono::high_resolution_clock::now();
    if ( _verbose )
      std::cout << std::fixed << "GaussNewton finished after " << iteration << " steps, time = "
                << std::chrono::duration<double, std::ratio<1> >( total_end - total_start ).count() << "s" << std::endl;

  }

protected:
  // E(x) = 0.5 f^T(x)f(x), J = Df, gradient is DE = J^Tf
  virtual RealType computeResidual( RealType energyOld, RealType energy, RealType tau, const VectorType &Direction ) const {
    // difference in energies
    //return energyOld - energy;

    // difference in positions (x_new = x_old + tau * dir)
    //return tau * Direction.norm();

    // norm of grad
    return Direction.norm();
  }

  // F(x) = 0.5 f^T(x)f(x),  then DF = J^Tf and D^2F \approx  J^TJ for J = Df
  // NOTE if we fix boundary nodes, J^T and rhs = DF have been masked before
  void computeNegativeDirection( const MatrixType &Jtrans, const VectorType &gradF, VectorType &solution ) const {

    // allocate and compute J^TJ, J = Df(x)
    auto t_start_hess = std::chrono::high_resolution_clock::now();
    MatrixType JtransJ = Jtrans * Jtrans.transpose();

    if ( this->_bdryMask.size() > 0 ) {
      applyMaskToSymmetricMatrix( this->_bdryMask, JtransJ );
    }
    else {
      // add eps * Id
      RealType eps = 1e-8;
      for ( int j = 0; j < JtransJ.rows(); j++ )
        JtransJ.coeffRef( j, j ) += eps;
    }
    auto t_end_hess = std::chrono::high_resolution_clock::now();
    status.additionalTimings["Hessian"] += std::chrono::duration<double, std::milli>( t_end_hess - t_start_hess ).count();

    // rhs = J^T*f(x) (wrong sign here!)
    // NOTE if we fix boundary nodes, J^T and rhs = DF have been masked before
    auto t_start_solve = std::chrono::high_resolution_clock::now();
    LinearSolver<ConfiguratorType>( ).solve( JtransJ, gradF, solution );
    auto t_end_solve = std::chrono::high_resolution_clock::now();
    status.additionalTimings["LinearSolver"] += std::chrono::duration<double, std::milli>( t_end_solve - t_start_solve ).count();

  }

};


#endif //GAUSSNEWTON_H
