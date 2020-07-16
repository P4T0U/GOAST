// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Header for a simple Newton's method
 * \author Heeren
 *
 * \todo Convert to new interface
 * \todo Documentation!
 */

#ifndef OPTIMIZATION_SIMPLENEWTON_H
#define OPTIMIZATION_SIMPLENEWTON_H

#include "optInterface.h"
#include "stepsizeControl.h"
#include <goast/Core/LinearSolver.h>
#include "optParameters.h"

//! Stepsize control for Newton
//! \author Heeren
template<typename ConfiguratorType>
class StepsizeControlForNewton : public StepsizeControlInterface<ConfiguratorType> {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::VectorType VectorType;

  const BaseOp<VectorType, VectorType> &_F;
  const MatrixType &_Jacobi;

public:
  StepsizeControlForNewton( const BaseOp<VectorType, VectorType> &F,
                            const MatrixType &Jacobi,
                            TIMESTEP_CONTROLLER timestepController = NEWTON_OPTIMAL,
                            RealType sigma = 0.1,
                            RealType beta = 0.9,
                            RealType startTau = 1.0,
                            RealType tauMin = 1e-12,
                            RealType tauMax = 4. ) : StepsizeControlInterface<ConfiguratorType>( sigma, beta,
                                                                                                 timestepController,
                                                                                                 startTau, tauMin,
                                                                                                 tauMax ), _F( F ),
                                                     _Jacobi( Jacobi ) {}

  //Returns the scalar objective function evaluated at CurrentPosition.
  RealType evaluateLinesearchFnc( const VectorType &Position ) const {
    VectorType temp( Position.size());
    _F.apply( Position, temp );
    if ( this->_bdryMask )
      applyMaskToVector( *this->_bdryMask, temp );
    return 0.5 * temp.squaredNorm();
  }

  //Returns the dot product of the energy derivative and the descent direction.
  RealType evaluateLinesearchGrad( const VectorType &Position, const VectorType &DescentDir ) const {
    VectorType outer = _Jacobi * DescentDir;
    VectorType grad( DescentDir.size());
    _F.apply( Position, grad );
    if ( this->_bdryMask )
      applyMaskToVector( *this->_bdryMask, grad );
    return outer.dot( grad );
  }

  RealType evaluateLinesearchFnc( const VectorType &CurrentPosition, const VectorType &DescentDir,
                                  RealType timestepWidth ) const {
    return evaluateLinesearchFnc( CurrentPosition + timestepWidth * DescentDir );
  }

  // Calculates the step size based on "Schaback, Werner - Numerische Mathematik, 4. Auflage, Seite 129". Is global convergent, as long as F fulfils certain regularity conditions.
  RealType getNewtonOptimalTimestepWidth( const VectorType &CurrentPosition, const VectorType &DescentDir ) const {

    const RealType DescentDirNormSqr = DescentDir.squaredNorm();

    if ( !(DescentDirNormSqr > 0.))
      return 0.;

    // initial guess for L: L = ||F(y)-F(x)-F'(x)(y-x)|| / ||y-x||^2, where x = CurrentPosition, y = x + DescentDir
    VectorType newPosition = CurrentPosition + DescentDir;

    VectorType pTmp( DescentDir.size());
    _F.apply( newPosition, pTmp );
    if ( this->_bdryMask )
      applyMaskToVector( *this->_bdryMask, pTmp );
    pTmp -= _Jacobi * DescentDir;
    VectorType Fx( DescentDir.size());
    _F.apply( CurrentPosition, Fx );
    if ( this->_bdryMask )
      applyMaskToVector( *this->_bdryMask, Fx );
    pTmp -= Fx;
    RealType L = pTmp.norm() / DescentDirNormSqr;

    // If F is locally linear and matches the gradient with numerical perfection take the full step
    if ( L < 1e-15 )
      return 1.0;

    RealType fNorm = Fx.norm();
    RealType tau = std::min( fNorm / (2. * L * DescentDirNormSqr), 1.0 );
    RealType temp1 = fNorm - std::sqrt( 2. * evaluateLinesearchFnc( CurrentPosition, DescentDir, tau ));
    RealType temp2 = tau * (fNorm - L * tau * DescentDirNormSqr);
    RealType tauCandidate;
    if ( temp1 < temp2 ) {
      do {
        L *= 2.;
        tau = std::min( fNorm / (2. * L * DescentDirNormSqr), 1.0 );
        temp1 = fNorm - sqrt( 2. * evaluateLinesearchFnc( CurrentPosition, DescentDir, tau ));
        temp2 = tau * (fNorm - L * tau * DescentDirNormSqr);
        // Prevent the while loop from getting stuck if DescentDir is not a descent direction.
        if ( tau < this->_tauMin ) {
          temp1 = temp2;
          tau = 0.;
        }
      } while ( temp1 < temp2 );
    } else {
      do {
        // Compute tauCandidate, temp1 and temp2 with L/2 instead of L to find out if a bigger step size is still feasible.
        tauCandidate = std::min( fNorm / (L * DescentDirNormSqr), 1.0 );
        if ( tauCandidate == 1.0 )
          break;
        temp1 = fNorm - std::sqrt( 2. * evaluateLinesearchFnc( CurrentPosition, DescentDir, tauCandidate ));
        temp2 = tauCandidate * (fNorm - 0.5 * L * tauCandidate * DescentDirNormSqr);

        if ( temp1 >= temp2 ) {
          L /= 2.;
          tau = tauCandidate;
        }
      } while ( temp1 >= temp2 );
    }

    // return
    return tau > this->_tauMin ? tau : 0.;

  }

  RealType getStepsize( const VectorType &CurrentPosition, const VectorType &CurrentGradient, const VectorType &descentDir, RealType tau_before, RealType currEnergy = -1. ) const {
    switch ( this->_timestepController ) {
      case NEWTON_OPTIMAL:
        return getNewtonOptimalTimestepWidth( CurrentPosition, descentDir );
      default:
        return StepsizeControlInterface<ConfiguratorType>::getStepsize( CurrentPosition, CurrentGradient, descentDir, tau_before, currEnergy );
    }
  }

};

//!==========================================================================================================
//! Newton method to find a root of a vector valued functional F: \R^n -> \R^n
//! \author Heeren
template<typename ConfiguratorType>
class NewtonMethod {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;

  const BaseOp<VectorType, VectorType> &_F;
  const BaseOp<VectorType, MatrixType> &_DF;
  mutable MatrixType _Jacobi;

  TIMESTEP_CONTROLLER _timestepController;
  std::unique_ptr<StepsizeControlInterface<ConfiguratorType>> _stepsizeControlPtr;

  LINEAR_SOLVER_TYPE _solverType;
  int _maxIterations;
  RealType _stopEpsilon;
  QUIET_MODE _quietMode;
  const std::vector<int> *_bdryMask;

public:
  mutable SolverStatus<ConfiguratorType> status;

  NewtonMethod( const BaseOp<VectorType, VectorType> &F,
                const BaseOp<VectorType, MatrixType> &DF,
                int MaxIterations,
                RealType StopEpsilon,
                TIMESTEP_CONTROLLER TimestepController,
                bool quiet,
                RealType sigma = 0.1,
                RealType tauMin = 1.e-6,
                RealType tauMax = 4. )
          : _F( F ), _DF( DF ), _timestepController( static_cast<TIMESTEP_CONTROLLER>(TimestepController)),
            _solverType( UMFPACK_LU_FACT ), _maxIterations( MaxIterations ), _stopEpsilon( StopEpsilon ),
            _quietMode(  quiet ? SUPERQUIET : SHOW_ALL ), _bdryMask( nullptr ) {

    _stepsizeControlPtr = std::unique_ptr<StepsizeControlInterface<ConfiguratorType>>(
            new StepsizeControlForNewton<ConfiguratorType>( _F, _Jacobi, _timestepController, sigma, 0.9, 1., tauMin,
                                                            tauMax ));
  }
  
  NewtonMethod( const BaseOp<VectorType, VectorType> &F,
                const BaseOp<VectorType, MatrixType> &DF,
                int MaxIterations = 1000,
                RealType StopEpsilon = 1e-8,
                TIMESTEP_CONTROLLER TimestepController = NEWTON_OPTIMAL,
                QUIET_MODE quietMode = SUPERQUIET,
                RealType sigma = 0.1,
                RealType tauMin = 1.e-6,
                RealType tauMax = 4. )
          : _F( F ), _DF( DF ), _timestepController( static_cast<TIMESTEP_CONTROLLER>(TimestepController)),
            _solverType( UMFPACK_LU_FACT ), _maxIterations( MaxIterations ), _stopEpsilon( StopEpsilon ),
            _quietMode( quietMode ), _bdryMask( nullptr ) {

    _stepsizeControlPtr = std::unique_ptr<StepsizeControlInterface<ConfiguratorType>>(
            new StepsizeControlForNewton<ConfiguratorType>( _F, _Jacobi, _timestepController, sigma, 0.9, 1., tauMin,
                                                            tauMax ));
  }

  NewtonMethod( const BaseOp<VectorType, VectorType> &F,
                const BaseOp<VectorType, MatrixType> &DF,
                const OptimizationParameters<ConfiguratorType> &optPars )
          : _F( F ), _DF( DF ), _timestepController( static_cast<TIMESTEP_CONTROLLER>(optPars.getNewtonTimeStepping())),
            _solverType( static_cast<LINEAR_SOLVER_TYPE>(optPars.getSolverType())),
            _maxIterations( optPars.getNewtonIterations()), _stopEpsilon( optPars.getStopEpsilon()),
            _quietMode( optPars.getQuietMode()), _bdryMask( nullptr ) {

    _stepsizeControlPtr = std::unique_ptr<StepsizeControlInterface<ConfiguratorType>>(
            new StepsizeControlForNewton<ConfiguratorType>( _F, _Jacobi, _timestepController, optPars.getSigma(),
                                                            optPars.getBeta(), optPars.getStartTau(),
                                                            optPars.getTauMin(), optPars.getTauMax()));
  }

  // residuum in iteration
  RealType computeErrorNorm( const VectorType &x_k, const VectorType &delta_x_k, const VectorType &F_x_k, RealType tau,
                             bool initial ) const {
    if ( initial )
      return 1.;
    else
      return delta_x_k.norm();
  }

  //
  void setSolver( LINEAR_SOLVER_TYPE solverType ) {
    _solverType = solverType;
  }

  // set boundary mask
  void setBoundaryMask( const std::vector<int> &Mask ) {
    _bdryMask = &Mask;
    _stepsizeControlPtr->setBoundaryMask( Mask );
  }

  // x^{k+1} = x^k + tau * d^k, where d^k solves D^2E[x^k] d^k = - DE[x^k]
  bool solve( const VectorType &x_0, VectorType &x_k ) const {
    status.Iteration = 0;
    status.totalTime = 0.;

    x_k = x_0;

    VectorType F_x_k( x_k.size());
    VectorType delta_x_k( x_k.size());
    _F.apply( x_k, F_x_k );
    if ( _bdryMask )
      applyMaskToVector( *_bdryMask, F_x_k );
    _Jacobi.resize( x_k.size(), x_k.size());

    RealType tau = _stepsizeControlPtr->getStartTau();
    RealType FNorm = computeErrorNorm( x_k, delta_x_k, F_x_k, tau, true );

    if ( _maxIterations == 0 )
      return FNorm <= _stopEpsilon;

    if ( _quietMode == SHOW_ALL ) {
      std::cout << "=========================================================================================" << std::endl;
      std::cout << "Start Newton method with " << _maxIterations << " iterations and eps = " << _stopEpsilon << "."
                << std::endl;
      writeOutput( x_k, 0, tau, FNorm, false );
      std::cout << "=========================================================================================" << std::endl;
    }

    int iterations = 0;
    while ( FNorm > _stopEpsilon && (iterations < _maxIterations) && tau > 0. ) {
      iterations++;

      auto t_start = std::chrono::high_resolution_clock::now();

      // Newton iteration given by x^{k+1} = x^k - tau D2F(x^k)^{-1}(DF(x^k))
      _DF.apply( x_k, _Jacobi );
      //TODO check whether mask also works for non-symmetric matrices!
      if ( _bdryMask )
        applyMaskToSymmetricMatrixAndVector( *_bdryMask, _Jacobi, F_x_k );

      //std::cout << _Jacobi << std::endl;
      VectorType rhs( F_x_k );
      rhs *= -1.;
      LinearSolver<ConfiguratorType>( _solverType ).solve( _Jacobi, rhs, delta_x_k );

      // get tau
      tau = _stepsizeControlPtr->getStepsize( x_k, F_x_k, delta_x_k, tau );

      if ( tau > 0 ) {
        // update position and descent direction
        x_k += tau * delta_x_k;
        _F.apply( x_k, F_x_k );
        if ( _bdryMask )
          applyMaskToVector( *_bdryMask, F_x_k );
        FNorm = computeErrorNorm( x_k, delta_x_k, F_x_k, tau, false );

        status.Iteration++;
      }

      auto t_end = std::chrono::high_resolution_clock::now();
      status.totalTime += std::chrono::duration<double, std::milli>( t_end - t_start ).count();

      if ( _quietMode == SHOW_ALL )
        writeOutput( x_k, iterations, tau, FNorm );
    } // end while

    if ( (_quietMode == SHOW_ALL) || (_quietMode == SHOW_TERMINATION_INFO) || ( (_quietMode == SHOW_ONLY_IF_FAILED) && (FNorm > _stopEpsilon) )  ) {
      std::cout << "=========================================================================================" << std::endl;
      std::cout << "Finished Newton's method after " << iterations << " steps (max. steps = " << _maxIterations << ", tol = " << _stopEpsilon << ")." << std::endl;
      std::cout << "Final stepsize = " << tau << ", final norm = " << FNorm << std::endl;
      writeOutput( x_k, iterations, tau, FNorm, false );
      std::cout << "=========================================================================================" << std::endl << std::endl;
    }

    return FNorm <= _stopEpsilon;

  }

protected:

  virtual void writeOutput( const VectorType &/*x_k*/, int iterations, RealType tau, RealType norm,
                            bool intermediate = true ) const {
    if ( intermediate )
      std::cout << std::scientific << "step = " << iterations << ", stepsize = " << tau << ", norm = " << norm << std::endl;
  }


};


//!==========================================================================================================
//! Newton method to optimize a functional F by looking for a root of DF
//! \author Heeren
template<typename ConfiguratorType>
class NewtonOptimizationMethod : public NewtonMethod<ConfiguratorType> {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;

  const BaseOp<VectorType, RealType> &_energy;

public:
  NewtonOptimizationMethod( const BaseOp<VectorType, RealType> &F,
                            const BaseOp<VectorType, VectorType> &DF,
                            const BaseOp<VectorType, MatrixType> &D2F,
                            const int MaxIterations = 1000,
                            const RealType StopEpsilon = 1e-8,
                            const TIMESTEP_CONTROLLER TimestepController = NEWTON_OPTIMAL,
                            const QUIET_MODE quietMode = SUPERQUIET,
                            const RealType sigma = 0.1,
                            const RealType tauMin = 1.e-6,
                            const RealType tauMax = 4. )
          : NewtonMethod<ConfiguratorType>( DF, D2F, MaxIterations, StopEpsilon, TimestepController, quietMode, sigma,
                                            tauMin, tauMax ), _energy( F ) { }
  
    NewtonOptimizationMethod( const BaseOp<VectorType, RealType> &F,
                              const BaseOp<VectorType, VectorType> &DF,
                              const BaseOp<VectorType, MatrixType> &D2F,
                              const OptimizationParameters<ConfiguratorType> &optPars )
          : NewtonMethod<ConfiguratorType>( DF, D2F, optPars ), _energy( F ) {  }


protected:
  virtual void writeOutput( const VectorType &x_k, int iterations, RealType tau, RealType norm,
                            bool intermediate = true ) const {
    RealType value;
    _energy.apply( x_k, value );
    if ( intermediate )
      std::cout << std::scientific << "step = " << iterations << ", stepsize = " << tau << ", norm = " << norm
                << ", energy = " << value << std::endl;
    else
      std::cout << std::scientific << "Energy = " << value << std::endl;
  }

};

#endif //OPTIMIZATION_SIMPLENEWTON_H
