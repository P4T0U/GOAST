// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Header for general stepsize control algorithms
 * \author Heeren
 */

#ifndef OPTIMIZATION_STEPSIZECONTROL_H
#define OPTIMIZATION_STEPSIZECONTROL_H

#include <goast/Core/Auxiliary.h>

enum TIMESTEP_CONTROLLER {
  CONST_TIMESTEP_CONTROL = -1,
  SIMPLE_TIMESTEP_CONTROL = 0,
  ARMIJO = 1,
  WOLFE = 2,
  NEWTON_OPTIMAL = 3,
  QUADRATIC = 4
};


//!==========================================================================================================
//! Pure virtual base class for simple, Armijo and Powell-Wolfe stepsize control based on inexact linesearch.
//! For an energy E, a position p and a descent direction d, the linesearch function is defined as f[t] := E[p + t*d].
//! In derived classes one has to provide:
//!  - the evaluation of the linesearch function at t=0 , i.e. evaluateLinesearchFnc()
//!  - the evaluation of the derivative of the linesearch function at t=0, i.e. evaluateLinesearchGrad()
//! \author Heeren
template<typename ConfiguratorType>
class StepsizeControlInterface {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;

  // The tolerence parameters.
  RealType _sigma;
  //! Additional parameter for getTimestepWidthWithPowellWolfeLineSearch.
  RealType _beta;

  mutable RealType _startTau;
  RealType _tauMin, _tauMax;
  bool _quadraticApprox;

  TIMESTEP_CONTROLLER _timestepController;

  const std::vector<int> *_bdryMask;

public:
  StepsizeControlInterface( RealType sigma, RealType beta, TIMESTEP_CONTROLLER timestepController,
                            RealType startTau = 1.0, RealType tauMin = 1e-12, RealType tauMax = 4. )
          : _sigma( sigma ), _beta( beta ), _startTau( startTau ), _tauMin( tauMin ), _tauMax( tauMax ),
            _quadraticApprox( true ), _bdryMask( NULL ), _timestepController( timestepController ) {}

  virtual ~StepsizeControlInterface() {}


  //Returns the scalar objective function evaluated at CurrentPosition.
  virtual RealType evaluateLinesearchFnc( const VectorType &CurrentPosition ) const = 0;

  //Returns the dot product of the energy derivative at position and the descent direction.
  virtual RealType evaluateLinesearchGrad( const VectorType &CurrentPosition, const VectorType &DescentDir ) const = 0;

  //Returns the scalar objective function evaluated at CurrentPosition+DescentDir*timestepWidth.
  virtual RealType evaluateLinesearchFnc( const VectorType &CurrentPosition, const VectorType &DescentDir, RealType timestepWidth ) const {
    return evaluateLinesearchFnc( CurrentPosition + timestepWidth * DescentDir );
  }

  virtual RealType evaluateLinesearchGrad( const VectorType &/*CurrentPosition*/, const VectorType &CurrentGradient, const VectorType &DescentDir ) const {
    return CurrentGradient.dot( DescentDir );
  }

  // set boundary mask
  void setBoundaryMask( const std::vector<int> &Mask ) {
    _bdryMask = &Mask;
  }

  void useQuadraticApproximation() {
    _quadraticApprox = true;
  }

  RealType getStartTau() const {
    return _startTau;
  }

  virtual RealType getStepsize( const VectorType &CurrentPosition, const VectorType &CurrentGrad, const VectorType &descentDir, RealType tau_before, RealType currEnergy = -1 ) const {
    switch ( _timestepController ) {
      case CONST_TIMESTEP_CONTROL:
        return getConstantTimestepWidth(  ); //
      case SIMPLE_TIMESTEP_CONTROL:
        return getTimestepWidthWithSimpleLineSearch( CurrentPosition, descentDir, tau_before, currEnergy ); //
      case ARMIJO:
        return getTimestepWidthWithArmijoLineSearch( CurrentPosition, CurrentGrad, descentDir, tau_before, currEnergy ); //
      case WOLFE:
        return getTimestepWidthWithPowellWolfeLineSearch( CurrentPosition, descentDir, currEnergy ); //
      default:
        throw BasicException( "StepsizeControlInterface::getStepsize: unknown stepsize method!" );
    }
  }

  //
  RealType getConstantTimestepWidth( ) const {
    RealType tau = 1e-5;
    return tau > _tauMin ? tau : 0.;
  }

  // Simple timestep width control, just ensures, that fnew < fnew
  RealType getTimestepWidthWithSimpleLineSearch( const VectorType &CurrentPosition, const VectorType &DescentDir, RealType OldTau, RealType CurrentEnergy ) const {

    RealType tau = std::min( std::max( 2 * OldTau, _tauMin ), _tauMax );

    RealType f = CurrentEnergy < 0. ? evaluateLinesearchFnc( CurrentPosition ) : CurrentEnergy;
    RealType fNew = evaluateLinesearchFnc( CurrentPosition, DescentDir, tau );

    if (std::isnan(fNew))
      fNew = std::numeric_limits<RealType>::infinity();

    while ((fNew >= f) && (tau >= _tauMin)) {
      tau = tau * 0.5;
      fNew = evaluateLinesearchFnc( CurrentPosition, DescentDir, tau );

      if (std::isnan(fNew))
        fNew = std::numeric_limits<RealType>::infinity();
    }

    return tau > _tauMin ? tau : 0.;
  }

  // Make use of local quadratic fitting to determine optimal stepsize
  RealType getTimestepWidthWithQuadraticLineSearch( const VectorType &CurrentPosition, const VectorType &CurrentGradient, const VectorType &DescentDir, RealType OldTau, RealType CurrentEnergy = -1. ) const {
    return zoomQuadratic( CurrentPosition, CurrentGradient, DescentDir, OldTau, CurrentEnergy );
  }

  // Armijo timestepping, i.e. E[p + tau*d] < E[p] + _sigma * tau * <DE[p], d>
  RealType getTimestepWidthWithArmijoLineSearch( const VectorType &CurrentPosition,
                                                 const VectorType &CurrentGradient,
                                                 const VectorType &DescentDir,
                                                 RealType OldTau,
                                                 RealType CurrentEnergy ) const {

    // Makes sure that tauMin <= startTau <= tauMax. For OldTau == 0, one descent direction component in the last step didn't lead to an energy reduction and the step is not done by setting tau to zero.
    RealType tau = std::min( std::max( OldTau, _tauMin ), _tauMax );

    const RealType f = CurrentEnergy < 0. ? evaluateLinesearchFnc( CurrentPosition ) : CurrentEnergy;
    RealType fNew = evaluateLinesearchFnc( CurrentPosition, DescentDir, tau );
    const RealType Df = evaluateLinesearchGrad( CurrentPosition, CurrentGradient, DescentDir );

    if ( std::isnan( fNew ))
      fNew = std::numeric_limits<RealType>::infinity();

    if ( Df >= 0 ) return 0;
    // use quadratic fit to initialize stepsize
    if ( _quadraticApprox ) {
      tau = zoomQuadratic( f, fNew, Df, tau );
      if ( tau < _tauMin )
        tau = std::min( std::max( OldTau, _tauMin ), _tauMax );
      else {
        tau = std::min( tau, _tauMax );
        fNew = evaluateLinesearchFnc( CurrentPosition, DescentDir, tau );
      }
      if (std::isnan(fNew))
        fNew = std::numeric_limits<RealType>::infinity();
    }

    RealType G = (fNew - f) / (tau * Df);

    //check
    if ( G >= _sigma && f >= fNew ) {
      //time step too small
      do {
        tau *= 2.;
        fNew = evaluateLinesearchFnc( CurrentPosition, DescentDir, tau );
        if (std::isnan(fNew))
          fNew = std::numeric_limits<RealType>::infinity();
        G = (fNew - f) / (tau * Df);
      } while ( G >= _sigma && f >= fNew && tau <= _tauMax );
      tau *= 0.5;
    }
    else {
      // time step too large
      do {
        if ( tau > _tauMin )
          tau *= 0.5;
        fNew = evaluateLinesearchFnc( CurrentPosition, DescentDir, tau );
        if (std::isnan(fNew))
          fNew = std::numeric_limits<RealType>::infinity();
        G = (fNew - f) / (Df * tau);
      } while (((G < _sigma || f < fNew)) && (tau > _tauMin));
    }

    return tau > _tauMin ? tau : 0.;
  }


  // CurrentPosition contains the currently optimal point,
  // DescentDir the linesearch direction,
  // Algorithm taken from: Nocedal, Wright: Numerical Optimization, Alg.3.5-3.6
  RealType getTimestepWidthWithPowellWolfeLineSearch( const VectorType &CurrentPosition,
                                                      const VectorType &DescentDir,
                                                      RealType CurrentEnergy = -1 ) const {

    RealType Tau = 1.; // important for Newton-based methods for superlinear convergence, since the first steplength is chosen which satisfies the Wolfe conditions
    RealType TauOld = 0.;
    RealType dE = -1. * DescentDir.squaredNorm();

    RealType energy = CurrentEnergy < 0. ? evaluateLinesearchFnc( CurrentPosition ) : CurrentEnergy;
    RealType energyBackup = energy;

    VectorType positionNew( CurrentPosition.size());

    bool stepSizeFound = false;
    do {
      positionNew = CurrentPosition;
      positionNew += Tau * DescentDir;
      RealType energyNew = evaluateLinesearchFnc( positionNew );

      if ( energyNew > std::min( energy + _sigma * Tau * dE, energyBackup )) {
        // Armijo condition violated...
        Tau = zoomWolfe( CurrentPosition, DescentDir, energy, dE, energyBackup, TauOld, Tau );
        stepSizeFound = true;

      } else {
        // Armijo condition fulfilled...
        RealType dENew = evaluateLinesearchGrad( positionNew, DescentDir );

        if ( std::abs( dENew ) <= -_beta * dE ) {
          // strong Wolfe condition (Armijo condition + curvature condition) fulfilled...
          energy = energyNew;
          stepSizeFound = true;
        } else {
          // curvature condition violated...
          if ( dENew >= 0. ) {
            // too long step
            Tau = zoomWolfe( CurrentPosition, DescentDir, energy, dE, energy, Tau, TauOld );
            stepSizeFound = true;
          } else {
            // too short step
            TauOld = Tau;
            energyBackup = energyNew;
            Tau *= 2.;
          }
        }
      }
    } while ( !stepSizeFound );

    return Tau > _tauMin ? Tau : 0.;
  }


protected:
  // auxiliary method for WOlfe condition
  RealType zoomWolfe( const VectorType &CurrentPosition,
                      const VectorType &DescentDir,
                      RealType &energy,
                      const RealType dE,
                      RealType ELo,
                      RealType TauLo,
                      RealType TauHi ) const {
    // either: TauLo satisfies the Armijo condition, but is too short for curvature condition, and TauHi > TauLo does not satisfy Armijo condition
    // or: both satisfy Armijo condition; TauLo is too long and TauHi too short for curvature condition
    RealType dENew, energyNew;
    VectorType positionNew( CurrentPosition.size() );
    RealType Tau = 0.;

    bool stepSizeFound = false;
    int counter = 0; // for non-termination due to rounding errors, we should fix the maximum number of iterations
    do {
      Tau = (TauLo + TauHi) / 2.;
      positionNew = CurrentPosition;
      positionNew += Tau * DescentDir;
      energyNew = evaluateLinesearchFnc( positionNew );

      if ( energyNew > std::min( energy + _sigma * Tau * dE, ELo ))
        // Armijo condition violated...
        TauHi = Tau;

      else {
        // Armijo condition fulfilled...
        dENew = evaluateLinesearchGrad( positionNew, DescentDir );

        if ( std::abs( dENew ) <= -_beta * dE ) {
          // strong Wolfe (Armijo + curvature) condition fulfilled...
          stepSizeFound = true;
        } else {
          // curvature condition violated...
          if ( dENew * (TauHi - TauLo) >= 0. )
            TauHi = TauLo;
          TauLo = Tau;
          ELo = energyNew;
        }
      }
    } while ( !stepSizeFound && (++counter < 30));

    if ( counter >= 30 ) {
      Tau = 0.;
    } else {
      energy = energyNew;
    }

    return Tau;
  }

  // Let f = E(pos), fNew = E(pos + tau*Dir), Df = DE[pos]*Dir = - |DE[pos]|^2 < 0, with tau being the old stepsize here.
  // Determine the quadratic function q such that q(0) = f, q(tau) = fNew and q'(0) = Df.
  // Set newTau such that q(newTau) is the unique minimum of q.
  RealType zoomQuadratic( const VectorType &CurrentPosition, const VectorType &CurrentGradient, const VectorType &DescentDir, RealType OldTau, RealType CurrentEnergy = -1. ) const {

    // Makes sure that tauMin <= startTau <= tauMax. For OldTau == 0, one descent direction component in the last step didn't lead to an energy reduction and the step is not done by setting tau to zero.
    RealType tau = std::min( std::max( OldTau, _tauMin ), _tauMax );

    RealType f  = CurrentEnergy < 0. ? evaluateLinesearchFnc( CurrentPosition ) : CurrentEnergy;
    RealType fNew = evaluateLinesearchFnc( CurrentPosition, DescentDir, tau );

    RealType Df = evaluateLinesearchGrad( CurrentPosition, CurrentGradient, DescentDir );
    if ( Df >= 0 ) return 0;

    return zoomQuadratic( f, fNew, Df, tau );
  }

  // see above
  RealType zoomQuadratic( const RealType &f, const RealType &fNew, const RealType &Df, const RealType &tau ) const {
    RealType denom = -2. * (fNew - f - tau * Df);
    return (denom > 0) ? 0. : Df * tau * tau / denom;
  }

};


//!==========================================================================================================
//! Class for simple, Armijo and Powell-Wolfe stepsize control
//! \author Heeren
template<typename ConfiguratorType>
class StepsizeControl : public StepsizeControlInterface<ConfiguratorType> {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;

  const BaseOp<VectorType, RealType> &_E;
  const BaseOp<VectorType, VectorType> &_DE;

public:
  StepsizeControl( const BaseOp<VectorType, RealType> &E,
                   const BaseOp<VectorType, VectorType> &DE,
                   TIMESTEP_CONTROLLER timestepController,
                   RealType sigma = 0.1,
                   RealType beta = 0.9,
                   RealType startTau = 1.0,
                   RealType tauMin = 1e-12,
                   RealType tauMax = 4. ) : StepsizeControlInterface<ConfiguratorType>( sigma, beta, timestepController,
                                                                                        startTau, tauMin, tauMax ),
                                            _E( E ), _DE( DE ) {}

  //Returns the scalar objective function evaluated at CurrentPosition.
  RealType evaluateLinesearchFnc( const VectorType &CurrentPosition ) const {
    return _E( CurrentPosition );
  }

  //Returns the dot product of the energy derivative and the descent direction.
  RealType evaluateLinesearchGrad( const VectorType &Position, const VectorType &DescentDir ) const {
    VectorType tmp( DescentDir );
    _DE.apply( Position, tmp );
    if ( this->_bdryMask )
      applyMaskToVector( *this->_bdryMask, tmp );
    return tmp.dot( DescentDir );
  }

};


#endif //REDUCEDBASIS_STEPSIZECONTROL_H
