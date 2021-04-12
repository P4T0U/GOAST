// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Header for quasi-Newton methods
 * \author Heeren
 *
 * \todo Convert to new interface
 * \todo Documentation!
 */

#ifndef OPTIMIZATION_QUASINEWTON_H
#define OPTIMIZATION_QUASINEWTON_H

#include "optInterface.h"
#include "optParameters.h"
#include "stepsizeControl.h"

/**
 * \brief This class implements a quasi-Newton method with BFGS-Update of the matrix for minimizing a scalar function.
 * \tparam ConfiguratorType Container with data types
 * \author Heeren
 */
template<typename ConfiguratorType>
class QuasiNewtonBFGS {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;

  const int _maxIterations;
  const RealType _stopEpsilon;
  const int _reset;
  QUIET_MODE _quietMode;

  // BFGS stored information
  mutable int _counter;
  mutable std::vector<VectorType> _y;
  mutable std::vector<VectorType> _dx;
  mutable std::vector<RealType> _dx_y;

  const BaseOp<VectorType, RealType> &_E;
  const BaseOp<VectorType, VectorType> &_DE;
  TIMESTEP_CONTROLLER _timeStepController;
  StepsizeControl<ConfiguratorType> _stepsizeControl;
  const std::vector<int> *_bdryMask;
  const BaseOp<VectorType, VectorType>* _inverseHessianOperator;

public:
  QuasiNewtonBFGS( const BaseOp<VectorType, RealType> &E,
                   const BaseOp<VectorType, VectorType> &DE,
                   int MaxIterations,
                   RealType StopEpsilon,
                   TIMESTEP_CONTROLLER TimestepController,
                   int Reset,
                   bool quiet,
                   RealType sigma = 0.1,
                   RealType tauMin = 1.e-6,
                   RealType tauMax = 4. ) :
          _maxIterations( MaxIterations ),
          _stopEpsilon( StopEpsilon ),
          _reset( Reset ),
          _quietMode( quiet ? SUPERQUIET : SHOW_ALL ),
          _counter( 0 ),
          _dx_y( Reset ),
          _E( E ),
          _DE( DE ),
          _timeStepController( static_cast<TIMESTEP_CONTROLLER>( TimestepController )),
          _stepsizeControl( E, DE, _timeStepController, sigma, 0.9, 1., tauMin, tauMax ),
          _bdryMask( nullptr ),
          _inverseHessianOperator(nullptr){}
          
  QuasiNewtonBFGS( const BaseOp<VectorType, RealType> &E,
                   const BaseOp<VectorType, VectorType> &DE,
                   int MaxIterations = 50,
                   RealType StopEpsilon = 1.e-6,
                   TIMESTEP_CONTROLLER TimestepController = ARMIJO,
                   int Reset = 50,
                   QUIET_MODE quietMode = SUPERQUIET,
                   RealType sigma = 0.1,
                   RealType tauMin = 1.e-6,
                   RealType tauMax = 4. ) :
          _maxIterations( MaxIterations ),
          _stopEpsilon( StopEpsilon ),
          _reset( Reset ),
          _quietMode( quietMode ),
          _counter( 0 ),
          _dx_y( Reset ),
          _E( E ),
          _DE( DE ),
          _timeStepController( static_cast<TIMESTEP_CONTROLLER>( TimestepController )),
          _stepsizeControl( E, DE, _timeStepController, sigma, 0.9, 1., tauMin, tauMax ),
          _bdryMask( nullptr ),
          _inverseHessianOperator(nullptr){}

  QuasiNewtonBFGS( const BaseOp<VectorType, RealType> &E,
                   const BaseOp<VectorType, VectorType> &DE,
                   const OptimizationParameters<ConfiguratorType> &optPars ) :
          _maxIterations( optPars.getBFGSIterations()),
          _stopEpsilon( optPars.getStopEpsilon()),
          _reset( optPars.getBFGSReset()),
          _quietMode( optPars.getQuietMode()),
          _counter( 0 ),
          _dx_y( _reset ),
          _E( E ),
          _DE( DE ),
          _timeStepController( static_cast<TIMESTEP_CONTROLLER>( optPars.getBFGSTimeStepping())),
          _stepsizeControl( E, DE, _timeStepController, optPars.getSigma(), optPars.getBeta(), optPars.getStartTau(),
                            optPars.getTauMin(), optPars.getTauMax()),
          _bdryMask( nullptr ),
          _inverseHessianOperator(nullptr){}


  void setBoundaryMask( const std::vector<int> &Mask ) {
    _bdryMask = &Mask;
    _stepsizeControl.setBoundaryMask( Mask );
  }

  void setInverseHessianOperator( const BaseOp<VectorType, VectorType>& InvHessianOp ) {
    _inverseHessianOperator = &InvHessianOp;
  }

  void solve( const VectorType &Arg, VectorType &Dest ) const {

    if ( _maxIterations == 0 )
      return;


    RealType energy;
    _E.apply( Arg, energy );
    if ( _quietMode == SHOW_ALL ) {
      std::cout << "=======================================================================" << std::endl;
      std::cout << "Start BFGS with " << _maxIterations << " iterations and eps = " << _stopEpsilon << "." << std::endl;
      std::cout << "Initial energy: " << energy << std::endl;
      std::cout << "=======================================================================" << std::endl;
    }

    Dest = Arg;

    VectorType Dx( Dest.size());
    VectorType f( Dest.size());
    VectorType Df( Dest.size());
    VectorType tmp( Dest.size());
    VectorType descentDir( Dest.size());

    // compute f
    _DE.apply( Dest, f );
    // Dirichlet boundary conditions?
    if ( _bdryMask )
      applyMaskToVector( *_bdryMask, f );

    RealType FNorm = 1e+15;
    RealType tau = 1.0;
    int iterations = 0;
    bool forcedReset = false;

    while ( FNorm > _stopEpsilon && (iterations < _maxIterations) && tau > 0. ) {
      auto t_start = std::chrono::high_resolution_clock::now();
      // Quasi-Newton-iteration given by x_{k+1} = x_k - tau_k B_k^{-1} f(x_k)
      // iteration number k (i.e. find x_{k+1} from current approximation x_k)

      // save x_k and f(x_k) for the computation of Dx = x_{k+1}-x_k and Df = f(x_{k+1})-f(x_k)
      Dx = Dest;
      Dx *= -1.;
      Df = f;
      Df *= -1.;

      applyInverse( f, descentDir, _counter );
      descentDir *= -1.;

      // get tau
      tau = _stepsizeControl.getStepsize( Dest, f, descentDir, tau, energy );

      if ((tau == 0) && !forcedReset ) {
        reset();
        tau = 1.;
        forcedReset = true;
        continue;
      }
      forcedReset = false;

      // update position
      Dest += tau * descentDir;

      _DE.apply( Dest, tmp );
      if ( _bdryMask )
        applyMaskToVector( *_bdryMask, tmp );
      FNorm = tmp.norm();

      Dx += Dest;
      Df += tmp;

      // update of B, B_{k+1} = B_k - \frac{B_k Dx Dx^T B_k}{Dx\cdot (B_k Dx)}+\frac{Df Df^T}{Df \cdot Dx}
      update( Dx, Df );
      iterations++;

      _E.apply( Dest, energy );
      f = tmp;
      auto t_end = std::chrono::high_resolution_clock::now();
      if ( _quietMode == SHOW_ALL )
        std::cout << std::scientific << "step = " << iterations << " , stepsize = " << tau
                  << ", energy = " << energy << ", error = " << FNorm
                  << std::fixed << ", time = " << std::chrono::duration<double, std::milli >( t_end - t_start ).count()
                  << "ms" << std::endl;

    } // end while

    if ( _quietMode != SUPERQUIET ) {
      _E.apply( Dest, energy );
      std::cout << "=======================================================================" << std::endl;
      std::cout << "Finished BFGS after " << iterations << " steps." << std::endl;
      std::cout << "Final stepsize = " << tau << ", energy = " << energy << ", error = " << FNorm << std::endl;
      std::cout << "=======================================================================" << std::endl;
    }
  }

protected:
  void update( const VectorType &DX, const VectorType &Y ) const {
    if ( _counter == _reset )
      reset();
    _counter++;

    _y.push_back(Y);
    _dx.push_back(DX);
    _dx_y[_counter - 1] = DX.dot( Y );
  }

  void reset() const {
    _y.clear();
    _dx.clear();
    _counter = 0;
  }

  // recursive function f_k( y ) = A_k * f_{k-1}(A_k^T y) + B_k y,
  // where A_k = I - \alpha_k dx_k y_k^T  and B_k = I - \alpha_k dx_k dx_k^T  with \alpha_k = 1 / (y_k^T dx_k)
  void applyInverse( const VectorType &Arg, VectorType &Dest, int iteration ) const {

    if( iteration == 0 ){
      Dest = Arg;
      if( _inverseHessianOperator )
        _inverseHessianOperator->apply( Arg, Dest );
      return;
    }

    // compute A_k^T y
    int prevIter = iteration - 1;
    RealType dotProd1 = ((_dx[prevIter]).dot( Arg )) / _dx_y[prevIter];
    VectorType prevArg = Arg - dotProd1 * _y[prevIter];

    // compute f_{k-1}(A_k^T y)
    applyInverse( prevArg, Dest, prevIter );

    // compute Dest =  A_k * f_{k-1}(A_k^T y) + B_k y
    RealType dotProdDiff = ( (_dx[prevIter]).dot( Arg ) - (_y[prevIter]).dot( Dest ) ) / _dx_y[prevIter];
    Dest += dotProdDiff * _dx[prevIter];
  }

};

#endif //OPTIMIZATION_QUASINEWTON_H
