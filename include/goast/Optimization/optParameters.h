// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Header for paremeter container for optimization methods
 * \author Heeren
 *
 * \todo Convert to new interface
 * \todo Documentation!
 */


#ifndef OPTIMIZATION_PARAMETERS_H
#define OPTIMIZATION_PARAMETERS_H

#include "optInterface.h"
#include "stepsizeControl.h"
#include <goast/Core/LinearSolver.h>
#include <goast/Core/ParameterParser.h>

//!
enum QUIET_MODE {
    SUPERQUIET = 0,
    SHOW_TERMINATION_INFO = 1,
    SHOW_ALL = 2,
    SHOW_ONLY_IF_FAILED = 3 
};    
    
//!     
template<typename ConfiguratorType>
class OptimizationParameters {

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;

  ParameterParser const *_pparser;
  mutable int _level;

  RealType _stopEpsilonCoarse, _stopEpsilonFine;

  int _gradItersFine, _gradItersCoarse;
  int _bfgsItersFine, _bfgsItersCoarse;
  int _newtonItersFine, _newtonItersCoarse;
  int _bfgsResetCoarse, _bfgsResetFine;

  TIMESTEP_CONTROLLER _gradTimeStepping, _bfgsTimeStepping, _newtonTimeStepping;
  RealType _startTau, _tauMin, _tauMax;
  RealType _sigma, _beta;

  int _solverType, _initializationScheme;
  bool _verboseFine, _verboseCoarse;
  QUIET_MODE _quietModeFine, _quietModeCoarse;

public:
  OptimizationParameters() : _pparser( NULL ), _level( 1 ), _stopEpsilonCoarse( 1e-8 ), _stopEpsilonFine( 1e-8 ),
                             _bfgsResetCoarse( 50 ), _bfgsResetFine( 50 ), _gradTimeStepping( ARMIJO ),
                             _bfgsTimeStepping( ARMIJO ), _newtonTimeStepping( NEWTON_OPTIMAL ), _startTau( 1. ),
                             _tauMin( 1e-8 ), _tauMax( 4. ), _sigma( 0.1 ), _beta( 0.9 ),
                             _solverType( UMFPACK_LU_FACT ), _initializationScheme(-1), 
                             _verboseFine( false ), _verboseCoarse( false ),
                             _quietModeFine(SUPERQUIET), _quietModeCoarse(SUPERQUIET){
    setZero();
  }

  OptimizationParameters( const ParameterParser &pparser ) : _pparser( &pparser ), _level( 1 ),
                                                             _stopEpsilonCoarse( 1e-8 ), _stopEpsilonFine( 1e-8 ),
                                                             _gradTimeStepping( ARMIJO ), _bfgsTimeStepping( ARMIJO ),
                                                             _newtonTimeStepping( NEWTON_OPTIMAL ), _startTau( 1. ),
                                                             _tauMin( 1e-8 ), _tauMax( 4. ), _sigma( 0.1 ),
                                                             _beta( 0.9 ), _solverType( UMFPACK_LU_FACT ), _initializationScheme(-1),
                                                             _verboseFine( false ), _verboseCoarse( false ),
                                                             _quietModeFine(SUPERQUIET), _quietModeCoarse(SUPERQUIET){
    setZero();

    _gradItersCoarse = pparser.getIntOrDefault( "gradItersCoarse", 0 );
    _gradItersFine = pparser.getIntOrDefault( "gradItersFine", 0 );

    _bfgsItersCoarse = pparser.getIntOrDefault( "bfgsItersCoarse", 0 );
    _bfgsItersFine = pparser.getIntOrDefault( "bfgsItersFine", 0 );

    _newtonItersCoarse = pparser.getIntOrDefault( "newtonItersCoarse", 0 );
    _newtonItersFine = pparser.getIntOrDefault( "newtonItersFine", 0 );

    _bfgsResetCoarse = pparser.getIntOrDefault( "bfgsResetCoarse", 50 );
    _bfgsResetFine = pparser.getIntOrDefault( "bfgsResetFine", 50 );

    _sigma = pparser.getDoubleOrDefault( "sigma", 0.1 );
    _tauMin = pparser.getDoubleOrDefault( "tauMin", 1e-8 );
    _tauMax = pparser.getDoubleOrDefault( "tauMax", 4. );

    _initializationScheme = pparser.getIntOrDefault("initializationScheme", -1 );
  }

  void setLevel( int Level ) const {
    _level = Level;
  }

  RealType getStopEpsilon() const {
    return (_level == 0) ? _stopEpsilonCoarse : _stopEpsilonFine;
  }

  int getGradientIterations() const {
    return (_level == 0) ? _gradItersCoarse : _gradItersFine;
  }

  int getBFGSIterations() const {
    return (_level == 0) ? _bfgsItersCoarse : _bfgsItersFine;
  }

  int getNewtonIterations() const {
    return (_level == 0) ? _newtonItersCoarse : _newtonItersFine;
  }

  int getBFGSReset() const {
    return (_level == 0) ? _bfgsResetCoarse : _bfgsResetFine;
  }

  int getSolverType() const {
    return _solverType;
  }

  void setGradientIterations( int Iters, int level = 1 ) {
    if ( level == 0 )
      _gradItersCoarse = Iters;
    else
      _gradItersFine = Iters;
  }

  void setBFGSIterations( int Iters, int level = 1 ) {
    if ( level == 0 )
      _bfgsItersCoarse = Iters;
    else
      _bfgsItersFine = Iters;
  }

  void setNewtonIterations( int Iters, int level = 1 ) {
    if ( level == 0 )
      _newtonItersCoarse = Iters;
    else
      _newtonItersFine = Iters;
  }

  void setEpsilon( RealType Eps, int level = 1 ){
    if ( level == 0 )
      _stopEpsilonCoarse = Eps;
    else
      _stopEpsilonFine = Eps;
  }

  void setQuiet( int level = 1 ) {
    if ( level == 0 )
      _quietModeCoarse = SUPERQUIET;
    else
      _quietModeFine = SUPERQUIET;
  }
  
  void setQuietMode( QUIET_MODE quietMode, int level = 1 ) {
    if ( level == 0 )
      _quietModeCoarse = quietMode;
    else
      _quietModeFine = quietMode;
  }
  
  QUIET_MODE getQuietMode(  ) const {
      return ( _level == 0 ) ? _quietModeCoarse : _quietModeFine;
  }

  void setSigma( RealType sigma ) {
    _sigma = sigma;
  }

  int getGradTimeStepping() const {
    return _gradTimeStepping;
  }

  int getBFGSTimeStepping() const {
    return _bfgsTimeStepping;
  }

  int getNewtonTimeStepping() const {
    return _newtonTimeStepping;
  }

  RealType getTauMin() const {
    return _tauMin;
  }

  RealType getTauMax() const {
    return _tauMax;
  }

  RealType getStartTau() const {
    return _startTau;
  }

  RealType getSigma() const {
    return _sigma;
  }

  RealType getBeta() const {
    return _beta;
  }

  int getInitializationScheme() const {
    return _initializationScheme;
  }

  void setInitializationScheme( int initScheme ) {
    _initializationScheme = initScheme;
  }

  void setSolverType(LINEAR_SOLVER_TYPE t) {
    _solverType = t;
  }

protected:
  void setZero() {
    _gradItersFine = 0;
    _gradItersCoarse = 0;
    _bfgsItersFine = 0;
    _bfgsItersCoarse = 0;
    _newtonItersFine = 0;
    _newtonItersCoarse = 0;
  }

};

#endif //REDUCEDBASIS_OPTPARAMETERS_H
