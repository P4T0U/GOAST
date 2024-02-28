// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DERIVATIVETESTS_HH
#define DERIVATIVETESTS_HH

//== INCLUDES =================================================================
#include <iostream>
#include <ctime>
#include <iomanip>
#include <fstream>
#include <random>

#include "Auxiliary.h"
#include "IO.h"


//! \brief Derivative test for functionals $F: \R^n \to \R$ with $\nabla F = DF^T \in \R^n$
//! \author Heeren
template<typename ConfiguratorType>
class ScalarValuedDerivativeTester {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;

  const BaseOp<VectorType, RealType> &_F;
  const BaseOp<VectorType, VectorType> &_DF;
  const RealType _stepSize;
  const int _numSteps;

public:
  ScalarValuedDerivativeTester( const BaseOp<VectorType, RealType> &F,
                                const BaseOp<VectorType, VectorType> &DF,
                                const RealType stepSize,
                                int numSteps = 50 )
          : _F( F ),
            _DF( DF ),
            _stepSize( stepSize ),
            _numSteps( numSteps ) {}

  // Plot f(t) := F[p + t*d] and g(t) = F[p] + t * <\nabla F[p], d>  for a test point p \in \R^n and a direction d \in \R^n
  // If \nabla F[p] is correct, g(t) should be tangential to f(t) at t=0;
  void plotSingleDirection( const VectorType &testPoint, const VectorType &testDirection,
                            const std::string saveName ) const {
    VectorType timeSteps( _numSteps ), energies( _numSteps ), derivs( _numSteps );

    VectorType gradient( testPoint.size());
    _DF.apply( testPoint, gradient );
    RealType deriv = gradient.dot( testDirection );

    RealType energy, energyShifted;
    _F.apply( testPoint, energy );
    for ( int i = 0; i < _numSteps; i++ ) {
      timeSteps[i] = ( i - ( _numSteps / 2 )) * _stepSize;

      VectorType shiftedPoint = testPoint + timeSteps[i] * testDirection;
      _F.apply( shiftedPoint, energyShifted );
      energies[i] = energyShifted;

      derivs[i] = energy + timeSteps[i] * deriv;
    }

    generatePNG( timeSteps, energies, derivs, saveName );
  }

  // test direction is d = e_i
  void plotSingleDirection( const VectorType &testPoint, int i, const std::string saveName ) const {
    VectorType testDirection( testPoint.size());
    testDirection.setZero();
    testDirection[i] = 1.;
    plotSingleDirection( testPoint, testDirection, saveName );
  }

  // Plot f(t) := F[p + t*d] and g(t) = F[p] + t * <\nabla F[p], d>  for a test point p \in \R^n 
  // and directions d = e_i \in \R^n, i = 1, ..., n, where e_i is the canonical basis vector of \R^n
  // If \nabla F[p] is correct, g(t) should be tangential to f(t) at t=0;
  void plotAllDirections( const VectorType &testPoint, const std::string saveNameStem ) const {
    VectorType timeSteps( _numSteps ), energies( _numSteps ), derivs( _numSteps );

    int numDofs = testPoint.size();
    VectorType gradient( testPoint.size()), testDirection( testPoint.size());
    _DF.apply( testPoint, gradient );

    RealType energy, energyShifted;
    _F.apply( testPoint, energy );
    for ( int j = 0; j < numDofs; j++ ) {
      std::ostringstream saveName;
      saveName << saveNameStem << j << ".png";
      testDirection.setZero();
      testDirection[j] = 1.;
      for ( int i = 0; i < _numSteps; i++ ) {
        timeSteps[i] = ( i - ( _numSteps / 2 )) * _stepSize;
        VectorType shiftedPoint = testPoint + timeSteps[i] * testDirection;
        _F.apply( shiftedPoint, energyShifted );
        energies[i] = energyShifted;
        derivs[i] = energy + timeSteps[i] * gradient[j];
      }

      // Shifting derivatives such that middle points (i.e. t=0) match
      RealType diff = energies[_numSteps / 2] - derivs[_numSteps / 2];
      for ( int i = 0; i < _numSteps; i++ )
        derivs[i] += diff;

      generatePNG( timeSteps, energies, derivs, saveName.str());
    }
  }

  // Plot f(t) := F[p + t*d_i] and g(t) = F[p] + t * <\nabla F[p], d_i>  for a test point p \in \R^n 
  // and random directions d_i = e_j \in \R^n, i = 1, ..., n, where e_j is the canonical basis vector of \R^n
  // If \nabla F[p] is correct, g(t) should be tangential to f(t) at t=0; 
  void plotRandomDirections( const VectorType &testPoint, int numOfRandomDir, const std::string saveNameStem ) const {

    VectorType timeSteps( _numSteps ), energies( _numSteps ), derivs( _numSteps );

    std::srand( std::time( 0 )); // use current time as seed for random generator
    VectorType randApproxDeriv( numOfRandomDir );
    std::vector<int> directionIndices( numOfRandomDir );

    int numDofs = testPoint.size();
    VectorType gradient( numDofs ), testDirection( numDofs );
    _DF.apply( testPoint, gradient );

    RealType energy, energyShifted;
    _F.apply( testPoint, energy );
    for ( int j = 0; j < numOfRandomDir; j++ ) {
      directionIndices[j] = std::rand() % numDofs;
      std::ostringstream saveName;
      saveName << saveNameStem << "_" << directionIndices[j] << ".png";
      testDirection.setZero();
      testDirection[directionIndices[j]] = 1.;
      for ( int i = 0; i < _numSteps; i++ ) {
        timeSteps[i] = ( i - ( _numSteps / 2 )) * _stepSize;
        VectorType shiftedPoint = testPoint + timeSteps[i] * testDirection;
        _F.apply( shiftedPoint, energyShifted );
        energies[i] = energyShifted;
        derivs[i] = energy + timeSteps[i] * gradient[directionIndices[j]];
      }

      // Shifting derivatives such that middle points (i.e. t=0) match
      RealType diff = energies[_numSteps / 2] - derivs[_numSteps / 2];
      for ( int i = 0; i < _numSteps; i++ )
        derivs[i] += diff;

      generatePNG( timeSteps, energies, derivs, saveName.str());
    }
  }

  // Compute \nabla_tau F[p] * d := (F[p + tau*d] - F[p]) / tau for a test point p and a test direction d
  RealType computeDiffQuotient( const VectorType &testPoint, const VectorType &testDirection, RealType tau,
                                bool centralDiff = false ) const {

    VectorType shiftedPoint = centralDiff ? ( testPoint - tau / 2. * testDirection )
                                          : testPoint; // - _stepSize/2. * testDirection;

    // evalualte J(m)
    RealType energy;
    _F.apply( shiftedPoint, energy );
    //evaluate J(m + h testDirection )
    RealType energyShifted;
    shiftedPoint += tau * testDirection; //  / 2.
    _F.apply( shiftedPoint, energyShifted );
    return ( energyShifted - energy ) / tau;
  }


  // Compute relative error of \nabla_tau F[p] * d (see above) and \nabla F[p] * d for a test point p and a test direction d
  void testSingleDirection( const VectorType &testPoint, const VectorType &testDirection, RealType tau = -1. ) const {
    //evaluate DJ(m)
    RealType stepSize = ( tau < 0 ) ? _stepSize : tau;
    VectorType derivative( testPoint.size());
    _DF.apply( testPoint, derivative );
    RealType gateaux = derivative.dot( testDirection );
    std::cerr << std::abs(( computeDiffQuotient( testPoint, testDirection, stepSize ) - gateaux ) / gateaux )
              << std::endl;
  }

  // Compute relative error of \nabla_tau F[p] * d (see above) and \nabla F[p] * d for a test point p and a test direction d
  void testSingleDirection( const VectorType &testPoint, int i, RealType tau = -1. ) const {
    VectorType testDirection( testPoint.size());
    testDirection.setZero();
    testDirection[i] = 1.;
    testSingleDirection( testPoint, testDirection, tau );
  }

  // Compute relative errors of \nabla_tau F[p] * d_i (see above) and \nabla F[p] * d_i for a test point p and all test direction d_i = e_i
  void testAllDirections( const VectorType &testPoint, const bool detailedOutput = false, RealType tau = -1. ) {
    int numDofs = testPoint.size();
    RealType stepSize = ( tau < 0 ) ? _stepSize : tau;
    //evaluate DJ(m)
    VectorType derivative( numDofs );
    _DF.apply( testPoint, derivative );
    // diffQuotient
    VectorType error( numDofs ), testDirection( numDofs ), diffQuot( numDofs );
    for ( unsigned i = 0; i < numDofs; ++i ) {
      testDirection.setZero();
      testDirection[i] = 1.;
      diffQuot[i] = computeDiffQuotient( testPoint, testDirection, stepSize );
      error[i] = std::abs(( diffQuot[i] - derivative[i] ) / derivative[i] );
    }

    if ( detailedOutput ) {
      printVector( diffQuot, 10 );
      std::cerr << std::endl;
      printVector( derivative, 10 );
      std::cerr << std::endl;
    }
    printVector( error, 10 );
  }

  // Compute relative errors of \nabla_tau F[p] * d_i (see above) and \nabla F[p] * d_i for a test point p and random test direction d_i = e_j
  void testRandomDirections( const VectorType &testPoint, unsigned numOfRandomDir, const bool detailedOutput = false,
                             RealType tau = -1. ) {

    int numDofs = testPoint.size();
    std::srand( std::time( 0 )); // use current time as seed for random generator
    RealType stepSize = ( tau < 0 ) ? _stepSize : tau;

    //evaluate DJ(m)
    VectorType derivative( numDofs ), testDirection( numDofs );
    _DF.apply( testPoint, derivative );

    // diffQuotient
    VectorType error( numOfRandomDir ), diffQuot( numOfRandomDir ), resDeriv( numOfRandomDir );
    std::vector<int> directionIndices( numOfRandomDir );
    for ( unsigned i = 0; i < numOfRandomDir; ++i ) {
      testDirection.setZero();
      directionIndices[i] = std::rand() % numDofs;
      testDirection[directionIndices[i]] = 1.;
      resDeriv[i] = derivative[directionIndices[i]];
      diffQuot[i] = computeDiffQuotient( testPoint, testDirection, stepSize );
      error[i] = std::abs(( diffQuot[i] - resDeriv[i] ) / resDeriv[i] );
    }
    if ( detailedOutput ) {
      printVector( directionIndices, 10 );
      printVector( diffQuot, 10 );
      printVector( resDeriv, 10 );
    }
    printVector( error, 10 );
  }

};


/**
 * \brief Derivative test for functionals $F: \R^n \to \R^m$ with $DF \in \R^{m,n}$.
 * \author Heeren
 * Here, $n$ is the number of degrees of freedom (dofs) and $m$ is the dimension of the range.
 * The dim. of the range can be set in the constructor, otherwise m = n.
 */
template<typename ConfiguratorType>
class VectorValuedDerivativeTester {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;

  const BaseOp<VectorType, VectorType> &_F;
  const BaseOp<VectorType, MatrixType> &_DF;

  const int _dimOfRange;
  RealType _stepSize;
  int _numSteps;

public:
  VectorValuedDerivativeTester( const BaseOp<VectorType, VectorType> &F,
                                const BaseOp<VectorType, MatrixType> &DF,
                                const RealType stepSize,
                                int dimOfRange = -1 )
          : _F( F ),
            _DF( DF ),
            _dimOfRange( dimOfRange ),
            _stepSize( stepSize ),
            _numSteps( 50 ) {}

  //
  void plotRandomDirections( const VectorType &testPoint, int numOfRandomDir, const std::string saveNameStem,
                             const bool testNonNegEntriesOfHessianOnly = true ) const {

    int numDofs = testPoint.size();
    int dimRange = _dimOfRange < 0 ? testPoint.size() : _dimOfRange;

    MatrixType JacobiMatrix( dimRange, numDofs );
    _DF.apply( testPoint, JacobiMatrix );
    MatrixType JacobiMatrixTrans = JacobiMatrix.transpose();

    VectorType timeSteps( _numSteps ), energies( _numSteps ), derivs( _numSteps );
    VectorType outerDirection( dimRange ), innerDirection( numDofs ), gradientWRTOuterDirection( numDofs );

    std::srand( std::time( 0 )); // use current time as seed for random generator
    //VectorType randApproxDeriv ( numOfRandomDir );    
    std::vector<int> directionOuterIndices( numOfRandomDir ), directionInnerIndices( numOfRandomDir );

    for ( int i = 0; i < numOfRandomDir; i++ ) {

      outerDirection.setZero();
      innerDirection.setZero();

      // outer index means, we test the ith component of  F = (f_1, ..., f_m), for 1 <= i <= m 
      directionOuterIndices[i] = std::rand() % dimRange;
      outerDirection[directionOuterIndices[i]] = 1.;
      // inner index means, we test the jth entry of the component f_i    
      directionInnerIndices[i] = std::rand() % numDofs;
      innerDirection[directionInnerIndices[i]] = 1.;

      // F = (f_1, ..., f_m)  =>  \nabla f_i = DF^T e_i (gradient of f_i is ith row of DF)
      RealType gateauxDerivative = JacobiMatrixTrans.coeff( directionInnerIndices[i], directionOuterIndices[i] );
      //gradientWRTOuterDirection  = JacobiMatrixTrans * outerDirection;
      //RealType gateauxDerivative = gradientWRTOuterDirection.dot( innerDirection );    

      std::ostringstream saveName;
      saveName << saveNameStem << "_" << directionOuterIndices[i] << "_" << directionInnerIndices[i] << ".png";

      if ( testNonNegEntriesOfHessianOnly && ( std::abs( gateauxDerivative ) < 1.e-10 )) {
        i--;
        continue;
      }

      VectorType tempVec( dimRange );
      _F.apply( testPoint, tempVec );
      RealType initialEnergy = tempVec.dot( outerDirection );

      for ( int j = 0; j < _numSteps; j++ ) {
        timeSteps[j] = ( j - ( _numSteps / 2 )) * _stepSize;
        VectorType shiftedPoint = testPoint + timeSteps[j] * innerDirection;
        _F.apply( shiftedPoint, tempVec );
        energies[j] = tempVec.dot( outerDirection );
        derivs[j] = initialEnergy + timeSteps[j] * gateauxDerivative;
      }

      generatePNG( timeSteps, energies, derivs, saveName.str());

    }

  }

  //
  void plotAllDirections( const VectorType &testPoint, const std::string saveNameStem ) const {

    int numDofs = testPoint.size();
    int dimRange = _dimOfRange < 0 ? testPoint.size() : _dimOfRange;

    MatrixType JacobiMatrix( dimRange, numDofs );
    _DF.apply( testPoint, JacobiMatrix );
    MatrixType JacobiMatrixTrans = JacobiMatrix.transpose();

    VectorType timeSteps( _numSteps ), energies( _numSteps ), derivs( _numSteps );
    VectorType outerDirection( dimRange ), innerDirection( numDofs ), gradientWRTOuterDirection( numDofs );

    for ( int k = 0; k < dimRange; k++ ) {

      outerDirection.setZero();
      outerDirection[k] = 1.;

      for ( int i = 0; i < numDofs; i++ ) {
        innerDirection.setZero();
        innerDirection[i] = 1.;

        // F = (f_1, ..., f_n)  =>  \nabla f_i = DF^T e_i
        gradientWRTOuterDirection = JacobiMatrixTrans * outerDirection;
        RealType gateauxDerivative = gradientWRTOuterDirection.dot( innerDirection );

        VectorType tempVec( dimRange );
        _F.apply( testPoint, tempVec );
        RealType initialEnergy = tempVec.dot( outerDirection );

        for ( int j = 0; j < _numSteps; j++ ) {
          timeSteps[j] = ( j - ( _numSteps / 2 )) * _stepSize;
          VectorType shiftedPoint = testPoint + timeSteps[j] * innerDirection;
          _F.apply( shiftedPoint, tempVec );
          energies[j] = tempVec.dot( outerDirection );
          derivs[j] = initialEnergy + timeSteps[j] * gateauxDerivative;
        }

        std::ostringstream saveName;
        saveName << saveNameStem << "_" << k << "_" << i << ".png";
        generatePNG( timeSteps, energies, derivs, saveName.str());
      }

    }

  }

  //
  void
  plotSingleDirection( const VectorType &testPoint, const VectorType &innerDirection, const VectorType &outerDirection,
                       const std::string saveName ) const {

    int numDofs = testPoint.size();
    int dimRange = _dimOfRange < 0 ? testPoint.size() : _dimOfRange;

    if ( outerDirection.size() != dimRange )
      throw BasicException( "VectorValuedDerivativeTester::plotSingleDirections: outerDir has wrong size!" );

    if ( innerDirection.size() != numDofs )
      throw BasicException( "VectorValuedDerivativeTester::plotSingleDirections: innerDir has wrong size!" );

    MatrixType JacobiMatrix( dimRange, numDofs );
    _DF.apply( testPoint, JacobiMatrix );
    MatrixType JacobiMatrixTrans = JacobiMatrix.transpose();

    VectorType timeSteps( _numSteps ), energies( _numSteps ), derivs( _numSteps );

    VectorType gradientWRTOuterDirection = JacobiMatrixTrans * outerDirection;
    RealType gateauxDerivative = gradientWRTOuterDirection.dot( innerDirection );

    VectorType tempVec( dimRange );
    _F.apply( testPoint, tempVec );
    RealType initialEnergy = tempVec.dot( outerDirection );

    for ( int j = 0; j < _numSteps; j++ ) {
      timeSteps[j] = ( j - ( _numSteps / 2 )) * _stepSize;
      VectorType shiftedPoint = testPoint + timeSteps[j] * innerDirection;
      _F.apply( shiftedPoint, tempVec );
      energies[j] = tempVec.dot( outerDirection );
      derivs[j] = initialEnergy + timeSteps[j] * gateauxDerivative;
    }

    generatePNG( timeSteps, energies, derivs, saveName );

  }

  //
  void plotSingleDirection( const VectorType &testPoint, int innerDirection, int outerDirection,
                            const std::string saveName ) const {

    int numDofs = testPoint.size();
    int dimRange = _dimOfRange < 0 ? testPoint.size() : _dimOfRange;

    if ( !( outerDirection < dimRange ))
      throw BasicException( "VectorValuedDerivativeTester::plotSingleDirections: outerDir too large!" );

    if ( !( innerDirection < numDofs ))
      throw BasicException( "VectorValuedDerivativeTester::plotSingleDirections: innerDir too large!" );

    VectorType innerDir( dimRange ), outerDir( numDofs );
    innerDir.setZero();
    innerDir[innerDirection] = 1.;
    outerDir.setZero();
    outerDir[outerDirection] = 1.;
    plotSingleDirection( testPoint, innerDir, outerDir, saveName );
  }

  //
  RealType computeDiffQuotient( const VectorType &testPoint, const VectorType &outerDir, const VectorType &innerDir,
                                RealType tau ) const {

    // evalualte DF(m)
    int numDofs = innerDir.size();
    int dimRange = _dimOfRange < 0 ? testPoint.size() : _dimOfRange;

    VectorType derivative( dimRange ), derivativeShifted( dimRange );
    _F.apply( testPoint, derivative );

    //evaluate DF(m + h testDirection )
    VectorType shiftedPoint = testPoint + tau * innerDir;
    _F.apply( shiftedPoint, derivativeShifted );
    return ( derivativeShifted.dot( outerDir ) - derivative.dot( outerDir )) / tau;
  }

  //
  void testSingleDirection( const VectorType &testPoint, const VectorType &outerDir, const VectorType &innerDir,
                            RealType tau = -1. ) const {

    int numDofs = testPoint.size();
    int dimRange = _dimOfRange < 0 ? testPoint.size() : _dimOfRange;

    RealType stepSize = ( tau < 0 ) ? _stepSize : tau;
    MatrixType JacobiMatrix( dimRange, numDofs );
    _DF.apply( testPoint, JacobiMatrix );
    MatrixType JacobiMatrixTrans = JacobiMatrix.transpose();
    VectorType gradientWRTOuterDirection = JacobiMatrixTrans * outerDir;
    RealType gateaux = gradientWRTOuterDirection.dot( innerDir );
    std::cerr << std::abs(( computeDiffQuotient( testPoint, outerDir, innerDir, stepSize ) - gateaux ) / gateaux )
              << std::endl;

  }

  //
  void testAllDirections( const VectorType &testPoint, bool detailedOutput = false, RealType tau = -1. ) const {

    int numDofs = testPoint.size();
    int dimRange = _dimOfRange < 0 ? testPoint.size() : _dimOfRange;
    RealType stepSize = ( tau < 0 ) ? _stepSize : tau;

    //evaluate DF and compute DF^T
    MatrixType JacobiMatrix( dimRange, numDofs );
    VectorType outerDir( dimRange ), innerDir( numDofs );
    _DF.apply( testPoint, JacobiMatrix );
    MatrixType JacobiMatrixTrans = JacobiMatrix.transpose();

    // run over outer directions

    for ( int i = 0; i < dimRange; i++ ) {
      VectorType error( numDofs ), diffQuot( numDofs ), resDeriv( numDofs );
      std::vector<std::tuple<int, int>> directionIndices( numDofs );
      outerDir.setZero();
      outerDir[i] = 1.;
      VectorType gradientWRTOuterDirection( numDofs );
      gradientWRTOuterDirection = JacobiMatrixTrans * outerDir;
      // run over inner directions
      for ( int j = 0; j < numDofs; j++ ) {
        innerDir.setZero();
        innerDir[j] = 1.;
        RealType gateaux = gradientWRTOuterDirection.dot( innerDir );
        resDeriv[j] = gateaux;
        diffQuot[j] = computeDiffQuotient( testPoint, outerDir, innerDir, stepSize );
        error[j] = std::abs(( diffQuot[j] - resDeriv[j] ) / resDeriv[j] );
        directionIndices[j] = std::make_tuple(i, j);
      }
      if ( detailedOutput ) {
        printVector( directionIndices, 10 );
        printVector( diffQuot, 10 );
        printVector( resDeriv, 10 );
      }
      printVector( error, 10 );
      std::cerr << std::endl;
    }
  }

  //
  void testRandomDirections( const VectorType &testPoint, unsigned numOfRandomDir, bool detailedOutput = false,
                             RealType tau = -1.,
                             bool testNonNegEntriesOfHessianOnly = true ) {

    int numDofs = testPoint.size();
    int dimRange = _dimOfRange < 0 ? testPoint.size() : _dimOfRange;
    std::srand( std::time( 0 )); // use current time as seed for random generator
    RealType stepSize = ( tau < 0 ) ? _stepSize : tau;

    MatrixType JacobiMatrix( dimRange, numDofs );
    VectorType outerDir( dimRange ), innerDir( numDofs ), error( numOfRandomDir );
    _DF.apply( testPoint, JacobiMatrix );
    MatrixType JacobiMatrixTrans = JacobiMatrix.transpose();
    RealType maxError = 0.;

    // diffQuotient
    VectorType randApproxDeriv( numOfRandomDir );

    VectorType diffQuot( numOfRandomDir ), resDeriv( numOfRandomDir );
    std::vector<std::string> directionIndices( numOfRandomDir );

    for ( unsigned i = 0; i < numOfRandomDir; ++i ) {
      int outerIdx = std::rand() % dimRange;
      int innerIdx = std::rand() % numDofs;
      outerDir.setZero();
      innerDir.setZero();
      outerDir[outerIdx] = 1.;
      innerDir[innerIdx] = 1.;
      VectorType gradientWRTOuterDirection( numDofs );
      gradientWRTOuterDirection = JacobiMatrixTrans * outerDir;
      RealType gateaux = gradientWRTOuterDirection.dot( innerDir );
      RealType diffQuotient = computeDiffQuotient( testPoint, outerDir, innerDir, stepSize );

      resDeriv[i] = gateaux;
      diffQuot[i] = diffQuotient;
      directionIndices[i] = std::to_string( innerIdx ) + "," + std::to_string( outerIdx );

      if ( testNonNegEntriesOfHessianOnly && std::abs( gateaux ) < 1.e-10 && std::abs( diffQuotient ) < 1.e-10 ) {
        i--;
        continue;
      }

      RealType denom = std::max( std::abs( diffQuotient ), std::abs( gateaux ));
      error[i] = denom > 1.e-10 ? std::abs(( diffQuotient - gateaux ) / denom ) : 0;
      if ( error[i] > maxError )
        maxError = error[i];
    }

    if ( detailedOutput ) {
      printVector( directionIndices, 10 );
      printVector( diffQuot, 10 );
      printVector( resDeriv, 10 );
    }

    printVector( error, 10 );
    std::cout << "max. error = " << maxError << std::endl << std::endl;
  }

  void setNumSteps( int Steps ) {
    _numSteps = Steps;
  }

  void setStepSize( RealType Stepsize ) {
    _stepSize = Stepsize;
  }

};


//! Derivative test for functionals $F: \R^n \to \R^{m,n}$, where $m$ is the dimension of the range
template<typename ConfiguratorType>
class TensorValuedDerivativeTester {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TensorType TensorType;

  const BaseOp<VectorType, MatrixType> &_F;
  const BaseOp<VectorType, TensorType> &_DF;
  const RealType _stepSize;
  const int _dimOfRange;
  const int _numSteps;

public:
  TensorValuedDerivativeTester( const BaseOp<VectorType, MatrixType> &F,
                                const BaseOp<VectorType, TensorType> &DF,
                                const RealType stepSize,
                                const int dimOfRange,
                                const int numSteps = 50 )
          : _F( F ),
            _DF( DF ),
            _stepSize( stepSize ),
            _dimOfRange( dimOfRange ),
            _numSteps( 50 ) {}


  void plotRandomDirections( const VectorType &testPoint, int numOfRandomDir, const std::string saveNameStem,
                             const bool testNonNegEntriesOfHessianOnly = true ) const {

    int numDofs = testPoint.size();

    TensorType Hessian( _dimOfRange, numDofs, numDofs );
    _DF.apply( testPoint, Hessian );

    VectorType timeSteps( _numSteps ), energies( _numSteps ), derivs( _numSteps );

    std::srand( std::time( 0 )); // use current time as seed for random generator

    for ( int i = 0; i < numOfRandomDir; i++ ) {

      int dim1 = std::rand() % _dimOfRange; //point at which hessian is checked
      int dim2 = std::rand() % numDofs;     //first derivative
      int dim3 = std::rand() % numDofs;     //second derivative

      std::ostringstream saveName;
      saveName << saveNameStem << "_" << dim1 << "_" << dim2 << "_" << dim3 << ".png";

      if ( testNonNegEntriesOfHessianOnly && ( std::abs( Hessian( dim1, dim2, dim3 )) < 1.e-10 )) {
        i--;
        continue;
      }

      RealType gateauxDerivative = Hessian( dim1, dim2, dim3 );

      MatrixType tempMat( _dimOfRange, numDofs );
      _F.apply( testPoint, tempMat );
      RealType initialEnergy = tempMat.coeffRef( dim1, dim2 );

      for ( int j = 0; j < _numSteps; j++ ) {
        timeSteps[j] = ( j - _numSteps / 2 ) * _stepSize;
        VectorType shiftedPoint = testPoint;
        shiftedPoint[dim3] += timeSteps[j];
        _F.apply( shiftedPoint, tempMat );
        energies[j] = tempMat.coeffRef( dim1, dim2 );
        derivs[j] = initialEnergy + timeSteps[j] * gateauxDerivative;
      }
      generatePNG( timeSteps, energies, derivs, saveName.str());

    }

  }

  void plotAllDirections( const VectorType &testPoint, const std::string saveNameStem ) const {
    int numDofs = testPoint.size();

    TensorType Hessian( _dimOfRange, numDofs, numDofs );
    _DF.apply( testPoint, Hessian );

    VectorType timeSteps( _numSteps ), energies( _numSteps ), derivs( _numSteps );

    for ( int dim1 = 0; dim1 < _dimOfRange; dim1++ ) { //point at which hessian is checked
      for ( int dim2 = 0; dim2 < numDofs; dim2++ ) { //first derivative
        for ( int dim3 = 0; dim3 < numDofs; dim3++ ) { //second derivative
          std::ostringstream saveName;
          saveName << saveNameStem << "_" << dim1 << "_" << dim2 << "_" << dim3 << ".png";

          RealType gateauxDerivative = Hessian( dim1, dim2, dim3 );

          MatrixType tempMat( _dimOfRange, numDofs );
          _F.apply( testPoint, tempMat );
          RealType initialEnergy = tempMat.coeffRef( dim1, dim2 );

          for ( int j = 0; j < _numSteps; j++ ) {
            timeSteps[j] = ( j - _numSteps / 2 ) * _stepSize;
            VectorType shiftedPoint = testPoint;
            shiftedPoint[dim3] += timeSteps[j];
            _F.apply( shiftedPoint, tempMat );
            energies[j] = tempMat.coeffRef( dim1, dim2 );
            derivs[j] = initialEnergy + timeSteps[j] * gateauxDerivative;
          }
          generatePNG( timeSteps, energies, derivs, saveName.str());

        }
      }
    }

  }


  void testAllDirections( const VectorType &testPoint, bool detailedOutput= false ) const {
    int numDofs = testPoint.size();

    TensorType Hessian( _dimOfRange, numDofs, numDofs );
    _DF.apply( testPoint, Hessian );

    VectorType outerDir( numDofs ), innerDir( numDofs );

    std::vector<RealType> error, diffQuot, resDeriv;
    std::vector<std::tuple<int,int,int>> directionIndices;

    for ( int dim1 = 0; dim1 < _dimOfRange; dim1++ ) { //point at which hessian is checked
      for ( int dim2 = 0; dim2 < numDofs; dim2++ ) { //first derivative
        outerDir.setZero();
        outerDir[dim2] = 1.;
        for ( int dim3 = 0; dim3 < numDofs; dim3++ ) { //second derivative
          innerDir.setZero();
          innerDir[dim3] = 1.;
          RealType gateauxDerivative = Hessian( dim1, dim2, dim3 );
          RealType diffQuotient = computeDiffQuotient(testPoint, dim1, outerDir, innerDir, _stepSize );

          diffQuot.push_back( diffQuotient );
          resDeriv.push_back( gateauxDerivative );
          error.push_back( std::abs(( diffQuotient - gateauxDerivative ) / gateauxDerivative ));
          directionIndices.emplace_back(  dim1, dim2, dim3 );
        }
      }
    }

    if ( detailedOutput ) {
      printVector( directionIndices, 10 );
      printVector( diffQuot, 10 );
      printVector( resDeriv, 10 );
    }

    printVector( error, 10 );

  }

  void testRandomDirections( const VectorType &testPoint, int numOfRandomDir, bool detailedOutput = false,
                             bool testNonNegEntriesOfHessianOnly = true ) const {
    int numDofs = testPoint.size();

    TensorType Hessian( _dimOfRange, numDofs, numDofs );
    _DF.apply( testPoint, Hessian );

    VectorType outerDir( numDofs ), innerDir( numDofs );

    std::vector<RealType> error, diffQuot, resDeriv;
    std::vector<std::tuple<int, int, int>> directionIndices;

    // Setup randomness
    std::random_device rd;
    std::mt19937 rng( rd());
    std::uniform_int_distribution<int> dist_range( 0, _dimOfRange - 1 );
    std::uniform_int_distribution<int> dist_domain( 0, numDofs - 1 );

    for ( unsigned i = 0; i < numOfRandomDir; ++i ) {
      int dim1 = dist_range( rng );
      int dim2 = dist_domain( rng );
      int dim3 = dist_domain( rng );

      outerDir.setZero();
      outerDir[dim2] = 1.;

      innerDir.setZero();
      innerDir[dim3] = 1.;

      RealType gateauxDerivative = Hessian( dim1, dim2, dim3 );
      RealType diffQuotient = computeDiffQuotient( testPoint, dim1, outerDir, innerDir, _stepSize );

      if ( testNonNegEntriesOfHessianOnly && std::abs( gateauxDerivative ) < 1.e-10 &&
           std::abs( diffQuotient ) < 1.e-10 ) {
        i--;
        continue;
      }

      diffQuot.push_back( diffQuotient );
      resDeriv.push_back( gateauxDerivative );
      error.push_back( std::abs(( diffQuotient - gateauxDerivative ) / gateauxDerivative ));
      directionIndices.emplace_back( dim1, dim2, dim3 );
    }


    if ( detailedOutput ) {
      printVector( directionIndices, 10 );
      printVector( diffQuot, 10 );
      printVector( resDeriv, 10 );
    }

    printVector( error, 10 );

  }

  RealType computeDiffQuotient( const VectorType &testPoint, int targetDim, const VectorType &outerDir, const VectorType &innerDir,
                                RealType tau ) const {

    // evalualte DF(m)
    int numDofs = innerDir.size();
    int dimRange = _dimOfRange < 0 ? testPoint.size() : _dimOfRange;

    MatrixType derivative, derivativeShifted;
    _F.apply( testPoint, derivative );

    //evaluate DF(m + h testDirection )
    VectorType shiftedPoint = testPoint + tau * innerDir;
    _F.apply( shiftedPoint, derivativeShifted );
    return (( derivativeShifted - derivative) * outerDir)[targetDim] / tau;
  }



  //! TODO
  void testSingleDirection( const VectorType &testPoint, const VectorType &dir1, const VectorType &dir2,
                            const VectorType &dir3, RealType tau = -1. ) const {

  }

  //! TODO

};

#endif
