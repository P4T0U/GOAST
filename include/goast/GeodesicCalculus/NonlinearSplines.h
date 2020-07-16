// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Nonlinear spline energies.
 * \author Heeren
 * 
 */
#ifndef NONLINEARSPLINES_HH
#define NONLINEARSPLINES_HH

#include <goast/Core.h>
#include <goast/Optimization.h>
#include "DiscreteGeodesicCalculus.h"

#include <iostream>

//===============================================================================================================================
//===============================================================================================================================

/**
 * \brief Nonlinear spline energy.
 * \author Heeren
 * 
 * This class implements the nonlinear spline energy \f[ F[s_0, \ldots, s_K] = 4K^3 \sum_{k=1}^{K-1} W[s_k, \tilde s_k]\, ,\f]
 * with the nonlinear constraint that \f$ (s_{k-1}, \tilde s_k, s_{k+1}) \f$ is a discrete geodesic for \f$ k = 1, \ldots, K-1 \f$.
 *
 * We shall refer to the constraint shapes \f$ \tilde s_1, \ldots, \tilde s_{K-1} \f$ as dual shapes.
 *
 * In practice, a subset of the shapes (including \f$ s_0 \f$ and \f$ s_K \f$!) are supposed to be fixed (so-called key frames or key shapes).
 *
 * This is determined by a vector \$ \{ 0,1 \}^{K+1} \$ given in the constructor specifying whether the kth shape is fixed (value = 1) or free (value = 0), \f$ k = 0, \ldots, K \f$.
 * 
 * \note See eq. (7) in \cite HeRuSc16
 */
template<typename ConfiguratorType>
class SplineEnergy :  public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;

    typedef typename ConfiguratorType::VectorType VectorType;

    const MeshTopologySaver& _topology;
    const DeformationBase<ConfiguratorType>& _W; 
    const OptimizationParameters<ConfiguratorType>& _optPars;
    int _numTimePoints;
    std::vector<int> _fixTimePoints;
    mutable std::vector<VectorType> _dualVariables;
    bool _quiet;    
    mutable VectorType _energies;    
    std::vector<int> _bdryMask;
    bool _periodicBC;

public:
    SplineEnergy( const MeshTopologySaver& Topology,
		  const DeformationBase<ConfiguratorType> &W, 
		  const std::vector<int>& FixTimePoints,
		  bool periodicBC,
		  const OptimizationParameters<ConfiguratorType>& optPars ) 
    : _topology(Topology),
      _W(W),
      _optPars(optPars),
      _numTimePoints( FixTimePoints.size() ), 
      _fixTimePoints( FixTimePoints ),
      _quiet( false ),
      _energies( _numTimePoints ),
      _periodicBC( periodicBC ){}
    
    void setQuietMode( bool quiet ){
      _quiet = quiet;
    }
    
    //! if there is a boundary to be fixed
    void setBoundaryMask( const std::vector<int>& mask ){
      _bdryMask =  mask;
    }
    
    const VectorType& getEnergies( ) const {
      return _energies;
    }
  
    // Arg contains concatenation of all K+1 shapes (fixed and free shapes!)
    void apply( const VectorType& Arg, RealType & Dest ) const {
      
      int numLocalDofs = 3 * _topology.getNumVertices();
      int numGlobalDofs = _numTimePoints * numLocalDofs;
      if( Arg.size() != numGlobalDofs )
	throw BasicException ( "SplineEnergy::apply(): argument has wrong size!!" );
      
      // bring into more convenient form
      std::vector<Eigen::Ref<const VectorType> > shapes;
      shapes.reserve(_numTimePoints);
      for( int k = 0; k < _numTimePoints; k++ )
        shapes.push_back( Arg.segment(k*numLocalDofs, numLocalDofs) );
      
      // initialize dual variables
      if( _dualVariables.size() == 0 )
	for( int k = 0; k < _numTimePoints; k++ )
	  _dualVariables.push_back( shapes[k] );
	
      _energies.setZero();
      
      // periodic boundary conditions?
      int minTildeObjectIdx = _periodicBC ? 0 : 1;
      int numOfActiveTimePoints = _periodicBC ? _numTimePoints - 1 : _numTimePoints;
	
#ifdef _OPENMP
#pragma omp parallel for
#endif	      
      for( int k = minTildeObjectIdx; k < _numTimePoints - 1; k++  ){
	
	//if( !_quiet ) std::cerr  << k << " of " << _numTimePoints - 1 << std::endl;
	
	// compute \tilde x_k (use old dual variable as initialization?)
	int prev = (k + numOfActiveTimePoints - 1) % numOfActiveTimePoints;
	int next = (k + 1) % numOfActiveTimePoints;
	// compute short geodesic
	convexInterpolation<ConfiguratorType>( _topology, shapes[prev], shapes[next], _bdryMask, _W, _optPars, 0.5, _dualVariables[k] );
	  
	// W[x_k, \tilde x_k]
	RealType temp;
	_W.applyEnergy( shapes[k], _dualVariables[k], temp );
        _energies[k]  = temp;
      }
      
      //TODO multilpy by correct factor 4K^3
      Dest = _energies.sum();
    }  
};

//! \brief Gradient of SplineEnergy
//! \author Heeren
template<typename ConfiguratorType>
class SplineGradient :  public BaseOp<typename ConfiguratorType::VectorType > {
  
protected:  
    typedef typename ConfiguratorType::RealType RealType;

    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;
    typedef std::vector<TripletType> TripletListType;

    const MeshTopologySaver& _topology;
    const DeformationBase<ConfiguratorType>& _W; 
    const OptimizationParameters<ConfiguratorType>& _optPars;
    int _numTimePoints;
    std::vector<int> _fixTimePoints;
    mutable std::vector<VectorType> _dualVariables, _scaledAdjoints;
    bool _quiet;     
    std::vector<int> _bdryMask;
    mutable TripletListType _constraintTriplets;
    bool _periodicBC, _gradientTest;

public:
    SplineGradient( const MeshTopologySaver& Topology,
		  const DeformationBase<ConfiguratorType> &W, 
		  const std::vector<int>& FixTimePoints,
		  bool periodicBC,
		  const OptimizationParameters<ConfiguratorType>& optPars ) 
    : _topology(Topology),
      _W(W),
      _optPars(optPars),
      _numTimePoints( FixTimePoints.size() ), 
      _fixTimePoints( FixTimePoints ),
      _scaledAdjoints( _numTimePoints ),
      _quiet( false ),
      _periodicBC( periodicBC ),
      _gradientTest(false){ }

    
    //! if there is a boundary to be fixed
    void setBoundaryMask( const std::vector<int>& mask ){
      _bdryMask =  mask;
    }
    
    void setQuietMode( bool quiet ){
      _quiet = quiet;
    }
    
    void setGradientTest( bool test ){
      _gradientTest = test;
    }
  
    // Arg contains concatenation of all K+1 shapes (fixed and free shapes!)
    void apply( const VectorType& Arg, VectorType& Dest ) const {  
      
      if( !_quiet ) std::cerr << std::endl << "------------------------------------------------" << std::endl;
      if( !_quiet ) std::cerr  << "Start spline gradient evaluation." << std::endl;
      auto t_start = std::chrono::high_resolution_clock::now();
      
      int numLocalDofs = 3 * _topology.getNumVertices();
      int numGlobalDofs = _numTimePoints * numLocalDofs;
      if( Arg.size() != numGlobalDofs )
	throw BasicException ( "SplineGradient::apply(): argument has wrong size!!" );
      
      if( Dest.size() != Arg.size() )
	Dest.resize( Arg.size() );
      Dest.setZero();
      
      // bring into more convenient form
      std::vector<Eigen::Ref<const VectorType> > shapes;
      shapes.reserve(_numTimePoints);
      std::vector<Eigen::Ref<VectorType> > partialGrads;
      partialGrads.reserve(_numTimePoints);
      for( int k = 0; k < _numTimePoints; k++ ){
        shapes.push_back( Arg.segment(k*numLocalDofs, numLocalDofs) );
	partialGrads.push_back( Dest.segment(k*numLocalDofs, numLocalDofs) );
	partialGrads[k].setZero();
      }
      
      // initialize dual variables
      if( _dualVariables.size() == 0 )
	for( int k = 0; k < _numTimePoints; k++ )
	  _dualVariables.push_back( shapes[k] );
	
      // periodic boundary conditions?
      int minTildeObjectIdx = _periodicBC ? 0 : 1;
      int numOfActiveTimePoints = _periodicBC ? _numTimePoints - 1 : _numTimePoints;
      
      if( !_quiet ) std::cerr << "Start to compute geodesic segments, adjoints and outer partial derivatives..." << std::endl;
      
#ifdef _OPENMP
#pragma omp parallel for
#endif	
      for( int k = minTildeObjectIdx; k < _numTimePoints - 1; k++ ){
	
	// define previous and next index
	int prev = (k + numOfActiveTimePoints - 1) % numOfActiveTimePoints;
	int next = (k + 1) % numOfActiveTimePoints;
        
	// compute \tilde x_k
	// if( !_fixTimePoints[k] ) _dualVariables[k] = shapes[k];
        convexInterpolation<ConfiguratorType>( _topology, shapes[prev], shapes[next], _bdryMask, _W, _optPars, 0.5, _dualVariables[k] ); 
	
	// Compute adjoint p_k
	computeAdjoint( shapes[prev], shapes[k], shapes[next], _dualVariables[k], _scaledAdjoints[k] );

	// gradient test?
	if( !_gradientTest && _fixTimePoints[k] )
	  continue;
        
	// Derivative w.r.t. geometries:
	// compute derivative of spline energy w.r.t. y_k (not accounting for the inner derivative \partial_{y_k} \tilde y_k here!)
        // \partial_{y_k} \E[\ldots] += \partial_{y_k} \W[y_k , \tilde y_k ]
	_W.applyAddUndefGradient( shapes[k], _dualVariables[k], partialGrads[k] );
      }
      
       if( !_quiet ) std::cerr << "Start to compute left mixed derivatives..." << std::endl;       
#ifdef _OPENMP
#pragma omp parallel for
#endif
      // partial_{y_{k-1}} \E[\ldots] +=  \partial_{y_{k-1}} \partial_{\tilde y} C^k[y_{k-1}, y_{k+1}, \tilde y_k] * p_k
      for ( int i = minTildeObjectIdx; i < _numTimePoints - 1; i++ ) {
	VectorType grad;
	// define previous index
	int prev = (i + numOfActiveTimePoints - 1) % numOfActiveTimePoints;        
	// gradient test or fix gradient?
	if( !_gradientTest && _fixTimePoints[prev] )
	  continue;
        addLeftMixedDerivative ( shapes[prev], _dualVariables[i], _scaledAdjoints[i], partialGrads[prev] );
      }

      if( !_quiet ) std::cerr << "Start to compute right mixed derivatives..." << std::endl;      
#ifdef _OPENMP
#pragma omp parallel for
#endif
      // partial_{y_{k+1}} \E[\ldots] +=  \partial_{y_{k+1}} \partial_{\tilde y} C^k[y_{k-1}, y_{k+1}, \tilde y_k]  * p_k
      for ( int i = minTildeObjectIdx; i < _numTimePoints - 1; i++ ) {
	VectorType grad;
	// define next index
	int next = (i + 1) % numOfActiveTimePoints;
	// gradient test or fix gradient?
	if( !_gradientTest && _fixTimePoints[next] )
	  continue;
        addRightMixedDerivative ( shapes[next], _dualVariables[i], _scaledAdjoints[i],  partialGrads[next]  );
      }
      
      // no masking in gradient test
      if( _gradientTest )
          return;     
	
#ifdef _OPENMP
#pragma omp parallel for
#endif
      // set gradient w.r.t. fixed data points zero    
      for ( int i = 0; i < _numTimePoints; i++ ){
	if( _fixTimePoints[i] )
	  partialGrads[i].setZero();
	else
          applyMaskToVector( _bdryMask, partialGrads[i] );	
      }
      
      auto t_end = std::chrono::high_resolution_clock::now();
      if( !_quiet ) std::cerr  << "Spline gradient evaluation done in " << std::chrono::duration<double, std::ratio<1> >(t_end - t_start).count() << " seconds." << std::endl << std::endl;
    }
    
protected:
  // compute adjoint p_k
  void computeAdjoint ( const VectorType& prevGeometry,
			const VectorType& currGeometry,
			const VectorType& nextGeometry,
                        const VectorType &dualGeometry,
                        VectorType &Adjoint ) const {     
      
      // compute rhs, i.e. \partial_{\tilde y_k} \E[\ldots] = \partial_{\tilde y_k} \W[\y_k, \tilde y_k]
      VectorType negativeRHS;
      _W.applyDefGradient( currGeometry, dualGeometry, negativeRHS );     
      applyMaskToVector( _bdryMask, negativeRHS );
      
      // compute matrix, i.e. \partial_{\tilde y} \partial_{\tilde y} C^k[y_{k-1}, y_{k+1},\tilde y_k]
      bool useLagrangeSetup( _bdryMask.size() == 0 );
      const int numLocalDofs = 3 * _topology.getNumVertices();
      const int dimMatrix = useLagrangeSetup ? numLocalDofs + 6 : numLocalDofs;
      
      // get triplets
      TripletListType tripletList;      
      tripletList.reserve( _W.numOfNonZeroHessianEntries() );
      _W.pushTripletsDefHessian ( prevGeometry, dualGeometry, tripletList, 0, 0, -1.);
      _W.pushTripletsUndefHessian ( dualGeometry, nextGeometry, tripletList, 0, 0, -1. );
      
      // if there is no boundary mask, we make use of Lagrange setting and push constraint triplets
      if( useLagrangeSetup ){
	  if( !_quiet ) std::cerr << "Account for Lagrange setup." << std::endl;
	  // push triplets to fix zeroth and first momentum
	  if( _constraintTriplets.size() == 0 ){
	    if( !_quiet ) std::cerr << "Push constraint triplets." << std::endl;
	    VectorType refGeom;
	    getGeometry( _topology.getGrid(), refGeom );
	    RigidBodyMotionsConstraintHandler<ConfiguratorType>(_topology, refGeom, 1).addConstraintHessian( _constraintTriplets );
	  }
	  tripletList.insert(tripletList.end(), _constraintTriplets.begin(), _constraintTriplets.end());
	  
	  // extend rhs
	  negativeRHS.conservativeResize( dimMatrix );
	  for( int i = 0; i < 6; i++ )
	    negativeRHS[numLocalDofs + i] = 0.;
      }
      
      // assemble matrix
      if( !_quiet ) std::cerr << "Assemble matrix." << std::endl; 
      MatrixType negativeHessian( dimMatrix, dimMatrix );
      negativeHessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );  
	
      // mask matrix and rhs
      if( !useLagrangeSetup ){
	  if( !_quiet ) std::cerr << "Mask system matrix and rhs." << std::endl;
	  applyMaskToSymmetricMatrixAndVector( _bdryMask, negativeHessian, negativeRHS );
      }
	
      // solve linear system
      if( !_quiet ) std::cerr << "Solve linear system for adjoint variable." << std::endl;
      LinearSolver<ConfiguratorType>().solve( negativeHessian, negativeRHS, Adjoint );
	
      // resize solution
      if( useLagrangeSetup )
	  Adjoint.conservativeResize(numLocalDofs);
  }
  
  // partial_{y_{k-1}} \E[\ldots] +=  \partial_{y_{k-1}} \partial_{\tilde y} C^k[y_{k-1}, y_{k+1}, \tilde y_k] * p_k
  void addLeftMixedDerivative ( const VectorType& prevGeometry, 
				const VectorType &dualGeometry, 
				const VectorType &Adjoint, 
				Eigen::Ref<VectorType> prevGradient ) const {   
      bool FirstDerivWRTDef = true;
      TripletListType tripletList;      
      tripletList.reserve( _W.numOfNonZeroHessianEntries() );
      _W.pushTripletsMixedHessian ( prevGeometry, dualGeometry, tripletList, 0, 0, FirstDerivWRTDef );	
      
      MatrixType  MixedHessian( 3 * _topology.getNumVertices(), 3 * _topology.getNumVertices() );
      MixedHessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );       
      applyMaskToMatrix( _bdryMask, MixedHessian );
      
      prevGradient += MixedHessian * Adjoint;  
  }
  
  // partial_{y_{k+1}} \E[\ldots] +=  \partial_{y_{k+1}} \partial_{\tilde y} C^k[y_{k-1}, y_{k+1}, \tilde y_k]  * p_k
  void addRightMixedDerivative ( const VectorType& nextGeometry, 
				 const VectorType &dualGeometry, 
				 const VectorType &Adjoint, 
				 Eigen::Ref<VectorType> nextGradient ) const {   
      
      bool FirstDerivWRTDef = false;
      TripletListType tripletList;      
      tripletList.reserve( _W.numOfNonZeroHessianEntries() );
      _W.pushTripletsMixedHessian ( dualGeometry, nextGeometry, tripletList, 0, 0, FirstDerivWRTDef );	
      
      MatrixType  MixedHessian( 3 * _topology.getNumVertices(), 3 * _topology.getNumVertices() );
      MixedHessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() ); 
      applyMaskToMatrix( _bdryMask, MixedHessian );
      
      nextGradient += MixedHessian * Adjoint;
  }
  
};

#endif
