// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Energies and optimization operator for discrete geodesic interpolation.
 * \author Heeren
 * 
 */
#ifndef GEODESICINTERPOLATION_HH
#define GEODESICINTERPOLATION_HH

#include <iostream>

#include <goast/Core.h>
#include <goast/Optimization.h>

#include <goast/NRIC/LeastSquaresReconstruction.h>

#include "ConvexGeodesicCombination.h"


//======================================================================================================================
// APPROXIMATIVE GEODESIC INTERPOLATIONS
//======================================================================================================================

/**
 * \brief Computes approximation of discrete geodesic by reconstruction from \f$ L\Theta \f$ - geodesic
 * \author Heeren
 *
 * For two given meshes \f$ s_A \f$ and \f$ s_B \f$, compute corresponding representations \f$ z_A \f$ and \f$ z_B \f$ in \f$ L \Theta \f$ - space,
 * i.e. by computing the vector of edge lengths and dihedral angles, perform linear interpolation of length \f$ K+1 \f$ between \f$ z_A \f$ and \f$ z_B \f$
 * and reconstruct shapes \f$ s_k \f$ from intermediate values \f$ z_k \f$ for \f$ k = 1, \ldots, K-1 \f$ by least squares optimization.
 *
 * Method as presented in the paper \cite FrBo11
 */
template<typename ConfiguratorType>
void computeLThetaInterpolation( const typename ConfiguratorType::VectorType& startGeom,
                                 const typename ConfiguratorType::VectorType& endGeom, 
                                 const MeshTopologySaver& Topology,
                                 const double bendingWeight,
                                 const int K, 
                                 const int maxIterations,
                                 const std::vector<int> * bdryMask,
                                 typename ConfiguratorType::VectorType& initialization ) {
            
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  
  RealType tau = 1. / K;
  int numV = Topology.getNumVertices();
  int totalNumDofs = 3*(K-1)*numV;
  bool hasInitialization( initialization.size() == totalNumDofs );
  if( !hasInitialization )
    initialization.resize( totalNumDofs );    
  
  // vector of lenghts and angles that are tried to be matched in the reconstruction
  VectorType startLengthsAngles, endLengthsAngles;
  computeLengthsAngles<ConfiguratorType>( Topology, startGeom, startLengthsAngles );
  computeLengthsAngles<ConfiguratorType>( Topology, endGeom, endLengthsAngles );
  
  // integration weights in reconstruction functionals 
  VectorType Weights;
  getReconstructionWeights<ConfiguratorType>( Topology, startGeom, Weights );
      
  // reconstruct in parallel
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for( int k = 1; k < K; k++ ){
    //std::cerr << "Start " << k << "th of " << K-1 << " reconstructions with iters = " << maxIterations << std::endl;
      
    // compute geodesic, i.e. linear interpolation, in LTheta space
    VectorType LengthsAngles = tau * (K-k) * startLengthsAngles + tau * k * endLengthsAngles;  
      
    // define reconstruction functionals
    RealType edgeLengthWeight = 1.;
    LinearReconstructionFunctional<ConfiguratorType> L( Topology, LengthsAngles, Weights, edgeLengthWeight, bendingWeight );
    LinearReconstructionDerivative<ConfiguratorType> DL( Topology, LengthsAngles, Weights, edgeLengthWeight, bendingWeight );
  
    // define initialization and apply Gauss-Newton
    VectorType init( startGeom ), destGeom;
    if( 2*k > K )
        init = endGeom;
    GaussNewtonAlgorithm<ConfiguratorType> GNOp( 2 * Topology.getNumEdges(), L, DL, maxIterations );
    GNOp.setQuiet();
    if( bdryMask ){
      GNOp.setBoundaryMask( *bdryMask );
      // use original content as initialization (at least at boundary)
      if( hasInitialization ){
          for( int i = 0; i < 3*numV; i++ )
              init[i] = initialization[(k-1)*3*numV + i];
      }          
    }
    GNOp.solve( init, destGeom );
  
    // write back 
    for( int i = 0; i < 3*numV; i++ )
      initialization[(k-1)*3*numV + i] = destGeom[i];
  }
    
}
  
/**
 * \brief Compute sequence of weighted 3-shape geodesics
 * \author Heeren, Perl
 *
 * Optimize \f$ \mapsto (K-k) \cdot W[S_0, S] + k \cdot W[S, S_K]\f$ for \f$ k = 1, ..., K-1 \f$ sequentially, using \f$ S_{k-1} \f$ as initialization.
 */
template<typename ConfiguratorType>    
void computeIterativeInterpolation( const typename ConfiguratorType::VectorType& startGeom,
                              const typename ConfiguratorType::VectorType& endGeom, 
                              const DeformationBase<ConfiguratorType>& W,
                              const MeshTopologySaver& Topology,
                              const int K, 
                              const OptimizationParameters<ConfiguratorType>& optPars,
                              const std::vector<int> * bdryMask,
                              typename ConfiguratorType::VectorType& initialization  ) {
    
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  
  RealType tau = 1. / K;
  int numV = Topology.getNumVertices();
  initialization.resize( 3*(K-1)*numV );
  VectorType activeShape = startGeom;      

  for( int k = 1; k < K; k++ ){
      RealType lambda = k * tau;
      ConvexInterpolationOp<ConfiguratorType> convexOp( Topology, startGeom, endGeom, W, lambda, optPars );
      if( bdryMask )
          convexOp.setBoundaryMask( *bdryMask );
      convexOp.execute( activeShape );   
      // write back 
      for( int i = 0; i < 3*numV; i++ )
            initialization[(k-1)*3*numV + i] = activeShape[i];
  }
}
  
  
//===============================================================================================================================
//===============================================================================================================================
// TIME-DISCRTE GEODESIC INTERPOLATION
//===============================================================================================================================
//===============================================================================================================================

/**
 * \brief Discrete path energy
 * \author Heeren
 *
 * For two given shapes \f$ s_A \f$ and \f$ s_B \f$, some integer \f$ K > 0 \f$ and some deformation functional \f$ W \f$,
 * this class realizes the discrete path energy \f[ E[s_1, \ldots, s_{K-1}] = K \sum_{k=1}^K W[s_{k-1}, s_k]\f]
 * with \f$ s_0 = s_A \f4 and \f$ s_K = s_B \f$.
 */
template<typename ConfiguratorType>
class DiscretePathEnergy : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType>,
                           public TimedClass<DiscretePathEnergy<ConfiguratorType>> {

protected:    
  typedef typename ConfiguratorType::RealType          RealType;

  typedef typename ConfiguratorType::VectorType        VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;
  typedef typename ConfiguratorType::TripletType       TripletType;

  using typename TimedClass<DiscretePathEnergy<ConfiguratorType>>::ScopeTimer;

  const DeformationBase<ConfiguratorType>& _W;
  const VectorType& _start, _end;
  int _K;
  
public:
  DiscretePathEnergy( const DeformationBase<ConfiguratorType>& W, int K, const VectorType& start, const VectorType& end ) : _W(W), _start(start), _end(end), _K(K) {}
    
  // Arg =  (s_1, ..., s_{K-1}) are intermediate shapes only, where s_k = (x_k, y_k, z_k)
  void apply ( const VectorType& Arg, RealType & Dest ) const {
    ScopeTimer timer("apply");
      
      int numOfFreeShapes = _K - 1;
      
      if( Arg.size()%numOfFreeShapes != 0 )
        throw BasicException("DiscretePathEnergy::apply: wrong number of dofs!");
      
      if( Arg.size()/numOfFreeShapes != _start.size() )
        throw BasicException("DiscretePathEnergy::apply: wrong size of dofs!");
      
      VectorType singleEnergies;
      evaluateSingleEnergies( Arg, singleEnergies );
      
      Dest = 0.;
      for( int k = 0; k < _K; k++ )
          Dest += _K * singleEnergies[k];
  
  }
  
  void evaluateSingleEnergies( const VectorType& Arg, VectorType& Energies ) const {
       int numOfFreeShapes = _K - 1;
      
      if( Arg.size()%numOfFreeShapes != 0 )
        throw BasicException("DiscretePathEnergy::evaluateSingleEnergies: wrong number of dofs!");
      
      if( Arg.size()/numOfFreeShapes != _start.size() )
        throw BasicException("DiscretePathEnergy::evaluateSingleEnergies: wrong size of dofs!");
      
      Energies.resize( _K );
      Energies.setZero();
          
      // bring into more convenient form
      std::vector< Eigen::Ref<const VectorType> > argRefs;
      argRefs.reserve(numOfFreeShapes);
      const int numLocalDofs = Arg.size() / numOfFreeShapes;
      
      for( int k = 0; k < numOfFreeShapes; k++ )
          argRefs.push_back( Arg.segment(k*numLocalDofs, numLocalDofs) );
      
      // compute path energy
      _W.applyEnergy ( _start, argRefs[0], Energies[0] );
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for( int k = 1; k < numOfFreeShapes; k++ )
          _W.applyEnergy ( argRefs[k-1], argRefs[k], Energies[k] );
      _W.applyEnergy ( argRefs[numOfFreeShapes-1], _end, Energies[numOfFreeShapes] );
  }
    
};

//! \brief Gradient of DiscretePathEnergy
//! \author Heeren
template<typename ConfiguratorType>
class DiscretePathEnergyGradient : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>,
                           public TimedClass<DiscretePathEnergyGradient<ConfiguratorType>> {

protected:    
  typedef typename ConfiguratorType::RealType          RealType;

  typedef typename ConfiguratorType::VectorType        VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;
  typedef typename ConfiguratorType::TripletType       TripletType;

  using typename TimedClass<DiscretePathEnergyGradient<ConfiguratorType>>::ScopeTimer;

  const DeformationBase<ConfiguratorType>& _W;
  const VectorType& _start, _end;
  int _K;
  std::vector<int> _mask;
  
public:
  DiscretePathEnergyGradient( const DeformationBase<ConfiguratorType>& W, int K, const VectorType& start, const VectorType& end ) : _W(W), _start(start), _end(end), _K(K) {}
  
  void setMask( const std::vector<int>& localMask ) { 
    fillPathMask( _K-1, _start.size(), localMask, _mask );
  }
  
  // Arg =  (s_1, ..., s_{K-1}) are intermediate shapes only, where s_k = (x_k, y_k, z_k)
  void apply ( const VectorType& Arg, VectorType& Dest ) const {
    ScopeTimer timer("apply");

      int numOfFreeShapes = _K - 1;
      
      if( Arg.size()%numOfFreeShapes != 0 )
          throw BasicException("DiscretePathEnergyGradient::apply: wrong number of dofs!");
      
      if( Arg.size()/numOfFreeShapes != _start.size() )
        throw BasicException("DiscretePathEnergyGradient::apply: wrong size of dofs!");
      
      if( Dest.size() != Arg.size() )
          Dest.resize( Arg.size() );      
      
      Dest.setZero();
      
      // bring into more convenient form
      std::vector<Eigen::Ref<const VectorType> > argRefs;
      argRefs.reserve(numOfFreeShapes);
      std::vector<Eigen::Ref<VectorType> > destRefs;
      destRefs.reserve(numOfFreeShapes);
      
      const int numLocalDofs = Arg.size() / numOfFreeShapes;
      for( int k = 0; k < numOfFreeShapes; k++ ){
          argRefs.push_back( Arg.segment(k*numLocalDofs, numLocalDofs) );
          destRefs.push_back( Dest.segment(k*numLocalDofs, numLocalDofs) );
      }
      
      // compute path energy gradient
      _W.applyAddDefGradient ( _start, argRefs[0], destRefs[0], _K );
      _W.applyAddUndefGradient ( argRefs[numOfFreeShapes-1], _end, destRefs[numOfFreeShapes-1], _K );
      
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for( int k = 0; k < numOfFreeShapes; k++ ){      
          if( k < numOfFreeShapes-1 )
            _W.applyAddUndefGradient ( argRefs[k], argRefs[k+1], destRefs[k], _K );
          if( k > 0 )
            _W.applyAddDefGradient ( argRefs[k-1], argRefs[k], destRefs[k], _K );
      }
      
      // mask? 
      if( _mask.size() > 0 )
          applyMaskToVector<VectorType>( _mask, Dest );

  }
    
};

//! \brief Hessian of DiscretePathEnergy
//! \author Heeren
template<typename ConfiguratorType>
class DiscretePathEnergyHessian : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType>,
                           public TimedClass<DiscretePathEnergyHessian<ConfiguratorType>> {

protected:    
  typedef typename ConfiguratorType::RealType          RealType;

  typedef typename ConfiguratorType::VectorType        VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;
  typedef typename ConfiguratorType::TripletType       TripletType;
  typedef std::vector<TripletType> TripletListType;

  using typename TimedClass<DiscretePathEnergyHessian<ConfiguratorType>>::ScopeTimer;

  const DeformationBase<ConfiguratorType>& _W;
  const VectorType& _start, _end;
  int _K;
  std::vector<int> _mask;
  
public:
  DiscretePathEnergyHessian( const DeformationBase<ConfiguratorType>& W, int K, const VectorType& start, const VectorType& end ) : _W(W), _start(start), _end(end), _K(K) {}
  
  void setMask( const std::vector<int>& localMask ) { 
    fillPathMask( _K-1, _start.size(), localMask, _mask );
  }
  
  // Arg =  (s_1, ..., s_{K-1}) are intermediate shapes only, where s_k = (x_k, y_k, z_k)
  void applyAdd ( const VectorType& Arg, MatrixType& Dest ) const {
    ScopeTimer timer("applyAdd");
      int numOfFreeShapes = _K - 1;
      if( Arg.size()%numOfFreeShapes != 0 )
        throw BasicException("DiscretePathEnergyHessian::applyAdd: wrong number of dofs!");      
      if( (Dest.cols() != Arg.size()) || (Dest.rows() != Arg.size()) )
          Dest.resize( Arg.size(), Arg.size() );    
      
      MatrixType Update;
      apply( Arg, Update );
      Dest += Update;
  }
  
  // Arg =  (s_1, ..., s_{K-1}) are intermediate shapes only, where s_k = (x_k, y_k, z_k)
  void apply( const VectorType& Arg, MatrixType& Dest ) const {
    ScopeTimer timer("apply");

      int numOfFreeShapes = _K - 1;
      
      if( Arg.size()%numOfFreeShapes != 0 )
        throw BasicException("DiscretePathEnergyHessian::apply: wrong number of dofs!");
      
      if( Arg.size()/numOfFreeShapes != _start.size() )
        throw BasicException("DiscretePathEnergyHessian::apply: wrong size of dofs!");
      
      if( (Dest.cols() != Arg.size()) || (Dest.rows() != Arg.size()) )
          Dest.resize( Arg.size(), Arg.size() ); 
      
      // fill triplet lists
      TripletListType tripletList;
      // we have approx. 3*(K-1) blocks
      tripletList.reserve( 3*numOfFreeShapes*_W.numOfNonZeroHessianEntries() );   
   
      pushTriplets( Arg, tripletList );
      
      Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
      
      // mask?
      if( _mask.size() > 0 )
          applyMaskToSymmetricMatrix<MatrixType>( _mask, Dest );
  }
  
  // fill triplets
  // Arg =  (s_1, ..., s_{K-1}) are intermediate shapes only, where s_k = (x_k, y_k, z_k)
  void pushTriplets( const VectorType& Arg, TripletListType& tripletList ) const {
          ScopeTimer("pushTriplets");
      int numOfFreeShapes = _K - 1;
      if( Arg.size()%numOfFreeShapes != 0 )
        throw BasicException("DiscretePathEnergyHessian::pushTriplets: wrong number of dofs!");           
      
      // bring into more convenient form
      std::vector<Eigen::Ref<const VectorType> > argRefs;
      argRefs.reserve(numOfFreeShapes);
      
      const int numLocalDofs = Arg.size() / numOfFreeShapes;
      for( int k = 0; k < numOfFreeShapes; k++ )
          argRefs.push_back( Arg.segment(k*numLocalDofs, numLocalDofs) );   
      
      //TODO parallelize by OpenMP, however, push_back of std::vector is not threadsafe!!!!
      // Use std::move for instance:
      // 
      // #include <iterator>
      // std::vector<T> a(100), b(100);
      // std::size_t n = a.size();
      // a.resize(a.size() + b.size());
      // std::move(b.begin(), b.end(), a.begin() + n);
      // OR: std::move(b.begin(), b.end(), std::back_inserter(a))
      
      // parts involving fixed shapes
      _W.pushTripletsDefHessian ( _start, argRefs[0], tripletList, 0, 0, _K );
      _W.pushTripletsUndefHessian ( argRefs[numOfFreeShapes-1], _end, tripletList, (numOfFreeShapes-1)*numLocalDofs, (numOfFreeShapes-1)*numLocalDofs, _K );
    
      // compute path energy Hessian on diagonal
      for( int k = 0; k < numOfFreeShapes-1; k++ ){      
          _W.pushTripletsDefHessian ( argRefs[k], argRefs[k+1], tripletList, (k+1)*numLocalDofs, (k+1)*numLocalDofs, _K );
          _W.pushTripletsUndefHessian ( argRefs[k], argRefs[k+1], tripletList, k*numLocalDofs, k*numLocalDofs, _K );
      }
      
      // compute path energy mixed Hessian
      for( int k = 0; k < numOfFreeShapes-1; k++ ){      
          _W.pushTripletsMixedHessian ( argRefs[k], argRefs[k+1], tripletList, k*numLocalDofs, (k+1)*numLocalDofs, false, _K );
          _W.pushTripletsMixedHessian ( argRefs[k], argRefs[k+1], tripletList, (k+1)*numLocalDofs, k*numLocalDofs, true, _K );
      }  
  }
      
};

//===============================================================================================================================
//===============================================================================================================================

/**
 * \brief Discrete geodesic interpolation operator to minimize DiscretePathEnergy
 * \brief Heeren
 *
 * Computes discrete geodesics of length \f$ K+1 \f$ by minimization of discrete path energy while fixing end points.
 * One can specify different initialization strategies in the execute() routines.
 *
 * \todo Add more optimization strategies, e.g. by sequential convex optim., red-black, cascadic or L-Theta reconstruction.
 */
template<typename ConfiguratorType>
class GeodesicInterpolation{
    
protected:    
  typedef typename ConfiguratorType::RealType          RealType;

  typedef typename ConfiguratorType::VectorType        VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;

  const MeshTopologySaver& _topology;
  const VectorType& _startGeom, _endGeom;
  const DeformationBase<ConfiguratorType>& _W;
  int _K;
  
  const OptimizationParameters<ConfiguratorType>& _optPars;
  bool _quiet;
  std::vector<int>* _fullPathMask;  
  std::vector<int> const* _singleShapeMask;
  
public:    
    GeodesicInterpolation( const MeshTopologySaver& Topology,
                                const VectorType& StartGeom, 
                                const VectorType& EndGeom,
                                const DeformationBase<ConfiguratorType>& W,
                                int lengthOfGeodesic, 
                                const OptimizationParameters<ConfiguratorType>& optPars,
                                bool quiet = true ) 
    :  _topology(Topology), _startGeom(StartGeom), _endGeom(EndGeom), _W(W), _K(lengthOfGeodesic-1), _optPars(optPars), _quiet(quiet), _fullPathMask(NULL), _singleShapeMask(NULL){ }
    
    ~GeodesicInterpolation(){
        if(_fullPathMask){
            if(!_quiet) std::cerr << "GeodesicInterpolation: Delete bdry mask." << std::endl;
            delete _fullPathMask;
        }
    }
    
    void setBoundaryMask( const std::vector<int>& localMask ) {
      if( localMask.size() > 0 ){
        _singleShapeMask = &localMask;
        _fullPathMask = new std::vector<int>;
        fillPathMask( _K-1, 3*_topology.getNumVertices(), localMask, *_fullPathMask );
      }
    }
    
    // initializationScheme = -1 means "There is already a meaningful initialization stored in Path".
    // initializationScheme: 0 = start shape, 1 = end shape, 2 = linear interpolation, 3 = iterative initialization, 4 = LTheta initialization
    void execute( VectorType& Path, int initializationScheme = -1, std::string saveNameStem = "" ) const {               
      
      int numOfFreeShapes = _K-1;
      bool saveResults( saveNameStem.size() > 0 );

      DiscretePathEnergy<ConfiguratorType>          E( _W, _K, _startGeom, _endGeom );
      DiscretePathEnergyGradient<ConfiguratorType> DE( _W, _K, _startGeom, _endGeom  );
      DiscretePathEnergyHessian<ConfiguratorType> D2E( _W, _K, _startGeom, _endGeom  );

      // initialization
      int numV = _topology.getNumVertices();
      int totalNumDofs = 3*numOfFreeShapes*numV;
      if( (initializationScheme == -1) && (Path.size() != totalNumDofs) )
          throw BasicException("GeodesicInterpolation::execute: size of argument is wrong!");
      if( initializationScheme != -1){
          if(!_quiet) std::cerr << "=============================================================" << std::endl;
          if(!_quiet) std::cerr << "Start initialization." << std::endl;  
          initializePath( Path, initializationScheme );   
          if(!_quiet) std::cerr << "Initialization done." << std::endl;  
          if(!_quiet) std::cerr << "=============================================================" << std::endl << std::endl;
          
          // saving?
          if( saveResults ){
            std::ostringstream solNameStem;
            solNameStem << saveNameStem << "_init";
            saveSolution( Path, solNameStem.str() );
          }
      }
      VectorType initialization( Path );

      RealType energy;
      VectorType grad;
      if(!_quiet){
        E.apply( initialization, energy );
        std::cerr << "Initial path energy = " << energy << std::endl;
        DE.apply( initialization, grad );
        if( _fullPathMask )
          applyMaskToVector<VectorType>( *_fullPathMask, grad );
        std::cerr << "Initial path energy gradient norm = " << grad.norm() << std::endl << std::endl;
      }
        
      // Optimization with gradient descent 
      if( _optPars.getGradientIterations( ) > 0 ){
        if(!_quiet)  std::cerr << "Start gradient descent... " << std::endl;
        GradientDescent<ConfiguratorType> GD( E, DE, _optPars );  
        if(_fullPathMask){
            if(!_quiet) std::cerr << "Set boundary mask... " << std::endl;
            GD.setBoundaryMask( *_fullPathMask );
        }
        GD.solve( initialization, Path );
        initialization = Path;
        
        // saving?
        if( saveResults ){
          std::ostringstream solNameStem;
          solNameStem << saveNameStem << "_grad";
          saveSolution( Path, solNameStem.str() );
        }
      }     
      
      // optmization with BFGS
      if( _optPars.getBFGSIterations()  > 0 ){          
        if(!_quiet) std::cerr << "Start Quasi-Newton... " << std::endl;
        QuasiNewtonBFGS<ConfiguratorType> QNM( E, DE, _optPars );
        if(_fullPathMask){
            if(!_quiet) std::cerr << "Set boundary mask... " << std::endl;
            QNM.setBoundaryMask( *_fullPathMask );
        }
        QNM.solve( initialization, Path );
        initialization = Path;
        
        // saving?
        if( saveResults ){
          std::ostringstream solNameStem;
          solNameStem << saveNameStem << "_bfgs";
          saveSolution( Path, solNameStem.str() );
        }
      }
      
      // optmization with Newton
      if( _optPars.getNewtonIterations( ) > 0 ){     
        if( _fullPathMask ){
          if(!_quiet) std::cerr << "Start Newton with boundary mask... " << std::endl;
          //NewtonOptimizationMethod<ConfiguratorType> NM( E, DE, D2E, _optPars );
          LineSearchNewton<ConfiguratorType> NM (E, DE, D2E, _optPars );
          if(_fullPathMask){
            if(!_quiet) std::cerr << "Set boundary mask... " << std::endl;
            NM.setBoundaryMask( *_fullPathMask );
          }
          NM.solve( initialization, Path );             
        }
        else{
          if(!_quiet) std::cerr << "Start Newton with Lagrange... " << std::endl;
   
          RigidBodyMotionsConstraintHandler<ConfiguratorType> constHandler( _topology, _startGeom, numOfFreeShapes );
          int numConstraints = numOfFreeShapes * 6;
        
          LagrangeFunctional<ConfiguratorType, DiscretePathEnergy<ConfiguratorType>, RigidBodyMotionsConstraintHandler<ConfiguratorType> > L( E, constHandler, totalNumDofs );
          LagrangeGradient<ConfiguratorType, DiscretePathEnergyGradient<ConfiguratorType>, RigidBodyMotionsConstraintHandler<ConfiguratorType> > DL(  DE, constHandler, totalNumDofs );
          LagrangeHessian<ConfiguratorType, DiscretePathEnergyHessian<ConfiguratorType>, RigidBodyMotionsConstraintHandler<ConfiguratorType> > D2L(  D2E, constHandler, totalNumDofs );
        
          initialization.conservativeResize( totalNumDofs + numConstraints );
          Path.conservativeResize( totalNumDofs + numConstraints );
        
          NewtonMethod<ConfiguratorType> NM( DL, D2L, _optPars );
          NM.solve( initialization, Path );
          Path.conservativeResize( totalNumDofs );
        }    
        
        // saving?
        if( saveResults ){
          std::ostringstream solNameStem;
          solNameStem << saveNameStem << "_newton";
          saveSolution( Path, solNameStem.str() );
        }
      }   
      
      if(!_quiet){
        E.apply( Path, energy );
        std::cerr << "Final path energy = " << energy << std::endl;
        DE.apply( Path, grad );
        if( _fullPathMask )
          applyMaskToVector<VectorType>( *_fullPathMask, grad );
        std::cerr << "Final path energy gradient norm = " << grad.norm() << std::endl << std::endl;
      }
    }
    
    // initializationScheme: 0 = start shape, 1 = end shape, 2 = linear interpolation, 3 = iterative initialization, 4 = LTheta initialization
    void execute( std::string saveNameStem, int initializationScheme = 2 ) const {      
        VectorType Path;
        execute( Path, initializationScheme, saveNameStem );    
    }    
    
    void checkIfGeodesicIsRelaxed( const VectorType& Path ) const {
        DiscretePathEnergy<ConfiguratorType>          E( _W, _K, _startGeom, _endGeom );
        DiscretePathEnergyGradient<ConfiguratorType> DE( _W, _K, _startGeom, _endGeom  );

      RealType energy;
        VectorType tempVec;
        E.apply( Path, energy );
        std::cerr << "Path energy is " << energy[0] << std::endl;
        E.evaluateSingleEnergies( Path, tempVec );
        std::cerr << "Intermediate deform. energies are " << std::endl;
        for( int i = 0; i < tempVec.size(); i++ )
                std::cerr << tempVec[i] << ", ";
        std::cerr << std::endl;
        DE.apply( Path, tempVec );
        std::cerr << "Path energy gradient norm is " << tempVec.norm() << std::endl << std::endl;
    }
    
protected:    
    //
    void saveSolution( const VectorType& Solution, std::string saveNameStem ) const {
      int numOfFreeShapes = _K-1;
      std::vector<Eigen::Ref<const VectorType> > argRefs;
      argRefs.reserve(numOfFreeShapes);
      
      const int numLocalDofs = Solution.size() / numOfFreeShapes;
      if( numLocalDofs == 0 )
          throw BasicException("GeodesicInterpolation::saveSolution: solution vector is empty!");
      
      for( int k = 0; k < numOfFreeShapes; k++ )
          argRefs.push_back( Solution.segment(k*numLocalDofs, numLocalDofs) ); 
      
      TriMesh auxMesh( _topology.getGrid() );
      for( int k = 0; k < numOfFreeShapes; k++ ){
        std::ostringstream saveName;
        saveName << saveNameStem << "_" << k+1 << ".ply";
        setGeometry( auxMesh, argRefs[k] );
        if ( !OpenMesh::IO::write_mesh( auxMesh, saveName.str() ) )
          throw BasicException("GeodesicInterpolation::saveSolution: could not write mesh!");
      }
    }
    
    //
    void initializeWithLTheta ( VectorType& initialization ) const {
        int numIterationsGaussNewton = 25;
        RealType bendingWeight = _W.getBendingWeight();;
        computeLThetaInterpolation<ConfiguratorType>( _startGeom, _endGeom, _topology, bendingWeight, _K, numIterationsGaussNewton, _singleShapeMask, initialization );
    }
    
    // initializationScheme: 0 = start shape, 1 = end shape, 2 = linear interpolation
    void initializePath( VectorType& initialization, int initializationScheme ) const {
        
      int numOfFreeShapes = _K-1;
      int numV = _topology.getNumVertices();
      int totalNumDofs = 3*numOfFreeShapes*numV;
      if( initialization.size() != totalNumDofs )
        initialization.resize( totalNumDofs );
      
      // initialize with start shape
      if( initializationScheme == 0 ){
        if(!_quiet) std::cerr << "Initialize with start shape." << std::endl;
        for( int k = 0; k < numOfFreeShapes; k++ )
          for( int i = 0; i < 3*numV; i++ )
            initialization[k*3*numV+i]  = _startGeom[i];
      }
      // initialize with end shape
      if( initializationScheme == 1 ){
        if(!_quiet) std::cerr << "Initialize with end shape." << std::endl;
        for( int k = 0; k < numOfFreeShapes; k++ )
          for( int i = 0; i < 3*numV; i++ )
            initialization[k*3*numV+i]  = _endGeom[i];
      }
      // linear interpolation
      if( initializationScheme == 2 ){
        if(!_quiet) std::cerr << "Initialize with linear interpolation of nodal values." << std::endl;
        RealType tau = 1. / _K;
        for( int k = 0; k < numOfFreeShapes; k++ )
          for( int i = 0; i < 3*numV; i++ )
            initialization[k*3*numV+i]  = _startGeom[i] + (k+1) * tau * (_endGeom[i] - _startGeom[i]);
      }    
      // iterative initialization
      if( initializationScheme == 3 ){
          if(!_quiet) std::cerr << "Initialize with iterative scheme." << std::endl;
          OptimizationParameters<ConfiguratorType> iterativeOptPars;
          iterativeOptPars.setGradientIterations( 25 );
          iterativeOptPars.setNewtonIterations( 10 );
          //iterativeOptPars.setVerbose();
          computeIterativeInterpolation<ConfiguratorType>( _startGeom, _endGeom, _W, _topology, _K, iterativeOptPars, _singleShapeMask, initialization );
      }
      
      // LTheta initialization
      if( initializationScheme == 4 ){
          if(!_quiet) std::cerr << "Initialize with LTheta scheme." << std::endl;
          initializeWithLTheta ( initialization );
      }

      if( initializationScheme == 5 ){
        if(!_quiet) std::cerr << "Initialize with start and shape." << std::endl;
        for( int k = 0; k < numOfFreeShapes; k++ )
          if ( k < numOfFreeShapes / 2 )
            for( int i = 0; i < 3*numV; i++ )
             initialization[k*3*numV+i]  = _startGeom[i];
          else
            for( int i = 0; i < 3*numV; i++ )
              initialization[k*3*numV+i]  = _endGeom[i];
      }
    }

};

#endif