// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Explicit classes for (weighted) 3-shape geodesic interpolation and extrapolation
 * \author Heeren
 * 
 */
#ifndef CONVEXGEODESICCOMBINATION_HH
#define CONVEXGEODESICCOMBINATION_HH

#include <goast/Core.h>
#include <goast/NRIC/LeastSquaresReconstruction.h>

#include <iostream>

//===============================================================================================================================
//===============================================================================================================================

//! \brief Wrapper function for Newton's method to solve for intermediate shape in 3-shape geodesic
//! \author Heeren
template<typename ConfiguratorType, typename GradientType, typename HessianType>
void solveNewtonSingleShape( const MeshTopologySaver& Topology, 
                             const GradientType& DF, 
                             const HessianType& D2F, 
                             const typename ConfiguratorType::VectorType& Init, 
                             typename ConfiguratorType::VectorType& Solution,
                             const OptimizationParameters<ConfiguratorType>& optParams,
                             const std::vector<int>* Mask = NULL,
                             bool quiet = true ) {        
        if( Mask ){
          if(!quiet) std::cerr << "Start Newton with boundary mask... " << std::endl;
          NewtonMethod<ConfiguratorType> NM( DF, D2F, optParams );
          NM.setBoundaryMask( *Mask );
          NM.solve( Init, Solution );             
        }
        else{
          if(!quiet) std::cerr << "Start Newton with Lagrange... " << std::endl;
   
          typename ConfiguratorType::VectorType refGeom, initialization(Init);
          getGeometry( Topology.getGrid(), refGeom );
          RigidBodyMotionsConstraintHandler<ConfiguratorType> constHandler( Topology, refGeom, 1 );
          int numConstraints = 6;
          int totalNumDofs = 3 * Topology.getNumVertices();

          LagrangeGradient<ConfiguratorType, GradientType, RigidBodyMotionsConstraintHandler<ConfiguratorType> > DL(  DF, constHandler, totalNumDofs );
          LagrangeHessian<ConfiguratorType, HessianType, RigidBodyMotionsConstraintHandler<ConfiguratorType> > D2L(  D2F, constHandler, totalNumDofs );
        
          initialization.conservativeResize( totalNumDofs + numConstraints );
          Solution.conservativeResize( totalNumDofs + numConstraints );
        
          NewtonMethod<ConfiguratorType> NM( DL, D2L, optParams );
          NM.solve( initialization, Solution );
          Solution.conservativeResize( totalNumDofs );
        }  
}

/**
 * \brief Wrapper function to perform gradient descent on least squares functional
 * \author Heeren
 *
 * For gradient class \f$ x \mapsto DF(x) \f$ define least squares energy \f$ E[x] = DF^T(x) DF(x)\f$
 * and perform gradient descent on \f$ E \f$.
 */
template<typename ConfiguratorType, typename GradientType, typename HessianType>
void solveLeastSquaresGradientDescentSingleShape( const MeshTopologySaver& Topology, 
                                                   const GradientType& DF, 
                                                   const HessianType& D2F, 
                                                   const typename ConfiguratorType::VectorType& Init, 
                                                   typename ConfiguratorType::VectorType& Solution,
                                                   const OptimizationParameters<ConfiguratorType>& optParams,
                                                   const std::vector<int>* Mask = NULL,
                                                   bool quiet = true ) {        
        if( Mask ){
          if(!quiet) std::cerr << "Start gradient descent of least squares functional F[x] = | E[x] |^2 with boundary mask... " << std::endl;
          SquaredFunctional<ConfiguratorType> SqE2E( DF );
          SquaredDerivative<ConfiguratorType> SqE2G( DF, D2F );
          
          GradientDescent<ConfiguratorType> GD( SqE2E, SqE2G, optParams );  
          GD.setBoundaryMask( *Mask );
          GD.solve( Init, Solution );             
        }
        else{
          throw BasicException("solveLeastSquaresGradientDescentSingleShape: Lagrange method has not been implemented yet");
        }  
}

/**
 * \brief Wrapper function to perform (approximative) Newton's method on least squares functional
 * \author Heeren
 *
 * Same procedure as in function solveLeastSquaresGradientDescentSingleShape above but using Newton's method with approximative Hessian.
 * In detail, \f$ E = DF^T Df \f$, then \f$ DE =  2 (D^2F)^T DF \f$ and \f$ D^2 E \approx 2 (D^2F)^T D^2F \f$
 */
template<typename ConfiguratorType, typename GradientType, typename HessianType>
void solveLeastSquaresLineSearchNewtonSingleShape( const MeshTopologySaver& Topology, 
                                                   const GradientType& DF, 
                                                   const HessianType& D2F, 
                                                   const typename ConfiguratorType::VectorType& Init, 
                                                   typename ConfiguratorType::VectorType& Solution,
                                                   const OptimizationParameters<ConfiguratorType>& optParams,
                                                   const std::vector<int>* Mask = NULL,
                                                   bool quiet = true ) {        
        if( Mask ){
          if(!quiet) std::cerr << "Start LineSearchNewton of least squares functional F[x] = | E[x] |^2 with boundary mask... " << std::endl;
          SquaredFunctional<ConfiguratorType> SqE2E( DF );
          SquaredDerivative<ConfiguratorType> SqE2G( DF, D2F );
          ReducedSquaredHessian<ConfiguratorType> SqE2H( DF, D2F  );
          
          LineSearchNewton<ConfiguratorType> NM ( SqE2E, SqE2G, SqE2H, optParams );
          NM.setBoundaryMask( *Mask );
          NM.solve( Init, Solution );             
        }
        else{
          throw BasicException("solveLeastSquaresLineSearchNewtonSingleShape: Lagrange method has not been implemented yet");
        }  
}

//===============================================================================================================================
//===============================================================================================================================
// CONVEX COMBINATION FOR INTERPOLATION
//===============================================================================================================================
//===============================================================================================================================

/**
 * \brief Weighted path energy for \f$ K=2 \f$
 * \author Heeren
 *
 * Realization of energy \f$ s \mapsto (1 - \lambda) * W[s_l, s] + \lambda W[s, s_r] \f$ for given \f$ s_l, s_r \f$ and \f$ 0 < \lambda < 1 \f$.
 * For \f$ \lambda = 0.5 \f$ this is the path energy for a 3-point discrete geodesic (i.e. \f$  K = 2 \f$) with free intermediate shape.
 */
template<typename ConfiguratorType>
class WeightedGeod3Energy : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {

protected:    
  typedef typename ConfiguratorType::RealType          RealType;

  typedef typename ConfiguratorType::VectorType        VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;
  typedef typename ConfiguratorType::TripletType       TripletType;

  const DeformationBase<ConfiguratorType>& _W;
  const VectorType& _start, _end;
  RealType _lambda;
  
public:
  WeightedGeod3Energy( const DeformationBase<ConfiguratorType>& W, const VectorType& start, const VectorType& end, RealType Lambda ) : _W(W), _start(start), _end(end), _lambda(Lambda) {}
    
  // argument is middle shape
  void apply ( const VectorType& Arg, RealType & Dest ) const {
      
      if( Arg.size() != _start.size() )
        throw BasicException("WeightedGeod3Energy::apply: wrong number of dofs!");

      Dest = 0.;
      _W.applyAddEnergy ( _start, Arg, Dest, 1. - _lambda );
      _W.applyAddEnergy ( Arg, _end, Dest,  _lambda );  
  }
    
};

//! \brief Gradient of WeightedGeod3Energy
//! \author Heeren
template<typename ConfiguratorType>
class WeightedGeod3Gradient : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {

protected:    
  typedef typename ConfiguratorType::RealType          RealType;

  typedef typename ConfiguratorType::VectorType        VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;
  typedef typename ConfiguratorType::TripletType       TripletType;

  const DeformationBase<ConfiguratorType>& _W;
  const VectorType& _start, _end;
  RealType _lambda;
  
public:
  WeightedGeod3Gradient( const DeformationBase<ConfiguratorType>& W, const VectorType& start, const VectorType& end, RealType Lambda ) : _W(W), _start(start), _end(end), _lambda(Lambda) {}
    
  // argument is middle shape
  void apply ( const VectorType& Arg, VectorType& Dest ) const {
      
      if( Arg.size() != _start.size() )
        throw BasicException("WeightedGeod3Gradient::apply: wrong number of dofs!");
      
      if( Dest.size() != _start.size() )
        Dest.resize( _start.size() );
      Dest.setZero();
      
      VectorType grad;   
      _W.applyDefGradient ( _start, Arg, grad );
      Dest += (1. - _lambda) * grad;
      _W.applyUndefGradient ( Arg, _end, grad );
      Dest += _lambda * grad;
  
  }
    
};

//! \brief Hessian of WeightedGeod3Energy
//! \author Heeren
template<typename ConfiguratorType>
class WeightedGeod3Hessian : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

protected:    
  typedef typename ConfiguratorType::RealType          RealType;

  typedef typename ConfiguratorType::VectorType        VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;
    typedef typename ConfiguratorType::TripletType       TripletType;
  typedef std::vector<TripletType> TripletListType;

  const DeformationBase<ConfiguratorType>& _W;
  const VectorType& _start, _end;
  RealType _lambda;
  
public:
  WeightedGeod3Hessian( const DeformationBase<ConfiguratorType>& W, const VectorType& start, const VectorType& end, RealType Lambda ) : _W(W), _start(start), _end(end), _lambda(Lambda) {}
    
  // argument is middle shape
  void apply ( const VectorType& Arg, MatrixType& Dest ) const {
      
      if( Arg.size() != _start.size() )
        throw BasicException("WeightedGeod3Hessian::apply: wrong number of dofs!");
      
      if( (Dest.rows() != _start.size()) || (Dest.cols() != _start.size()) )
        Dest.resize( _start.size(), _start.size() );
      Dest.setZero();
      
      // fill triplet lists
      TripletListType tripletList;      
      // we have approx. 3*(K-1) blocks
      tripletList.reserve( _W.numOfNonZeroHessianEntries() );   

      pushTriplets( Arg, tripletList );
      Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend() );  
  }
  
  void pushTriplets( const VectorType& Arg, TripletListType& tripletList ) const {    
      if( Arg.size() != _start.size() )
        throw BasicException("WeightedGeod3Hessian::pushTriplets: wrong number of dofs!");
        
      _W.pushTripletsDefHessian ( _start, Arg, tripletList, 0, 0, 1. - _lambda );
      _W.pushTripletsUndefHessian ( Arg, _end, tripletList, 0, 0, _lambda );
  }
    
};

//! \brief Interpolation operator to optimize WeightedGeod3Energy
//! \author Heeren
template<typename ConfiguratorType>
class ConvexInterpolationOp{
    
protected:    
  typedef typename ConfiguratorType::RealType          RealType;

  typedef typename ConfiguratorType::VectorType        VectorType;

  const MeshTopologySaver& _topology;
  const VectorType& _startGeom, _endGeom;
  const DeformationBase<ConfiguratorType>& _W;
  RealType  _lambda;
  
  const OptimizationParameters<ConfiguratorType>& _optPars;
  bool _quiet;
  std::vector<int> const* _mask;
  
  
public:    
    ConvexInterpolationOp( const MeshTopologySaver& Topology, 
                           const VectorType& StartGeom, 
                           const VectorType& EndGeom,
                           const DeformationBase<ConfiguratorType>& W,
                           RealType Lambda,
                           const OptimizationParameters<ConfiguratorType>& optPars,
                           bool quiet = true ) 
    :  _topology(Topology), _startGeom(StartGeom), _endGeom(EndGeom), _W(W), _lambda(Lambda), _optPars(optPars), _quiet(quiet), _mask(NULL){}
    
    void setBoundaryMask( const std::vector<int>& Mask ) {
      _mask = &Mask;
    }
    
    //
    void execute( VectorType& MiddleShape ) const {               
      
      WeightedGeod3Energy<ConfiguratorType>    E( _W, _startGeom, _endGeom, _lambda );
      WeightedGeod3Gradient<ConfiguratorType> DE( _W, _startGeom, _endGeom, _lambda );
      WeightedGeod3Hessian<ConfiguratorType> D2E( _W, _startGeom, _endGeom, _lambda );

      // initialization     
      if( MiddleShape.size() != _startGeom.size() )
         MiddleShape = ( _lambda < 0.5 ) ? _startGeom : _endGeom;
      
      RealType energy;
      VectorType grad;
      if(!_quiet) {
        E.apply( MiddleShape, energy );
        std::cerr << "Initial energy = " << energy << std::endl;
        DE.apply( MiddleShape, grad );
        std::cerr << "Initial gradient norm = " << grad.norm() << std::endl << std::endl;
      }
      VectorType initialization( MiddleShape );
        
      // Optimization with gradient descent 
      if( _optPars.getGradientIterations() > 0 ){
        if(!_quiet) std::cerr << "Start gradient descent... " << std::endl;
        GradientDescent<ConfiguratorType> GD( E, DE, _optPars );  
        if(_mask){
            if(!_quiet)  std::cerr << "Set boundary mask... " << std::endl;
            GD.setBoundaryMask( *_mask );
        }
        GD.solve( initialization, MiddleShape );
        initialization = MiddleShape;
      }     
      
      // optmization with BFGS
      if( _optPars.getBFGSIterations() > 0 ){          
        if(!_quiet) std::cerr << "Start Quasi-Newton... " << std::endl;
        QuasiNewtonBFGS<ConfiguratorType> QNM( E, DE, _optPars );
        if(_mask){
            if(!_quiet) std::cerr << "Set boundary mask... " << std::endl;
            QNM.setBoundaryMask( *_mask );
        }
        QNM.solve( initialization, MiddleShape );
        initialization = MiddleShape;
      }
      
      // optmization with Newton
      if( _optPars.getNewtonIterations() > 0 ){   
        typedef WeightedGeod3Gradient<ConfiguratorType> GradType;
        typedef WeightedGeod3Hessian<ConfiguratorType> HessType;
        if(!_quiet) std::cerr << "Start Newton... " << std::endl;
        solveNewtonSingleShape<ConfiguratorType, GradType, HessType >( _topology, DE, D2E, initialization, MiddleShape, _optPars, _mask );       
      }   
      
      if(!_quiet){
        E.apply( MiddleShape, energy );
        std::cerr << "Final energy = " << energy << std::endl;
        DE.apply( MiddleShape, grad );
        std::cerr << "Final gradient norm = " << grad.norm() << std::endl << std::endl;
      }
    }

};



//===============================================================================================================================
//===============================================================================================================================
// CONVEX COMBINATION FOR EXTRAPOLATION
//===============================================================================================================================
//===============================================================================================================================

/**
 * \brief Vector-valued, weighted discrete Exp2 operator
 * \author Heeren
 *
 * Gradient of \f$ s \mapsto (1 - \lambda) * W[s_l, s_v] + \lambda W[s_v, s] \f$ for given \f$ s_l, s_v \f$ and \f$ 0 < \lambda < 1\ f$.
 * For \f$ \lambda = 0.5 \f$ this is the gradient of the path energy for a 3-point discrete geodesic (i.e. \f$  K = 2 \f$).
 */
template <typename ConfiguratorType >
class WeightedExp2Energy : public BaseOp< typename ConfiguratorType::VectorType > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;

  const DeformationBase<ConfiguratorType>& _W;
  const VectorType& _shape0, _shape1;
  RealType _lambda;
  VectorType _constPart;
  int _numDofs;

public:
    WeightedExp2Energy ( const DeformationBase<ConfiguratorType>& W,
                         const VectorType& shape0,
                         const VectorType& shape1,
                         RealType Lambda ) :
      _W(W),
      _shape0(shape0),
      _shape1(shape1),
      _lambda( Lambda ),
      _numDofs(shape0.size()){
        // initialization
        calcConstPartOfEnergy();
      }


  //! The vertex positions of S_2 are given as argument.
  void apply( const VectorType& shape2, VectorType& Dest ) const {

    if( shape2.size() != _numDofs )
      throw BasicException ( "WeightedExp2Energy::apply(): arg has wrong size!" );
    if( Dest.size() != _numDofs )
      Dest.resize( _numDofs );
    
    // add constant partial
    Dest = _constPart;
    // add other part
    _W.applyAddUndefGradient( _shape1, shape2, Dest, _lambda );

  }

protected:
  // pre-cimpute constant part of energy W[S_0, S_1]
  void calcConstPartOfEnergy() {
    _constPart.resize( _numDofs );
    _W.applyDefGradient( _shape0, _shape1, _constPart );
    _constPart *= 1. - _lambda;
  }

};

//! \brief Matrix-valeud Derivative of WeightedExp2Energy
//! \author Heeren
template <typename ConfiguratorType >
class WeightedExp2Gradient : public BaseOp< typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType       TripletType;
  typedef std::vector<TripletType> TripletListType;

  const DeformationBase<ConfiguratorType>& _W;
  const VectorType& _shape1;
  RealType _lambda;
  int _numDofs;

public:
    WeightedExp2Gradient ( const DeformationBase<ConfiguratorType>& W,
                 const VectorType& shape0,
                 const VectorType& shape1,
                 RealType Lambda ) :
      _W(W),
      _shape1(shape1),
      _lambda( Lambda ),
      _numDofs(shape0.size()){  }
      
  //! The vertex positions of S_2 are given as argument.
  void apply( const VectorType& shape2, MatrixType& Dest ) const {

    if( shape2.size() != _numDofs )
      throw BasicException ( "WeightedExp2Gradient::apply(): arg has wrong size!" );
    if( (Dest.rows() != _numDofs) || (Dest.cols() != _numDofs) )
      Dest.resize( _numDofs, _numDofs );

    // fill triplet lists
    TripletListType tripletList;
    tripletList.reserve( _W.numOfNonZeroHessianEntries() );      
         
    pushTriplets( shape2, tripletList );
    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
  
  void pushTriplets( const VectorType& Arg, TripletListType& tripletList ) const {
    if( Arg.size() != _numDofs )
      throw BasicException ( "WeightedExp2Gradient::pushTriplets(): arg has wrong size!" );
    _W.pushTripletsMixedHessian( _shape1, Arg, tripletList, 0, 0, false, _lambda );
  }


};

//! \brief Extrapolation operator to compute root of WeightedExp2Energy
//! \author Heeren
template<typename ConfiguratorType>
class ConvexExtrapolationOp{
    
protected:    
  typedef typename ConfiguratorType::RealType          RealType;

  typedef typename ConfiguratorType::VectorType        VectorType;

  const MeshTopologySaver& _topology;
  const VectorType& _startGeom, _varGeom;
  const DeformationBase<ConfiguratorType>& _W;
  RealType  _lambda;
  
  const OptimizationParameters<ConfiguratorType>& _optPars;
  bool _quiet;
  std::vector<int> const* _mask;
  
  
public:    
    ConvexExtrapolationOp( const MeshTopologySaver& Topology, 
                           const VectorType& StartGeom, 
                           const VectorType& VariationGeom,
                           const DeformationBase<ConfiguratorType>& W,
                           RealType Lambda,
                           const OptimizationParameters<ConfiguratorType>& optPars,
                           bool quiet = true ) 
    :  _topology(Topology), _startGeom(StartGeom), _varGeom(VariationGeom), _W(W), _lambda(Lambda), _optPars(optPars), _quiet(quiet), _mask(NULL){}
    
    void setBoundaryMask( const std::vector<int>& Mask ) {
      _mask = &Mask;
    }
    
    //
    void execute( VectorType& nextShape ) const {               
      
      WeightedExp2Energy<ConfiguratorType>    E( _W, _startGeom, _varGeom, _lambda );
      WeightedExp2Gradient<ConfiguratorType> DE( _W, _startGeom, _varGeom, _lambda );

      // initialization     
      if( nextShape.size() != _varGeom.size() )
         nextShape = _varGeom;
      
      VectorType grad;
      if(!_quiet) {
        E.apply( nextShape, grad );
        std::cerr << "Initial gradient norm = " << grad.norm() << std::endl << std::endl;
      }
      VectorType initialization( nextShape );
      
      // optmization with Newton
      if( _optPars.getNewtonIterations() > 0 ){   
        typedef WeightedExp2Energy<ConfiguratorType> GradType;
        typedef WeightedExp2Gradient<ConfiguratorType> HessType;
        if(!_quiet) std::cerr << "Start Newton... " << std::endl;
        solveNewtonSingleShape<ConfiguratorType, GradType, HessType >( _topology, E, DE, initialization, nextShape, _optPars, _mask );       
      }   
      
      if(!_quiet){
        E.apply( nextShape, grad );
        std::cerr << "Final gradient norm = " << grad.norm() << std::endl << std::endl;
      }
    }

};

#endif
