// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef OPTIMIZATION_LAGRANGE_H
#define OPTIMIZATION_LAGRANGE_H

#include <vector>

#include "Objectives.h"

//!==========================================================================================================
//!==========================================================================================================
//!==========================================================================================================

//!==========================================================================================================
//!
template<typename ConfiguratorType, typename EnergyType, typename ConstraintHandlerType>
class LagrangeFunctional : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType>{
    
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef  typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;
  
  const EnergyType& _E;
  const ConstraintHandlerType& _constraintHandler;
  int _numDofs;
  
public:
  LagrangeFunctional( const EnergyType& E, const ConstraintHandlerType& constraintHandler, int numDofs ) 
  : _E(E), _constraintHandler(constraintHandler), _numDofs( numDofs ) {}
  
  void apply( const VectorType& Arg, RealType & Dest ) const {
    
    if( Arg.size() != _numDofs + _constraintHandler.getNumOfTotalConstraints() )
      throw BasicException("LagrangeFunctional::apply: argument has wrong size!");  
 
    Eigen::Ref<const VectorType> geomRef = Arg.segment(0, _numDofs );
    _E.apply( geomRef, Dest );

    Eigen::Ref<const VectorType> multRef = Arg.segment(_numDofs, _constraintHandler.getNumOfTotalConstraints() );    
    _constraintHandler.addConstraintEnergy( geomRef, Dest, multRef );

  }
  
};

//!
template<typename ConfiguratorType, typename GradientType, typename ConstraintHandlerType>
class LagrangeGradient : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>{
    
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef  typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;
  
  const GradientType& _DE;
  const ConstraintHandlerType& _constraintHandler;
  int _numDofs;
  
public:
  LagrangeGradient( const GradientType& DE, const ConstraintHandlerType& constraintHandler, int numDofs ) 
  : _DE(DE), _constraintHandler(constraintHandler),_numDofs(numDofs) {}
  
  void apply( const VectorType& Arg, VectorType& Dest ) const {
      
    int argSize = _numDofs + _constraintHandler.getNumOfTotalConstraints();
    if( Arg.size() != argSize ){
        std::cerr << "Size arg in grad  = " << Arg.size() << std::endl;
      throw BasicException("LagrangeGradient::apply: argument has wrong size!");
    }
    
    if( Dest.size() != argSize )
        Dest.resize( argSize );
    Dest.setZero();
      
    Eigen::Ref<const VectorType> geomRef = Arg.segment(0, _numDofs );
    VectorType gradRef( _numDofs );
    _DE.apply( geomRef, gradRef );
    Dest.segment( 0, _numDofs ) = gradRef;
    
    Eigen::Ref<const VectorType> multRef = Arg.segment(_numDofs, _constraintHandler.getNumOfTotalConstraints() );
    _constraintHandler.addConstraintGrad( geomRef, Dest, multRef );    
  }

};

//!
template<typename ConfiguratorType, typename HessianType, typename ConstraintHandlerType>
class LagrangeHessian : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType>{
    
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef  typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;
  
  const HessianType& _D2E;
  const ConstraintHandlerType& _constraintHandler;
  int _numDofs;
  
public:
  LagrangeHessian( const HessianType& D2E, const ConstraintHandlerType& constraintHandler, int numDofs ) 
  : _D2E(D2E), _constraintHandler(constraintHandler),_numDofs(numDofs) {}
  
  void apply( const VectorType& Arg, MatrixType& Dest ) const {
      
    int argSize = _numDofs + _constraintHandler.getNumOfTotalConstraints();  
    if( Arg.size() != argSize ){
        std::cerr << "Size arg in Hessian  = " << Arg.size() << ", where dim = " << argSize << std::endl;
      throw BasicException("LagrangeHessian::apply: argument has wrong size!");
    }
    
    if( (Dest.rows() != argSize) || (Dest.cols() != argSize) )
      Dest.resize( argSize, argSize );
    Dest.setZero();
      
    TripletListType tripletList;      
    Eigen::Ref<const VectorType> geomRef = Arg.segment(0, _numDofs );
    _D2E.pushTriplets( geomRef, tripletList );
    VectorType multiplierDummy;
    _constraintHandler.addConstraintHessian( geomRef, tripletList, multiplierDummy );
    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
  
};

#endif