// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef __BASEOPINTERFACE_H
#define __BASEOPINTERFACE_H

//== INCLUDES =================================================================
#include "Auxiliary.h"

//==========================================================================================================
// ABSTRACT OPERATOR INTERFACES
//==========================================================================================================
template <typename _DomainType, typename _RangeType = _DomainType>
class BaseOp {
  typedef Eigen::Triplet<double> _TripletType;
  typedef std::vector<_TripletType> _TripletListType;
  typedef typename _DomainType::Scalar _RealType;

public:
  typedef _DomainType DomainType;
  typedef _RangeType  RangeType;

  BaseOp() { }

  // Destroy polymorphic BaseOps correctly, important!
  virtual ~BaseOp () {}

  virtual void apply ( const DomainType &Arg, RangeType &Dest ) const = 0;

  virtual void applyAdd ( const DomainType &Arg, RangeType &Dest ) const {
    RangeType Update;
    apply( Arg, Update );
    Dest += Update;    
  }

  virtual RangeType operator() ( const DomainType &Arg ) const {
    RangeType Update;
    apply( Arg, Update );
    return Update;
  }
  
  virtual void applyTransposed ( const DomainType &Arg, RangeType &Dest ) const {
    throw std::logic_error("BaseOp::applyTransposed(): Unimplemented function! Should be provided in derived class!");
  }

  virtual void apply( const DomainType &Arg, RangeType &Dest, int hessSize, int hessOffset ) const {
    throw std::logic_error("BaseOp::applyTransposed(): Unimplemented function! Should be provided in derived class!");
  }

  virtual void pushTriplets ( const DomainType &Arg, _TripletListType &Dest ) const {
    throw std::logic_error("BaseOp::pushTriplets(): Unimplemented function! Should be provided in derived class!");
  }

  virtual void pushTriplets ( const DomainType &Arg, _TripletListType &Dest, const DomainType &Lambda ) const {
    throw std::logic_error("BaseOp::pushTriplets(): Unimplemented function! Should be provided in derived class!");
  }

  virtual void pushTriplets ( const DomainType &Arg, _TripletListType &Dest, _RealType factor,
          int rowOffset, int colOffset) const {
    throw std::logic_error("BaseOp::pushTriplets(): Unimplemented function! Should be provided in derived class!");
  }

  virtual void pushTriplets( const DomainType &Arg, _TripletListType &Dest, const DomainType &Lambda,
                             int hessOffset ) const {
    throw std::logic_error("BaseOp::pushTriplets(): Unimplemented function! Should be provided in derived class!");
  }

  virtual void pushTriplets( const DomainType &Arg, _TripletListType &Dest, int hessOffset) const {
    throw std::logic_error("BaseOp::pushTriplets(): Unimplemented function! Should be provided in derived class!");
  }

  virtual void setTriplets( const DomainType &Arg, std::vector<_TripletListType> &Dest, int hessOffset = 0 ) const {
    throw std::logic_error("BaseOp::setTriplets(): Unimplemented function! Should be provided in derived class!");
  }



  virtual int getTargetDimension() const {
    throw std::logic_error("BaseOp::getTargetDimension(): Unimplemented function! "
                           "Should be provided in derived class!");
  }

  virtual int getNNZ() const {
    throw std::logic_error("BaseOp::getNNZ(): Unimplemented function! Should be provided in derived class!");
  }

};
 

//==========================================================================================================
// WRAPPER CLASSES FOR GAUSS-NEWTON FUNCTIONALS
//==========================================================================================================
/**
 * Square of a vector-valued functional: $F[x] = \alpha f[x]^T f[x]$
 * \tparam ConfiguratorType Underlying types for scales, vectors, matrices, etc.
 * \author Sassen
 */
template<typename ConfiguratorType>
class SquaredFunctional
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;

  const RealType _alpha;
  const BaseOp<typename ConfiguratorType::VectorType> &_F;


public:
  /**
   * Construct giving functional
   * \param F functional, which will be squared
   * \param alpha scalar weight for resulting functional, default is 1/2
   */
  SquaredFunctional( const BaseOp<typename ConfiguratorType::VectorType> &F, RealType alpha = 0.5 )
          : _F( F ), _alpha( alpha ) {}


  /**
   * Apply squared functional
   * \param Arg vector to evaluate at
   * \param Dest value
   */
  void apply( const VectorType &Arg, RealType &Dest ) const {
    VectorType linearPart;
    _F.apply( Arg, linearPart );
    Dest[0] = _alpha * linearPart.squaredNorm(); // 0.5 *
  }

};

/**
 * Derivative of squared vector-value functional: DF[x] = 2 \alpha Df[x]^T f[x]
 * \tparam ConfiguratorType Underlying types for scales, vectors, matrices, etc.
 * \author Sassen
 */
template<typename ConfiguratorType>
class SquaredDerivative : public BaseOp<typename ConfiguratorType::VectorType> {

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;

  const BaseOp<typename ConfiguratorType::VectorType> &_F;
  const BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> &_DF;
  const RealType _alpha;

public:
  /**
   * Construct giving functional and its derivative
   * \param F functional, which will be squared
   * \param DF derivative of functional
   * \param alpha scalar weight for resulting functional, default is 1/2
   */
  SquaredDerivative( const BaseOp<typename ConfiguratorType::VectorType> &F,
                     const BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> &DF,
                     RealType alpha = 0.5 )
          : _F( F ),
            _DF( DF ), _alpha( alpha ) {}

  /**
   * Apply derivative of squared functional
   * \param Arg vector to evaluate at
   * \param Dest value
   */
  void apply( const VectorType &Arg, VectorType &Dest ) const {
    VectorType linearPart;
    _F.apply( Arg, linearPart );

    MatrixType JacobianOfLinearPart;
    _DF.apply( Arg, JacobianOfLinearPart );

    Dest = _alpha * 2 * JacobianOfLinearPart.transpose() * linearPart;
  }

};

/**
 * Approximated hessian of squared vector-value functional, neglecting second-derivatives of f:
 * D^2F[x] \approx 2 \alpha Df[x]^T Df[x]
 * \tparam ConfiguratorType Underlying types for scales, vectors, matrices, etc.
 * \author Sassen
 */
template<typename ConfiguratorType>
class ReducedSquaredHessian
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  const BaseOp<typename ConfiguratorType::VectorType> &_F;
  const BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> &_DF;
  const RealType _alpha;

public:
  /**
   * Construct giving functional and its derivative
   * \param F functional, which will be squared
   * \param DF derivative of functional
   * \param alpha scalar weight for resulting functional, default is 1/2
   */
  ReducedSquaredHessian( const BaseOp<typename ConfiguratorType::VectorType> &F,
                         const BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> &DF,
                         RealType alpha = 0.5 )
          : _F( F ), _DF( DF ), _alpha( alpha ) {}

  /**
   * Apply approximated hessian of squared functional
   * \param Arg vector to evaluate at
   * \param Dest value
   */
  void apply( const VectorType &Arg, MatrixType &Dest ) const {
    MatrixType JacobianOfLinearPart;
    _DF.apply( Arg, JacobianOfLinearPart );

    Dest = _alpha * 2 * JacobianOfLinearPart.transpose() * JacobianOfLinearPart;
  }

  void pushTriplets( const VectorType &Arg, TripletListType &Dest ) const override {
    MatrixType JacobianOfLinearPart;
    _DF.apply( Arg, JacobianOfLinearPart );

    MatrixType JacobianSquared( Arg.size(), Arg.size());
    JacobianSquared = _alpha * 2 * (JacobianOfLinearPart.transpose() * JacobianOfLinearPart);

    Dest.reserve( Dest.size() + JacobianSquared.nonZeros());

    for ( int k = 0; k < JacobianSquared.outerSize(); ++k )
      for ( typename MatrixType::InnerIterator it( JacobianSquared, k ); it; ++it )
        Dest.emplace_back( it.row(), it.col(), it.value());

  }

};

/**
 * Hessian of squared vector-value functional, neglecting second-derivatives of f:
 * D^2F[x] = 2 \alpha (Df[x]^T Df[x] +  D^2f[x] f[x])
 * \tparam ConfiguratorType Underlying types for scales, vectors, matrices, etc.
 * \author Sassen
 */
template<typename ConfiguratorType>
class SquaredHessian
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TensorType TensorType;

  const BaseOp<typename ConfiguratorType::VectorType> &_F;
  const BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> &_DF;
  const BaseOp<typename ConfiguratorType::VectorType, TensorType> &_D2F;
  const RealType _alpha;

public:
  /**
   * Construct giving functional and its derivative
   * \param F functional, which will be squared
   * \param DF derivative of functional
   * \param alpha scalar weight for resulting functional, default is 1/2
   */
  SquaredHessian( const BaseOp<typename ConfiguratorType::VectorType> &F,
                  const BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> &DF,
                  const BaseOp<typename ConfiguratorType::VectorType, TensorType> &D2F,
                  RealType alpha = 0.5 )
          : _F( F ), _DF( DF ), _D2F( D2F ), _alpha( alpha ) {}

  /**
   * Apply approximated hessian of squared functional
   * \param Arg vector to evaluate at
   * \param Dest value
   */
  void apply( const VectorType &Arg, MatrixType &Dest ) const {
    VectorType linearPart;
    _F.apply( Arg, linearPart );

    MatrixType JacobianOfLinearPart;
    _DF.apply( Arg, JacobianOfLinearPart );

    Dest = _alpha * 2 * JacobianOfLinearPart.transpose() * JacobianOfLinearPart;

    TensorType HessianOfLinearPart;
    _D2F.apply( Arg, HessianOfLinearPart );

    HessianOfLinearPart.applyVector( linearPart, JacobianOfLinearPart );

    Dest += _alpha * 2 * JacobianOfLinearPart;
  }

};
#endif