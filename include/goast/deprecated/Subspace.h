// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by sassen on 22.06.18.
//

#ifndef NRIC_SUBSPACE_H
#define NRIC_SUBSPACE_H


#include <chrono>
#include <ctime>
#include <goast/Core.h>

//! \brief Distance of vector to affine subspace given by orthonormal basis
//! \author Sassen
//!
//! Given basis u_0,...,u_J, reference point \bar{z}, and weights of scalar product M this implements
//! dist^2(z, {\bar{z} + U\omega}) = (z - \bar{z})^T (Id - UU^T M)^T \, M\, (Id - UU^T M) (z - \bar{z})
//! for a vector z.
template<typename ConfiguratorType>
class SubspaceDistanceFunctional
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;

  const VectorType &_referencePoint;
  const FullMatrixType &_basis;
  const VectorType &_masses;
  FullMatrixType _distOp;

  const int _numDOFs, _dimension;

public:
  SubspaceDistanceFunctional( const VectorType &masses,
                              const VectorType &referencePoint,
                              const FullMatrixType &Basis )
          : _referencePoint( referencePoint ), _basis( Basis ), _masses ( masses ), _dimension( Basis.cols() ),
            _numDOFs( referencePoint.size() )  {
    if ( referencePoint.size() != Basis.rows())
      throw BasicException( "SubspaceDistanceFunctional: size of basis and reference point don't match!" );

    if ( referencePoint.size() != masses.size())
      throw BasicException( "SubspaceDistanceFunctional: size of basis and scalar product weights/masses don't match" );

    // Id - U U^T M =: A
//    _distOp = FullMatrixType::Identity(_numDOFs, _numDOFs) - _basis * _basis.transpose() * _masses.asDiagonal();
    // A^T M A
//    _distOp = _distOp.transpose() * _masses.asDiagonal() * _distOp;


  }

  void apply( const VectorType &Arg, RealType &Dest ) const override {
    if ( Arg.size() != _numDOFs )
      throw BasicException( "SubspaceDistanceFunctional::apply: wrong size of dofs!" );

    Dest.setZero();

    VectorType projectionCoeff = _basis.transpose() * _masses.asDiagonal() * (Arg - _referencePoint);
    VectorType diff = (Arg - _referencePoint) - _basis * projectionCoeff;

    Dest[0] = diff.transpose() * _masses.asDiagonal() * diff;
  }
};

template<typename ConfiguratorType>
class SubspaceDistanceGradient
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;

  const VectorType &_referencePoint;
  const FullMatrixType &_basis;
  const VectorType &_masses;
  FullMatrixType _distOp;

  const int _numDOFs, _dimension;

public:
  SubspaceDistanceGradient( const VectorType &masses,
                            const VectorType &referencePoint,
                            const FullMatrixType &Basis )
          : _referencePoint( referencePoint ), _basis( Basis ), _masses ( masses ), _dimension( Basis.cols() ),
            _numDOFs( referencePoint.size() )  {
    if ( referencePoint.size() != Basis.rows())
      throw BasicException( "SubspaceDistanceGradient: size of basis and reference point don't match!" );

    if ( referencePoint.size() != masses.size())
      throw BasicException( "SubspaceDistanceGradient: size of basis and scalar product weights/masses don't match" );

    // Id - U U^T M =: A
//    _distOp = FullMatrixType::Identity(_numDOFs, _numDOFs) - _basis * _basis.transpose() * _masses.asDiagonal();
    // A^T M A
//    _distOp = _distOp.transpose() * _masses.asDiagonal() * _distOp;


  }

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    if ( Arg.size() != _numDOFs )
      throw BasicException( "SubspaceDistanceGradient::apply: wrong size of dofs!" );

    if ( Dest.size() != _numDOFs )
      Dest.resize( _numDOFs );

    Dest.setZero();

    VectorType projectionCoeff = _basis.transpose() * _masses.asDiagonal() * (Arg - _referencePoint);
    VectorType diff = (Arg - _referencePoint) - _basis * projectionCoeff;

    Dest = 2 * _masses.asDiagonal() * diff;
  }
};

template<typename ConfiguratorType>
class SubspaceDistanceHessian
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;

  const VectorType &_referencePoint;
  const FullMatrixType &_basis;
  const VectorType &_masses;

  const int _numDOFs, _dimension;

public:
  SubspaceDistanceHessian( const VectorType &masses,
                            const VectorType &referencePoint,
                            const FullMatrixType &Basis )
          : _referencePoint( referencePoint ), _basis( Basis ), _masses ( masses ), _dimension( Basis.cols() ),
            _numDOFs( referencePoint.size() )  {
    if ( referencePoint.size() != Basis.rows())
      throw BasicException( "SubspaceDistanceHessian: size of basis and reference point don't match!" );

    if ( referencePoint.size() != masses.size())
      throw BasicException( "SubspaceDistanceHessian: size of basis and scalar product weights/masses don't match" );

    // Id - U U^T M =: A
//    _distOp = FullMatrixType::Identity(_numDOFs, _numDOFs) - _basis * _basis.transpose() * _masses.asDiagonal();
    // A^T M A
//    _distOp = _distOp.transpose() * _masses.asDiagonal() * _distOp;

  }

  void apply( const VectorType &Arg, MatrixType &Dest ) const override {
    if ( Arg.size() != _numDOFs )
      throw BasicException( "SubspaceDistanceHessian::apply: wrong size of dofs!" );

    if ( Dest.cols() != _numDOFs || Dest.rows() != _numDOFs)
      Dest.resize( _numDOFs, _numDOFs );

    MatrixType sparseMasses;
    sparseMasses.resize( _numDOFs, _numDOFs );
    sparseMasses.setIdentity();
    sparseMasses.diagonal() = _masses.asDiagonal();

    MatrixType sparseBasis = _basis.sparseView();

    Dest.setZero();

    Dest = -1. * sparseMasses * sparseBasis * sparseBasis.transpose() * sparseMasses;
    Dest += sparseMasses;
    Dest *= 2;
  }
};

#endif //NRIC_SUBSPACE_H
