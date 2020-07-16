// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Operators determining if a set of lenghts and angles is admissible or not
 * \author Sassen
 *
 * Based on Yue Wang, B. Liu, and Y. Tong. Linear surface reconstruction from discrete fundamental forms on triangle
 * meshes. Comput. Graph. Forum, 31:2277–2287, 2012.
 *
 * \note This currently relies directly on the Axis-Angle and Quaternion implementations from Eigen
 */

#ifndef NRIC_ADMISSIBILITY_H
#define NRIC_ADMISSIBILITY_H

#include <utility>

#include <goast/Core/Topology.h>

#include "TriangleInequality.h"
#include "IntegrabilityAngle.h"
#include "IntegrabilityEuler.h"
#include "IntegrabilityQuaternions.h"

/**
 * \brief Operator combining discrete integrability and triangle inequality map
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 */
template<typename ConfiguratorType>
class AdmissibilityOp
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;

  const MeshTopologySaver &_topology;
  const TriangleInequalityTripleOp<ConfiguratorType> _triOp;
  const QuaternionIntegrabilityOp<ConfiguratorType> _intOp;

  const int _triDim;
  const int _intDim;

public:
  explicit AdmissibilityOp( const MeshTopologySaver &topology ) :
          _topology( topology ),
          _triOp( topology ),
          _intOp( topology ),
          _triDim( _triOp.getTargetDimension()),
          _intDim( _intOp.getTargetDimension()) {}

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    if ( Arg.size() < _topology.getNumEdges())
      throw std::length_error( "AdmissibilityOp::apply(): Arg too small!" );

    Dest.resize( _intDim + _triDim );
    Dest.setZero();

    VectorType intPart, triPart;

    _intOp.apply( Arg, intPart );
    _triOp.apply( Arg, triPart );

    Dest.segment( 0, _intDim ) = intPart;
    Dest.segment( _intDim, _triDim ) = triPart;

  }

  int getTargetDimension() const override {
    return _intDim + _triDim;
  }
};

/**
 * \brief Gradient of operator combining gradients of simplified discrete integrability and triangle inequality map
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 * \sa AdmissibilityOp
 */
template<typename ConfiguratorType>
class AdmissibilityGrad
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver &_topology;

  TriangleInequalityTripleGradient<ConfiguratorType> _triGrad;
  QuaternionIntegrabilityGradient<ConfiguratorType> _intGrad;

  const int _triDim;
  const int _intDim;
public:
  explicit AdmissibilityGrad( const MeshTopologySaver &topology ) :
          _topology( topology ), _triGrad( topology ),
          _intGrad( topology ),
          _triDim( _triGrad.getTargetDimension()),
          _intDim( _intGrad.getTargetDimension()) {}

  void apply( const VectorType &Arg, SparseMatrixType &Dest ) const override {
    if ( Arg.size() != 2 * _topology.getNumEdges())
      throw std::length_error( "AdmissibilityGrad::apply(): Arg wrong size!" );

    if (( Dest.cols() != 2 * _topology.getNumEdges()) || ( Dest.rows() != _intDim + _triDim ))
      Dest.resize( _intDim + _triDim, 2 * _topology.getNumEdges());

    Dest.setZero();

    TripletListType tripletList;
    pushTriplets( Arg, tripletList );

    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend());
  }

  void pushTriplets( const VectorType &Arg, TripletListType &tripletList ) const {
    TripletListType intTriplets;
    _intGrad.pushTriplets( Arg, intTriplets );

    TripletListType triTriplets;
    _triGrad.pushTriplets( Arg, triTriplets, 1., _intDim );

    tripletList.insert( tripletList.end(), intTriplets.begin(), intTriplets.end());
    tripletList.insert( tripletList.end(), triTriplets.begin(), triTriplets.end());
  }

  int getTargetDimension() const override {
    return _intDim + _triDim;
  }

  int getNNZ() const override {
    return _intGrad.getNNZ() + _triGrad.getNNZ();
  }
};

#endif //ADMISSIBILITY_H
