// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Triangle inequaliy operators for edge lengths (and angles)
 * \author Sassen
 */

#ifndef NRIC_TRIANGLEINEQUALITY_H
#define NRIC_TRIANGLEINEQUALITY_H

#include <goast/Core/Topology.h>

/**
 * \brief Operator evaluating triangle inequality condition on edge lengths
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 *
 * For each triangle, this operator computes \f$a + b + c - 2\max(a,b,c)\f$, where a, b, and c are the three edge
 * lengths of the triangle and the result is stored in a vector. Fulfilling the triangle inequality for all triangles
 * is then equivalent to all entries of this vector being positive.
 */
template<typename ConfiguratorType>
class TriangleInequalityOp
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;

  const MeshTopologySaver &_topology;

public:
  explicit TriangleInequalityOp( const MeshTopologySaver &topology ) : _topology( topology ) {}

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    if ( Arg.size() < _topology.getNumEdges())
      throw std::length_error( "TriangleInequalityOp::apply(): Arg too small!" );

    Dest.resize( _topology.getNumFaces());
    Dest.setZero();

    for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); faceIdx++ ) {
      RealType a = Arg[_topology.getEdgeOfTriangle( faceIdx, 0 )];
      RealType b = Arg[_topology.getEdgeOfTriangle( faceIdx, 1 )];
      RealType c = Arg[_topology.getEdgeOfTriangle( faceIdx, 2 )];

      Dest[faceIdx] = a + b + c - 2 * std::max( { a, b, c } );
    }
  }
};

/**
 * \brief Derivative of triangle inequality operator
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 *
 * For each triangle, this operator computes the derivative of \f$a + b + c - 2\max(a,b,c)\f$, where a, b, and c are
 * the three edge lengths of the triangle and the result is stored in a vector. As the maximum is not continuously
 * differentiable at configurations where the two longest edge lengths are equal, at those points we take the
 * subgradient given by the first (in order of local indexing) edge length with maximal length.
 */
template<typename ConfiguratorType>
class TriangleInequalityGradient
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver &_topology;

public:
  explicit TriangleInequalityGradient( const MeshTopologySaver &topology ) : _topology( topology ) {}

  void apply( const VectorType &Arg, SparseMatrixType &Dest ) const override {
    if ( Arg.size() != 2 * _topology.getNumEdges())
      throw std::length_error( "TriangleInequalityOp::apply(): Arg wrong size!" );

    if (( Dest.cols() != _topology.getNumEdges()) || ( Dest.rows() != _topology.getNumFaces()))
      Dest.resize( _topology.getNumFaces(), 2 * _topology.getNumEdges());

    Dest.setZero();

    TripletListType tripletList;
    tripletList.reserve( 4 * _topology.getNumFaces());
    pushTriplets( Arg, tripletList );

    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend());
  }

  void pushTriplets( const VectorType &Arg, TripletListType &tripletList ) const override {
    for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); faceIdx++ ) {
      RealType a = Arg[_topology.getEdgeOfTriangle( faceIdx, 0 )];
      RealType b = Arg[_topology.getEdgeOfTriangle( faceIdx, 1 )];
      RealType c = Arg[_topology.getEdgeOfTriangle( faceIdx, 2 )];

      tripletList.emplace_back( faceIdx, _topology.getEdgeOfTriangle( faceIdx, 0 ), 1 );
      tripletList.emplace_back( faceIdx, _topology.getEdgeOfTriangle( faceIdx, 1 ), 1 );
      tripletList.emplace_back( faceIdx, _topology.getEdgeOfTriangle( faceIdx, 2 ), 1 );

      RealType max_length = std::max( { a, b, c } );

      if ( a == max_length ) {
        tripletList.emplace_back( faceIdx, _topology.getEdgeOfTriangle( faceIdx, 0 ), -2 );
      }
      else if ( b == max_length ) {
        tripletList.emplace_back( faceIdx, _topology.getEdgeOfTriangle( faceIdx, 1 ), -2 );
      }
      else {
        tripletList.emplace_back( faceIdx, _topology.getEdgeOfTriangle( faceIdx, 2 ), -2 );
      }
    }
  }
};

/**
 * \brief Operator evaluating triangle inequality condition on edge lengths
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 *
 * For each triangle, this operator computes \f$(a+b-c, a-b+c, -a+b+c)\f$, where a, b, and c are the three edge lengths
 * of the triangle and the results are combined in a single vector. Fulfilling the triangle inequality for all triangles
 * is then equivalent to all entries of this vector being positive.
 */
template<typename ConfiguratorType>
class TriangleInequalityTripleOp
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;

  const MeshTopologySaver &_topology;

public:
  explicit TriangleInequalityTripleOp( const MeshTopologySaver &topology ) : _topology( topology ) {}

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    if ( Arg.size() < _topology.getNumEdges())
      throw std::length_error( "TriangleInequalityOp::apply(): Arg too small!" );

    Dest.resize( 3 * _topology.getNumFaces());
    Dest.setZero();

    for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); faceIdx++ ) {
      RealType a = Arg[_topology.getEdgeOfTriangle( faceIdx, 0 )];
      RealType b = Arg[_topology.getEdgeOfTriangle( faceIdx, 1 )];
      RealType c = Arg[_topology.getEdgeOfTriangle( faceIdx, 2 )];

      Dest[3 * faceIdx] = a + b - c;
      Dest[3 * faceIdx + 1] = a - b + c;
      Dest[3 * faceIdx + 2] = -a + b + c;
    }
  }

  int getTargetDimension() const override {
    return 3 * _topology.getNumFaces();
  }
};

/**
 * \brief Derivative of triple triangle inequality operator
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 */
template<typename ConfiguratorType>
class TriangleInequalityTripleGradient
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver &_topology;

public:
  explicit TriangleInequalityTripleGradient( const MeshTopologySaver &topology ) : _topology( topology ) {}

  void apply( const VectorType &Arg, SparseMatrixType &Dest ) const override {
    if ( Arg.size() != 2 * _topology.getNumEdges())
      throw std::length_error( "TriangleInequalityOp::apply(): Arg wrong size!" );

    if (( Dest.cols() != 2 * _topology.getNumEdges()) || ( Dest.rows() != 3 * _topology.getNumFaces()))
      Dest.resize( 3 * _topology.getNumFaces(), 2 * _topology.getNumEdges());

    Dest.setZero();

    TripletListType tripletList;
    tripletList.reserve( 9 * _topology.getNumFaces());
    pushTriplets( Arg, tripletList );

    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend());
  }

  void pushTriplets( const VectorType &Arg, TripletListType &Dest ) const override {
    pushTriplets( Arg, Dest, 1., 0, 0 );
  }

  void pushTriplets( const VectorType &Arg, TripletListType &tripletList, RealType factor, int rowOffset,
                     int colOffset ) const override {
    for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); faceIdx++ ) {
      int aIdx = _topology.getEdgeOfTriangle( faceIdx, 0 );
      int bIdx = _topology.getEdgeOfTriangle( faceIdx, 1 );
      int cIdx = _topology.getEdgeOfTriangle( faceIdx, 2 );

      tripletList.emplace_back( 3 * faceIdx + rowOffset, aIdx + colOffset, factor );
      tripletList.emplace_back( 3 * faceIdx + rowOffset, bIdx + colOffset, factor );
      tripletList.emplace_back( 3 * faceIdx + rowOffset, cIdx + colOffset, -factor );

      tripletList.emplace_back( 3 * faceIdx + 1 + rowOffset, aIdx + colOffset, factor );
      tripletList.emplace_back( 3 * faceIdx + 1 + rowOffset, bIdx + colOffset, -factor );
      tripletList.emplace_back( 3 * faceIdx + 1 + rowOffset, cIdx + colOffset, factor );

      tripletList.emplace_back( 3 * faceIdx + 2 + rowOffset, aIdx + colOffset, -factor );
      tripletList.emplace_back( 3 * faceIdx + 2 + rowOffset, bIdx + colOffset, factor );
      tripletList.emplace_back( 3 * faceIdx + 2 + rowOffset, cIdx + colOffset, factor );
    }
  }

  int getTargetDimension() const override {
    return 3 * _topology.getNumFaces();
  }

  int getNNZ() const override {
    return 9 * _topology.getNumFaces();
  }
};

#endif //NRIC_TRIANGLEINEQUALITY_H
