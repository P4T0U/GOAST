// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Operators from elementary trigonometry on lengths and angles
 * \author Sassen
 */

#ifndef NRIC_TRIGONOMETRYOPERATORS_H
#define NRIC_TRIGONOMETRYOPERATORS_H

#include <goast/Core/Auxiliary.h>
#include <goast/Core/Topology.h>

#include <algorithm>
#include <array>

/**
 * \brief Operator determining triangle areas from edge lengths
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 *
 * For each triangle, this operator computes its area from its edge lengths via Heron's formula.
 * If set of edge lengths does not give a valid triangle, the corresponding entry is NaN.
 */
template<typename ConfiguratorType>
class TriangleAreaOp
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;

  const MeshTopologySaver &_topology;

public:
  explicit TriangleAreaOp( const MeshTopologySaver &topology ) : _topology( topology ) {}

  /**
   * \brief Evaluate operator
   * \param Arg Edges lengths (and dihedral angles) as vector of size at least \f$|E|\f$, only those first entries are considered
   * \param Dest Triangle areas as vector of size \f$|F|\f$, NaN entries for invalid triangles
   */
  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    if ( Arg.size() < _topology.getNumEdges())
      throw std::length_error( "TriangleAreaOp::apply(): Arg too small!" );

    Dest.resize( _topology.getNumFaces());
    Dest.setZero();

    for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); faceIdx++ ) {
      // Sort lengths such that evaluation is numerically stable
      // (cf. https://people.eecs.berkeley.edu/~wkahan/Triangle.pdf)
      std::array<RealType, 3> lengths = { Arg[_topology.getEdgeOfTriangle( faceIdx, 0 )],
                                          Arg[_topology.getEdgeOfTriangle( faceIdx, 1 )],
                                          Arg[_topology.getEdgeOfTriangle( faceIdx, 2 )] };
      std::sort( lengths.begin(), lengths.end(), std::greater<RealType>());

      Dest[faceIdx] = sqrt((lengths[0] + (lengths[1] + lengths[2])) * (lengths[2] - (lengths[0] - lengths[1])) *
                           (lengths[2] + (lengths[0] - lengths[1])) * (lengths[0] + (lengths[1] - lengths[2]))) / 4.;
    }
  }
};

/**
 * \brief Operator determining interior angles from edge lengths
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 *
 * For each triangle, this operator computes its interior angles from its edge lengths via Heron's formula.
 * If set of edge lengths does not give a valid triangle, the corresponding entry is NaN.
 */
template<typename ConfiguratorType>
class InteriorAngleOp
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;

  const MeshTopologySaver &_topology;

public:
  explicit InteriorAngleOp( const MeshTopologySaver &topology ) : _topology( topology ) {}

  /**
   * \brief Evaluate operator
   * \param Arg Edges lengths (and dihedral angles) as vector of size at least \f$|E|\f$, only those first entries are considered
   * \param Dest Interior angles as vector of size \f$3*|F|\f$, order such that local angle i is opposite of edge i; NaN entries for invalid triangles
   */
  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    if ( Arg.size() < _topology.getNumEdges())
      throw std::length_error( "InteriorAngleOp::apply(): Arg too small!" );

    Dest.resize( 3 * _topology.getNumFaces());
    Dest.setZero();

#ifdef GOAST_WITH_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); faceIdx++ ) {
      // Numerically stable implementation by appropriate ordering
      // cf. https://people.eecs.berkeley.edu/~wkahan/Triangle.pdf
      std::array<RealType, 3> lengths = { Arg[_topology.getEdgeOfTriangle( faceIdx, 0 )],
                                          Arg[_topology.getEdgeOfTriangle( faceIdx, 1 )],
                                          Arg[_topology.getEdgeOfTriangle( faceIdx, 2 )] };

      std::array<RealType, 3> sqLengths = { lengths[0] * lengths[0], lengths[1] * lengths[1], lengths[2] * lengths[2] };

      int a, b;
      RealType mu;

      // Opposite of c
      for ( int c : { 0, 1, 2 } ) {
        if ( lengths[(c + 1) % 3] > lengths[(c + 2) % 3] ) {
          a = (c + 1) % 3;
          b = (c + 2) % 3;
        }
        else {
          a = (c + 2) % 3;
          b = (c + 1) % 3;
        }

        if ( lengths[c] > lengths[b] )
          mu = lengths[b] - (lengths[a] - lengths[c]);
        else
          mu = lengths[c] - (lengths[a] - lengths[b]);

        Dest[3 * faceIdx + c] = 2 * atan( sqrt(((lengths[a] - lengths[b]) + lengths[c]) * mu /
                                               ((lengths[a] + (lengths[b] + lengths[c])) *
                                                ((lengths[a] - lengths[c]) + lengths[b]))));

//        Dest[3 * faceIdx + c] = std::acos((sqLengths[a] + sqLengths[b] - sqLengths[c]) / (2 * lengths[a] * lengths[b]));
      }
    }
  }

  int getTargetDimension() const override {
    return 3 * _topology.getNumFaces();
  }

};

/**
 * \brief Jacobian of operator determining interior angles from edge lengths
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 * \see InteriorAngleOp
 */
template<typename ConfiguratorType>
class InteriorAngleGradient
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver &_topology;

  TriangleAreaOp<ConfiguratorType> _areaOp;

public:
  explicit InteriorAngleGradient( const MeshTopologySaver &topology ) : _topology( topology ), _areaOp( _topology ) {}

  void apply( const VectorType &Arg, MatrixType &Dest ) const override {
    if ( Arg.size() < _topology.getNumEdges())
      throw std::length_error( "InteriorAngleGradient::apply(): Arg too small!" );

    Dest.resize( 3 * _topology.getNumFaces(), 2 * _topology.getNumEdges());
    Dest.reserve( 9 * _topology.getNumFaces());
    Dest.setZero();

    TripletListType tripletList;
    pushTriplets( Arg, tripletList );

    Dest.setFromTriplets( tripletList.begin(), tripletList.end());

  }

  void pushTriplets( const VectorType &Arg, TripletListType &Dest, RealType factor = 1., int rowOffset = 0,
                     int colOffset = 0 ) const {
    if ( Arg.size() < _topology.getNumEdges())
      throw std::length_error( "InteriorAngleGradient::pushTriplets(): Arg too small!" );

    Dest.reserve( 9 * _topology.getNumFaces());

    VectorType triangleAreas;
    _areaOp.apply( Arg, triangleAreas );

//#pragma omp parallel for schedule(static)
    for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); faceIdx++ ) {
      // Numerically stable implementation by appropriate ordering
      // cf. https://people.eecs.berkeley.edu/~wkahan/Triangle.pdf
      std::array<RealType, 3> lengths = { Arg[_topology.getEdgeOfTriangle( faceIdx, 0 )],
                                          Arg[_topology.getEdgeOfTriangle( faceIdx, 1 )],
                                          Arg[_topology.getEdgeOfTriangle( faceIdx, 2 )] };

      std::array<RealType, 3> sqLengths = { lengths[0] * lengths[0], lengths[1] * lengths[1], lengths[2] * lengths[2] };

      std::array<int, 3> e = { _topology.getEdgeOfTriangle( faceIdx, 0 ),
                               _topology.getEdgeOfTriangle( faceIdx, 1 ),
                               _topology.getEdgeOfTriangle( faceIdx, 2 ) };

      int a, b;

      // Opposite of c
      for ( int c : { 0, 1, 2 } ) {
        a = (c + 1) % 3;
        b = (c + 2) % 3;

        Dest.emplace_back( 3 * faceIdx + c + rowOffset,
                           e[a] + colOffset,
                           factor * (-sqLengths[a] + sqLengths[b] - sqLengths[c]) /
                           (4 * lengths[a] * triangleAreas[faceIdx]));

        Dest.emplace_back( 3 * faceIdx + c + rowOffset,
                           e[b] + colOffset,
                           factor * (sqLengths[a] - sqLengths[b] - sqLengths[c]) /
                           (4 * lengths[b] * triangleAreas[faceIdx]));

        Dest.emplace_back( 3 * faceIdx + c + rowOffset,
                           e[c] + colOffset,
                           factor * 2 * lengths[c] / (4 * triangleAreas[faceIdx]));

      }
    }
  }

  int getTargetDimension() const override {
    return 3 * _topology.getNumFaces();
  }

  int getNNZ() const override {
    return 9 * _topology.getNumFaces();
  }
};

/**
 * \brief Hessian of operator determining interior angles from edge lengths
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 * \see InteriorAngleOp
 */
template<typename ConfiguratorType>
class InteriorAngleHessian
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::TensorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;
  typedef typename ConfiguratorType::TensorType TensorType;

  const MeshTopologySaver &_topology;

  TriangleAreaOp<ConfiguratorType> _areaOp;

  const int _numEdges, _numFaces;

public:
  explicit InteriorAngleHessian( const MeshTopologySaver &topology ) : _topology( topology ), _areaOp( _topology ),
                                                                       _numFaces( _topology.getNumFaces()),
                                                                       _numEdges( _topology.getNumEdges()) {}

  void apply( const VectorType &Arg, TensorType &Dest ) const {
    apply( Arg, Dest, -1, 0 );
  }

  void apply( const VectorType &Arg, TensorType &Dest, int hessSize, int hessOffset ) const {
    if ( Arg.size() < _topology.getNumEdges())
      throw std::length_error( "InteriorAngleHessian::apply(): Arg too small!" );

    Dest.resize( 3 * _numFaces );

    std::vector<TripletListType> vertexTripletLists;
    vertexTripletLists.resize( 3 * _numFaces );

    setTriplets( Arg, vertexTripletLists, hessOffset );


#ifdef GOAST_WITH_OPENMP
#pragma omp parallel for
#endif
    for ( int faceIdx = 0; faceIdx < _numFaces; faceIdx++ ) {
      for ( int t : { 0, 1, 2 } ) {
        if ( hessSize == -1 )
          Dest[3 * faceIdx + t].resize( 2 * _numEdges, 2 * _numEdges );
        else
          Dest[3 * faceIdx + t].resize( hessSize, hessSize );

        Dest[3 * faceIdx + t].reserve( vertexTripletLists[3 * faceIdx + t].size());

        Dest[3 * faceIdx + t].setFromTriplets( vertexTripletLists[3 * faceIdx + t].begin(),
                                               vertexTripletLists[3 * faceIdx + t].end());
      }
    }

  }

  void pushTriplets( const VectorType &Arg, TripletListType &Dest, int hessOffset = 0 ) const {
    if ( Arg.size() < _topology.getNumEdges())
      throw std::length_error( "InteriorAngleHessian::pushTriplets(): Arg too small!" );

    std::vector<TripletListType> vertexTripletLists;
    vertexTripletLists.resize( 3 * _numFaces );

    setTriplets( Arg, vertexTripletLists, hessOffset );

    int totalNumTriplets = 0;
    for ( const auto &vTL : vertexTripletLists ) {
      totalNumTriplets += vTL.size();
    }

    Dest.reserve( totalNumTriplets );

    for ( const auto &vTL : vertexTripletLists ) {
      for ( const auto &trip : vTL )
        Dest.emplace_back( trip.row(), trip.col(), trip.value());
    }

  }

  void pushTriplets( const VectorType &Arg, TripletListType &Dest,
                     const VectorType &Lambda, int hessOffset = 0 ) const {
    if ( Arg.size() < _topology.getNumEdges())
      throw std::length_error( "InteriorAngleHessian::pushTriplets(): Arg too small!" );
    if ( Lambda.size() != 3 * _numFaces )
      throw std::length_error( "InteriorAngleHessian::pushTriplets(): Lambda too small!" );

    std::vector<TripletListType> vertexTripletLists;
    vertexTripletLists.resize( 3 * _numFaces );

    setTriplets( Arg, vertexTripletLists, hessOffset );

    int totalNumTriplets = 0;
    for ( const auto &vTL : vertexTripletLists ) {
      totalNumTriplets += vTL.size();
    }

    Dest.reserve( totalNumTriplets );

    for ( int i = 0; i < 3 * _numFaces; i++ ) {
      for ( const auto &trip : vertexTripletLists[i] )
        Dest.emplace_back( trip.row(), trip.col(), Lambda[i] * trip.value());
    }

  }

  void setTriplets( const VectorType &Arg, std::vector<std::vector<TripletType>> &Dest, int hessOffset = 0 ) const {
    if ( Arg.size() < 2 * _topology.getNumEdges())
      throw std::length_error( "SimplifiedIntegrabilityGradient::apply(): Arg too small!" );

    Dest.resize( 3 * _numFaces );
    for ( auto &vTL : Dest )
      vTL.clear();

    VectorType triangleAreas;
    _areaOp.apply( Arg, triangleAreas );

//#pragma omp parallel for schedule(static)
    for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); faceIdx++ ) {
      // Numerically stable implementation by appropriate ordering
      // cf. https://people.eecs.berkeley.edu/~wkahan/Triangle.pdf
      std::array<RealType, 3> lengths = { Arg[_topology.getEdgeOfTriangle( faceIdx, 0 )],
                                          Arg[_topology.getEdgeOfTriangle( faceIdx, 1 )],
                                          Arg[_topology.getEdgeOfTriangle( faceIdx, 2 )] };

      std::array<RealType, 3> sqLengths = { lengths[0] * lengths[0], lengths[1] * lengths[1], lengths[2] * lengths[2] };
      std::array<RealType, 3> quLengths = { sqLengths[0] * sqLengths[0], sqLengths[1] * sqLengths[1],
                                            sqLengths[2] * sqLengths[2] };
      std::array<RealType, 3> siLengths = { sqLengths[0] * quLengths[0], sqLengths[1] * quLengths[1],
                                            sqLengths[2] * quLengths[2] };

      std::array<int, 3> e = { _topology.getEdgeOfTriangle( faceIdx, 0 ),
                               _topology.getEdgeOfTriangle( faceIdx, 1 ),
                               _topology.getEdgeOfTriangle( faceIdx, 2 ) };

      int a, b;

      // Opposite of c
      for ( int c : { 0, 1, 2 } ) {
        a = (c + 1) % 3;
        b = (c + 2) % 3;

        RealType value;
        RealType denom = 64 * triangleAreas[faceIdx] * triangleAreas[faceIdx] * triangleAreas[faceIdx];

        // Diagonal elements
        value = siLengths[a] -
                3 * quLengths[a] * (sqLengths[b] - sqLengths[c]) +
                sqLengths[a] * (3 * quLengths[b] + 2 * sqLengths[b] * sqLengths[c] - 5 * quLengths[c]) -
                (sqLengths[b] - sqLengths[c]) * (sqLengths[b] - sqLengths[c]) * (sqLengths[b] - sqLengths[c]);
        value /= -sqLengths[a] * denom;
        Dest[3 * faceIdx + c].emplace_back( e[a] + hessOffset, e[a] + hessOffset, value );

        value = siLengths[b] -
                3 * quLengths[b] * (sqLengths[a] - sqLengths[c]) +
                sqLengths[b] * (3 * quLengths[a] + 2 * sqLengths[a] * sqLengths[c] - 5 * quLengths[c]) -
                (sqLengths[a] - sqLengths[c]) * (sqLengths[a] - sqLengths[c]) * (sqLengths[a] - sqLengths[c]);
        value /= -sqLengths[b] * denom;
        Dest[3 * faceIdx + c].emplace_back( e[b] + hessOffset, e[b] + hessOffset, value );

        value = -(2 * (sqLengths[a] - sqLengths[b]) * (sqLengths[a] - sqLengths[b]) - 2 * sqLengths[c] * sqLengths[c]);
        value /= denom;
        Dest[3 * faceIdx + c].emplace_back( e[c] + hessOffset, e[c] + hessOffset, value );

        // a, b
        value = 8 * lengths[a] * lengths[b] * sqLengths[c] / denom;
        Dest[3 * faceIdx + c].emplace_back( e[a] + hessOffset, e[b] + hessOffset, value );
        Dest[3 * faceIdx + c].emplace_back( e[b] + hessOffset, e[a] + hessOffset, value );

        // a, c
        value = 4 * lengths[a] * lengths[c] * (sqLengths[a] - sqLengths[b] - sqLengths[c]) / denom;
        Dest[3 * faceIdx + c].emplace_back( e[a] + hessOffset, e[c] + hessOffset, value );
        Dest[3 * faceIdx + c].emplace_back( e[c] + hessOffset, e[a] + hessOffset, value );

        // b, c
        value = 4 * lengths[b] * lengths[c] * (sqLengths[b] - sqLengths[a] - sqLengths[c]) / denom;
        Dest[3 * faceIdx + c].emplace_back( e[b] + hessOffset, e[c] + hessOffset, value );
        Dest[3 * faceIdx + c].emplace_back( e[c] + hessOffset, e[b] + hessOffset, value );
      }
    }

  }


  int getTargetDimension() const override {
    return 3 * _topology.getNumFaces();
  }

//  int getNNZ() const override {
//    return 9 * _topology.getNumFaces();
//  }
};

#endif //NRIC_TRIGONOMETRYOPERATORS_H
