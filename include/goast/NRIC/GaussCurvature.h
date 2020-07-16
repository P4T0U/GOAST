// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Header for the angle-defect Gauss curvature on NRIC
 * \author Sassen
 */

#ifndef NRIC_GAUSSCURVATURE_H
#define NRIC_GAUSSCURVATURE_H

#include <goast/Core.h>

#include <algorithm>
#include <array>

#include "TrigonometryOperators.h"

/**
 * \brief Operator determining angle-defect Gauss curvature
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 *
 */
template<typename ConfiguratorType>
class GaussCurvatureOp
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;

  const MeshTopologySaver &_topology;
  InteriorAngleOp<ConfiguratorType> _gammaOp;


  std::vector<std::vector<int>> associatedAngles;

public:
  int _numInteriorVertices;
  std::vector<int> interiorVertices;

  explicit GaussCurvatureOp( const MeshTopologySaver &topology ) : _topology( topology ), _gammaOp( _topology ) {
    _numInteriorVertices = 0;
    for ( int vertexIdx = 0; vertexIdx < _topology.getNumVertices(); vertexIdx++ ) {

      bool interiorVertex = true;
      for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); edgeIdx++ ) {
        if ( _topology.getAdjacentNodeOfEdge( edgeIdx, 0 ) == vertexIdx ||
             _topology.getAdjacentNodeOfEdge( edgeIdx, 1 ) == vertexIdx ) {
          int f0 = _topology.getOppositeNodeOfEdge( edgeIdx, 0 );
          int f1 = _topology.getOppositeNodeOfEdge( edgeIdx, 1 );

          if ( f0 == -1 || f1 == -1 ) {
            interiorVertex = false;
            break;
          }
        }
      }

      if ( interiorVertex ) {
        _numInteriorVertices++;
        interiorVertices.push_back( vertexIdx );
      }
    }

    associatedAngles.resize( _topology.getNumVertices());

    for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); faceIdx++ ) {
      std::array<int, 3> v = { _topology.getNodeOfTriangle( faceIdx, 0 ),
                               _topology.getNodeOfTriangle( faceIdx, 1 ),
                               _topology.getNodeOfTriangle( faceIdx, 2 ) };

      for ( int i : { 0, 1, 2 } ) {
        associatedAngles[v[i]].push_back( 3 * faceIdx + i );
      }
    }
  }

  /**
   * \brief Evaluate operator
   * \param Arg Edges lengths (and dihedral angles) as vector of size at least \f$|E|\f$, only those first entries are considered
   * \param Dest Vertex-wise Gauss curvature i.e. angle defect
   */
  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    if ( Arg.size() < _topology.getNumEdges())
      throw std::length_error( "InteriorAngleOp::apply(): Arg too small!" );

    Dest.resize( _numInteriorVertices );
    Dest.setConstant( 2 * M_PI );

    VectorType gamma;
    _gammaOp.apply( Arg, gamma );

#pragma omp parallel for schedule(static)
    for ( int ivIdx = 0; ivIdx < _numInteriorVertices; ivIdx++ ) {
      int vertexIdx = interiorVertices[ivIdx];
      for ( const auto &angleIdx : associatedAngles[vertexIdx] )
        Dest[ivIdx] -= gamma[angleIdx];
    }

  }

  int getTargetDimension() const override {
    return _numInteriorVertices;
  }
};

/**
 * \brief Jacobian of angle-defect Gauss curvature
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 * \see GaussCurvatureOp
 */
template<typename ConfiguratorType>
class GaussCurvatureGradient
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
  InteriorAngleOp<ConfiguratorType> _gammaOp;
  InteriorAngleGradient<ConfiguratorType> _gammaGrad;

  std::vector<std::vector<int>> associatedAngles;
  std::vector<int> associatedVertex;

  std::vector<int> interiorVerticesIndices;

  int _numNonzeros;
public:
  int _numInteriorVertices;
  std::vector<int> interiorVertices;

  /**
   * \brief Construct operator
   * \param topology class containing information on the topology (i.e. connectivity) of the mesh.
   *
   * \todo Move computation of the needed topological information to a separate topology class
   */
  explicit GaussCurvatureGradient( const MeshTopologySaver &topology ) : _topology( topology ), _gammaOp( _topology ),
                                                                         _gammaGrad( _topology ) {
    _numInteriorVertices = 0;
    interiorVerticesIndices.resize( _topology.getNumVertices(), -1 );
    for ( int vertexIdx = 0; vertexIdx < _topology.getNumVertices(); vertexIdx++ ) {

      bool interiorVertex = true;
      for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); edgeIdx++ ) {
        if ( _topology.getAdjacentNodeOfEdge( edgeIdx, 0 ) == vertexIdx ||
             _topology.getAdjacentNodeOfEdge( edgeIdx, 1 ) == vertexIdx ) {
          int f0 = _topology.getOppositeNodeOfEdge( edgeIdx, 0 );
          int f1 = _topology.getOppositeNodeOfEdge( edgeIdx, 1 );

          if ( f0 == -1 || f1 == -1 ) {
            interiorVertex = false;
            break;
          }
        }
      }

      if ( interiorVertex ) {
        interiorVerticesIndices[vertexIdx] = _numInteriorVertices;
        _numInteriorVertices++;
        interiorVertices.push_back( vertexIdx );
      }
    }

    associatedAngles.resize( _topology.getNumVertices());
    associatedVertex.resize( 3 * _topology.getNumFaces());


    for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); faceIdx++ ) {
      std::array<int, 3> v = { _topology.getNodeOfTriangle( faceIdx, 0 ),
                               _topology.getNodeOfTriangle( faceIdx, 1 ),
                               _topology.getNodeOfTriangle( faceIdx, 2 ) };

      for ( int i : { 0, 1, 2 } ) {
        associatedAngles[v[i]].push_back( 3 * faceIdx + i );
        associatedVertex[3 * faceIdx + i] = v[i];
      }
    }


    _numNonzeros = 0;
    for ( int ivIdx = 0; ivIdx < _numInteriorVertices; ivIdx++ ) {
      int vertexIdx = interiorVertices[ivIdx];
      _numNonzeros += 2 * associatedAngles[vertexIdx].size();

    }

  }

  /**
   * \brief Evaluate Jacobian
   * \param Arg Edges lengths (and dihedral angles) as vector of size at least \f$|E|\f$, only those first entries are considered
   * \param Dest Jacobian of vertex-wise Gauss curvature
   */
  void apply( const VectorType &Arg, MatrixType &Dest ) const override {
    if ( Arg.size() < _topology.getNumEdges())
      throw std::length_error( "InteriorAngleOp::apply(): Arg too small!" );

    Dest.resize( _numInteriorVertices, 2 * _topology.getNumEdges());
    Dest.reserve( 9 * _topology.getNumFaces());
    Dest.setZero();

    TripletListType tripletList;
    pushTriplets( Arg, tripletList );

    Dest.setFromTriplets( tripletList.begin(), tripletList.end());
  }

  /**
   * \brief Evaluate derivative to triplet list
   * \tparam transposed determines whether to return transposed Jacobian
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[out] Dest Jacobian as triplet list
   * \param[in] factor scaling of the Jacobian
   * \param[in] rowOffset Row offset of the triplets
   * \param[in] colOffset Column offset of the triplets
   */
  void pushTriplets( const VectorType &Arg, TripletListType &Dest, RealType factor = 1., int rowOffset = 0,
                     int colOffset = 0 ) const {
    if ( Arg.size() < _topology.getNumEdges())
      throw std::length_error( "InteriorAngleGradient::pushTriplets(): Arg too small!" );

    Dest.reserve( 9 * _topology.getNumFaces());

    TripletListType angleTriplets;
    _gammaGrad.pushTriplets( Arg, angleTriplets, factor );

    for ( const auto &triplet : angleTriplets ) {
      const auto &assocVertexIdx = associatedVertex[triplet.row()];
      const auto &assocIvIdx = interiorVerticesIndices[assocVertexIdx];
      if ( assocIvIdx >= 0 )
        Dest.emplace_back( assocIvIdx, triplet.col(), -triplet.value());
    }

  }

  int getTargetDimension() const override {
    return _numInteriorVertices;
  }

  int getNNZ() const override {
    return _numNonzeros;
  }
};

/**
 * \brief Hessian of angle-defect Gauss curvature
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 * \see GaussCurvatureOp
 */
template<typename ConfiguratorType>
class GaussCurvatureHessian
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
  InteriorAngleHessian<ConfiguratorType> _gammaHess;

  std::vector<std::vector<int>> associatedAngles;
  std::vector<int> associatedVertex;

  std::vector<int> interiorVerticesIndices;

  const int _numEdges, _numFaces;

public:
  int _numInteriorVertices;
  std::vector<int> interiorVertices;

  /**
   * \brief Construct operator
   * \param topology class containing information on the topology (i.e. connectivity) of the mesh.
   *
   * \todo Move computation of the needed topological information to a separate topology class
   */
  explicit GaussCurvatureHessian( const MeshTopologySaver &topology ) : _topology( topology ), _gammaHess( _topology ),
                                                                        _numFaces( _topology.getNumFaces()),
                                                                        _numEdges( _topology.getNumEdges()) {
    _numInteriorVertices = 0;
    interiorVerticesIndices.resize( _topology.getNumVertices(), -1 );
    for ( int vertexIdx = 0; vertexIdx < _topology.getNumVertices(); vertexIdx++ ) {

      bool interiorVertex = true;
      for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); edgeIdx++ ) {
        if ( _topology.getAdjacentNodeOfEdge( edgeIdx, 0 ) == vertexIdx ||
             _topology.getAdjacentNodeOfEdge( edgeIdx, 1 ) == vertexIdx ) {
          int f0 = _topology.getOppositeNodeOfEdge( edgeIdx, 0 );
          int f1 = _topology.getOppositeNodeOfEdge( edgeIdx, 1 );

          if ( f0 == -1 || f1 == -1 ) {
            interiorVertex = false;
            break;
          }
        }
      }

      if ( interiorVertex ) {
        interiorVerticesIndices[vertexIdx] = _numInteriorVertices;
        _numInteriorVertices++;
        interiorVertices.push_back( vertexIdx );
      }
    }

    associatedAngles.resize( _topology.getNumVertices());
    associatedVertex.resize( 3 * _topology.getNumFaces());


    for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); faceIdx++ ) {
      std::array<int, 3> v = { _topology.getNodeOfTriangle( faceIdx, 0 ),
                               _topology.getNodeOfTriangle( faceIdx, 1 ),
                               _topology.getNodeOfTriangle( faceIdx, 2 ) };

      for ( int i : { 0, 1, 2 } ) {
        associatedAngles[v[i]].push_back( 3 * faceIdx + i );
        associatedVertex[3 * faceIdx + i] = v[i];
      }
    }
  }

  /**
   * \brief Evaluate Hessian
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[out] Dest Hessian of vertex-wise Gauss curvature as vector of sparse matrices
   * \param[in] hessSize (optional) modified size of the components of the Hessian
   * \param[in] hessOffset (optional) symmetric offset for the entires of the individual components of the Hessian
   */
  void apply( const VectorType &Arg, std::vector<MatrixType> &Dest, int hessSize = -1, int hessOffset = 0 ) const {
    if ( Arg.size() < _topology.getNumEdges())
      throw std::length_error( "GaussCurvatureHessian::apply(): Arg too small!" );

    Dest.resize( _numInteriorVertices );

    std::vector<TripletListType> vertexTripletLists;
    vertexTripletLists.resize( _numInteriorVertices );

    setTriplets( Arg, vertexTripletLists, hessOffset );

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int ivIdx = 0; ivIdx < _numInteriorVertices; ivIdx++ ) {
      if ( hessSize == -1 )
        Dest[ivIdx].resize( 2 * _numEdges, 2 * _numEdges );
      else
        Dest[ivIdx].resize( hessSize, hessSize );

      Dest[ivIdx].reserve( vertexTripletLists[ivIdx].size());
      Dest[ivIdx].setFromTriplets( vertexTripletLists[ivIdx].begin(), vertexTripletLists[ivIdx].end());
    }

  }

  /**
   * \brief Evaluate Hessian
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[out] Dest Hessian of discrete integrability map as GenericTensor
   */
  void apply( const VectorType &Arg, TensorType &Dest ) const {
    apply( Arg, Dest, -1, 0 );
  }

  /**
   * \brief Evaluate Hessian
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[out] Dest Hessian of discrete integrability map as GenericTensor
   * \param[in] hessSize (optional) modified size of the components of the Hessian
   * \param[in] hessOffset (optional) symmetric offset for the entires of the individual components of the Hessian
   */
  void apply( const VectorType &Arg, TensorType &Dest, int hessSize, int hessOffset ) const {
    if ( Arg.size() < _topology.getNumEdges())
      throw std::length_error( "GaussCurvatureHessian::apply(): Arg too small!" );

    Dest.resize( _numInteriorVertices );

    std::vector<TripletListType> vertexTripletLists;
    vertexTripletLists.resize( _numInteriorVertices );

    setTriplets( Arg, vertexTripletLists, hessOffset );

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int ivIdx = 0; ivIdx < _numInteriorVertices; ivIdx++ ) {
      if ( hessSize == -1 )
        Dest[ivIdx].resize( 2 * _numEdges, 2 * _numEdges );
      else
        Dest[ivIdx].resize( hessSize, hessSize );

      Dest[ivIdx].reserve( vertexTripletLists[ivIdx].size());
      Dest[ivIdx].setFromTriplets( vertexTripletLists[ivIdx].begin(), vertexTripletLists[ivIdx].end());
    }

  }

  /**
   * \brief Evaluate sum of Hessian components to triplet list
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[out] Dest Entries of sum of Hessian components as triplet list
   * \param[in] hessSize (optional) modified size of the components of the Hessian
   */
  void pushTriplets( const VectorType &Arg, TripletListType &Dest, int hessOffset = 0 ) const {
    if ( Arg.size() < _topology.getNumEdges())
      throw std::length_error( "GaussCurvatureHessian::pushTriplets(): Arg too small!" );

    VectorType Lambda( _numInteriorVertices );
    Lambda.setConstant( 1. );

    pushTriplets( Arg, Dest, Lambda, hessOffset );
  }

  /**
   * \brief Evaluate weighted sum of Hessian components to triplet list
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[out] Dest Entries of weighted sum of Hessian components as triplet list
   * \param[in] Lambda Weights in the sum
   * \param[in] hessSize (optional) modified size of the components of the Hessian
   */
  void pushTriplets( const VectorType &Arg, TripletListType &Dest, const VectorType &Lambda,
                     int hessOffset = 0 ) const {
    if ( Arg.size() < _topology.getNumEdges())
      throw std::length_error( "GaussCurvatureHessian::pushTriplets(): Arg too small!" );
    if ( Lambda.size() != _numInteriorVertices )
      throw std::length_error( "GaussCurvatureHessian::pushTriplets(): Lambda too small!" );

    std::vector<TripletListType> vertexTripletLists;
    vertexTripletLists.resize( _numInteriorVertices );

    setTriplets( Arg, vertexTripletLists, hessOffset );

    int totalNumTriplets = 0;
    for ( const auto &vTL : vertexTripletLists ) {
      totalNumTriplets += vTL.size();
    }

    Dest.reserve( totalNumTriplets );

    for ( int i = 0; i < _numInteriorVertices; i++ ) {
      for ( const auto &trip : vertexTripletLists[i] )
        Dest.emplace_back( trip.row(), trip.col(), Lambda[i] * trip.value());
    }
  }

  /**
   * \brief Evaluate Hessian to triplet lists
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[out] Dest Entries of Hessian as vector of triplet lists
   * \param[in] hessSize (optional) modified size of the components of the Hessian
   */
  void setTriplets( const VectorType &Arg, std::vector<TripletListType> &Dest, int hessOffset = 0 ) const {
    if ( Arg.size() < _topology.getNumEdges())
      throw std::length_error( "GaussCurvatureHessian::apply(): Arg too small!" );

    Dest.resize( _numInteriorVertices );
    for ( auto &vTL : Dest )
      vTL.clear();

    std::vector<TripletListType> angleTriplets;
    _gammaHess.setTriplets( Arg, angleTriplets, hessOffset );

    for ( int ivIdx = 0; ivIdx < _numInteriorVertices; ivIdx++ ) {
      int vertexIdx = interiorVertices[ivIdx];

      int numTriplets = 0;
      for ( const auto &angleIdx : associatedAngles[vertexIdx] )
        numTriplets += angleTriplets[angleIdx].size();

      Dest[ivIdx].reserve( numTriplets );

      for ( const auto &angleIdx : associatedAngles[vertexIdx] ) {
        for ( const auto &triplet : angleTriplets[angleIdx] ) {
          Dest[ivIdx].emplace_back( triplet.row(), triplet.col(), -triplet.value());
        }
      }
    }

  }

  int getTargetDimension() const override {
    return _numInteriorVertices;
  }
};

#endif //REDUCEDBASIS_GAUSSCURVATURE_H
