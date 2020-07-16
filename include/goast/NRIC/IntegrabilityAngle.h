// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Operators yielding the integrability violation of given lengths and angles as angle per vertex
 * \author Sassen
 *
 * Based on
 * Wang, Y., Liu, B., & Tong, Y. (2012). Linear surface reconstruction from discrete fundamental forms on triangle
 * meshes. Computer Graphics Forum, 31(8), 2277–2287.
 *
 * \note This currently relies directly on the dense tensors from Eigen
 * \todo Remove Eigen's tensors and cleanup!
 */

#ifndef NRIC_TRACEINTEGRABILITY_H
#define NRIC_TRACEINTEGRABILITY_H

#include <goast/Core/Auxiliary.h>
#include <goast/Core/LocalMeshGeometry.h>
#include <goast/Core/Topology.h>

#include "TriangleInequality.h"
#include "TransitionRotations.h"

/**
 * \brief Operator evaluating the integrability violation on edge lengths and dihedral angles as angle per vertex
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 *
 * This operator computes the discrete integrability map, i.e. the equation of the discrete integrability conditions,
 * for each vertex, by chaining the transition rotations around the adjacent edges and returns the angle of the
 * resulting roation for each interior vertex.
 *
 * \todo Move computation of the needed topological information to a topology class?
 */
template<typename ConfiguratorType>
class AngleIntegrabilityOp
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::RotationType RotationType;

  const MeshTopologySaver &_topology;
  const int _numVertices;
  const int _numEdges;
  std::vector<bool> boundaryVertex;
  std::vector<std::vector<std::pair<int, bool>>> orientedEdgeRings;
  std::vector<std::vector<std::tuple<int, int, int>>> assocEdgeLengths;


public:
  int _numInteriorVertices;
  std::vector<int> interiorVertices;

  /**
   * \brief Construct operator
   * \param topology class containing information on the topology (i.e. connectivity) of the mesh.
   */
  explicit AngleIntegrabilityOp( const MeshTopologySaver &topology ) : _topology( topology ),
                                                                       _numVertices( _topology.getNumVertices()),
                                                                       _numEdges( _topology.getNumEdges()) {
    boundaryVertex.resize( _numVertices, false );
    orientedEdgeRings.resize( _numVertices, {} );
    assocEdgeLengths.resize( _numVertices, {} );

    interiorVertices.reserve( _numVertices );
    _numInteriorVertices = 0;


    // Compute which vertices are boundary vertices, and for those which are not compute the oriented 1-ring of adjacent
    // edges
    for ( int vertexIdx = 0; vertexIdx < _topology.getNumVertices(); vertexIdx++ ) {
      // Determine edges connected to vertex
      std::vector<int> edges;
      for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); edgeIdx++ ) {
        if ( _topology.getAdjacentNodeOfEdge( edgeIdx, 0 ) == vertexIdx ||
             _topology.getAdjacentNodeOfEdge( edgeIdx, 1 ) == vertexIdx ) {
          edges.push_back( edgeIdx );
        }
      }

      if ( edges.empty()) {
        boundaryVertex[vertexIdx] = true;
        continue;
      }

      // order the edges and determine if it is a boundary edge
      std::vector<std::pair<int, bool>> orderedEdges;
      orderedEdges.emplace_back( edges[0], false );

      while ( orderedEdges.size() < edges.size()) { // while we have not sorted in all edges.
        int activeEdge = orderedEdges.back().first; // last sorted in edge
        bool activeEdgeFlipped = orderedEdges.back().second;
        int f1 = _topology.getAdjacentTriangleOfEdge( activeEdge, activeEdgeFlipped );
        int f2 = _topology.getAdjacentTriangleOfEdge( activeEdge, !activeEdgeFlipped );

        // if edge has only one neighboring triangle it is a boundary edge and thus we have boundary vertex
        if ( f1 == -1 || f2 == -1 ) {
          boundaryVertex[vertexIdx] = true;
          break;
        }

        // possible next edges
        std::array<int, 6> candidateEdges = {
                _topology.getEdgeOfTriangle( f2, 0 ),
                _topology.getEdgeOfTriangle( f2, 1 ),
                _topology.getEdgeOfTriangle( f2, 2 )
        };

        // find a candidate edges which we haven't considered so far and is adjacent to our vertex
        for ( auto e : candidateEdges ) {
          if ( std::find_if( orderedEdges.begin(), orderedEdges.end(),
                             [e]( const std::pair<int, bool> &s ) { return s.first == e; } ) == orderedEdges.end() &&
               std::find( edges.begin(), edges.end(), e ) != edges.end()) {
            // found the next edge, now look determine the orientation of the edge, i.e. if we have to invert the
            // transition rotation
            int e0a = _topology.getAdjacentTriangleOfEdge( e, 0 );
            int e0b = _topology.getAdjacentTriangleOfEdge( e, 1 );
            orderedEdges.emplace_back( e, e0a != f2 );
            break;
          }
        }
      }

      if ( boundaryVertex[vertexIdx] )
        continue;

      interiorVertices.push_back( vertexIdx );
      _numInteriorVertices++;

      orientedEdgeRings[vertexIdx] = orderedEdges;

      // Determine for each edge around the vertex which edge lengths determine the corresponding interior angle at v
      std::vector<std::tuple<int, int, int>> lengthIdx;
      for ( auto oe : orderedEdges ) {
        int edgeIdx = oe.first;
        bool orientation = oe.second;

        int f1 = _topology.getAdjacentTriangleOfEdge( edgeIdx, orientation );
        int f2 = _topology.getAdjacentTriangleOfEdge( edgeIdx, !orientation );

        std::array<int, 3> v2 = { _topology.getNodeOfTriangle( f2, 0 ),
                                  _topology.getNodeOfTriangle( f2, 1 ),
                                  _topology.getNodeOfTriangle( f2, 2 ) };

        std::array<int, 3> e2 = { _topology.getEdgeOfTriangle( f2, 0 ),
                                  _topology.getEdgeOfTriangle( f2, 1 ),
                                  _topology.getEdgeOfTriangle( f2, 2 ) };

        // Determine which local indices the edge has and the local indices of the angle between the two
        long vIdx_f2 = std::distance( v2.begin(), std::find( v2.begin(), v2.end(), vertexIdx ));
        long eIdx_f2 = std::distance( e2.begin(), std::find( e2.begin(), e2.end(), edgeIdx ));
        int li_0 = ( vIdx_f2 == 0 ) ? e2[1] : ( vIdx_f2 == 1 ) ? e2[2] : e2[0];
        int li_1 = ( vIdx_f2 == 0 ) ? e2[2] : ( vIdx_f2 == 1 ) ? e2[0] : e2[1];

        lengthIdx.emplace_back( li_1, li_0, e2[vIdx_f2] );
      }

      assocEdgeLengths[vertexIdx] = lengthIdx;
    }

  }


  /**
   * \brief Evaluate operator
   * \param Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param Dest Trace of chained rotation matrices for each vertex as vector
   */
  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    if ( Arg.size() != 2 * _numEdges )
      throw std::length_error( "AngleIntegrabilityOp::apply(): Arg too small!" );

    Dest.resize( _numInteriorVertices );
    Dest.setZero();

#if !defined(NDEBUG)
    if ( Arg.segment( 0, _numEdges ).minCoeff() < 0 ) {
      std::cerr << "WARNING(AngleIntegrabilityOp::apply):  negative edge lengths!" << std::endl;
    }

    TriangleInequalityTripleOp<ConfiguratorType> triqOp( _topology );
    VectorType triangleInequalityValues;
    triqOp.apply( Arg, triangleInequalityValues );

    if ( triangleInequalityValues.minCoeff() <= 0 ) {
      std::cerr << "WARNING(AngleIntegrabilityOp::apply): triangle inequality violated!" << std::endl;
    }
#endif

    LocalTransitionRotationOp<ConfiguratorType> rotationGen;

    for ( int ivIdx = 0; ivIdx < _numInteriorVertices; ivIdx++ ) {
      int vertexIdx = interiorVertices[ivIdx];
      RotationType accRotation = RotationType::Identity();
      for ( int i = 0; i < orientedEdgeRings[vertexIdx].size(); i++ ) {
        int edgeIdx = orientedEdgeRings[vertexIdx][i].first;
        accRotation = accRotation * rotationGen( Arg[_numEdges + edgeIdx],
                                                 Arg[std::get<0>( assocEdgeLengths[vertexIdx][i] )],
                                                 Arg[std::get<1>( assocEdgeLengths[vertexIdx][i] )],
                                                 Arg[std::get<2>( assocEdgeLengths[vertexIdx][i] )] );
      }

      Dest[ivIdx] = accRotation( 0, 0 ) + accRotation( 1, 1 ) + accRotation( 2, 2 );
    }

  }

  int getTargetDimension() const {
    return _numInteriorVertices;
  }
};

/**
 * \brief Jacobian of the simplified discrete integrability map on edge lengths and dihedral angles
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 * \see TraceIntegrabilityOp
 */
template<typename ConfiguratorType>
class AngleIntegrabilityGradient
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::RotationType RotationType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  template<int R>
  using TensorType = Eigen::Tensor<RealType, R>;

  const MeshTopologySaver &_topology;
  const int _numVertices;
  const int _numEdges;
  int _numNonzero;
  int _numTriplets;
  std::vector<bool> boundaryVertex;
  std::vector<std::vector<std::pair<int, bool>>> orientedEdgeRings;
  std::vector<std::vector<int>> interiorAngleIdx;
  std::vector<std::vector<std::tuple<int, int, int>>> assocEdgeLengths;

  int _numInteriorVertices;
  std::vector<int> interiorVertices;

public:
  /**
   * \brief Construct operator
   * \param topology class containing information on the topology (i.e. connectivity) of the mesh.
   */
  explicit AngleIntegrabilityGradient( const MeshTopologySaver &topology ) :
          _topology( topology ),
          _numVertices( _topology.getNumVertices()),
          _numEdges( _topology.getNumEdges()) {
    boundaryVertex.resize( _numVertices, false );
    orientedEdgeRings.resize( _numVertices, {} );
    interiorAngleIdx.resize( _numVertices, {} );
    assocEdgeLengths.resize( _numVertices, {} );

    _numNonzero = 0;
    _numTriplets = 0;

    interiorVertices.reserve( _numVertices );
    _numInteriorVertices = 0;

    // Compute which vertices are boundary vertices, and for those which are not compute the oriented 1-ring of adjacent
    // edges
    for ( int vertexIdx = 0; vertexIdx < _topology.getNumVertices(); vertexIdx++ ) {
      // Determine edges connected to vertex
      std::vector<int> edges;
      for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); edgeIdx++ ) {
        if ( _topology.getAdjacentNodeOfEdge( edgeIdx, 0 ) == vertexIdx ||
             _topology.getAdjacentNodeOfEdge( edgeIdx, 1 ) == vertexIdx ) {
          edges.push_back( edgeIdx );
        }
      }

      if ( edges.empty()) {
        boundaryVertex[vertexIdx] = true;
        continue;
      }

      // order the edges and determine if it is a boundary edge
      std::vector<std::pair<int, bool>> orderedEdges;
      orderedEdges.emplace_back( edges[0], false );

      while ( orderedEdges.size() < edges.size()) { // while we have not sorted in all edges.
        int activeEdge = orderedEdges.back().first; // last sorted in edge
        bool activeEdgeFlipped = orderedEdges.back().second;
        int f1 = _topology.getAdjacentTriangleOfEdge( activeEdge, activeEdgeFlipped );
        int f2 = _topology.getAdjacentTriangleOfEdge( activeEdge, !activeEdgeFlipped );

        // if edge has only one neighboring triangle it is a boundary edge and thus we have boundary vertex
        if ( f1 == -1 || f2 == -1 ) {
          boundaryVertex[vertexIdx] = true;
          break;
        }

        // possible next edges
        std::array<int, 6> candidateEdges = {
                _topology.getEdgeOfTriangle( f2, 0 ),
                _topology.getEdgeOfTriangle( f2, 1 ),
                _topology.getEdgeOfTriangle( f2, 2 )
        };

        // find a candidate edges which we haven considered so far and is adjacent to our vertex
        for ( auto e : candidateEdges ) {
          if ( std::find_if( orderedEdges.begin(), orderedEdges.end(),
                             [e]( const std::pair<int, bool> &s ) { return s.first == e; } ) == orderedEdges.end() &&
               std::find( edges.begin(), edges.end(), e ) != edges.end()) {
            // found the next edge, now look determine the orientation of the edge, i.e. if we have to invert the
            // transition rotation
            int e0a = _topology.getAdjacentTriangleOfEdge( e, 0 );
            int e0b = _topology.getAdjacentTriangleOfEdge( e, 1 );
            orderedEdges.emplace_back( e, e0a != f2 );
            break;
          }
        }
      }

      if ( boundaryVertex[vertexIdx] )
        continue;

      interiorVertices.push_back( vertexIdx );
      _numInteriorVertices++;

      orientedEdgeRings[vertexIdx] = orderedEdges;

      _numNonzero += 3 * orderedEdges.size();
      _numTriplets += 4 * orderedEdges.size();

      std::vector<int> angleIdx;
      std::vector<std::tuple<int, int, int>> lengthIdx;
      for ( auto oe : orderedEdges ) {
        int edgeIdx = oe.first;
        bool orientation = oe.second;

        int f1 = _topology.getAdjacentTriangleOfEdge( edgeIdx, orientation );
        int f2 = _topology.getAdjacentTriangleOfEdge( edgeIdx, !orientation );

        std::array<int, 3> v2 = { _topology.getNodeOfTriangle( f2, 0 ),
                                  _topology.getNodeOfTriangle( f2, 1 ),
                                  _topology.getNodeOfTriangle( f2, 2 ) };

        std::array<int, 3> e2 = { _topology.getEdgeOfTriangle( f2, 0 ),
                                  _topology.getEdgeOfTriangle( f2, 1 ),
                                  _topology.getEdgeOfTriangle( f2, 2 ) };

        // Determine which local indices the edge has and the local indices of the angle between the two
        long vIdx_f2 = std::distance( v2.begin(), std::find( v2.begin(), v2.end(), vertexIdx ));
        long eIdx_f2 = std::distance( e2.begin(), std::find( e2.begin(), e2.end(), edgeIdx ));
        int li_0 = ( vIdx_f2 == 0 ) ? e2[1] : ( vIdx_f2 == 1 ) ? e2[2] : e2[0];
        int li_1 = ( vIdx_f2 == 0 ) ? e2[2] : ( vIdx_f2 == 1 ) ? e2[0] : e2[1];

        angleIdx.emplace_back( 3 * f2 + vIdx_f2 );
        lengthIdx.emplace_back( li_1, li_0, e2[vIdx_f2] );

      }

      interiorAngleIdx[vertexIdx] = angleIdx;
      assocEdgeLengths[vertexIdx] = lengthIdx;
    }

  }


  /**
   * \brief Evaluate derivative
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[out] Dest Jacobian of discrete integrability map
   */
  void apply( const VectorType &Arg, MatrixType &Dest ) const override {
    if ( Arg.size() != 2 * _topology.getNumEdges())
      throw std::length_error( "AngleIntegrabilityGradient::apply(): Arg too small!" );

    Dest.resize( _numInteriorVertices, 2 * _numEdges );
    Dest.reserve( _numNonzero );
    Dest.setZero();

    TripletListType tripletList;
    pushTriplets( Arg, tripletList );

    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend());
  }

  /**
   * \brief Evaluate derivative to triplet list
   * \tparam transposed determines whether to return transposed Jacobian
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[out] Dest Jacobian as triplet list
   */
  template<bool transposed = false>
  void pushTriplets( const VectorType &Arg,
                     TripletListType &Dest,
                     RealType factor = 1.,
                     int rowOffset = 0,
                     int colOffset = 0 ) const {
    if ( Arg.size() != 2 * _topology.getNumEdges())
      throw std::length_error( "AngleIntegrabilityGradient::pushTriplets(): Arg too small!" );

    LocalTransitionRotationOp<ConfiguratorType> rotationGen;
    LocalTransitionRotationGradient<ConfiguratorType> rotationGradGen;
    Dest.reserve( _numTriplets );

#ifdef GOAST_WITH_OPENMP
#pragma omp parallel for
#endif
    for ( int ivIdx = 0; ivIdx < _numInteriorVertices; ivIdx++ ) {
      int vertexIdx = interiorVertices[ivIdx];
      TripletListType localTripletList;
      localTripletList.reserve( 4 * orientedEdgeRings[vertexIdx].size());

      // Computing local gradients of the transition rotations
      std::vector<TensorType<3>> transitionGradients( orientedEdgeRings[vertexIdx].size());
      for ( int i = 0; i < orientedEdgeRings[vertexIdx].size(); i++ ) {
        int edgeIdx = orientedEdgeRings[vertexIdx][i].first;
        auto ael = assocEdgeLengths[vertexIdx][i];
        transitionGradients[i] = rotationGradGen( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )],
                                                  Arg[std::get<1>( ael )], Arg[std::get<2>( ael )] );
      }

      // Summands of the product rule
      for ( int i = 0; i < orientedEdgeRings[vertexIdx].size(); i++ ) {
        // Multiplication of rotations left of the partial derivative
        RotationType preRotation = RotationType::Identity();
        for ( int j = 0; j < i; j++ ) {
          int edgeIdx = orientedEdgeRings[vertexIdx][j].first;
          auto ael = assocEdgeLengths[vertexIdx][j];

          preRotation = preRotation * rotationGen( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )],
                                                   Arg[std::get<1>( ael )], Arg[std::get<2>( ael )] );
        }

        // Multiplication of rotations right of the partial derivative
        RotationType postRotation = RotationType::Identity();
        for ( int j = i + 1; j < orientedEdgeRings[vertexIdx].size(); j++ ) {
          int edgeIdx = orientedEdgeRings[vertexIdx][j].first;
          auto ael = assocEdgeLengths[vertexIdx][j];

          postRotation = postRotation * rotationGen( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )],
                                                     Arg[std::get<1>( ael )], Arg[std::get<2>( ael )] );
        }

        int edgeIdx = orientedEdgeRings[vertexIdx][i].first;
        auto ael = assocEdgeLengths[vertexIdx][i];

        // Consider matrices as 2-tensors
        Eigen::TensorMap<TensorType<2>> preTensor( preRotation.data(), 3, 3 );
        Eigen::TensorMap<TensorType<2>> postTensor( postRotation.data(), 3, 3 );

        // The contractions we perform are in fact just ordinary matrix multiplications
        Eigen::array<Eigen::IndexPair<int>, 1> product_dims = { Eigen::IndexPair<int>( 1, 0 ) };

        // Compute the partial derivatives by multiplying preRotation * transitionGradient * postRotation
        TensorType<2> partialDeriv( 3, 3 );
        RealType trace;

        partialDeriv = preTensor
                .contract( transitionGradients[i].chip( 0, 2 ), product_dims )
                .contract( postTensor, product_dims );
        trace = partialDeriv( 0, 0 ) + partialDeriv( 1, 1 ) + partialDeriv( 2, 2 );
        if ( transposed )
          localTripletList.emplace_back( _numEdges + edgeIdx + rowOffset, ivIdx + colOffset, factor * trace );
        else
          localTripletList.emplace_back( ivIdx + rowOffset, _numEdges + edgeIdx + colOffset, factor * trace );

        partialDeriv = preTensor
                .contract( transitionGradients[i].chip( 1, 2 ), product_dims )
                .contract( postTensor, product_dims );
        trace = partialDeriv( 0, 0 ) + partialDeriv( 1, 1 ) + partialDeriv( 2, 2 );
        if ( transposed )
          localTripletList.emplace_back( std::get<0>( ael ) + rowOffset, ivIdx + colOffset, factor * trace );
        else
          localTripletList.emplace_back( ivIdx + rowOffset, std::get<0>( ael ) + colOffset, factor * trace );

        partialDeriv = preTensor
                .contract( transitionGradients[i].chip( 2, 2 ), product_dims )
                .contract( postTensor, product_dims );
        trace = partialDeriv( 0, 0 ) + partialDeriv( 1, 1 ) + partialDeriv( 2, 2 );

        if ( transposed )
          localTripletList.emplace_back( std::get<1>( ael ) + rowOffset, ivIdx + colOffset, factor * trace );
        else
          localTripletList.emplace_back( ivIdx + rowOffset, std::get<1>( ael ) + colOffset, factor * trace );

        partialDeriv = preTensor
                .contract( transitionGradients[i].chip( 3, 2 ), product_dims )
                .contract( postTensor, product_dims );
        trace = partialDeriv( 0, 0 ) + partialDeriv( 1, 1 ) + partialDeriv( 2, 2 );

        if ( transposed )
          localTripletList.emplace_back( std::get<2>( ael ) + rowOffset, ivIdx + colOffset, factor * trace );
        else
          localTripletList.emplace_back( ivIdx + rowOffset, std::get<2>( ael ) + colOffset, factor * trace );

      }

#ifdef GOAST_WITH_OPENMP
#pragma omp critical
#endif
      Dest.insert( Dest.end(), localTripletList.begin(), localTripletList.end());
    }
  }

  int getTargetDimension() const override {
    return _numInteriorVertices;
  }

  int getNNZ() const override {
    return _numNonzero;
  }
};

/**
 * \brief Hessian of the simplified discrete integrability map on edge lengths and dihedral angles
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 * \see TraceIntegrabilityOp
 */
template<typename ConfiguratorType>
class AngleIntegrabilityHessian
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::TensorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
  typedef typename ConfiguratorType::RotationType RotationType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::TensorType SparseTensorType;
  typedef std::vector<TripletType> TripletListType;

  template<int R>
  using TensorType = Eigen::Tensor<RealType, R>;

  const MeshTopologySaver &_topology;
  const int _numVertices;
  const int _numEdges;
  int _numNonzero;
  int _numTriplets;
  std::vector<bool> boundaryVertex;
  std::vector<std::vector<std::pair<int, bool>>> orientedEdgeRings;
  std::vector<std::vector<int>> interiorAngleIdx;
  std::vector<std::vector<std::array<int, 3>>> assocEdgeLengths;
  std::vector<std::vector<std::array<int, 4>>> assocDOF;

  int _numInteriorVertices;
  std::vector<int> interiorVertices;

public:
  /**
   * \brief Construct operator
   * \param topology class containing information on the topology (i.e. connectivity) of the mesh.
   */
  explicit AngleIntegrabilityHessian( const MeshTopologySaver &topology ) :
          _topology( topology ),
          _numVertices( _topology.getNumVertices()),
          _numEdges( _topology.getNumEdges()) {
    boundaryVertex.resize( _numVertices, false );
    orientedEdgeRings.resize( _numVertices, {} );
    interiorAngleIdx.resize( _numVertices, {} );
    assocEdgeLengths.resize( _numVertices, {} );
    assocDOF.resize( _numVertices, {} );

    _numNonzero = 0;
    _numTriplets = 0;

    interiorVertices.reserve( _numVertices );
    _numInteriorVertices = 0;

    // Compute which vertices are boundary vertices, and for those which are not compute the oriented 1-ring of adjacent
    // edges
    for ( int vertexIdx = 0; vertexIdx < _topology.getNumVertices(); vertexIdx++ ) {
      // Determine edges connected to vertex
      std::vector<int> edges;
      for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); edgeIdx++ ) {
        if ( _topology.getAdjacentNodeOfEdge( edgeIdx, 0 ) == vertexIdx ||
             _topology.getAdjacentNodeOfEdge( edgeIdx, 1 ) == vertexIdx ) {
          edges.push_back( edgeIdx );
        }
      }

      if ( edges.empty()) {
        boundaryVertex[vertexIdx] = true;
        continue;
      }

      // order the edges and determine if it is a boundary edge
      std::vector<std::pair<int, bool>> orderedEdges;
      orderedEdges.emplace_back( edges[0], false );

      while ( orderedEdges.size() < edges.size()) { // while we have not sorted in all edges.
        int activeEdge = orderedEdges.back().first; // last sorted in edge
        bool activeEdgeFlipped = orderedEdges.back().second;
        int f1 = _topology.getAdjacentTriangleOfEdge( activeEdge, activeEdgeFlipped );
        int f2 = _topology.getAdjacentTriangleOfEdge( activeEdge, !activeEdgeFlipped );

        // if edge has only one neighboring triangle it is a boundary edge and thus we have boundary vertex
        if ( f1 == -1 || f2 == -1 ) {
          boundaryVertex[vertexIdx] = true;
          break;
        }

        // possible next edges
        std::array<int, 6> candidateEdges = {
                _topology.getEdgeOfTriangle( f2, 0 ),
                _topology.getEdgeOfTriangle( f2, 1 ),
                _topology.getEdgeOfTriangle( f2, 2 )
        };

        // find a candidate edges which we haven considered so far and is adjacent to our vertex
        for ( auto e : candidateEdges ) {
          if ( std::find_if( orderedEdges.begin(), orderedEdges.end(),
                             [e]( const std::pair<int, bool> &s ) { return s.first == e; } ) == orderedEdges.end() &&
               std::find( edges.begin(), edges.end(), e ) != edges.end()) {
            // found the next edge, now look determine the orientation of the edge, i.e. if we have to invert the
            // transition rotation
            int e0a = _topology.getAdjacentTriangleOfEdge( e, 0 );
            int e0b = _topology.getAdjacentTriangleOfEdge( e, 1 );
            orderedEdges.emplace_back( e, e0a != f2 );
            break;
          }
        }
      }

      if ( boundaryVertex[vertexIdx] )
        continue;

      interiorVertices.push_back( vertexIdx );
      _numInteriorVertices++;

      orientedEdgeRings[vertexIdx] = orderedEdges;

      _numNonzero += 9 * orderedEdges.size() * orderedEdges.size();
      _numTriplets += 2 * 16 * orderedEdges.size() * ( orderedEdges.size() + 1 ) / 2;

      std::vector<int> angleIdx;
      std::vector<std::array<int, 3>> lengthIdx;
      std::vector<std::array<int, 4>> dof;
      for ( auto oe : orderedEdges ) {
        int edgeIdx = oe.first;
        bool orientation = oe.second;

        int f1 = _topology.getAdjacentTriangleOfEdge( edgeIdx, orientation );
        int f2 = _topology.getAdjacentTriangleOfEdge( edgeIdx, !orientation );

        std::array<int, 3> v2 = { _topology.getNodeOfTriangle( f2, 0 ),
                                  _topology.getNodeOfTriangle( f2, 1 ),
                                  _topology.getNodeOfTriangle( f2, 2 ) };

        std::array<int, 3> e2 = { _topology.getEdgeOfTriangle( f2, 0 ),
                                  _topology.getEdgeOfTriangle( f2, 1 ),
                                  _topology.getEdgeOfTriangle( f2, 2 ) };

        // Determine which local indices the edge has and the local indices of the angle between the two
        long vIdx_f2 = std::distance( v2.begin(), std::find( v2.begin(), v2.end(), vertexIdx ));
        long eIdx_f2 = std::distance( e2.begin(), std::find( e2.begin(), e2.end(), edgeIdx ));
        int li_0 = ( vIdx_f2 == 0 ) ? e2[1] : ( vIdx_f2 == 1 ) ? e2[2] : e2[0];
        int li_1 = ( vIdx_f2 == 0 ) ? e2[2] : ( vIdx_f2 == 1 ) ? e2[0] : e2[1];

        angleIdx.emplace_back( 3 * f2 + vIdx_f2 );
        lengthIdx.push_back( { li_1, li_0, e2[vIdx_f2] } );
        dof.push_back( { _numEdges + edgeIdx, li_1, li_0, e2[vIdx_f2] } );

      }

      interiorAngleIdx[vertexIdx] = angleIdx;
      assocEdgeLengths[vertexIdx] = lengthIdx;
      assocDOF[vertexIdx] = dof;
    }

  }


  /**
   * \brief Evaluate Hessian
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[out] Dest Hessian of discrete integrability map as GenericTensor
   * \param[in] hessSize (optional) modified size of the components of the Hessian
   * \param[in] hessOffset (optional) symmetric offset for the entires of the individual components of the Hessian
   */
  void apply( const VectorType &Arg, SparseTensorType &Dest, int hessSize = -1, int hessOffset = 0 ) const {
    apply( Arg, static_cast<std::vector<MatrixType>>(Dest), -1, 0 );
  }

  /**
   * \brief Evaluate Hessian
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[out] Dest Hessian of discrete integrability map as GenericTensor
   */
  void apply( const VectorType &Arg, SparseTensorType &Dest ) const {
    apply( Arg, static_cast<std::vector<MatrixType>>(Dest), -1, 0 );
  }

  /**
   * \brief Evaluate Hessian
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[out] Dest Hessian of discrete integrability map
   */
  void apply( const VectorType &Arg, std::vector<MatrixType> &Dest,
              int hessSize = -1, int hessOffset = 0 ) const {
    if ( Arg.size() != 2 * _topology.getNumEdges())
      throw std::length_error( "AngleIntegrabilityHessian::apply(): Arg too small!" );

    Dest.resize( _numInteriorVertices );

    LocalTransitionRotationOp<ConfiguratorType> rotationGen;
    LocalTransitionRotationGradient<ConfiguratorType> rotationGradGen;
    LocalTransitionRotationHessian<ConfiguratorType> rotationGradHess;


#ifdef GOAST_WITH_OPENMP
#pragma omp parallel for
#endif
    for ( int ivIdx = 0; ivIdx < _numInteriorVertices; ivIdx++ ) {
      int vertexIdx = interiorVertices[ivIdx];
      MatrixType vertexHessian;
      if ( hessSize == -1 )
        vertexHessian.resize( 2 * _numEdges, 2 * _numEdges );
      else
        vertexHessian.resize( hessSize, hessSize );
      vertexHessian.reserve( 9 * orientedEdgeRings[vertexIdx].size() * orientedEdgeRings[vertexIdx].size());
      vertexHessian.setZero();

      TripletListType vertexTripletList;
      vertexTripletList.reserve(
              2 * 16 * orientedEdgeRings[vertexIdx].size() * ( orientedEdgeRings[vertexIdx].size() + 1 ) / 2 );

      // Computing local gradients of the transition rotations
      std::vector<TensorType<3>> transitionGradients( orientedEdgeRings[vertexIdx].size());
      std::vector<TensorType<4>> transitionHessians( orientedEdgeRings[vertexIdx].size());
      for ( int i = 0; i < orientedEdgeRings[vertexIdx].size(); i++ ) {
        int edgeIdx = orientedEdgeRings[vertexIdx][i].first;
        auto ael = assocEdgeLengths[vertexIdx][i];
        transitionGradients[i] = rotationGradGen( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )],
                                                  Arg[std::get<1>( ael )], Arg[std::get<2>( ael )] );
        transitionHessians[i] = rotationGradHess( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )],
                                                  Arg[std::get<1>( ael )], Arg[std::get<2>( ael )] );
      }

      // Summands of the first product rule
      for ( int i = 0; i < orientedEdgeRings[vertexIdx].size(); i++ ) {
        // Summands of the second product rule
        for ( int j = 0; j <= i; j++ ) {
//          std::cout << ivIdx << " - " << i << "," << j << std::endl;
          int edgeIdx_i = orientedEdgeRings[vertexIdx][i].first;
          auto ael_i = assocEdgeLengths[vertexIdx][i];
          auto dof_i = assocDOF[vertexIdx][i];

          int edgeIdx_j = orientedEdgeRings[vertexIdx][j].first;
          auto ael_j = assocEdgeLengths[vertexIdx][j];
          auto dof_j = assocDOF[vertexIdx][j];

          // Multiplication of rotations left of the partial derivatives
          RotationType preRotation = RotationType::Identity();
          for ( int k = 0; k < j; k++ ) {
            int edgeIdx = orientedEdgeRings[vertexIdx][k].first;
            auto ael = assocEdgeLengths[vertexIdx][k];

            preRotation = preRotation * rotationGen( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )],
                                                     Arg[std::get<1>( ael )], Arg[std::get<2>( ael )] );
          }

          // Multiplication of rotations right of the partial derivatives
          RotationType postRotation = RotationType::Identity();
          for ( int k = i + 1; k < orientedEdgeRings[vertexIdx].size(); k++ ) {
            int edgeIdx = orientedEdgeRings[vertexIdx][k].first;
            auto ael = assocEdgeLengths[vertexIdx][k];

            postRotation = postRotation * rotationGen( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )],
                                                       Arg[std::get<1>( ael )], Arg[std::get<2>( ael )] );
          }

          if ( i == j ) {  // Second partial derivatives
            // Consider matrices as 2-tensors
            Eigen::TensorMap<TensorType<2>> preTensor( preRotation.data(), 3, 3 );
            Eigen::TensorMap<TensorType<2>> postTensor( postRotation.data(), 3, 3 );

            // The contractions we perform are in fact just ordinary matrix multiplications
            Eigen::array<Eigen::IndexPair<int>, 1> product_dims = { Eigen::IndexPair<int>( 1, 0 ) };

            // Compute the partial derivatives by multiplying preRotation * localDeriv * midRotation * localDeriv * postRotation
            TensorType<2> partialDeriv( 3, 3 );


            for ( int l : { 0, 1, 2, 3 } ) {
              for ( int h : { 0, 1, 2, 3 } ) {
                partialDeriv = preTensor
                        .contract( transitionHessians[i].chip( l, 3 ).chip( h, 2 ), product_dims )
                        .contract( postTensor, product_dims );

                RealType trace = partialDeriv( 0, 0 ) + partialDeriv( 1, 1 ) + partialDeriv( 2, 2 );

                vertexTripletList.push_back( TripletType( dof_i[h] + hessOffset, dof_i[l] + hessOffset, trace ));
              }
            }
          }
          else { // Mixed first partial derivatives
            // Multiplication of rotations between the partial derivatives
            RotationType midRotation = RotationType::Identity();
            for ( int k = j + 1; k < i; k++ ) {
              int edgeIdx = orientedEdgeRings[vertexIdx][k].first;
              auto ael = assocEdgeLengths[vertexIdx][k];

              midRotation = midRotation * rotationGen( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )],
                                                       Arg[std::get<1>( ael )], Arg[std::get<2>( ael )] );
            }

            // Consider matrices as 2-tensors
            Eigen::TensorMap<TensorType<2>> preTensor( preRotation.data(), 3, 3 );
            Eigen::TensorMap<TensorType<2>> midTensor( midRotation.data(), 3, 3 );
            Eigen::TensorMap<TensorType<2>> postTensor( postRotation.data(), 3, 3 );

            // The contractions we perform are in fact just ordinary matrix multiplications
            Eigen::array<Eigen::IndexPair<int>, 1> product_dims = { Eigen::IndexPair<int>( 1, 0 ) };

            // Compute the partial derivatives by multiplying preRotation * localDeriv * midRotation * localDeriv * postRotation
            TensorType<2> partialDeriv( 3, 3 );

            for ( int l : { 0, 1, 2, 3 } ) {
              for ( int h : { 0, 1, 2, 3 } ) {
                partialDeriv = preTensor
                        .contract( transitionGradients[j].chip( l, 2 ), product_dims )
                        .contract( midTensor, product_dims )
                        .contract( transitionGradients[i].chip( h, 2 ), product_dims )
                        .contract( postTensor, product_dims );

                RealType trace = partialDeriv( 0, 0 ) + partialDeriv( 1, 1 ) + partialDeriv( 2, 2 );

                vertexTripletList.push_back( TripletType( dof_j[l] + hessOffset, dof_i[h] + hessOffset, trace ));
                vertexTripletList.push_back( TripletType( dof_i[h] + hessOffset, dof_j[l] + hessOffset, trace ));
              }
            }
          }

        }
      }

      vertexHessian.setFromTriplets( vertexTripletList.cbegin(), vertexTripletList.cend());
//      std::cout << "Diagonal (" << ivIdx << "): ";
//      for (int kk = 0; kk < 2*_numEdges;kk++)
//        std::cout << vertexHessian.coeffRef(kk,kk) << " ";
//      std::cout << std::endl;
      Dest[ivIdx] = vertexHessian;
    }

  }

  /**
   * \brief Evaluate sum of Hessian components to triplet list
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[in] weights Vector containing the weights for the different components of the Hessian
   * \param[out] Dest Hessian as triplet list
   *
   * Computes \f$\sum_i \lambda_i Hess f_i \in R^{2m\times2m}\f$
   */
  void pushTriplets( const VectorType &Arg, const VectorType &weights, TripletListType &Dest ) const {
    if ( Arg.size() != 2 * _topology.getNumEdges())
      throw std::length_error( "AngleIntegrabilityHessian::apply(): Arg too small!" );

    if ( weights.size() != _numInteriorVertices )
      throw std::length_error( "AngleIntegrabilityHessian::apply(): Weights too small!" );


    Dest.reserve( _numTriplets );

    LocalTransitionRotationOp<ConfiguratorType> rotationGen;
    LocalTransitionRotationGradient<ConfiguratorType> rotationGradGen;
    LocalTransitionRotationHessian<ConfiguratorType> rotationGradHess;


#ifdef GOAST_WITH_OPENMP
#pragma omp parallel for
#endif
    for ( int ivIdx = 0; ivIdx < _numInteriorVertices; ivIdx++ ) {
      int vertexIdx = interiorVertices[ivIdx];
      TripletListType vertexTripletList;
      vertexTripletList.reserve(
              2 * 16 * orientedEdgeRings[vertexIdx].size() * ( orientedEdgeRings[vertexIdx].size() + 1 ) / 2 );

      // Computing local gradients of the transition rotations
      std::vector<TensorType<3>> transitionGradients( orientedEdgeRings[vertexIdx].size());
      std::vector<TensorType<4>> transitionHessians( orientedEdgeRings[vertexIdx].size());
      for ( int i = 0; i < orientedEdgeRings[vertexIdx].size(); i++ ) {
        int edgeIdx = orientedEdgeRings[vertexIdx][i].first;
        auto ael = assocEdgeLengths[vertexIdx][i];
        transitionGradients[i] = rotationGradGen( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )],
                                                  Arg[std::get<1>( ael )], Arg[std::get<2>( ael )] );
        transitionHessians[i] = rotationGradHess( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )],
                                                  Arg[std::get<1>( ael )], Arg[std::get<2>( ael )] );
      }

      // Summands of the first product rule
      for ( int i = 0; i < orientedEdgeRings[vertexIdx].size(); i++ ) {
        // Summands of the second product rule
        for ( int j = 0; j <= i; j++ ) {
          int edgeIdx_i = orientedEdgeRings[vertexIdx][i].first;
          auto ael_i = assocEdgeLengths[vertexIdx][i];
          auto dof_i = assocDOF[vertexIdx][i];

          int edgeIdx_j = orientedEdgeRings[vertexIdx][j].first;
          auto ael_j = assocEdgeLengths[vertexIdx][j];
          auto dof_j = assocDOF[vertexIdx][j];

          // Multiplication of rotations left of the partial derivatives
          RotationType preRotation = RotationType::Identity();
          for ( int k = 0; k < j; k++ ) {
            int edgeIdx = orientedEdgeRings[vertexIdx][k].first;
            auto ael = assocEdgeLengths[vertexIdx][k];

            preRotation = preRotation * rotationGen( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )],
                                                     Arg[std::get<1>( ael )], Arg[std::get<2>( ael )] );
          }

          // Multiplication of rotations right of the partial derivatives
          RotationType postRotation = RotationType::Identity();
          for ( int k = i + 1; k < orientedEdgeRings[vertexIdx].size(); k++ ) {
            int edgeIdx = orientedEdgeRings[vertexIdx][k].first;
            auto ael = assocEdgeLengths[vertexIdx][k];

            postRotation = postRotation * rotationGen( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )],
                                                       Arg[std::get<1>( ael )], Arg[std::get<2>( ael )] );
          }

          if ( i == j ) {  // Second partial derivatives
            // Consider matrices as 2-tensors
            Eigen::TensorMap<TensorType<2>> preTensor( preRotation.data(), 3, 3 );
            Eigen::TensorMap<TensorType<2>> postTensor( postRotation.data(), 3, 3 );

            // The contractions we perform are in fact just ordinary matrix multiplications
            Eigen::array<Eigen::IndexPair<int>, 1> product_dims = { Eigen::IndexPair<int>( 1, 0 ) };

            // Compute the partial derivatives by multiplying preRotation * localDeriv * midRotation * localDeriv * postRotation
            TensorType<2> partialDeriv( 3, 3 );

            for ( int l : { 0, 1, 2, 3 } ) {
              for ( int h : { 0, 1, 2, 3 } ) {
                partialDeriv = preTensor
                        .contract( transitionHessians[i].chip( l, 3 ).chip( h, 2 ), product_dims )
                        .contract( postTensor, product_dims );

                RealType trace = partialDeriv( 0, 0 ) + partialDeriv( 1, 1 ) + partialDeriv( 2, 2 );

                vertexTripletList.push_back( TripletType( dof_i[h], dof_i[l], trace * weights[ivIdx] ));
              }
            }
          }
          else { // Mixed first partial derivatives
            // Multiplication of rotations between the partial derivatives
            RotationType midRotation = RotationType::Identity();
            for ( int k = j + 1; k < i; k++ ) {
              int edgeIdx = orientedEdgeRings[vertexIdx][k].first;
              auto ael = assocEdgeLengths[vertexIdx][k];

              midRotation = midRotation * rotationGen( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )],
                                                       Arg[std::get<1>( ael )], Arg[std::get<2>( ael )] );
            }

            // Consider matrices as 2-tensors
            Eigen::TensorMap<TensorType<2>> preTensor( preRotation.data(), 3, 3 );
            Eigen::TensorMap<TensorType<2>> midTensor( midRotation.data(), 3, 3 );
            Eigen::TensorMap<TensorType<2>> postTensor( postRotation.data(), 3, 3 );

            // The contractions we perform are in fact just ordinary matrix multiplications
            Eigen::array<Eigen::IndexPair<int>, 1> product_dims = { Eigen::IndexPair<int>( 1, 0 ) };

            // Compute the partial derivatives by multiplying preRotation * localDeriv * midRotation * localDeriv * postRotation
            TensorType<2> partialDeriv( 3, 3 );

            for ( int l : { 0, 1, 2, 3 } ) {
              for ( int h : { 0, 1, 2, 3 } ) {
                partialDeriv = preTensor
                        .contract( transitionGradients[j].chip( l, 2 ), product_dims )
                        .contract( midTensor, product_dims )
                        .contract( transitionGradients[i].chip( h, 2 ), product_dims )
                        .contract( postTensor, product_dims );

                RealType trace = partialDeriv( 0, 0 ) + partialDeriv( 1, 1 ) + partialDeriv( 2, 2 );

                vertexTripletList.push_back( TripletType( dof_j[l], dof_i[h], trace * weights[ivIdx] ));
                vertexTripletList.push_back( TripletType( dof_i[h], dof_j[l], trace * weights[ivIdx] ));
              }
            }
          }

        }
      }
#ifdef GOAST_WITH_OPENMP
#pragma omp critical
#endif
      Dest.insert( Dest.end(), vertexTripletList.begin(), vertexTripletList.end());
    }

  }

  int getTargetDimension() const {
    return _numInteriorVertices;
  }

  int getNNZ() const {
    return _numNonzero;
  }

  void applyAdd( const VectorType &Arg, std::vector<MatrixType> &Dest ) const {
    throw std::logic_error( "applyAdd: Not implemented!" );
  }
};

#endif //NRIC_TRACEINTEGRABILITY_H
