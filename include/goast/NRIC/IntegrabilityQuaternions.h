// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Operators yielding the integrability violation of given lengths and angles as (vector part of) a quaternion per vertex
 * \author Sassen
 *
 * Based on
 * Wang, Y., Liu, B., & Tong, Y. (2012). Linear surface reconstruction from discrete fundamental forms on triangle
 * meshes. Computer Graphics Forum, 31(8), 2277–2287.
 * and
 * Sassen, J., Heeren, B., Hildebrandt, K., & Rumpf, M. (2020). Geometric optimization using nonlinear
 * rotation-invariant coordinates. Computer Aided Geometric Design, 77, 101829.
 *
 * \note This currently relies directly on the Axis-Angle and Quaternion implementations from Eigen
 */

#ifndef NRIC_INTEGRABILITYQUATERNIONS_H
#define NRIC_INTEGRABILITYQUATERNIONS_H

#include <goast/Core/Auxiliary.h>
#include <goast/Core/LocalMeshGeometry.h>
#include <goast/Core/Topology.h>

#include "TriangleInequality.h"
#include "TransitionQuaternions.h"

/**
 * \brief Operator evaluating discrete integrability conditions using quaternions on edge lengths and dihedral angles
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 *
 * This operator computes the modified discrete integrability map, i.e. the equation of the discrete integrability
 * conditions, for each vertex, by chaining the transition quaternion around the adjacent edges.
 *
 * \todo Move computation of the needed topological information to a topology class
 */
template<typename ConfiguratorType>
class QuaternionIntegrabilityOp
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;
  typedef typename ConfiguratorType::QuaternionType QuaternionType;
  typedef QuaternionType RotationType;

  const MeshTopologySaver &_topology;
  const int _numVertices;
  const int _numEdges;
  std::vector<bool> boundaryVertex;
  std::vector<std::vector<std::pair<int, bool>>> orientedEdgeRings;
  std::vector<std::vector<int>> interiorAngleIdx;
  std::vector<std::vector<std::tuple<int, int, int>>> assocEdgeLengths;

  const bool _onlyVec;

public:
  int _numInteriorVertices;
  std::vector<int> interiorVertices;

  /**
   * \brief Construct operator
   * \param topology class containing information on the topology (i.e. connectivity) of the mesh.
   * \param onlyVec whether to ouput only the vector part of quaternions, default is true
   */
  explicit QuaternionIntegrabilityOp( const MeshTopologySaver &topology, bool onlyVec = true )
          : _topology( topology ),
            _numVertices( _topology.getNumVertices()),
            _numEdges( _topology.getNumEdges()),
            _onlyVec( onlyVec ) {
    boundaryVertex.resize( _numVertices, false );
    orientedEdgeRings.resize( _numVertices, {} );
    interiorAngleIdx.resize( _numVertices, {} );
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
        int li_0 = (vIdx_f2 == 0) ? e2[1] : (vIdx_f2 == 1) ? e2[2] : e2[0];
        int li_1 = (vIdx_f2 == 0) ? e2[2] : (vIdx_f2 == 1) ? e2[0] : e2[1];

        angleIdx.emplace_back( 3 * f2 + vIdx_f2 );
        lengthIdx.emplace_back( li_1, li_0, e2[vIdx_f2] );

      }

      interiorAngleIdx[vertexIdx] = angleIdx;
      assocEdgeLengths[vertexIdx] = lengthIdx;
    }

//    interiorVertices.pop_back();
//    _numInteriorVertices--;

  }


  /**
   * \brief Evaluate operator
   * \param Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param Dest Euler angles of chained rotation matrices for each vertex as vector
   */
  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    if ( Arg.size() != 2 * _numEdges )
      throw std::length_error( "QuaternionIntegrabilityOp::apply(): Arg too small!" );

    Dest.resize( _onlyVec ? 3 * _numInteriorVertices : 4 * _numInteriorVertices );
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

    LocalTransitionQuaternionOp<ConfiguratorType> rotationGen;
#ifdef GOAST_WITH_OPENMP
#pragma omp parallel for
#endif
    for ( int ivIdx = 0; ivIdx < _numInteriorVertices; ivIdx++ ) {
      const int &vertexIdx = interiorVertices[ivIdx];

      RotationType accRotation = RotationType::Identity();
      std::vector<RotationType, Eigen::aligned_allocator<RotationType>> transitionRotations( orientedEdgeRings[vertexIdx].size() );
      for ( int i = 0; i < orientedEdgeRings[vertexIdx].size(); i++ ) {
        int edgeIdx = orientedEdgeRings[vertexIdx][i].first;

        rotationGen.apply( Arg[_numEdges + edgeIdx],
                           Arg[std::get<0>( assocEdgeLengths[vertexIdx][i] )],
                           Arg[std::get<1>( assocEdgeLengths[vertexIdx][i] )],
                           Arg[std::get<2>( assocEdgeLengths[vertexIdx][i] )], transitionRotations[i] );

        accRotation = accRotation * transitionRotations[i];
      }

      if ( _onlyVec ) {
        Dest[3 * ivIdx] = accRotation.x();
        Dest[3 * ivIdx + 1] = accRotation.y();
        Dest[3 * ivIdx + 2] = accRotation.z();
      }
      else {
        Dest[4 * ivIdx] = accRotation.w();
        Dest[4 * ivIdx + 1] = accRotation.x();
        Dest[4 * ivIdx + 2] = accRotation.y();
        Dest[4 * ivIdx + 3] = accRotation.z();
      }
    }

    for ( int j = 0; j < Dest.size(); j++ )
      if ( std::isnan( Dest[j] ))
        Dest[j] = std::numeric_limits<RealType>::infinity();

  }

  int getTargetDimension() const {
    return _onlyVec ? 3 * _numInteriorVertices : 4 * _numInteriorVertices;
  }
};

/**
 * \brief Jacobian of the modified discrete integrability map on edge lengths and dihedral angles
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 * \see QuaternionIntegrabilityOp
 */
template<typename ConfiguratorType>
class QuaternionIntegrabilityGradient
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;
  typedef typename ConfiguratorType::QuaternionType QuaternionType;
  typedef QuaternionType RotationType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;


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

  const bool _onlyVec;
  const int _vertexDOF;

public:
  /**
   * \brief Construct operator
   * \param topology class containing information on the topology (i.e. connectivity) of the mesh.
   */
  explicit QuaternionIntegrabilityGradient( const MeshTopologySaver &topology, bool onlyVec = true  ) :
          _topology( topology ),
          _numVertices( _topology.getNumVertices()),
          _numEdges( _topology.getNumEdges()), _onlyVec ( onlyVec ), _vertexDOF( onlyVec ? 3 : 4) {
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

      _numNonzero += _vertexDOF * 3 * orderedEdges.size();
      _numTriplets += _vertexDOF * 4 * orderedEdges.size();

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
        int li_0 = (vIdx_f2 == 0) ? e2[1] : (vIdx_f2 == 1) ? e2[2] : e2[0];
        int li_1 = (vIdx_f2 == 0) ? e2[2] : (vIdx_f2 == 1) ? e2[0] : e2[1];

        angleIdx.emplace_back( 3 * f2 + vIdx_f2 );
        lengthIdx.emplace_back( li_1, li_0, e2[vIdx_f2] );

      }

      interiorAngleIdx[vertexIdx] = angleIdx;
      assocEdgeLengths[vertexIdx] = lengthIdx;
    }

//    std::cerr << "Remove vertex " << interiorVertices.back() << std::endl;
//    interiorVertices.pop_back();
//    _numInteriorVertices--;

//    std::cerr << "Remove vertex " << interiorVertices[0] << std::endl;
//    interiorVertices.erase(interiorVertices.begin());
//    _numInteriorVertices--;

  }


  /**
   * \brief Evaluate derivative
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[out] Dest Jacobian of discrete integrability map
   */
  void apply( const VectorType &Arg, MatrixType &Dest ) const override {
    if ( Arg.size() != 2 * _topology.getNumEdges())
      throw std::length_error( "QuaternionIntegrabilityGradient::apply(): Arg too small!" );

    Dest.resize( _vertexDOF * _numInteriorVertices, 2 * _numEdges );
    Dest.reserve( 3 * _numNonzero );
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
   * \param[in] factor scaling of the Jacobian
   * \param[in] rowOffset Row offset of the triplets
   * \param[in] colOffset Column offset of the triplets
   */
  //  template<bool transposed = false>
  void pushTriplets( const VectorType &Arg, TripletListType &Dest,
                     RealType factor = 1., int rowOffset = 0, int colOffset = 0 ) const {
    const bool transposed = false;
    if ( Arg.size() != 2 * _topology.getNumEdges())
      throw std::length_error( "QuaternionIntegrabilityGradient::pushTriplets(): Arg too small!" );

    LocalTransitionQuaternionOp<ConfiguratorType> rotationGen;
    LocalTransitionQuaternionGradient<ConfiguratorType> rotationGradGen;
    Dest.reserve( _numTriplets );

#ifdef GOAST_WITH_OPENMP
#pragma omp parallel for
#endif
    for ( int ivIdx = 0; ivIdx < _numInteriorVertices; ivIdx++ ) {
      int vertexIdx = interiorVertices[ivIdx];
      TripletListType localTripletList;
      localTripletList.reserve( _vertexDOF * 4 * orientedEdgeRings[vertexIdx].size());

      // Computing local gradients of the transition rotations
      RotationType accRotation = RotationType::Identity();
      std::vector<RotationType, Eigen::aligned_allocator<RotationType>> transitionRotations( orientedEdgeRings[vertexIdx].size());
      std::vector<std::array<RotationType, 4>, Eigen::aligned_allocator<std::array<RotationType, 4>>> transitionGradients( orientedEdgeRings[vertexIdx].size());
      for ( int i = 0; i < orientedEdgeRings[vertexIdx].size(); i++ ) {
        int edgeIdx = orientedEdgeRings[vertexIdx][i].first;
        auto ael = assocEdgeLengths[vertexIdx][i];
        rotationGen.apply( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )], Arg[std::get<1>( ael )],
                           Arg[std::get<2>( ael )], transitionRotations[i] );

        rotationGradGen.apply( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )], Arg[std::get<1>( ael )],
                               Arg[std::get<2>( ael )], transitionGradients[i] );

        accRotation = accRotation * transitionRotations[i];
      }

      RotationType preRotation = RotationType::Identity();
      RotationType postRotation = accRotation;

      // Summands of the product rule
      for ( int i = 0; i < orientedEdgeRings[vertexIdx].size(); i++ ) {
        // Multiplication of rotations right of the partial derivative
        postRotation = transitionRotations[i].conjugate() * postRotation;

        int edgeIdx = orientedEdgeRings[vertexIdx][i].first;
        auto ael = assocEdgeLengths[vertexIdx][i];

        // Compute the partial derivatives by multiplying preRotation * transitionGradient * postRotation
        RotationType partialDeriv;

        // Partial derivative w.r.t. dihedral angle
        partialDeriv = preRotation * transitionGradients[i][0] * postRotation;
        pushToList<transposed>( ivIdx, _numEdges + edgeIdx, partialDeriv,
                                localTripletList, factor, rowOffset, colOffset );

        // Partial derivative w.r.t. first edge length
        partialDeriv = preRotation * transitionGradients[i][1] * postRotation;
        pushToList<transposed>( ivIdx, std::get<0>( ael ), partialDeriv,
                                localTripletList, factor, rowOffset, colOffset );

        // Partial derivative w.r.t. second edge length
        partialDeriv = preRotation * transitionGradients[i][2] * postRotation;
        pushToList<transposed>( ivIdx, std::get<1>( ael ), partialDeriv,
                                localTripletList, factor, rowOffset, colOffset );

        // Partial derivative w.r.t. third edge length
        partialDeriv = preRotation * transitionGradients[i][3] * postRotation;
        pushToList<transposed>( ivIdx, std::get<2>( ael ), partialDeriv,
                                localTripletList, factor, rowOffset, colOffset );

        // Multiplication of rotations left of the partial derivative
        preRotation = preRotation * transitionRotations[i];
      }

#ifdef GOAST_WITH_OPENMP
#pragma omp critical
#endif
      Dest.insert( Dest.end(), localTripletList.begin(), localTripletList.end());
    }
  }

  int getTargetDimension() const override {
    return _vertexDOF * _numInteriorVertices;
  }

  int getNNZ() const override {
    return _numNonzero;
  }

private:
  template<bool transposed = false>
  void pushToList( const int ivIdx, const int dofIdx, const QuaternionType &value, TripletListType &tripletList,
                   const RealType factor = 1., const int rowOffset = 0, const int colOffset = 0 ) const {
    if ( transposed ) {
      if ( _onlyVec ) {
        tripletList.emplace_back( dofIdx + rowOffset, 3 * ivIdx + colOffset, factor * value.x());
        tripletList.emplace_back( dofIdx + rowOffset, 3 * ivIdx + 1 + colOffset, factor * value.y());
        tripletList.emplace_back( dofIdx + rowOffset, 3 * ivIdx + 2 + colOffset, factor * value.z());
      }
      else {
        tripletList.emplace_back( dofIdx + rowOffset, 4 * ivIdx + colOffset, factor * value.w());
        tripletList.emplace_back( dofIdx + rowOffset, 4 * ivIdx + 1 + colOffset, factor * value.x());
        tripletList.emplace_back( dofIdx + rowOffset, 4 * ivIdx + 2 + colOffset, factor * value.y());
        tripletList.emplace_back( dofIdx + rowOffset, 4 * ivIdx + 3 + colOffset, factor * value.z());
      }
    }
    else {
      if ( _onlyVec ) {
        tripletList.emplace_back( 3 * ivIdx + rowOffset, dofIdx + colOffset, factor * value.x());
        tripletList.emplace_back( 3 * ivIdx + 1 + rowOffset, dofIdx + colOffset, factor * value.y());
        tripletList.emplace_back( 3 * ivIdx + 2 + rowOffset, dofIdx + colOffset, factor * value.z());
      }
      else {
        tripletList.emplace_back( 4 * ivIdx + rowOffset, dofIdx + colOffset, factor * value.w());
        tripletList.emplace_back( 4 * ivIdx + 1 + rowOffset, dofIdx + colOffset, factor * value.x());
        tripletList.emplace_back( 4 * ivIdx + 2 + rowOffset, dofIdx + colOffset, factor * value.y());
        tripletList.emplace_back( 4 * ivIdx + 3 + rowOffset, dofIdx + colOffset, factor * value.z());
      }
    }
  }
};

/**
 * \brief Hessian of the modified discrete integrability map on edge lengths and dihedral angles
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 * \see QuaternionIntegrabilityOp
 */
template<typename ConfiguratorType>
class QuaternionIntegrabilityHessian
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::TensorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;
  typedef typename ConfiguratorType::QuaternionType QuaternionType;
  typedef QuaternionType RotationType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;
  typedef typename ConfiguratorType::TensorType TensorType;

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
  std::vector<int> numAssocDOF;
  std::vector<int> tripletsStartIdx;
  std::vector<std::vector<int>> combinedAssocDOF;
  std::vector<std::vector<int>> caDOF_localToGlobal;

  std::vector<int> interiorVertices;

  const bool _onlyVec;
  const int _vertexDOF;

public:
  int _numInteriorVertices;
  mutable std::map<std::string, RealType> timings;

  /**
   * \brief Construct operator
   * \param topology class containing information on the topology (i.e. connectivity) of the mesh.
   */
  explicit QuaternionIntegrabilityHessian( const MeshTopologySaver &topology, bool onlyVec = true ) :
          _topology( topology ),
          _numVertices( _topology.getNumVertices()),
          _numEdges( _topology.getNumEdges()), _onlyVec ( onlyVec ), _vertexDOF( onlyVec ? 3 : 4) {
    boundaryVertex.resize( _numVertices, false );
    orientedEdgeRings.resize( _numVertices, {} );
    interiorAngleIdx.resize( _numVertices, {} );
    assocEdgeLengths.resize( _numVertices, {} );
    assocDOF.resize( _numVertices, {} );
    numAssocDOF.resize( _numVertices, 0 );
    combinedAssocDOF.resize( _numVertices, {} );
    caDOF_localToGlobal.resize( _numVertices, {} );
    tripletsStartIdx.resize( _numVertices, 0 );

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

      //! \todo Figure out the right values!
      _numNonzero += 9 * orderedEdges.size() * orderedEdges.size();

      std::vector<int> angleIdx;
      std::vector<std::array<int, 3>> lengthIdx;
      std::vector<std::array<int, 4>> dof;
      std::vector<int> combinedDof;
      std::vector<int> cdof_localToGlobal;
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
        int li_0 = (vIdx_f2 == 0) ? e2[1] : (vIdx_f2 == 1) ? e2[2] : e2[0];
        int li_1 = (vIdx_f2 == 0) ? e2[2] : (vIdx_f2 == 1) ? e2[0] : e2[1];

        angleIdx.emplace_back( 3 * f2 + vIdx_f2 );
        lengthIdx.push_back( { li_1, li_0, e2[vIdx_f2] } );
        dof.push_back( { _numEdges + edgeIdx, li_1, li_0, e2[vIdx_f2] } );

        for ( auto &d : { _numEdges + edgeIdx, li_1, li_0, e2[vIdx_f2] } ) {
          if ( std::find( combinedDof.begin(), combinedDof.end(), d ) == combinedDof.end()) {
            combinedDof.push_back( d );
          }
        }
      }


      interiorAngleIdx[vertexIdx] = angleIdx;
      assocEdgeLengths[vertexIdx] = lengthIdx;
      assocDOF[vertexIdx] = dof;
      combinedAssocDOF[vertexIdx] = combinedDof;

      numAssocDOF[vertexIdx] = combinedDof.size();

      tripletsStartIdx[vertexIdx] = _numTriplets;
      _numTriplets += _vertexDOF * numAssocDOF[vertexIdx] * numAssocDOF[vertexIdx];
    }

//    interiorVertices.pop_back();
//    _numInteriorVertices--;

  }


  /**
   * \brief Evaluate Hessian
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[out] Dest Hessian of discrete integrability map as GenericTensor
   * \param[in] hessSize (optional) modified size of the components of the Hessian
   * \param[in] hessOffset (optional) symmetric offset for the entires of the individual components of the Hessian
   */
  void apply( const VectorType &Arg, TensorType &Dest, int hessSize = -1, int hessOffset = 0 ) const {
    if ( Arg.size() != 2 * _topology.getNumEdges())
      throw std::length_error( "QuaternionIntegrabilityHessian::apply(): Arg too small!" );

    if ( hessSize == -1 )
      Dest.resize( _vertexDOF * _numInteriorVertices, 2 * _numEdges, 2 * _numEdges );
    else
      Dest.resize( _vertexDOF * _numInteriorVertices, hessSize, hessSize );

    std::vector<TripletListType> vertexTripletLists;
    vertexTripletLists.resize( _vertexDOF * _numInteriorVertices );

    setTriplets( Arg, vertexTripletLists, hessOffset );

    std::vector<int> subIndices;
    if (_onlyVec) {
      subIndices = {0,1,2};
    }
    else {
      subIndices = {0,1,2,3};
    }



#ifdef GOAST_WITH_OPENMP
#pragma omp parallel for
#endif
    for ( int ivIdx = 0; ivIdx < _numInteriorVertices; ivIdx++ ) {

      for ( int t : subIndices ) {
        Dest[_vertexDOF * ivIdx + t].reserve( vertexTripletLists[_vertexDOF * ivIdx + t].size());
        Dest[_vertexDOF * ivIdx + t].setFromTriplets( vertexTripletLists[_vertexDOF * ivIdx + t].begin(),
                                                      vertexTripletLists[_vertexDOF * ivIdx + t].end());
      }
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
   * \param[out] Dest Hessian of discrete integrability map as vector of sparse matrices
   * \param[in] hessSize (optional) modified size of the components of the Hessian
   * \param[in] hessOffset (optional) symmetric offset for the entires of the individual components of the Hessian
   */
  void apply( const VectorType &Arg, std::vector<SparseMatrixType> &Dest,
              int hessSize = -1, int hessOffset = 0 ) const {
    if ( Arg.size() != 2 * _topology.getNumEdges())
      throw std::length_error( "QuaternionIntegrabilityHessian::apply(): Arg too small!" );

    Dest.resize( _vertexDOF * _numInteriorVertices );

    std::vector<TripletListType> vertexTripletLists;
    vertexTripletLists.resize( _vertexDOF * _numInteriorVertices );

    setTriplets( Arg, vertexTripletLists, hessOffset );


    std::vector<int> subIndices;
    if (_onlyVec) {
      subIndices = {0,1,2};
    }
    else {
      subIndices = {0,1,2,3};
    }

#ifdef GOAST_WITH_OPENMP
#pragma omp parallel for
#endif
    for ( int ivIdx = 0; ivIdx < _numInteriorVertices; ivIdx++ ) {
      for ( int t : subIndices ) {
        if ( hessSize == -1 )
          Dest[_vertexDOF * ivIdx + t].resize( 2 * _numEdges, 2 * _numEdges );
        else
          Dest[_vertexDOF * ivIdx + t].resize( hessSize, hessSize );

        Dest[_vertexDOF * ivIdx + t].reserve( vertexTripletLists[_vertexDOF * ivIdx + t].size());

        Dest[_vertexDOF * ivIdx + t].setFromTriplets( vertexTripletLists[_vertexDOF * ivIdx + t].begin(),
                                                      vertexTripletLists[_vertexDOF * ivIdx + t].end());
      }
    }

  }


  /**
   * \brief Evaluate sum of Hessian components to triplet lists
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[out] Dest Entries of sum of Hessian components as triplet list
   * \param[in] hessSize (optional) modified size of the components of the Hessian
   */
  void pushTriplets( const VectorType &Arg, TripletListType &Dest, int hessOffset ) const {
    if ( Arg.size() != 2 * _topology.getNumEdges())
      throw std::length_error( "QuaternionIntegrabilityHessian::apply(): Arg too small!" );

    std::vector<TripletListType> vertexTripletLists;
    vertexTripletLists.resize( _vertexDOF * _numInteriorVertices );

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


  /**
   * \brief Evaluate weighted sum of Hessian components to triplet list
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[out] Dest Entries of weighted sum of Hessian components as triplet list
   * \param[in] Lambda Weights in the sum
   * \param[in] hessSize (optional) modified size of the components of the Hessian
   */
  void pushTriplets( const VectorType &Arg, TripletListType &Dest, const VectorType &Lambda,
                     int hessOffset ) const {
    if ( Arg.size() != 2 * _topology.getNumEdges())
      throw std::length_error( "QuaternionIntegrabilityHessian::pushTriplets(): Arg too small! " );
    if ( Lambda.size() != _vertexDOF * _numInteriorVertices )
      throw std::length_error( "QuaternionIntegrabilityHessian::pushTriplets(): Lambda too small!" );

    int start = Dest.size();

    Dest.resize( Dest.size() + _numTriplets, TripletType (0,0,0.));

    std::vector<int> subIndices;
    if (_onlyVec) {
      subIndices = {0,1,2};
    }
    else {
      subIndices = {0,1,2,3};
    }

    timings["computePartialSecondDerivatives"] = 0.;
    timings["setValues"] = 0.;
    timings["Local Derivatives"] = 0;
    timings["Product Rule"] = 0;


#ifdef GOAST_WITH_OPENMP
#pragma omp parallel for default(shared)
#endif
    for ( int ivIdx = 0; ivIdx < _numInteriorVertices; ivIdx++ ) {
      int vertexIdx = interiorVertices[ivIdx];

      const auto &ADOF = combinedAssocDOF[vertexIdx];
      const auto &numADOF = numAssocDOF[vertexIdx];

      // Create space for  second derivatives
      std::vector<std::vector<RotationType, Eigen::aligned_allocator<RotationType>>> partialSecondDerivatives;
      partialSecondDerivatives.resize( numAssocDOF[vertexIdx] );
      for ( auto &pSD : partialSecondDerivatives ) {
        pSD.resize( numAssocDOF[vertexIdx], QuaternionType( 0, 0, 0, 0 ));
      }

      auto t_start = std::chrono::high_resolution_clock::now();
      computePartialSecondDerivatives( Arg, vertexIdx, partialSecondDerivatives );
      auto t_end = std::chrono::high_resolution_clock::now();

      t_start = std::chrono::high_resolution_clock::now();
      for ( int i = 0; i < numADOF; i++ ) {
        for ( int j = 0; j < numADOF; j++ ) {
          const QuaternionType &value = partialSecondDerivatives[i][j];
          const int startIdx = start + tripletsStartIdx[vertexIdx];
          const int localIdx = numADOF*i+j;
          if ( _onlyVec ) {
            Dest[startIdx + 3 * localIdx] = TripletType ( ADOF[i] + hessOffset, ADOF[j] + hessOffset, Lambda[3*ivIdx] * value.x());
            Dest[startIdx + 3 * localIdx + 1] = TripletType( ADOF[i] + hessOffset, ADOF[j] + hessOffset, Lambda[3*ivIdx + 1] * value.y());
            Dest[startIdx + 3 * localIdx + 2] = TripletType( ADOF[i] + hessOffset, ADOF[j] + hessOffset, Lambda[3*ivIdx + 2] * value.z());
          }
          else {
            Dest[startIdx + 4 * localIdx] = TripletType( ADOF[i] + hessOffset, ADOF[j] + hessOffset, Lambda[4*ivIdx] * value.w());
            Dest[startIdx + 4 * localIdx + 1] = TripletType( ADOF[i] + hessOffset, ADOF[j] + hessOffset, Lambda[4*ivIdx + 1] * value.x());
            Dest[startIdx + 4 * localIdx + 2] = TripletType( ADOF[i] + hessOffset, ADOF[j] + hessOffset, Lambda[4*ivIdx + 2] * value.y());
            Dest[startIdx + 4 * localIdx + 3] = TripletType( ADOF[i] + hessOffset, ADOF[j] + hessOffset, Lambda[4*ivIdx + 3] * value.z());
          }

        }
      }
      t_end = std::chrono::high_resolution_clock::now();

    }
  }


  /**
   * \brief Evaluate sum of Hessian components to triplet list
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[out] Dest Entries of sum of Hessian components as triplet list
   */
  void pushTriplets( const VectorType &Arg, TripletListType &Dest ) const override {
    pushTriplets( Arg, Dest, 0 );
  }

  /**
   * \brief Evaluate weighted sum of Hessian components to triplet list
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[out] Dest Entries of weighted sum of Hessian components as triplet lists
   * \param[in] Lambda Weights in the sum
   */
  void pushTriplets( const VectorType &Arg, TripletListType &Dest, const VectorType &Lambda ) const override {
    pushTriplets( Arg, Dest, Lambda, 0 );
  }

  /**
   * \brief Evaluate Hessian to triplet lists
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[out] Dest Entries of Hessian as vector of triplet lists
   * \param[in] hessSize (optional) modified size of the components of the Hessian
   */
  void setTriplets( const VectorType &Arg, std::vector<std::vector<TripletType>> &Dest, int hessOffset = 0 ) const {
    if ( Arg.size() != 2 * _topology.getNumEdges())
      throw std::length_error( "QuaternionIntegrabilityHessian::apply(): Arg too small!" );

    Dest.resize( _vertexDOF * _numInteriorVertices );
    for ( auto &vTL : Dest )
      vTL.clear();

    std::vector<int> subIndices;
    if (_onlyVec) {
      subIndices = {0,1,2};
    }
    else {
      subIndices = {0,1,2,3};
    }

    timings["computePartialSecondDerivatives"] = 0.;
    timings["pushToLists"] = 0.;
    timings["Local Derivatives"] = 0;
    timings["Product Rule"] = 0;


#ifdef GOAST_WITH_OPENMP
#pragma omp parallel for
#endif
    for ( int ivIdx = 0; ivIdx < _numInteriorVertices; ivIdx++ ) {
      int vertexIdx = interiorVertices[ivIdx];

      const auto &ADOF = combinedAssocDOF[vertexIdx];
      const auto &numADOF = numAssocDOF[vertexIdx];

      for ( int t : subIndices ) {
        Dest[_vertexDOF * ivIdx + t].reserve( numADOF * numADOF );
      }

      // Create space for  second derivatives
      std::vector<std::vector<RotationType, Eigen::aligned_allocator<RotationType>>> partialSecondDerivatives;
      partialSecondDerivatives.resize( numAssocDOF[vertexIdx] );
      for ( auto &pSD : partialSecondDerivatives ) {
        pSD.resize( numAssocDOF[vertexIdx], QuaternionType( 0, 0, 0, 0 ));
      }


      auto t_start = std::chrono::high_resolution_clock::now();
      computePartialSecondDerivatives( Arg, vertexIdx, partialSecondDerivatives );
      auto t_end = std::chrono::high_resolution_clock::now();
//#pragma omp critical
//      timings["computePartialSecondDerivatives"] += std::chrono::duration<double, std::milli >(t_end - t_start).count();

      t_start = std::chrono::high_resolution_clock::now();
      for ( int i = 0; i < numADOF; i++ )
        for ( int j = 0; j < numADOF; j++ )
          pushToLists(ivIdx, ADOF[i], ADOF[j], partialSecondDerivatives[i][j], Dest, 1., hessOffset);
      t_end = std::chrono::high_resolution_clock::now();

//#pragma omp critical
//      timings["pushToLists"] += std::chrono::duration<double, std::milli >(t_end - t_start).count();

    }
  }


  int getTargetDimension() const {
    return _vertexDOF * _numInteriorVertices;
  }

  int getNNZ() const {
    return _numNonzero;
  }

protected:
  // const RealType factor = 1.,
  void pushToLists( const int ivIdx, const int row, const int col, const QuaternionType &value,
                    std::vector<TripletListType> &tripletLists, const RealType factor = 1., const int hessOffset = 0 ) const {
    if ( _onlyVec ) {
      tripletLists[3 * ivIdx].emplace_back( row + hessOffset, col + hessOffset, factor * value.x());
      tripletLists[3 * ivIdx + 1].emplace_back( row + hessOffset, col + hessOffset, factor * value.y());
      tripletLists[3 * ivIdx + 2].emplace_back( row + hessOffset, col + hessOffset, factor * value.z());
    }
    else {
      tripletLists[4 * ivIdx].emplace_back( row + hessOffset, col + hessOffset, factor * value.w());
      tripletLists[4 * ivIdx + 1].emplace_back( row + hessOffset, col + hessOffset, factor * value.x());
      tripletLists[4 * ivIdx + 2].emplace_back( row + hessOffset, col + hessOffset, factor * value.y());
      tripletLists[4 * ivIdx + 3].emplace_back( row + hessOffset, col + hessOffset, factor * value.z());
    }
  }

  void computePartialSecondDerivatives( const VectorType &Arg, const int vertexIdx,
                                        std::vector<std::vector<RotationType, Eigen::aligned_allocator<RotationType>>> &partialSecondDerivatives ) const {
    const auto &ADOF = combinedAssocDOF[vertexIdx];
    const auto &numADOF = numAssocDOF[vertexIdx];

    auto t_start = std::chrono::high_resolution_clock::now();
    LocalTransitionQuaternionOp<ConfiguratorType> rotationGen;
    LocalTransitionQuaternionGradient<ConfiguratorType> rotationGradGen;
    LocalTransitionQuaternionHessian<ConfiguratorType> rotationGradHess;

    // Computing local gradients of the transition rotations
    RotationType accRotation = RotationType::Identity();
    std::vector<RotationType, Eigen::aligned_allocator<RotationType>> transitionRotations( orientedEdgeRings[vertexIdx].size());
    std::vector<std::array<RotationType, 4>, Eigen::aligned_allocator<std::array<RotationType, 4>>> transitionGradients( orientedEdgeRings[vertexIdx].size());
    std::vector<std::array<std::array<RotationType, 4>, 4>, Eigen::aligned_allocator<std::array<std::array<RotationType, 4>, 4>> > transitionHessians( orientedEdgeRings[vertexIdx].size());
    for ( int i = 0; i < orientedEdgeRings[vertexIdx].size(); i++ ) {
      int edgeIdx = orientedEdgeRings[vertexIdx][i].first;
      auto ael = assocEdgeLengths[vertexIdx][i];
      rotationGen.apply( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )], Arg[std::get<1>( ael )],
                         Arg[std::get<2>( ael )], transitionRotations[i] );
      rotationGradGen.apply( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )], Arg[std::get<1>( ael )],
                             Arg[std::get<2>( ael )], transitionGradients[i] );
      rotationGradHess.apply( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )], Arg[std::get<1>( ael )],
                              Arg[std::get<2>( ael )], transitionHessians[i] );
      accRotation = accRotation * transitionRotations[i];
    }
    auto t_end = std::chrono::high_resolution_clock::now();
//#pragma omp critical
//    timings["Local Derivatives"] += std::chrono::duration<double, std::milli >(t_end - t_start).count();

    t_start = std::chrono::high_resolution_clock::now();

    RotationType postRotation = accRotation;

    // Summands of the first product rule
    for ( int i = 0; i < orientedEdgeRings[vertexIdx].size(); i++ ) {
      // Multiplication of rotations right of the partial derivative
      postRotation = transitionRotations[i].conjugate() * postRotation;

      RotationType preRotation = RotationType::Identity();

      // Summands of the second product rule
      for ( int j = 0; j <= i; j++ ) {
        int edgeIdx_i = orientedEdgeRings[vertexIdx][i].first;
        auto ael_i = assocEdgeLengths[vertexIdx][i];
        auto dof_i = assocDOF[vertexIdx][i];

        int edgeIdx_j = orientedEdgeRings[vertexIdx][j].first;
        auto ael_j = assocEdgeLengths[vertexIdx][j];
        auto dof_j = assocDOF[vertexIdx][j];

        if ( i == j ) {  // Second partial derivatives
          // Compute the partial derivatives by multiplying preRotation * localDeriv * midRotation * localDeriv * postRotation
          RotationType partialDeriv;

          for ( int l : { 0, 1, 2, 3 } ) {
            int localADOF_l = std::distance( ADOF.begin(), std::find( ADOF.begin(), ADOF.end(), dof_i[l] ));
            for ( int h : { 0, 1, 2, 3 } ) {
              int localADOF_h = std::distance( ADOF.begin(), std::find( ADOF.begin(), ADOF.end(), dof_i[h] ));

              partialDeriv = preRotation * transitionHessians[i][h][l] * postRotation;

              partialSecondDerivatives[localADOF_h][localADOF_l].coeffs() += partialDeriv.coeffs();

            }
          }
        }
        else { // Mixed first partial derivatives
          // Multiplication of rotations between the partial derivatives
          RotationType midRotation = RotationType::Identity();
          for ( int k = j + 1; k < i; k++ ) {
            midRotation = midRotation * transitionRotations[k];
          }

          // Compute the partial derivatives by multiplying preRotation * localDeriv * midRotation * localDeriv * postRotation
          RotationType partialDeriv;

          for ( int l : { 0, 1, 2, 3 } ) {
            for ( int h : { 0, 1, 2, 3 } ) {
              int localADOF_l = std::distance( ADOF.begin(), std::find( ADOF.begin(), ADOF.end(), dof_j[l] ));
              int localADOF_h = std::distance( ADOF.begin(), std::find( ADOF.begin(), ADOF.end(), dof_i[h] ));

              partialDeriv = preRotation * transitionGradients[j][l] * midRotation * transitionGradients[i][h] *
                             postRotation;

              partialSecondDerivatives[localADOF_h][localADOF_l].coeffs() += partialDeriv.coeffs();
              partialSecondDerivatives[localADOF_l][localADOF_h].coeffs() += partialDeriv.coeffs();
            }
          }
        }
        // Multiplication of rotations left of the partial derivatives
        preRotation = preRotation * transitionRotations[j];
      }
    }

    t_end = std::chrono::high_resolution_clock::now();

//#pragma omp critical
//    timings["Product Rule"] += std::chrono::duration<double, std::milli >(t_end - t_start).count();

  }

};

#endif //NRIC_INTEGRABILITYQUATERNIONS_H
