// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Operators yielding the integrability violation of given lengths and angles as Euler angles per vertex
 * \author Sassen
 *
 * Based on
 * Wang, Y., Liu, B., & Tong, Y. (2012). Linear surface reconstruction from discrete fundamental forms on triangle
 * meshes. Computer Graphics Forum, 31(8), 2277–2287.
 * and
 * Sassen, J. (2019). Discrete Gauß-Codazzi Equations for Efficient Triangle Mesh Processing (M.Sc. thesis). University
 * of Bonn.
 *
 * \note This currently relies directly on the Axis-Angle and Quaternion implementations from Eigen
 */

#ifndef NRIC_INTEGRABILITYEULER_H
#define NRIC_INTEGRABILITYEULER_H


/**
 * \brief Operator evaluating discrete integrability conditions using Euler angles on edge lengths and dihedral angles
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 *
 * This operator computes the modified discrete integrability map, i.e. the equation of the discrete integrability
 * conditions, for each vertex, by chaining the transition rotations around the adjacent edges and computing the Euler
 * angles of the resulting rotation.
 *
 * \todo Move computation of the needed topological information to a topology class
 */
template<typename ConfiguratorType>
class EulerIntegrabilityOp
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
  std::vector<std::vector<int>> interiorAngleIdx;
  std::vector<std::vector<std::tuple<int, int, int>>> assocEdgeLengths;

public:
  int _numInteriorVertices;
  std::vector<int> interiorVertices;

  /**
   * \brief Construct operator
   * \param topology class containing information on the topology (i.e. connectivity) of the mesh.
   */
  explicit EulerIntegrabilityOp( const MeshTopologySaver &topology ) : EulerIntegrabilityOp( topology, {} ) {}

  /**
   * \brief Construct operator and ignore certain vertices
   * \param topology Class containing information on the topology (i.e. connectivity) of the mesh.
   * \param ignoredVertices Vector containing the indices of the vertices to be ignored
   */
  EulerIntegrabilityOp( const MeshTopologySaver &topology, const std::vector<int> &ignoredVertices ) :
          _topology( topology ), _numVertices( _topology.getNumVertices()), _numEdges( _topology.getNumEdges()) {
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

      if ( boundaryVertex[vertexIdx] ||
           std::find( ignoredVertices.begin(), ignoredVertices.end(), vertexIdx ) != ignoredVertices.end())
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
   * \brief Evaluate operator
   * \param Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param Dest Euler angles of chained rotation matrices for each vertex as vector
   */
  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    if ( Arg.size() != 2 * _numEdges )
      throw std::length_error( "SimplifiedIntegrabilityOp::apply(): Arg too small!" );

    Dest.resize( 3 * _numInteriorVertices );
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
#ifdef GOAST_WITH_OPENMP
#pragma omp parallel for
#endif
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

      Dest[3 * ivIdx] = -std::asin( accRotation( 2, 0 )); // -ArcSin[31]
      Dest[3 * ivIdx + 1] = std::atan2( accRotation( 2, 1 ), accRotation( 2, 2 )); // atan2[32, 33]
      Dest[3 * ivIdx + 2] = std::atan2( accRotation( 1, 0 ), accRotation( 0, 0 )); // atan2[21, 11]
    }

    for ( int j = 0; j < Dest.size(); j++ )
      if ( std::isnan( Dest[j] ))
        Dest[j] = std::numeric_limits<RealType>::infinity();

  }

  int getTargetDimension() const {
    return 3 * _numInteriorVertices;
  }
};

/**
 * \brief Jacobian of the modified discrete integrability map on edge lengths and dihedral angles
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 * \see EulerIntegrabilityOp
 */
template<typename ConfiguratorType>
class EulerIntegrabilityGradient
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::RotationType RotationType;
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

public:
  /**
   * \brief Construct operator
   * \param topology class containing information on the topology (i.e. connectivity) of the mesh.
   */
  explicit EulerIntegrabilityGradient( const MeshTopologySaver &topology ) : EulerIntegrabilityGradient( topology,
                                                                                                         {} ) {}

  /**
   * \brief Construct operator and ignore certain vertices
   * \param topology Class containing information on the topology (i.e. connectivity) of the mesh.
   * \param ignoredVertices Vector containing the indices of the vertices to be ignored
   */
  EulerIntegrabilityGradient( const MeshTopologySaver &topology, const std::vector<int> &ignoredVertices ) :
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

      if ( boundaryVertex[vertexIdx] ||
           std::find( ignoredVertices.begin(), ignoredVertices.end(), vertexIdx ) != ignoredVertices.end())
        continue;


      interiorVertices.push_back( vertexIdx );
      _numInteriorVertices++;


      orientedEdgeRings[vertexIdx] = orderedEdges;

      _numNonzero += 3 * 3 * orderedEdges.size();
      _numTriplets += 3 * 4 * orderedEdges.size();

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

//    interiorVertices.pop_back();
//    _numInteriorVertices--;

  }


  /**
   * \brief Evaluate derivative
   * \param[in] Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param[out] Dest Jacobian of discrete integrability map
   */
  void apply( const VectorType &Arg, MatrixType &Dest ) const override {
    if ( Arg.size() != 2 * _topology.getNumEdges())
      throw std::length_error( "SimplifiedIntegrabilityGradient::apply(): Arg too small!" );

    Dest.resize( 3 * _numInteriorVertices, 2 * _numEdges );
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
      throw std::length_error( "SimplifiedIntegrabilityGradient::pushTriplets(): Arg too small!" );

    LocalTransitionRotationOp<ConfiguratorType> rotationGen;
    LocalTransitionRotationGradient<ConfiguratorType> rotationGradGen;
    Dest.reserve( _numTriplets );

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int ivIdx = 0; ivIdx < _numInteriorVertices; ivIdx++ ) {
      int vertexIdx = interiorVertices[ivIdx];
      TripletListType localTripletList;
      localTripletList.reserve( 3 * 4 * orientedEdgeRings[vertexIdx].size());

      // Computing local gradients of the transition rotations
      RotationType accRotation = RotationType::Identity();
      std::vector<RotationType> transitionRotations( orientedEdgeRings[vertexIdx].size());
      std::vector<std::array<RotationType, 4>> transitionGradients( orientedEdgeRings[vertexIdx].size());
      for ( int i = 0; i < orientedEdgeRings[vertexIdx].size(); i++ ) {
        int edgeIdx = orientedEdgeRings[vertexIdx][i].first;
        auto ael = assocEdgeLengths[vertexIdx][i];
        transitionRotations[i] = rotationGen( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )],
                                              Arg[std::get<1>( ael )], Arg[std::get<2>( ael )] );
        rotationGradGen.apply( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )], Arg[std::get<1>( ael )],
                               Arg[std::get<2>( ael )], transitionGradients[i] );
        accRotation = accRotation * transitionRotations[i];
      }

      RotationType preRotation = RotationType::Identity();
      RotationType postRotation = accRotation;

      // Summands of the product rule
      for ( int i = 0; i < orientedEdgeRings[vertexIdx].size(); i++ ) {
        // Multiplication of rotations right of the partial derivative
        postRotation = transitionRotations[i].transpose() * postRotation;

        int edgeIdx = orientedEdgeRings[vertexIdx][i].first;
        auto ael = assocEdgeLengths[vertexIdx][i];

        // Compute the partial derivatives by multiplying preRotation * transitionGradient * postRotation
        RotationType partialDeriv;
        RealType value;
        RealType Dasin = -1. / std::sqrt( 1 - accRotation( 2, 0 ) * accRotation( 2, 0 )); // derivative of arcsin term
        RealType denom1 = ( accRotation( 2, 1 ) * accRotation( 2, 1 ) + accRotation( 2, 2 ) * accRotation( 2, 2 ));
        RealType denom2 = ( accRotation( 1, 0 ) * accRotation( 1, 0 ) + accRotation( 0, 0 ) * accRotation( 0, 0 ));

        // Partial derivative w.r.t. dihedral angle
        partialDeriv = preRotation * transitionGradients[i][0] * postRotation;

        // Chain rule for Euler angles
        /// \todo Turn this into a better readble loop
        // -ArcSin[20] -- -arcsin'(20) * 20'
        value = Dasin * partialDeriv( 2, 0 );
        pushToList<transposed>( 3 * ivIdx, _numEdges + edgeIdx, value, localTripletList, factor, rowOffset, colOffset );

        // atan2[21, 22] -- d_1 atan2(21,22) * 21' + d_2 atan2(21,22) * 22' // x=22, y=21
        value = accRotation( 2, 2 ) / denom1 * partialDeriv( 2, 1 ) -
                accRotation( 2, 1 ) / denom1 * partialDeriv( 2, 2 );
        pushToList<transposed>( 3 * ivIdx + 1, _numEdges + edgeIdx, value, localTripletList, factor, rowOffset,
                                colOffset );

        // atan2[10, 00] -- d_1 atan2(10,00) * 10' + d_2 atan2(10,00) * 00' // x=00, y=10
        value = accRotation( 0, 0 ) / denom2 * partialDeriv( 1, 0 ) -
                accRotation( 1, 0 ) / denom2 * partialDeriv( 0, 0 );
        pushToList<transposed>( 3 * ivIdx + 2, _numEdges + edgeIdx, value, localTripletList, factor, rowOffset,
                                colOffset );

        // Partial derivative w.r.t. first edge length
        partialDeriv = preRotation * transitionGradients[i][1] * postRotation;

        // -ArcSin[20] -- -arcsin'(20) * 20'
        value = Dasin * partialDeriv( 2, 0 );
        pushToList<transposed>( 3 * ivIdx, std::get<0>( ael ), value, localTripletList, factor, rowOffset, colOffset );

        // atan2[21, 22] -- d_1 atan2(21,22) * 21' + d_2 atan2(21,22) * 22' // x=22, y=21
        value = accRotation( 2, 2 ) / denom1 * partialDeriv( 2, 1 ) -
                accRotation( 2, 1 ) / denom1 * partialDeriv( 2, 2 );
        pushToList<transposed>( 3 * ivIdx + 1, std::get<0>( ael ), value, localTripletList, factor, rowOffset,
                                colOffset );

        // atan2[10, 00] -- d_1 atan2(10,00) * 10' + d_2 atan2(10,00) * 00' // x=00, y=10
        value = accRotation( 0, 0 ) / denom2 * partialDeriv( 1, 0 ) -
                accRotation( 1, 0 ) / denom2 * partialDeriv( 0, 0 );
        pushToList<transposed>( 3 * ivIdx + 2, std::get<0>( ael ), value, localTripletList, factor, rowOffset,
                                colOffset );

        // Partial derivative w.r.t. second edge length
        partialDeriv = preRotation * transitionGradients[i][2] * postRotation;

        // -ArcSin[20] -- -arcsin'(20) * 20'
        value = Dasin * partialDeriv( 2, 0 );
        pushToList<transposed>( 3 * ivIdx, std::get<1>( ael ), value, localTripletList, factor, rowOffset, colOffset );

        // atan2[21, 22] -- d_1 atan2(21,22) * 21' + d_2 atan2(21,22) * 22' // x=22, y=21
        value = accRotation( 2, 2 ) / denom1 * partialDeriv( 2, 1 ) -
                accRotation( 2, 1 ) / denom1 * partialDeriv( 2, 2 );
        pushToList<transposed>( 3 * ivIdx + 1, std::get<1>( ael ), value, localTripletList, factor, rowOffset,
                                colOffset );

        // atan2[10, 00] -- d_1 atan2(10,00) * 10' + d_2 atan2(10,00) * 00' // x=00, y=10
        value = accRotation( 0, 0 ) / denom2 * partialDeriv( 1, 0 ) -
                accRotation( 1, 0 ) / denom2 * partialDeriv( 0, 0 );
        pushToList<transposed>( 3 * ivIdx + 2, std::get<1>( ael ), value, localTripletList, factor, rowOffset,
                                colOffset );

        // Partial derivative w.r.t. third edge length
        partialDeriv = preRotation * transitionGradients[i][3] * postRotation;

        // -ArcSin[20] -- -arcsin'(20) * 20'
        value = Dasin * partialDeriv( 2, 0 );
        pushToList<transposed>( 3 * ivIdx, std::get<2>( ael ), value, localTripletList, factor, rowOffset, colOffset );

        // atan2[21, 22] -- d_1 atan2(21,22) * 21' + d_2 atan2(21,22) * 22' // x=22, y=21
        value = accRotation( 2, 2 ) / denom1 * partialDeriv( 2, 1 ) -
                accRotation( 2, 1 ) / denom1 * partialDeriv( 2, 2 );
        pushToList<transposed>( 3 * ivIdx + 1, std::get<2>( ael ), value, localTripletList, factor, rowOffset,
                                colOffset );

        // atan2[10, 00] -- d_1 atan2(10,00) * 10' + d_2 atan2(10,00) * 00' // x=00, y=10
        value = accRotation( 0, 0 ) / denom2 * partialDeriv( 1, 0 ) -
                accRotation( 1, 0 ) / denom2 * partialDeriv( 0, 0 );
        pushToList<transposed>( 3 * ivIdx + 2, std::get<2>( ael ), value, localTripletList, factor, rowOffset,
                                colOffset );

        // Multiplication of rotations left of the partial derivative
        preRotation = preRotation * transitionRotations[i];

      }

#ifdef _OPENMP
#pragma omp critical
#endif
      Dest.insert( Dest.end(), localTripletList.begin(), localTripletList.end());
    }
  }

  int getTargetDimension() const override {
    return 3 * _numInteriorVertices;
  }

  int getNNZ() const override {
    return _numNonzero;
  }

private:
  template<bool transposed = false>
  void pushToList( const int row, const int col, const RealType value, TripletListType &tripletList,
                   const RealType factor = 1., const int rowOffset = 0, const int colOffset = 0 ) const {
    if ( transposed )
      tripletList.emplace_back( col + rowOffset, row + colOffset, factor * value );
    else
      tripletList.emplace_back( row + rowOffset, col + colOffset, factor * value );
  }
};

/**
 * \brief Hessian of the modified discrete integrability map on edge lengths and dihedral angles
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 * \see EulerIntegrabilityOp
 */
template<typename ConfiguratorType>
class EulerIntegrabilityHessian
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::TensorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
  typedef typename ConfiguratorType::RotationType RotationType;
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
  std::vector<std::vector<int>> combinedAssocDOF;
  std::vector<std::vector<int>> caDOF_localToGlobal;

  std::vector<int> interiorVertices;

public:
  int _numInteriorVertices;
  mutable std::map<std::string, RealType> timings;

/**
 * \brief Construct operator
 * \param topology class containing information on the topology (i.e. connectivity) of the mesh.
 */
  EulerIntegrabilityHessian( const MeshTopologySaver &topology ) : EulerIntegrabilityHessian( topology, {} ) {}

/**
 * \brief Construct operator and ignore certain vertices
 * \param topology Class containing information on the topology (i.e. connectivity) of the mesh.
 * \param ignoredVertices Vector containing the indices of the vertices to be ignored
 */
  EulerIntegrabilityHessian( const MeshTopologySaver &topology, const std::vector<int> &ignoredVertices ) :
          _topology( topology ),
          _numVertices( _topology.getNumVertices()),
          _numEdges( _topology.getNumEdges()) {
    boundaryVertex.resize( _numVertices, false );
    orientedEdgeRings.resize( _numVertices, {} );
    interiorAngleIdx.resize( _numVertices, {} );
    assocEdgeLengths.resize( _numVertices, {} );
    assocDOF.resize( _numVertices, {} );
    numAssocDOF.resize( _numVertices, 0 );
    combinedAssocDOF.resize( _numVertices, {} );
    caDOF_localToGlobal.resize( _numVertices, {} );

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

      if ( boundaryVertex[vertexIdx] ||
           std::find( ignoredVertices.begin(), ignoredVertices.end(), vertexIdx ) != ignoredVertices.end())
        continue;

      interiorVertices.push_back( vertexIdx );
      _numInteriorVertices++;

      orientedEdgeRings[vertexIdx] = orderedEdges;

      _numNonzero += 9 * orderedEdges.size() * orderedEdges.size();
      _numTriplets += 2 * 16 * orderedEdges.size() * ( orderedEdges.size() + 1 ) / 2;

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
        int li_0 = ( vIdx_f2 == 0 ) ? e2[1] : ( vIdx_f2 == 1 ) ? e2[2] : e2[0];
        int li_1 = ( vIdx_f2 == 0 ) ? e2[2] : ( vIdx_f2 == 1 ) ? e2[0] : e2[1];

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
      throw std::length_error( "EulerIntegrabilityHessian::apply(): Arg too small!" );

    if ( hessSize == -1 )
      Dest.resize( 3 * _numInteriorVertices, 2 * _numEdges, 2 * _numEdges );
    else
      Dest.resize( 3 * _numInteriorVertices, hessSize, hessSize );

    std::vector<TripletListType> vertexTripletLists;
    vertexTripletLists.resize( 3 * _numInteriorVertices );

    setTriplets( Arg, vertexTripletLists, hessOffset );


#ifdef GOAST_WITH_OPENMP
#pragma omp parallel for
#endif
    for ( int ivIdx = 0; ivIdx < _numInteriorVertices; ivIdx++ ) {
      for ( int t : { 0, 1, 2 } ) {
        Dest[3 * ivIdx + t].reserve( vertexTripletLists[3 * ivIdx + t].size());
        Dest[3 * ivIdx + t].setFromTriplets( vertexTripletLists[3 * ivIdx + t].begin(),
                                             vertexTripletLists[3 * ivIdx + t].end());
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
  void apply( const VectorType &Arg, std::vector<MatrixType> &Dest,
              int hessSize = -1, int hessOffset = 0 ) const {
    if ( Arg.size() != 2 * _topology.getNumEdges())
      throw std::length_error( "EulerIntegrabilityHessian::apply(): Arg too small!" );

    Dest.resize( 3 * _numInteriorVertices );

    std::vector<TripletListType> vertexTripletLists;
    vertexTripletLists.resize( 3 * _numInteriorVertices );

    setTriplets( Arg, vertexTripletLists, hessOffset );


#ifdef GOAST_WITH_OPENMP
#pragma omp parallel for
#endif
    for ( int ivIdx = 0; ivIdx < _numInteriorVertices; ivIdx++ ) {
      for ( int t : { 0, 1, 2 } ) {
        if ( hessSize == -1 )
          Dest[3 * ivIdx + t].resize( 2 * _numEdges, 2 * _numEdges );
        else
          Dest[3 * ivIdx + t].resize( hessSize, hessSize );

        Dest[3 * ivIdx + t].reserve( vertexTripletLists[3 * ivIdx + t].size());

        Dest[3 * ivIdx + t].setFromTriplets( vertexTripletLists[3 * ivIdx + t].begin(),
                                             vertexTripletLists[3 * ivIdx + t].end());
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
      throw std::length_error( "EulerIntegrabilityHessian::pushTriplets(): Arg too small!" );

    std::vector<TripletListType> vertexTripletLists;
    vertexTripletLists.resize( 3 * _numInteriorVertices );

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
      throw std::length_error( "EulerIntegrabilityHessian::pushTriplets(): Arg too small! " );
    if ( Lambda.size() != 3 * _numInteriorVertices )
      throw std::length_error( "EulerIntegrabilityHessian::pushTriplets(): Lambda too small!" );

    std::vector<TripletListType> vertexTripletLists;
    vertexTripletLists.resize( 3 * _numInteriorVertices );

    auto t_start = std::chrono::high_resolution_clock::now();
    setTriplets( Arg, vertexTripletLists, hessOffset );
    auto t_end = std::chrono::high_resolution_clock::now();
    timings["setTriplets"] = std::chrono::duration<double, std::milli>( t_end - t_start ).count();

    int totalNumTriplets = 0;
    for ( const auto &vTL : vertexTripletLists ) {
      totalNumTriplets += vTL.size();
    }

    Dest.reserve( Dest.size() + totalNumTriplets );

    t_start = std::chrono::high_resolution_clock::now();
    for ( int i = 0; i < 3 * _numInteriorVertices; i++ ) {
      for ( const auto &trip : vertexTripletLists[i] )
        Dest.emplace_back( trip.row(), trip.col(), Lambda[i] * trip.value());
    }
    t_end = std::chrono::high_resolution_clock::now();
    timings["combineTriplets"] = std::chrono::duration<double, std::milli>( t_end - t_start ).count();

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
      throw std::length_error( "EulerIntegrabilityHessian::setTriplets(): Arg too small!" );

    Dest.resize( 3 * _numInteriorVertices );
    for ( auto &vTL : Dest )
      vTL.clear();

    LocalTransitionRotationOp<ConfiguratorType> rotationGen;
    LocalTransitionRotationGradient<ConfiguratorType> rotationGradGen;
    LocalTransitionRotationHessian<ConfiguratorType> rotationGradHess;

    timings["Local Derivatives"] = 0.;
    timings["Chain Rule"] = 0.;
    timings["Product Rule"] = 0.;


#ifdef GOAST_WITH_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for ( int ivIdx = 0; ivIdx < _numInteriorVertices; ivIdx++ ) {
      int vertexIdx = interiorVertices[ivIdx];

      const auto &ADOF = combinedAssocDOF[vertexIdx];
      const auto &numADOF = numAssocDOF[vertexIdx];

      for ( int t : { 0, 1, 2 } ) {
        Dest[3 * ivIdx + t].reserve( numADOF * numADOF );
      }

      auto t_start = std::chrono::high_resolution_clock::now();
      // Create space for first and second derivatives
      std::vector<RotationType> partialFirstDerivatives;
      partialFirstDerivatives.resize( numAssocDOF[vertexIdx], RotationType::Zero());

      std::vector<std::vector<RotationType>> partialSecondDerivatives;
      partialSecondDerivatives.resize( numAssocDOF[vertexIdx] );
      for ( auto &pSD : partialSecondDerivatives ) {
        pSD.resize( numAssocDOF[vertexIdx], RotationType::Zero());
      }


      // Computing local gradients of the transition rotations
      RotationType accRotation = RotationType::Identity();
      std::vector<RotationType> transitionRotations( orientedEdgeRings[vertexIdx].size());
      std::vector<std::array<RotationType, 4>> transitionGradients( orientedEdgeRings[vertexIdx].size());
      std::vector<std::array<std::array<RotationType, 4>, 4>> transitionHessians( orientedEdgeRings[vertexIdx].size());
      for ( int i = 0; i < orientedEdgeRings[vertexIdx].size(); i++ ) {
        int edgeIdx = orientedEdgeRings[vertexIdx][i].first;
        auto ael = assocEdgeLengths[vertexIdx][i];
        transitionRotations[i] = rotationGen( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )],
                                              Arg[std::get<1>( ael )], Arg[std::get<2>( ael )] );
        rotationGradGen.apply( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )], Arg[std::get<1>( ael )],
                               Arg[std::get<2>( ael )], transitionGradients[i] );
        rotationGradHess.apply( Arg[_numEdges + edgeIdx], Arg[std::get<0>( ael )], Arg[std::get<1>( ael )],
                                Arg[std::get<2>( ael )], transitionHessians[i] );
        accRotation = accRotation * transitionRotations[i];
      }
      auto t_end = std::chrono::high_resolution_clock::now();
//#pragma omp critical
//      timings["Local Derivatives"] += std::chrono::duration<double, std::milli>( t_end - t_start ).count();


      t_start = std::chrono::high_resolution_clock::now();
      // Summands of the first product rule
      for ( int i = 0; i < orientedEdgeRings[vertexIdx].size(); i++ ) {
        // Multiplication of rotations right of the partial derivatives
        RotationType postRotation = RotationType::Identity();
        for ( int k = i + 1; k < orientedEdgeRings[vertexIdx].size(); k++ ) {
          postRotation = postRotation * transitionRotations[k];
        }

        RotationType preRotation = RotationType::Identity();

        // Summands of the second product rule
        for ( int j = 0; j <= i; j++ ) {
//          std::cout << ivIdx << " - " << i << "," << j << std::endl;
          int edgeIdx_i = orientedEdgeRings[vertexIdx][i].first;
          auto ael_i = assocEdgeLengths[vertexIdx][i];
          auto dof_i = assocDOF[vertexIdx][i];

          int edgeIdx_j = orientedEdgeRings[vertexIdx][j].first;
          auto ael_j = assocEdgeLengths[vertexIdx][j];
          auto dof_j = assocDOF[vertexIdx][j];

          if ( i == j ) {  // Second partial derivatives
            // Compute the partial derivatives by multiplying preRotation * localDeriv * midRotation * localDeriv * postRotation
            RotationType partialDeriv( 3, 3 );

            for ( int l : { 0, 1, 2, 3 } ) {
              int localADOF_l = std::distance( ADOF.begin(), std::find( ADOF.begin(), ADOF.end(), dof_i[l] ));
              partialFirstDerivatives[localADOF_l].noalias() += preRotation * transitionGradients[i][l] * postRotation;

              for ( int h : { 0, 1, 2, 3 } ) {
                int localADOF_h = std::distance( ADOF.begin(), std::find( ADOF.begin(), ADOF.end(), dof_i[h] ));

                partialSecondDerivatives[localADOF_h][localADOF_l].noalias() +=
                        preRotation * transitionHessians[i][h][l] * postRotation;
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
            RotationType partialDeriv( 3, 3 );

            for ( int l : { 0, 1, 2, 3 } ) {
              for ( int h : { 0, 1, 2, 3 } ) {
                int localADOF_l = std::distance( ADOF.begin(), std::find( ADOF.begin(), ADOF.end(), dof_j[l] ));
                int localADOF_h = std::distance( ADOF.begin(), std::find( ADOF.begin(), ADOF.end(), dof_i[h] ));

                partialDeriv.noalias() =
                        preRotation * transitionGradients[j][l] * midRotation * transitionGradients[i][h] *
                        postRotation;

                partialSecondDerivatives[localADOF_h][localADOF_l] += partialDeriv;
                partialSecondDerivatives[localADOF_l][localADOF_h] += partialDeriv;
              }
            }
          }
          // Multiplication of rotations left of the partial derivatives
          preRotation = preRotation * transitionRotations[j];
        }
      }
      t_end = std::chrono::high_resolution_clock::now();

//#pragma omp critical
//      timings["Product Rule"] += std::chrono::duration<double, std::milli>( t_end - t_start ).count();

      t_start = std::chrono::high_resolution_clock::now();
      // Chain rule for Euler angles
      for ( int i = 0; i < numADOF; i++ ) {
        for ( int j = 0; j < numADOF; j++ ) {
          const auto &pFD_i = partialFirstDerivatives[i];
          const auto &pFD_j = partialFirstDerivatives[j];
          const auto &pSD = partialSecondDerivatives[i][j];

          RealType denom = std::sqrt( 1 - accRotation( 2, 0 ) * accRotation( 2, 0 ));
          denom = denom * denom * denom;
          RealType nom = accRotation( 2, 0 ) * pFD_i( 2, 0 ) * pFD_j( 2, 0 ) +
                         ( 1 - accRotation( 2, 0 ) * accRotation( 2, 0 )) * pSD( 2, 0 );

          Dest[3 * ivIdx + 0].emplace_back( ADOF[i] + hessOffset, ADOF[j] + hessOffset, -nom / denom );

          denom = accRotation( 2, 1 ) * accRotation( 2, 1 ) + accRotation( 2, 2 ) * accRotation( 2, 2 );

          RealType value = ( pFD_i( 2, 1 ) * pFD_j( 2, 2 ) - pFD_j( 2, 1 ) * pFD_i( 2, 2 )) / denom;
          value += ( pSD( 2, 1 ) * accRotation( 2, 2 ) - accRotation( 2, 1 ) * pSD( 2, 2 )) / denom;
          value -= ( 2 * accRotation( 2, 1 ) * pFD_j( 2, 1 ) + 2 * accRotation( 2, 2 ) * pFD_j( 2, 2 )) *
                   ( pFD_i( 2, 1 ) * accRotation( 2, 2 ) - accRotation( 2, 1 ) * pFD_i( 2, 2 )) / ( denom * denom );

          Dest[3 * ivIdx + 1].emplace_back( ADOF[i] + hessOffset, ADOF[j] + hessOffset, value );

          denom = accRotation( 1, 0 ) * accRotation( 1, 0 ) + accRotation( 0, 0 ) * accRotation( 0, 0 );

          value = ( pFD_i( 1, 0 ) * pFD_j( 0, 0 ) - pFD_j( 1, 0 ) * pFD_i( 0, 0 )) / denom;
          value += ( pSD( 1, 0 ) * accRotation( 0, 0 ) - accRotation( 1, 0 ) * pSD( 0, 0 )) / denom;
          value -= ( 2 * accRotation( 1, 0 ) * pFD_j( 1, 0 ) + 2 * accRotation( 0, 0 ) * pFD_j( 0, 0 )) *
                   ( pFD_i( 1, 0 ) * accRotation( 0, 0 ) - accRotation( 1, 0 ) * pFD_i( 0, 0 )) / ( denom * denom );

          Dest[3 * ivIdx + 2].emplace_back( ADOF[i] + hessOffset, ADOF[j] + hessOffset, value );
        }
      }
      t_end = std::chrono::high_resolution_clock::now();

//#pragma omp critical
//      timings["Chain Rule"] += std::chrono::duration<double, std::milli>( t_end - t_start ).count();

    }

  }

  int getTargetDimension() const {
    return 3 * _numInteriorVertices;
  }

  int getNNZ() const {
    return _numNonzero;
  }
};


#endif //NRIC_INTEGRABILITYEULER_H
