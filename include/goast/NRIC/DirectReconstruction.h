// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Header containing classes for the frame-based direct reconstruction
 * \author Sassen
 *
 * Based on
 * Wang, Y., Liu, B., & Tong, Y. (2012). Linear surface reconstruction from discrete fundamental forms on triangle
 * meshes. Computer Graphics Forum, 31(8), 2277–2287.
 * and
 * Sassen, J. (2019). Discrete Gauß-Codazzi Equations for Efficient Triangle Mesh Processing (M.Sc. thesis). University
 * of Bonn.
 *
 * \note This header currently directly uses Eigen types.
 */

#ifndef NRIC_DIRECTRECONSTRUCTION_H
#define NRIC_DIRECTRECONSTRUCTION_H

#include <goast/Core/Auxiliary.h>
#include <goast/Core/Topology.h>
#include "TransitionRotations.h"

#include <queue>

/**
 * \brief Operator reconstructing nodal positions from lengths and angles by traversing a dual spanning tree
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 */
template<typename ConfiguratorType>
class DirectReconstruction
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::RotationType RotationType;
  typedef typename ConfiguratorType::FrameType FrameType;

  const MeshTopologySaver &_topology;
  const int _numVertices;
  const int _numEdges;
  const int _numFaces;

  // Starting information for reconstruction
  const int _initVertexIndex;
  const VecType &_initVertexPosition;

  const int _initFrameIndex;
  const FrameType &_initFrame;

  // Operators and memory for interior angles, frames and transition rotations
  InteriorAngleOp<ConfiguratorType> interiorAngleOp;
  TransitionRotationsOp<ConfiguratorType> transitionRotationsOp;

  mutable VectorType _interiorAngles;
  mutable std::vector<FrameType> _frames;

  // Ordering of reconstruction
  std::vector<int> dualOrdering;
  std::vector<int> edgeOrdering;
  std::vector<int> vertexReconstructionOrdering;

  // Additional stuff useful for visualizations
  std::vector<int> faceOrdering;
  std::vector<int> vertexOrdering;
  std::vector<int> faceLevels;

public:
  /**
   * \brief Construct operator and construct dual spanning tree via breadth-first search
   * \param Topology class containing information on the topology (i.e. connectivity) of the mesh.
   * \param initVertexIndex index of the vertex from which the reconstruction will start
   * \param initVertexPosition prescribed position of the initial vertex
   * \param initFrameIndex index of the face from which the reconstruction will start
   * \param initFrame prescribed frame (=orientation) of the initial face
   */
  DirectReconstruction( const MeshTopologySaver &Topology, const int initVertexIndex,
                        const VecType &initVertexPosition, const int initFrameIndex, const FrameType &initFrame )
          : _topology( Topology ),
            _numVertices( Topology.getNumVertices()), _numEdges( Topology.getNumEdges()),
            _numFaces( Topology.getNumFaces()),
            _initVertexIndex( initVertexIndex ),
            _initVertexPosition( initVertexPosition ),
            _initFrameIndex( initFrameIndex ),
            _initFrame( initFrame ),
            transitionRotationsOp( _topology ), interiorAngleOp( _topology ) {
    _frames.resize( _numFaces, FrameType::Zero());
    _frames[_initFrameIndex] = _initFrame;

    buildDualSpanningTree();
  }

  /**
   * \brief Create spanning tree of the dual graph via breadth-first search
   */
  void buildDualSpanningTree() {
    /**
     * We treverse the dual graph breadth-first (i.e. iterate over the faces going to adjacent ones) starting at the
     * given initial node and create a dual spanning tree. The resulting edges of the primal graph (which correspond to
     * edges of the dual graph) are saved in breadth-first ordering s.t. they can be used for reconstruction. Along the
     * way we also save if and which vertex should be reconstructed in during each step
     */

    // Clear existing order
    dualOrdering.clear();
    edgeOrdering.clear();
    faceOrdering.clear();
    vertexReconstructionOrdering.clear();
    vertexOrdering.clear();
    faceLevels.clear();

    // Keep track of already visited faces and vertices
    std::vector<bool> faceVisited( _numFaces, false );
    std::vector<bool> vertexVisited( _numVertices, false );

    // .. and in which level of the tree they are
    faceLevels.resize( _numFaces, -1 );

    // Queue of the next faces to visit
    std::queue<int> faceQueue;
    faceQueue.emplace( _initFrameIndex );

    faceOrdering.push_back( _initFrameIndex );

    // Mark initial face and its vertices as visited
    faceVisited[_initFrameIndex] = true;
    for ( auto localIdx : { 0, 1, 2 } ) {
      vertexVisited[_topology.getNodeOfTriangle( _initFrameIndex, localIdx )] = true;
      vertexOrdering.push_back( _topology.getNodeOfTriangle( _initFrameIndex, localIdx ));
    }

    faceLevels[_initFrameIndex] = 0;


    while ( !faceQueue.empty()) {
      // Get the next dual vertex = primal face
      const int currentFaceIdx = faceQueue.front();

      // Determine its edges
      std::array<int, 3> edges = { _topology.getEdgeOfTriangle( currentFaceIdx, 0 ),
                                   _topology.getEdgeOfTriangle( currentFaceIdx, 1 ),
                                   _topology.getEdgeOfTriangle( currentFaceIdx, 2 ) };

      for ( const auto &edgeIdx : edges ) {
        // Determine the orientation of the edge, i.e. on which side the current face is, this information will be
        // needed later on to determine if we have to invert the transition rotation of this edge.
        int f1 = _topology.getAdjacentTriangleOfEdge( edgeIdx, 0 );
        int f2 = _topology.getAdjacentTriangleOfEdge( edgeIdx, 1 );

        if ( f1 == -1 || f2 == -1 )
          continue;

        int orientation = (f1 == currentFaceIdx) ? 1 : -1;
        int thisFaceIdx = (f1 == currentFaceIdx) ? f1 : f2;
        int othFaceIdx = (f1 == currentFaceIdx) ? f2 : f1;
        int othVertexIdx = _topology.getOppositeNodeOfEdge( edgeIdx, (f1 == currentFaceIdx));

        if ( !faceVisited[othFaceIdx] ) {
          // If this dual edge does not form a cycle, i.e. visits a new dual vertex, add it to the spanning tree and
          // mark the face ready to visit
          // Add one to the index s.t. we don't get in trouble if it is zero
          //! \todo Clean this up with proper data structures
          dualOrdering.push_back( orientation * (edgeIdx + 1));
          edgeOrdering.push_back( edgeIdx ); // To make calling the transition rotations operator easier

          faceVisited[othFaceIdx] = true;
          faceQueue.emplace( othFaceIdx );

          // Additional output information useful for visualizations
          faceOrdering.push_back( othFaceIdx );
          faceLevels[othFaceIdx] = faceLevels[thisFaceIdx] + 1;

          //! \todo Figure out how to avoid visiting unnecessary faces (so far problems with cactus!)
          if ( !vertexVisited[othVertexIdx] ) {
            // If we get to see a new vertex, save that we want to reconstruct its position
            vertexVisited[othVertexIdx] = true;
            // Additional output information useful for visualizations
            vertexReconstructionOrdering.push_back( othVertexIdx );
            vertexOrdering.push_back( othVertexIdx );
          }
          else {
            vertexReconstructionOrdering.push_back( -1 );
          }
        }
      }
      // Remove the currently visited dual vertex from the queue
      faceQueue.pop();
    }
  }

  /**
   * \brief Create minimal spanning tree of the dual graph
   * \param Weights weights of the dual edges
   */
  void buildDualMST( const VectorType &Weights ) {
    dualOrdering.clear();
    edgeOrdering.clear();
    faceOrdering.clear();
    vertexReconstructionOrdering.clear();
    vertexOrdering.clear();
    faceLevels.clear();

    std::function<bool( std::pair<int, double>, std::pair<int, double> )> compare = []( std::pair<int, double> a,
                                                                                        std::pair<int, double> b ) {
      return a.second > b.second;
    };

    // Keep track of already visited faces and vertices
    std::vector<bool> faceVisited( _numFaces, false );
    std::vector<bool> vertexVisited( _numVertices, false );

    // .. and in which level of the tree they are
    faceLevels.resize( _numFaces, -1 );

    // Queue of the next faces to visit
    std::vector<std::pair<int, double>> edgeQueue;
    faceOrdering.push_back( _initFrameIndex );

    // Mark initial face and its vertices as visited
    faceVisited[_initFrameIndex] = true;
    for ( auto localIdx : { 0, 1, 2 } ) {
      vertexVisited[_topology.getNodeOfTriangle( _initFrameIndex, localIdx )] = true;
      vertexOrdering.push_back( _topology.getNodeOfTriangle( _initFrameIndex, localIdx ));

      const int edgeIdx = _topology.getEdgeOfTriangle( _initFrameIndex, localIdx );
      int f1 = _topology.getAdjacentTriangleOfEdge( edgeIdx, 0 );
      edgeQueue.emplace_back(((f1 == _initFrameIndex) ? 1 : -1) * (edgeIdx + 1), Weights[edgeIdx] );

    }

    faceLevels[_initFrameIndex] = 0;

    std::make_heap( edgeQueue.begin(), edgeQueue.end(), compare );

    while ( !edgeQueue.empty()) {
      // Get the next candidate dual edge
      const int currentEdge = edgeQueue.front().first;
      const int currentEdgeIdx = std::abs( edgeQueue.front().first ) - 1;
      const bool currentOrientation = std::signbit( edgeQueue.front().first );
      const double currentWeight = edgeQueue.front().second;

      // Remove edge from heap
      std::pop_heap( edgeQueue.begin(), edgeQueue.end(), compare );
      edgeQueue.pop_back();

      int f1 = _topology.getAdjacentTriangleOfEdge( currentEdgeIdx, 0 );
      int f2 = _topology.getAdjacentTriangleOfEdge( currentEdgeIdx, 1 );

      if ( f1 == -1 || f2 == -1 )
        continue;

      int thisFaceIdx = !currentOrientation ? f1 : f2;
      int othFaceIdx = !currentOrientation ? f2 : f1;

      int othVertexIdx = _topology.getOppositeNodeOfEdge( currentEdgeIdx, !currentOrientation ? 1 : 0 );


      if ( !faceVisited[othFaceIdx] ) {
        // If this dual edge does not form a cycle, i.e. visits a new dual vertex, add it to the spanning tree and
        // mark the face ready to visit
        // Add one to the index s.t. we don't get in trouble if it is zero
        //! \todo Clean this up with proper data structures
        dualOrdering.push_back( currentEdge );
        edgeOrdering.push_back( currentEdgeIdx ); // To make calling the transition rotations operator easier

        faceOrdering.push_back( othFaceIdx );
        faceVisited[othFaceIdx] = true;
        faceLevels[othFaceIdx] = faceLevels[thisFaceIdx] + 1;

        //! \todo Figure out how to avoid visiting unnecessary faces (so far problems with cactus!)
        if ( !vertexVisited[othVertexIdx] ) {
          // If we get to see a new vertex, save that we want to reconstruct its position
          vertexReconstructionOrdering.push_back( othVertexIdx );
          vertexOrdering.push_back( othVertexIdx );
          vertexVisited[othVertexIdx] = true;
        }
        else {
          vertexReconstructionOrdering.push_back( -1 );
        }

        // Extend heap by edges of other face
        std::array<int, 3> edges = { _topology.getEdgeOfTriangle( othFaceIdx, 0 ),
                                     _topology.getEdgeOfTriangle( othFaceIdx, 1 ),
                                     _topology.getEdgeOfTriangle( othFaceIdx, 2 ) };
        for ( const auto &edgeIdx : edges ) {
          if ( edgeIdx != currentEdgeIdx ) {
            int nf1 = _topology.getAdjacentTriangleOfEdge( edgeIdx, 0 );

            edgeQueue.emplace_back(((nf1 == othFaceIdx) ? 1 : -1) * (edgeIdx + 1), Weights[edgeIdx] );
            std::push_heap( edgeQueue.begin(), edgeQueue.end(), compare );
          }
        }
      }
    }

  }

  /**
   * \brief Create shortest path tree of the dual graph
   * \param Weights weights of the dual edges
   */
  void buildDualSPT( const VectorType &Weights ) {
    dualOrdering.clear();
    edgeOrdering.clear();
    faceOrdering.clear();
    vertexReconstructionOrdering.clear();
    vertexOrdering.clear();
    faceLevels.clear();

    std::function<bool( std::pair<int, double>, std::pair<int, double> )> compare = []( std::pair<int, double> a,
                                                                                        std::pair<int, double> b ) {
      return a.second > b.second;
    };

    // Keep track of already visited faces and vertices
    std::vector<bool> faceVisited( _numFaces, false );
    std::vector<bool> vertexVisited( _numVertices, false );

    // .. and in which level of the tree they are
    faceLevels.resize( _numFaces, -1 );

    // Queue of the next faces to visit
    std::vector<std::pair<int, double>> edgeQueue;
    faceOrdering.push_back( _initFrameIndex );

    // Mark initial face and its vertices as visited
    faceVisited[_initFrameIndex] = true;
    for ( auto localIdx : { 0, 1, 2 } ) {
      vertexVisited[_topology.getNodeOfTriangle( _initFrameIndex, localIdx )] = true;
      vertexOrdering.push_back( _topology.getNodeOfTriangle( _initFrameIndex, localIdx ));

      const int edgeIdx = _topology.getEdgeOfTriangle( _initFrameIndex, localIdx );
      int f1 = _topology.getAdjacentTriangleOfEdge( edgeIdx, 0 );
      edgeQueue.emplace_back(((f1 == _initFrameIndex) ? 1 : -1) * (edgeIdx + 1), Weights[edgeIdx] );

    }

    faceLevels[_initFrameIndex] = 0;

    std::make_heap( edgeQueue.begin(), edgeQueue.end(), compare );

    while ( !edgeQueue.empty()) {
      // Get the next candidate dual edge
      const int currentEdge = edgeQueue.front().first;
      const int currentEdgeIdx = std::abs( edgeQueue.front().first ) - 1;
      const bool currentOrientation = std::signbit( edgeQueue.front().first );
      const double currentWeight = edgeQueue.front().second;

      // Remove edge from heap
      std::pop_heap( edgeQueue.begin(), edgeQueue.end(), compare );
      edgeQueue.pop_back();

      int f1 = _topology.getAdjacentTriangleOfEdge( currentEdgeIdx, 0 );
      int f2 = _topology.getAdjacentTriangleOfEdge( currentEdgeIdx, 1 );

      int thisFaceIdx = !currentOrientation ? f1 : f2;
      int othFaceIdx = !currentOrientation ? f2 : f1;

      int othVertexIdx = _topology.getOppositeNodeOfEdge( currentEdgeIdx, !currentOrientation ? 1 : 0 );

      if ( !faceVisited[othFaceIdx] ) {
        // If this dual edge does not form a cycle, i.e. visits a new dual vertex, add it to the spanning tree and
        // mark the face ready to visit
        // Add one to the index s.t. we don't get in trouble if it is zero
        //! \todo Clean this up with proper data structures
        dualOrdering.push_back( currentEdge );
        edgeOrdering.push_back( currentEdgeIdx ); // To make calling the transition rotations operator easier

        faceOrdering.push_back( othFaceIdx );
        faceVisited[othFaceIdx] = true;
        faceLevels[othFaceIdx] = faceLevels[thisFaceIdx] + 1;

        //! \todo Figure out how to avoid visiting unnecessary faces (so far problems with cactus!)
        if ( !vertexVisited[othVertexIdx] ) {
          // If we get to see a new vertex, save that we want to reconstruct its position
          vertexReconstructionOrdering.push_back( othVertexIdx );
          vertexOrdering.push_back( othVertexIdx );
          vertexVisited[othVertexIdx] = true;
        }
        else {
          vertexReconstructionOrdering.push_back( -1 );
        }

        // Extend heap by edges of other face
        std::array<int, 3> edges = { _topology.getEdgeOfTriangle( othFaceIdx, 0 ),
                                     _topology.getEdgeOfTriangle( othFaceIdx, 1 ),
                                     _topology.getEdgeOfTriangle( othFaceIdx, 2 ) };
        for ( const auto &edgeIdx : edges ) {
          if ( edgeIdx != currentEdgeIdx ) {
            int nf1 = _topology.getAdjacentTriangleOfEdge( edgeIdx, 0 );

            edgeQueue.emplace_back(((nf1 == othFaceIdx) ? 1 : -1) * (edgeIdx + 1), currentWeight + Weights[edgeIdx] );
            std::push_heap( edgeQueue.begin(), edgeQueue.end(), compare );
          }
        }
      }
    }
  }

  /**
   * \brief Update the position of a given vertex from the given frame
   * \param Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param Dest Nodal positions of resulting immersion
   * \param faceIdx Index of the face which is used to construct the nodal position
   * \param vertexIdx Index of the vertex which position will be constructed
   *
   * Implements the construction of edges vectors described in Section 3.4.2 of the author's Master's thesis
   */
  void updateVertex( const VectorType &Arg, VectorType &Dest, const int faceIdx, const int vertexIdx ) const {
    std::array<int, 3> edges = { _topology.getEdgeOfTriangle( faceIdx, 0 ),
                                 _topology.getEdgeOfTriangle( faceIdx, 1 ),
                                 _topology.getEdgeOfTriangle( faceIdx, 2 ) };

    std::array<int, 3> nodes = { _topology.getNodeOfTriangle( faceIdx, 0 ),
                                 _topology.getNodeOfTriangle( faceIdx, 1 ),
                                 _topology.getNodeOfTriangle( faceIdx, 2 ) };

    long localEdgeIdx = std::distance( nodes.begin(), std::find( nodes.begin(), nodes.end(), vertexIdx ));
    long localRecEdgeIdx = (localEdgeIdx + 1) % 3;

    int recEdgeIdx = edges[localRecEdgeIdx];
    int v1 = _topology.getAdjacentNodeOfEdge( recEdgeIdx, 0 );
    int v2 = _topology.getAdjacentNodeOfEdge( recEdgeIdx, 1 );

    long localNodeIdx0 = std::distance( nodes.begin(), std::find( nodes.begin(), nodes.end(), v1 ));
    long localNodeIdx1 = std::distance( nodes.begin(), std::find( nodes.begin(), nodes.end(), v2 ));

    VecType E;

    if ( localRecEdgeIdx == 0 )
      E = VecType( _frames[faceIdx]( 0, 0 ), _frames[faceIdx]( 1, 0 ), _frames[faceIdx]( 2, 0 ));
    else {
      int angleId = (localRecEdgeIdx == 1) ? 3 * faceIdx + 2 : 3 * faceIdx + 1;
      RealType interiorAngle = _interiorAngles[angleId];
      Eigen::Matrix<RealType, 3, 1> rotVec;
      if ( localRecEdgeIdx == 1 ) {
        rotVec = Eigen::Matrix<RealType, 3, 1>( std::cos( interiorAngle ), std::sin( interiorAngle ), 0 );
      }
      else {
        rotVec = Eigen::Matrix<RealType, 3, 1>( std::cos( interiorAngle ), -std::sin( interiorAngle ), 0 );
      }

      Eigen::Matrix<RealType, 3, 1> eigE = _frames[faceIdx] * rotVec;
      E = VecType( eigE( 0 ), eigE( 1 ), eigE( 2 ));
    }

    if ( localNodeIdx0 > localNodeIdx1 && localRecEdgeIdx == 2 )
      E *= -1;

    if ( localNodeIdx0 < localNodeIdx1 && localRecEdgeIdx != 2 )
      E *= -1;

    if ( v1 == vertexIdx ) {
      VecType V1, V2;
      getXYZCoord<VectorType, VecType>( Dest, V2, v2 );
      V1 = V2 - E * Arg[recEdgeIdx];

      setXYZCoord<VectorType, VecType>( Dest, V1, v1 );
    }
    else if ( v2 == vertexIdx ) {
      VecType V1, V2;
      getXYZCoord<VectorType, VecType>( Dest, V1, v1 );

      V2 = V1 + E * Arg[recEdgeIdx];

      setXYZCoord<VectorType, VecType>( Dest, V2, v2 );
    }
  }

  /**
   * \brief Apply construction of immersion to given lengths and angles
   * \param Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param Dest Nodal positions of resulting immersion
   */
  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    Dest.resize( 3 * _numVertices );
    Dest.setZero();

    // Compute interior angles
    interiorAngleOp.apply( Arg, _interiorAngles );

    // Set NaNs to zero
    if ( _interiorAngles.hasNaN())
      for ( int i = 0; i < _interiorAngles.size(); i++ )
        if ( std::isnan( _interiorAngles[i] ))
          _interiorAngles[i] = 0;

    // Compute (necessary) transition rotations
    std::vector<RotationType> transitionRotations;
    transitionRotationsOp.apply( Arg, _interiorAngles, transitionRotations, edgeOrdering );

    // Seed (or initial) frame and vertices
    _frames[_initFrameIndex] = _initFrame;
    setXYZCoord<VectorType, VecType>( Dest, _initVertexPosition, _initVertexIndex );

    // Compute vertex positions of initial face
    std::array<int, 3> nodes = { _topology.getNodeOfTriangle( _initFrameIndex, 0 ),
                                 _topology.getNodeOfTriangle( _initFrameIndex, 1 ),
                                 _topology.getNodeOfTriangle( _initFrameIndex, 2 ) };
    long localVertexIdx = std::distance( nodes.begin(), std::find( nodes.begin(), nodes.end(), _initVertexIndex ));

    updateVertex( Arg, Dest, _initFrameIndex, nodes[(localVertexIdx + 1) % 3] );
    updateVertex( Arg, Dest, _initFrameIndex, nodes[(localVertexIdx + 2) % 3] );

    // Iterate over spanning tree
    for ( int i = 0; i < dualOrdering.size(); i++ ) {
      // Extract id of edge and from which side we are iterating over it
      const int edgeIdx = std::abs( dualOrdering[i] ) - 1;
      const bool orientation = std::signbit( dualOrdering[i] );

      // Get IDs of adjacent triangles
      const int f1 = _topology.getAdjacentTriangleOfEdge( edgeIdx, 0 );
      const int f2 = _topology.getAdjacentTriangleOfEdge( edgeIdx, 1 );

      // Determine new frame
      if ( orientation )
        _frames[f1] = _frames[f2] * transitionRotations[edgeIdx].transpose();
      else
        _frames[f2] = _frames[f1] * transitionRotations[edgeIdx];

      // If vertex should be updated in this iteration then do it
      if ( vertexReconstructionOrdering[i] >= 0 ) {
        updateVertex( Arg, Dest, orientation ? f1 : f2, vertexReconstructionOrdering[i] );
      }
    }
  }

  /**
   * \brief Compute minimal spanning tree and apply construction of immersion to given lengths and angles
   * \param Arg Edges lengths and dihedral angles as vector of size \f$2*|E|\f$
   * \param Dest Nodal positions of resulting immersion
   * \param Weights weights of the dual edges for the MST
   */
  void apply( const VectorType &Arg, VectorType &Dest, const VectorType &Weights ) {
    buildDualMST( Weights );
    apply( Arg, Dest );
  }

  const std::vector<FrameType> &getFrames() const {
    return _frames;
  }

  const std::vector<int> &getFaceOrdering() const {
    return faceOrdering;
  }

  const std::vector<int> &getVertexOrdering() const {
    return vertexOrdering;
  }

  const std::vector<int> &getFaceLevels() const {
    return faceLevels;
  }

  const std::vector<int> &getVertexReconstructionOrdering() const {
    return vertexReconstructionOrdering;
  }
};

/**
 * \brief Extract a specific standard discrete frame from the given mesh
 * \tparam VectorType Vector data type used
 * \tparam FrameType Frame data type used (i.e. type of 3x3 matrix)
 * \param Geometry nodal positions of the mesh as vector
 * \param Topology of the mesh
 * \param faceIdx index of the face which frame will be returned
 * \return extracted frame
 */
template<typename RealType, typename VectorType, typename FrameType>
FrameType extractFrame( const VectorType &Geometry, const MeshTopologySaver &Topology, int faceIdx ) {
  //! \todo get rid of explicit dependence on Eigen
  typedef Eigen::Matrix<RealType, 3, 1> VecType;

  FrameType localFrame = FrameType::Identity();

  std::array<int, 3> nodes = { Topology.getNodeOfTriangle( faceIdx, 0 ),
                               Topology.getNodeOfTriangle( faceIdx, 1 ),
                               Topology.getNodeOfTriangle( faceIdx, 2 ) };

  int e0 = Topology.getEdgeOfTriangle( faceIdx, 0 );
  int p0 = Topology.getAdjacentNodeOfEdge( e0, 0 );
  int p1 = Topology.getAdjacentNodeOfEdge( e0, 1 );

  long localNodeIdx0 = std::distance( nodes.begin(), std::find( nodes.begin(), nodes.end(), p0 ));
  long localNodeIdx1 = std::distance( nodes.begin(), std::find( nodes.begin(), nodes.end(), p1 ));

  VecType P0, P1;
  getXYZCoord<VectorType, VecType>( Geometry, P0, p0 );
  getXYZCoord<VectorType, VecType>( Geometry, P1, p1 );

  // set up deformed vertices and edges
  VecType Pi, Pj, Pk;
  getXYZCoord<VectorType, VecType>( Geometry, Pi, nodes[0] );
  getXYZCoord<VectorType, VecType>( Geometry, Pj, nodes[1] );
  getXYZCoord<VectorType, VecType>( Geometry, Pk, nodes[2] );

  // First edge as first frame vector
  if ( localNodeIdx1 > localNodeIdx0 )
    localFrame.col( 0 ) = P0 - P1;
  else
    localFrame.col( 0 ) = P1 - P0;
  localFrame.col( 0 ).normalize();

  // Normal as third frame vector
  localFrame.col( 2 ) = (Pk - Pj).cross( Pi - Pk );
  localFrame.col( 2 ).normalize();

  // Choosing the second frame vector such that the first two form a orthonormal basis of the tangent space
  localFrame.col( 1 ) = localFrame.col( 0 ).cross( localFrame.col( 2 ));
  localFrame.col( 1 ).normalize();

  return localFrame;
}


/**
 * \brief Extract all standard discrete frames from a given mesh
 * \tparam VectorType Vector data type used
 * \tparam FrameType Frame data type used (i.e. type of 3x3 matrix)
 * \param Geometry nodal positions of the mesh as vector
 * \param Topology of the mesh
 * \return Vector containing the frames
 */
template<typename RealType, typename VectorType, typename FrameType>
std::vector<FrameType> extractFrames( const VectorType &Geometry, const MeshTopologySaver &Topology ) {
  const int numFaces = Topology.getNumFaces();

  std::vector<FrameType> Frames( numFaces );

  for ( int faceIdx = 0; faceIdx < numFaces; faceIdx++ ) {
    Frames[faceIdx] = extractFrame<RealType, VectorType, FrameType>( Geometry, Topology, faceIdx );
  }

  return Frames;
}

#endif //NRIC_DIRECTRECONSTRUCTION_H
