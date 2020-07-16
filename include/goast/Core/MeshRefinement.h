// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef MESHREFINEMENT_HH
#define MESHREFINEMENT_HH


//== INCLUDES =================================================================
#include "Auxiliary.h"
#include "Topology.h"



//=================================================================================

/**
 * \brief Subdivide triangle mesh topologically by means of Loop subdivision
 * \author Heeren
 *
 * No geometric update of old vertices, all new vertices are inserted as midpoints of (old) edges.
 * The index map of a new vertex v, which has been defined as midpoint of edge e = (i,j), returns the pair of the two indices (i,j).
 * The index map of an old vertex v is simply the global index of that vertex.
 */
template<typename ConfiguratorType>
void subdivideMeshTopologically( const TriMesh &oldMesh, TriMesh &newMesh, std::vector<std::pair<int, int> > &indexMap,
                                 std::map<int, std::pair<int, bool> > &edgeMap ) {

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;

  MeshTopologySaver Topol( oldMesh );
  VectorType Geometry;
  getGeometry( oldMesh, Geometry );

  int numV = Topol.getNumVertices();
  int numF = Topol.getNumFaces();
  indexMap.clear();

  //
  std::vector<std::vector<int> > idxNewVerts;
  //
  std::vector<bool> faceHasBeenSubdivided( numF );
  for ( int i = 0; i < numF; i++ )
    faceHasBeenSubdivided[i] = false;

  // add old vertices
  for ( int i = 0; i < numV; i++ ) {
    TriMesh::Point newVertex;
    getXYZCoord<VectorType, TriMesh::Point>( Geometry, newVertex, i );
    newMesh.add_vertex( newVertex );
    indexMap.push_back( std::make_pair( i, i ));
  }

  // add new vertices and faces
  for ( int currFaceIdx = 0; currFaceIdx < numF; currFaceIdx++ ) {
    //std::cerr << "-----------------------------------" << std::endl;
    //std::cerr << "face idx " << currFaceIdx << std::endl;
    int currNumV = newMesh.n_vertices();
    std::vector<int> newIdx( 3 );

    // get old vertices of this face
    std::vector<TriMesh::Point> pos( 3 );
    std::vector<int> idx( 3 );
    //std::cerr << "node indices = ";
    for ( int j = 0; j < 3; j++ ) {
      idx[j] = Topol.getNodeOfTriangle( currFaceIdx, j );
      //std::cerr << idx[j] << " ";
      getXYZCoord<VectorType, TriMesh::Point>( Geometry, pos[j], idx[j] );
    }
    //std::cerr << std::endl;

    // add new vertices on three edges (if necessary)
    for ( int j = 0; j < 3; j++ ) {
      int neighboringFaceIdx = Topol.getNeighbourOfTriangle( currFaceIdx, j );
      //std::cerr << "neighboring face " << j << " has index " << neighboringFaceIdx << std::endl;
      // edge j is subdivided now
      if ((neighboringFaceIdx < 0) || !faceHasBeenSubdivided[neighboringFaceIdx] ) {
        //std::cerr << "subdivide local edge " << j << std::endl;
        newIdx[j] = newMesh.n_vertices();
        TriMesh::Point newVertex;
        for ( int k = 0; k < 3; ++k )
          newVertex[k] = (pos[(j + 1) % 3][k] + pos[(j + 2) % 3][k]) / 2.;
        newMesh.add_vertex( newVertex );
        indexMap.emplace_back( idx[(j + 1) % 3], idx[(j + 2) % 3] );
        //std::cerr << "new index " << j << " is " << newIdx[j] << std::endl;
      }
        // edge j has already been subdivided 
      else {            // find local neighboring index of this triangle from neighboring triangle
        int localIdx = -1;
        for ( int k = 0; k < 3; k++ )
          if ( Topol.getNeighbourOfTriangle( neighboringFaceIdx, k ) == currFaceIdx )
            localIdx = k;
        if ( localIdx < 0 )
          throw BasicException( "subdivideMesh: local index is negative - should not happen!" );
        if ( !(neighboringFaceIdx < idxNewVerts.size()))
          throw BasicException( "subdivideMesh: size of vector of new vertices is wrong - should not happen!" );
        newIdx[j] = idxNewVerts[neighboringFaceIdx][localIdx];
      }
    }

    // add four faces, where (0,1,2) are old indices and (i,ii,iii) new indices
    //   2
    //   |\
    //   | \
    //   |  \
    // ii|___\i
    //   |\  |\
    //   | \ | \
    //   |__\|__\
    //   0  iii  1
    //
    // (0, i, iii)
    //std::cerr << newMesh.vertex_handle(idx[0]) << " " << newMesh.vertex_handle(newIdx[2]) << " " <<  newMesh.vertex_handle(newIdx[1]) << std::endl;
    TriMesh::FaceHandle fh;
    fh = newMesh.add_face( newMesh.vertex_handle( idx[0] ), newMesh.vertex_handle( newIdx[2] ),
                           newMesh.vertex_handle( newIdx[1] ));
    int locIdx = 1;
    for (TriMesh::ConstFaceEdgeIter cfe_it = newMesh.cfe_iter( fh ); cfe_it.is_valid(); ++cfe_it) {
      if ( locIdx == 0 ) {
        edgeMap[cfe_it->idx()] = std::make_pair( Topol.getEdgeOfTriangle( currFaceIdx, 0 ), true );
      }
      else if ( locIdx == 1 ) {
        edgeMap[cfe_it->idx()] = std::make_pair( Topol.getEdgeOfTriangle( currFaceIdx, 1 ), false );
      }
      else if ( locIdx == 2 ) {
        edgeMap[cfe_it->idx()] = std::make_pair( Topol.getEdgeOfTriangle( currFaceIdx, 2 ), false );
      }

      locIdx = (locIdx + 1) % 3;
    }


    // (i,1,ii)
    //std::cerr << newMesh.vertex_handle(newIdx[2]) << " " <<  newMesh.vertex_handle(idx[1]) << " " <<  newMesh.vertex_handle(newIdx[0]) << std::endl;
    fh = newMesh.add_face( newMesh.vertex_handle( newIdx[2] ), newMesh.vertex_handle( idx[1] ),
                      newMesh.vertex_handle( newIdx[0] ));
    locIdx = 1;
    for (TriMesh::ConstFaceEdgeIter cfe_it = newMesh.cfe_iter( fh ); cfe_it.is_valid(); ++cfe_it) {
      if ( locIdx == 0 ) {
        edgeMap[cfe_it->idx()] = std::make_pair( Topol.getEdgeOfTriangle( currFaceIdx, 0 ), false );
      }
      else if ( locIdx == 1 ) {
        edgeMap[cfe_it->idx()] = std::make_pair( Topol.getEdgeOfTriangle( currFaceIdx, 1 ), true );
      }
      else if ( locIdx == 2 ) {
        edgeMap[cfe_it->idx()] = std::make_pair( Topol.getEdgeOfTriangle( currFaceIdx, 2 ), false );
      }

      locIdx = (locIdx + 1) % 3;
    }
    // (ii, 2, iii)
    //std::cerr << newMesh.vertex_handle(newIdx[0]) << " " << newMesh.vertex_handle(idx[2]) << " " <<  newMesh.vertex_handle(newIdx[1]) << std::endl;
    fh = newMesh.add_face( newMesh.vertex_handle( newIdx[0] ), newMesh.vertex_handle( idx[2] ),
                      newMesh.vertex_handle( newIdx[1] ));
    locIdx = 1;
    for (TriMesh::ConstFaceEdgeIter cfe_it = newMesh.cfe_iter( fh ); cfe_it.is_valid(); ++cfe_it) {
      if ( locIdx == 0 ) {
        edgeMap[cfe_it->idx()] = std::make_pair( Topol.getEdgeOfTriangle( currFaceIdx, 1 ), false );
      }
      else if ( locIdx == 1 ) {
        edgeMap[cfe_it->idx()] = std::make_pair( Topol.getEdgeOfTriangle( currFaceIdx, 2 ), true );
      }
      else if ( locIdx == 2 ) {
        edgeMap[cfe_it->idx()] = std::make_pair( Topol.getEdgeOfTriangle( currFaceIdx, 0 ), false );
      }

      locIdx = (locIdx + 1) % 3;
    }
    // (i,ii,iii)
    //std::cerr << newMesh.vertex_handle(newIdx[0]) << " " << newMesh.vertex_handle(newIdx[1]) << " " << newMesh.vertex_handle(newIdx[2]) << std::endl;
    fh = newMesh.add_face( newMesh.vertex_handle( newIdx[0] ), newMesh.vertex_handle( newIdx[1] ),
                      newMesh.vertex_handle( newIdx[2] ));

    //
    faceHasBeenSubdivided[currFaceIdx] = true;
    idxNewVerts.push_back( newIdx );
  }

}

//! \brief See documentation of previous class.
//! \author Heeren
template<typename ConfiguratorType>
void subdivideMeshTopologically( const TriMesh &oldMesh, TriMesh &newMesh ) {
  std::vector<std::pair<int, int> > indexMap;
  std::map<int, std::pair<int, bool> >  edgeMap;
  subdivideMeshTopologically<ConfiguratorType>( oldMesh, newMesh, indexMap, edgeMap );
}

/**
 * \brief Prolongate one-dimensional data vector
 * \author Heeren
 *
 * Assume a mesh with n vertices has been subdivided to a mesh with N vertices (by one of the classes above).
 * Then the coarse n-dim. vector \f$X\f$ is prolongated to the fine N-dim. vector \f$Y\f$.
 * If the i-th vertex was already part of the coarse mesh, we set \f$Y[i] = X[i]\f$.
 * If the i-th vertex has been inserted as midpoint of edge $e = (k,l)$, we set \f$Y[i] = 0.5 \cdot (X[k] + X[l])\f$.
 */
template<typename VectorType>
void prolongateVectorToSubdividedGeometry( const std::vector< std::pair<int,int> >& indexMap, const VectorType& CoarseVector, VectorType& FineVector ){
    int newSize = FineVector.size();
    FineVector = CoarseVector;
    FineVector.conservativeResize( newSize );
    
    int numOldVerts = CoarseVector.size();
    int numNewVerts = newSize - numOldVerts;
    for( int i = 0; i < numNewVerts; i++ ){
        if( !( numOldVerts + i < indexMap.size() ) )
            throw BasicException("prolongateVectorToSubdividedGeometry: index map out of bounds!");
        int idx1 = indexMap.at(numOldVerts + i).first;
        int idx2 = indexMap.at(numOldVerts + i).second;
        FineVector[numOldVerts + i] = (CoarseVector[idx1] + CoarseVector[idx2]) / 2.;
    }    
} 

#endif