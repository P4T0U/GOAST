// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
//=============================================================================
//
//  Geometric functions and derivatives related to a flap of triangles
//
//=============================================================================

#ifndef LOCALMESHGEOMETRY_HH
#define LOCALMESHGEOMETRY_HH

#include "Auxiliary.h"
#include "SmallVecMat.h"
#include "TriangleGeometry.h"
#include "Topology.h"

//! \brief Compute squared edge length \f$(l_e)^2\f$ of one particular edge $e$
//! \author Heeren
template<typename ConfiguratorType>
double getEdgeLengthSqr( const MeshTopologySaver& Topology, int edgeIdx, const typename ConfiguratorType::VectorType& Geometry ) {
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  int pi(Topology.getAdjacentNodeOfEdge(edgeIdx, 0)),
      pj(Topology.getAdjacentNodeOfEdge(edgeIdx, 1));
  VecType Pi, Pj;
  getXYZCoord<VectorType, VecType>(Geometry, Pi, pi);
  getXYZCoord<VectorType, VecType>(Geometry, Pj, pj);
  Pi -= Pj;
  return Pi.normSqr();
}   

//! \brief Compute squared edge lengths \f$((l_e)^2)_e\f$ of all edges
//! \author Heeren
template<typename ConfiguratorType>
void getEdgeLengthsSqr( const MeshTopologySaver& Topology, const typename ConfiguratorType::VectorType& Geometry, typename ConfiguratorType::VectorType& SqrEdgeLengths ) {
  SqrEdgeLengths.resize( Topology.getNumEdges() );
  for( int i = 0; i < Topology.getNumEdges(); i++ )
      SqrEdgeLengths[i] = getEdgeLengthSqr<ConfiguratorType>( Topology, i, Geometry );
} 

//! \brief Compute edge length \f$l_e\f$ of one particular edge \f$e\f$
//! \author Heeren
template<typename ConfiguratorType>
double getEdgeLength( const MeshTopologySaver& Topology, int edgeIdx, const typename ConfiguratorType::VectorType& Geometry ) {
  return std::sqrt( getEdgeLengthSqr<ConfiguratorType>(Topology, edgeIdx, Geometry) );
}   

//! \brief Compute edge lengths \f$(l_e)_e\f$ of all edges
//! \author Heeren
template<typename ConfiguratorType>
void getEdgeLengths( const MeshTopologySaver& Topology, const typename ConfiguratorType::VectorType& Geometry, typename ConfiguratorType::VectorType& EdgeLengths ) {
  EdgeLengths.resize( Topology.getNumEdges() );
  for( int i = 0; i < Topology.getNumEdges(); i++ )
      EdgeLengths[i] = getEdgeLength<ConfiguratorType>( Topology, i, Geometry );
} 

/**
 * \brief Compute dual edge length \f$d_e\f$ of one particular edge \f$e\f$
 * \author Heeren
 * \tparam ConfiguratorType Container with data types
 *
 * If \f$e = (i,j)\f$ is edge with neighboring vertices \f$k\f$ and \f$l\f$, then \f$d_e = \| x_k - x_l \|\f$.
 */
template<typename ConfiguratorType>
double getDualEdgeLength( const MeshTopologySaver& Topology, int edgeIdx, const typename ConfiguratorType::VectorType& Geometry );

template<typename ConfiguratorType>
double getDualEdgeLength(const MeshTopologySaver &Topology, int edgeIdx, const typename ConfiguratorType::VectorType &Geometry) {
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    int pi(Topology.getOppositeNodeOfEdge(edgeIdx, 0)),
            pj(Topology.getOppositeNodeOfEdge(edgeIdx, 1));
    VecType Pi, Pj;
    getXYZCoord<VectorType, VecType>(Geometry, Pi, pi);
    getXYZCoord<VectorType, VecType>(Geometry, Pj, pj);
    Pi -= Pj;
    return Pi.norm();
}

//! \brief Compute dual edge lengths \f$(d_e)_e\f$ of all dual edges
//! \author Heeren
template<typename ConfiguratorType>
void getDualEdgeLengths( const MeshTopologySaver& Topology, const typename ConfiguratorType::VectorType& Geometry, typename ConfiguratorType::VectorType& dualEdgeLengths ) {
  dualEdgeLengths.resize( Topology.getNumEdges() );
  for( int i = 0; i < Topology.getNumEdges(); i++ )
      dualEdgeLengths[i] = getDualEdgeLength<ConfiguratorType>( Topology, i, Geometry );
}

//! \brief Compute normal and area of one particular face
//! \author Heeren
template<typename ConfiguratorType>    
double getNormalAndArea( const MeshTopologySaver& Topology, int faceIdx, const typename ConfiguratorType::VectorType& Geometry, typename ConfiguratorType::VecType& normal  ) {
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  int pi(Topology.getNodeOfTriangle(faceIdx, 0)),
      pj(Topology.getNodeOfTriangle(faceIdx, 1)),
      pk(Topology.getNodeOfTriangle(faceIdx, 2));
  //! get the deformed values
  VecType Pi, Pj, Pk;
  getXYZCoord<VectorType, VecType>(Geometry, Pi, pi);
  getXYZCoord<VectorType, VecType>(Geometry, Pj, pj);
  getXYZCoord<VectorType, VecType>(Geometry, Pk, pk);
                        
  normal.makeCrossProduct( Pk-Pj, Pi-Pk);
  double norm =  normal.norm();
  normal /= norm;
  return norm / 2.;
}

//! \brief Compute face areas of all faces.
//! \author Heeren
template<typename ConfiguratorType>    
void getFaceAreas( const MeshTopologySaver& Topology,
                   const typename ConfiguratorType::VectorType& Geometry,
                   typename ConfiguratorType::VectorType& FaceAreas ) {
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  
  FaceAreas.resize( Topology.getNumFaces() );
  VecType Pi, Pj, Pk;
  
  for( int faceIdx = 0; faceIdx < Topology.getNumFaces(); faceIdx++ ){
  
    int pi(Topology.getNodeOfTriangle(faceIdx, 0)),
        pj(Topology.getNodeOfTriangle(faceIdx, 1)),
        pk(Topology.getNodeOfTriangle(faceIdx, 2));
  
    getXYZCoord<VectorType, VecType>(Geometry, Pi, pi);
    getXYZCoord<VectorType, VecType>(Geometry, Pj, pj);
    getXYZCoord<VectorType, VecType>(Geometry, Pk, pk);
    FaceAreas[faceIdx] = getArea( Pi, Pj, Pk );
  }
}

//! \brief Compute (cosine of) dihedral angles of all edges.
//! \author Heeren
//! \tparam ConfiguratorType Container with data types
template<typename ConfiguratorType>
void getDihedralAngles( const MeshTopologySaver& Topology, 
                        const typename ConfiguratorType::VectorType& ReferenceGeometry, 
                        typename ConfiguratorType::VectorType& DihedralAngles, 
                        bool computeCosine = false) {      
#ifdef DEBUGMODE
      if( hasNanEntries(ReferenceGeometry) )
          throw BasicException("getDihedralAngles(): NaN error in reference geometry!");
#endif   
      
    if( ReferenceGeometry.size() != 3 * Topology.getNumVertices() )  
        throw BasicException("getDihedralAngles(): geometry has wrong size!");
      
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
      
    DihedralAngles.resize(Topology.getNumEdges() ); 
    DihedralAngles.setZero();
      
    // run over all edges
    for ( int edgeIdx = 0; edgeIdx < Topology.getNumEdges(); ++edgeIdx ){
      
      if( !(Topology.isEdgeValid(edgeIdx)) )
	continue;

      int pi( Topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( Topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
          pk( Topology.getOppositeNodeOfEdge(edgeIdx,0) ),
          pl( Topology.getOppositeNodeOfEdge(edgeIdx,1) );

      // no bending at boundary edges
      if( std::min( pl, pk) < 0 )
        continue;

      // set up vertices and edges
      VecType Pi, Pj, Pk, Pl, temp;

      // get reference geometry
      getXYZCoord<VectorType, VecType>( ReferenceGeometry, Pi, pi);
      getXYZCoord<VectorType, VecType>( ReferenceGeometry, Pj, pj);
      getXYZCoord<VectorType, VecType>( ReferenceGeometry, Pk, pk);
      getXYZCoord<VectorType, VecType>( ReferenceGeometry, Pl, pl);

      // get dihedral angle
      DihedralAngles[edgeIdx] = computeCosine ? getCosOfDihedralAngle(Pi, Pj, Pk, Pl) : getDihedralAngle( Pi, Pj, Pk, Pl );  
      
#ifdef DEBUGMODE
      if( std::isnan(DihedralAngles[edgeIdx]) ) ){
          std::cerr << "NaN in simple bending functional in edge " << edgeIdx << "! " << std::endl;
          std::cerr << "theta = " << DihedralAngles[edgeIdx] << std::endl;            
          throw BasicException("getDihedralAngles(): NaN Error!");
      }
#endif
    }
}

//! \brief Compute dihedral angles and corresponding gradients for all edges.
//! \author Heeren
//! \tparam ConfiguratorType Container with data types
template<typename ConfiguratorType>
void getDihedralAnglesGradients( const MeshTopologySaver& Topology, 
                                 const typename ConfiguratorType::VectorType& ReferenceGeometry, 
                                 typename ConfiguratorType::VectorType& DihedralAngles,
                                 std::vector<typename ConfiguratorType::VecType>& gradThetaPi,
                                 std::vector<typename ConfiguratorType::VecType>& gradThetaPj,
                                 std::vector<typename ConfiguratorType::VecType>& gradThetaPk,
                                 std::vector<typename ConfiguratorType::VecType>& gradThetaPl ) {      
#ifdef DEBUGMODE
      if( hasNanEntries(ReferenceGeometry) )
          throw BasicException("getDihedralAnglesGradients(): NaN error in reference geometry!");
#endif   
      
    if( ReferenceGeometry.size() != 3 * Topology.getNumVertices() )  
        throw BasicException("getDihedralAnglesGradients(): geometry has wrong size!");
      
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
      
    gradThetaPi.resize( Topology.getNumEdges() ); 
    gradThetaPj.resize( Topology.getNumEdges() ); 
    gradThetaPk.resize( Topology.getNumEdges() ); 
    gradThetaPl.resize( Topology.getNumEdges() ); 
    
    DihedralAngles.resize( Topology.getNumEdges() ); 
    DihedralAngles.setZero();
      
    // run over all edges
    for ( int edgeIdx = 0; edgeIdx < Topology.getNumEdges(); ++edgeIdx ){
      
      if( !(Topology.isEdgeValid(edgeIdx)) )
	continue;

      int pi( Topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( Topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
          pk( Topology.getOppositeNodeOfEdge(edgeIdx,0) ),
          pl( Topology.getOppositeNodeOfEdge(edgeIdx,1) );

      // no bending at boundary edges
      if( std::min( pl, pk) < 0 )
        continue;

      // set up vertices
      VecType Pi, Pj, Pk, Pl;
      getXYZCoord<VectorType, VecType>( ReferenceGeometry, Pi, pi);
      getXYZCoord<VectorType, VecType>( ReferenceGeometry, Pj, pj);
      getXYZCoord<VectorType, VecType>( ReferenceGeometry, Pk, pk);
      getXYZCoord<VectorType, VecType>( ReferenceGeometry, Pl, pl);
      
      // dihedral angle
      DihedralAngles[edgeIdx] = getDihedralAngle( Pi, Pj, Pk, Pl ); 
  
      // compute gradients
      getThetaGradK(Pi, Pj, Pk, gradThetaPk[edgeIdx]);
      getThetaGradK(Pj, Pi, Pl, gradThetaPl[edgeIdx]);
      getThetaGradI(Pi, Pj, Pk, Pl, gradThetaPi[edgeIdx]);
      getThetaGradJ(Pi, Pj, Pk, Pl, gradThetaPj[edgeIdx]);
    }
}

/**
 * \brief Compute area, area gradients and edges for all faces.
 * \author Heeren
 * \tparam ConfiguratorType Container with data types
 */
template<typename ConfiguratorType>
void getAreaGradients( const MeshTopologySaver& Topology, 
                        const typename ConfiguratorType::VectorType& ReferenceGeometry, 
                        typename ConfiguratorType::VectorType& FaceAreas,
                        std::vector<typename ConfiguratorType::VecType>& gradAreaPi,
                        std::vector<typename ConfiguratorType::VecType>& gradAreaPj,
                        std::vector<typename ConfiguratorType::VecType>& gradAreaPk,
                        std::vector<typename ConfiguratorType::VecType>& edges    ) {
    
#ifdef DEBUGMODE
      if( hasNanEntries(ReferenceGeometry) )
          throw BasicException("getAreaGradients(): NaN error in reference geometry!");
#endif   
      
    if( ReferenceGeometry.size() != 3 * Topology.getNumVertices() )  
        throw BasicException("getAreaGradients(): geometry has wrong size!");
      
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
      
    gradAreaPi.resize( Topology.getNumFaces() ); 
    gradAreaPj.resize( Topology.getNumFaces() ); 
    gradAreaPk.resize( Topology.getNumFaces() ); 
    edges.resize( 3 * Topology.getNumFaces() ); 
    
    FaceAreas.resize( Topology.getNumFaces() ); 
    FaceAreas.setZero();
    
    std::vector<VecType> nodes(3);
    VecType normal;
      
    // run over all faces
    for ( int faceIdx = 0; faceIdx < Topology.getNumFaces(); ++faceIdx ){
      
      //! get reference positions
      for (int j = 0; j < 3; j++)
	getXYZCoord<VectorType, VecType>(ReferenceGeometry, nodes[j], Topology.getNodeOfTriangle(faceIdx, j));
      
      // compute volume
      normal.makeCrossProduct(nodes[2] - nodes[1], nodes[0] - nodes[2]);
      FaceAreas[faceIdx] = normal.norm() / 2.;
      normal /= 4. * FaceAreas[faceIdx];

      // compute gradients
      gradAreaPk[faceIdx].makeCrossProduct( normal, nodes[1] - nodes[0]); 
      gradAreaPi[faceIdx].makeCrossProduct( normal, nodes[2] - nodes[1]);
      gradAreaPj[faceIdx].makeCrossProduct( normal, nodes[0] - nodes[2]);
      
      // compute edges
      for (int j = 0; j < 3; j++)
        edges[3*faceIdx + j] = nodes[(j + 2) % 3] - nodes[(j + 1) % 3];
   }  
}


//======================================================================================================================================

/**
 * \brief Compute edge lengths and dihedral angles
 * \author Heeren
 * \tparam ConfiguratorType Container with data types
 *
 * Consider a mesh (with n vertices and m edges) given as topology information and geometry
 * (where geometry refers to the nodal positions vector in 3D Euclidean space),
 * this method computes the m-dim. vectors of all edge lengths, dihedral angles and face volumes
 */
template<typename ConfiguratorType>
void computeLengthsAngles( const MeshTopologySaver &Topology,
                           const typename ConfiguratorType::VectorType &Geometry,
                           typename ConfiguratorType::VectorType &LengthsAngles ) {

  if ( Geometry.size() != 3 * Topology.getNumVertices())
    throw BasicException( "computeLengthsAngles(): input geometry has wrong size!!!" );

  LengthsAngles.resize( 2 * Topology.getNumEdges());

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::RealType RealType;

  for ( int edgeIdx = 0; edgeIdx < Topology.getNumEdges(); ++edgeIdx ) {

    int pi( Topology.getAdjacentNodeOfEdge( edgeIdx, 0 )),
            pj( Topology.getAdjacentNodeOfEdge( edgeIdx, 1 )),
            pk( Topology.getOppositeNodeOfEdge( edgeIdx, 0 )),
            pl( Topology.getOppositeNodeOfEdge( edgeIdx, 1 ));

    // set up vertices
    VecType Pi, Pj;
    getXYZCoord<VectorType, VecType>( Geometry, Pi, pi );
    getXYZCoord<VectorType, VecType>( Geometry, Pj, pj );

    // compute squared length of edge
    LengthsAngles[edgeIdx] = std::sqrt( dotProduct( Pj - Pi, Pj - Pi ));

    // compute dihedral angle  (no bending at boundary edges)
    if ( std::min( pl, pk ) == -1 )
      continue;

    VecType Pk, Pl;
    getXYZCoord<VectorType, VecType>( Geometry, Pk, pk );
    getXYZCoord<VectorType, VecType>( Geometry, Pl, pl );
    LengthsAngles[Topology.getNumEdges() + edgeIdx] = getDihedralAngle( Pi, Pj, Pk, Pl );
  }

}

#endif
