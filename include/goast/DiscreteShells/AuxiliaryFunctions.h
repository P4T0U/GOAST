// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Helper functions to compute fixed geometric quantities on a reference mesh.
 * \author Heeren
 *
 */
 #ifndef AUXILIARYFUNCTIONS_HH
#define AUXILIARYFUNCTIONS_HH


//== INCLUDES =================================================================
#include <goast/Core/Auxiliary.h>
#include <goast/Core/LocalMeshGeometry.h>
#include <goast/Core/Topology.h>

//#define DEBUGMODE

/**
 * \brief Compute dihedral angles, edge areas and squared edge lengths of a reference shape
 * \author Heeren
 *
 * For given mesh (given as topology plus geometry) with m edges,
 * this function computes the vectors of all dihedral angles, all edge areas and all squared edge lengths,
 * i.e. three vectors living in \f$ \R^m \f$.
 * If edge e has adjacent faces $l$ and $r$ with face areas \f$ a_l \f$ and  \f$ a_r \f$,
 * then the edge area \f$ a_e \f$ is given by \f$ a_e = (a_l + a_r) / 3 \f$.
 * \attention We omit the factor 1/3 here!
 * If computeCosine = true, we actually compute the cosine of the dihedral angles instead.
 */
template<typename ConfiguratorType>
void computeReferenceQuantitiesBending( const MeshTopologySaver& Topology, 
                                 const typename ConfiguratorType::VectorType& ReferenceGeometry, 
                                 typename ConfiguratorType::VectorType& DihedralAngles, 
                                 typename ConfiguratorType::VectorType& EdgeAreas,
                                 typename ConfiguratorType::VectorType& SqrEdgeLengths,
                                 bool computeCosine = false) {      
#ifdef DEBUGMODE
      if( hasNanEntries(ReferenceGeometry) )
          throw BasicException("computeReferenceQuantitiesBending(): NaN error in reference geometry!");
#endif   
      
    if( ReferenceGeometry.size() != 3 * Topology.getNumVertices() )  
        throw BasicException("computeReferenceQuantitiesBending(): geometry has wrong size!");
      
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
      
    DihedralAngles.resize(Topology.getNumEdges() ); 
    DihedralAngles.setZero();
    EdgeAreas.resize(Topology.getNumEdges() ); 
    EdgeAreas.setZero();
    SqrEdgeLengths.resize(Topology.getNumEdges() ); 
    SqrEdgeLengths.setZero();
      
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

      // compute volume and length of edge 
      // A_e = h_e * l_e = (|T_1| + |T_2|)/3 if T_1 and T_2 share edge e
      // Furthermore, |T| = 0.5 * |e_1 x e_2|, if e_1 and e_2 are edges of T
      temp.makeCrossProduct( Pk-Pj, Pi-Pk );
      EdgeAreas[edgeIdx] = temp.norm() / 2.;
      temp.makeCrossProduct( Pl-Pi, Pj-Pl );
      EdgeAreas[edgeIdx] += temp.norm() / 2.;
      SqrEdgeLengths[edgeIdx] = dotProduct(Pj-Pi,Pj-Pi);       
      
#ifdef DEBUGMODE
      if( std::isnan(DihedralAngles[edgeIdx]) || std::isnan(SqrEdgeLengths[edgeIdx]) || std::isnan(EdgeAreas[edgeIdx]) ){
          std::cerr << "NaN in simple bending functional in edge " << edgeIdx << "! " << std::endl;
          std::cerr << "theta = " << DihedralAngles[edgeIdx] << std::endl;            
          std::cerr << "elengthSqr = " << SqrEdgeLengths[edgeIdx] << std::endl;
          std::cerr << "vol = " << EdgeAreas[edgeIdx] << std::endl;  
          throw BasicException("computeReferenceQuantitiesBending(): NaN Error!");
      }
#endif
    }
}

/**
 * \brief Compute integration weights for Discrete Shells length functional
 * \author Heeren
 *
 * Compute for all edges e the integration weight \f$ w_e = a_e / l_e^2 \f$,
 * where \f$ l_e \f$ is the edge length and \f$ a_e \f$ is associated edge area.
 * If edge e has adjacent faces l and r with face areas \f$ a_l \f$ and \f$ a_r \f$,
 * then the edge area \f$ a_e \f$ is given by \f$ a_e = (a_l + a_r) / 3 \f$.
 */
template<typename ConfiguratorType>
void computeIntegrationWeightsLengthFunctional( const MeshTopologySaver& Topology, 
                                                const typename ConfiguratorType::VectorType& ReferenceGeometry, 
                                                typename ConfiguratorType::VectorType& IntegrationWeights ) {      
 
      
    if( ReferenceGeometry.size() != 3 * Topology.getNumVertices() )  
        throw BasicException("computeIntegrationWeightsLengthFunctional(): geometry has wrong size!");
      
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
      
    IntegrationWeights.resize(Topology.getNumEdges() ); 
    IntegrationWeights.setZero();
      
    // run over all edges
    for ( int edgeIdx = 0; edgeIdx < Topology.getNumEdges(); ++edgeIdx ){
      
      if( !(Topology.isEdgeValid(edgeIdx)) )
	continue;

      int pi( Topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( Topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
          pk( Topology.getOppositeNodeOfEdge(edgeIdx,0) ),
          pl( Topology.getOppositeNodeOfEdge(edgeIdx,1) );

      // set up vertices and edges
      VecType Pi, Pj, Pk, Pl, temp;
      getXYZCoord<VectorType, VecType>( ReferenceGeometry, Pi, pi);
      getXYZCoord<VectorType, VecType>( ReferenceGeometry, Pj, pj);      

      // compute volume and length of edge 
      // A_e = h_e * l_e = (|T_1| + |T_2|)/3 if T_1 and T_2 share edge e
      // Furthermore, |T| = 0.5 * |e_1 x e_2|, if e_1 and e_2 are edges of T
      if( !(pk < 0) ){
        getXYZCoord<VectorType, VecType>( ReferenceGeometry, Pk, pk);
        temp.makeCrossProduct( Pk-Pj, Pi-Pk );
        IntegrationWeights[edgeIdx] += temp.norm() / 6.;
      }
      if( !(pl < 0) ){
        getXYZCoord<VectorType, VecType>( ReferenceGeometry, Pl, pl);
        temp.makeCrossProduct( Pl-Pi, Pj-Pl );
        IntegrationWeights[edgeIdx] += temp.norm() / 6.;
      }
      IntegrationWeights[edgeIdx] /= dotProduct(Pj-Pi,Pj-Pi);       

    }
}

/**
 * \brief For all faces, compute face areas and (weighted) dot products of edges
 * \author Heeren
 *
 * For each face $f$, with area \f$ a_f \f$ and edges \f$ e_0 \f$, \f$ e_1 \f$ and  \f$ e_2 \f$,
 * compute three (weighted) dot products \f$ \langle e_{i+2}, e_{i+1} \rangle / a_f \f$, for \f$ i = 0,1,2 \f$.
 */
template<typename ConfiguratorType>    
void getFaceAreasAndFactors( const MeshTopologySaver& Topology, 
                             const typename ConfiguratorType::VectorType& Geometry, 
                             typename ConfiguratorType::VectorType& FaceAreas, 
                             typename ConfiguratorType::VectorType& Factors ) {
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  
  FaceAreas.resize( Topology.getNumFaces() );
  Factors.resize( 3 * Topology.getNumFaces() );
  VecType Pi, Pj, Pk;
  std::vector<VecType> nodes(3), edges(3);
  
  for( int faceIdx = 0; faceIdx < Topology.getNumFaces(); faceIdx++ ){
  
    // compute nodes
    for( int j = 0; j < 3; j++ )
        getXYZCoord<VectorType, VecType>(Geometry, nodes[j], Topology.getNodeOfTriangle(faceIdx, j) );
        
    // compute edges
    for (int j = 0; j < 3; j++)
        edges[j] = nodes[(j + 2) % 3] - nodes[(j + 1) % 3];
      
    // compute face area
    FaceAreas[faceIdx] = getArea( nodes[0], nodes[1], nodes[2] );
    
    // compute factors
    for (int j = 0; j < 3; j++)
         Factors[3*faceIdx + j] = dotProduct(edges[(j + 2) % 3], edges[(j + 1) % 3]) / FaceAreas[faceIdx];

  }
}

/**
 * \brief For all faces, compute (squared) face areas and dot products of edges
 * \author Heeren
 *
 * For each face \f$ f \f$, with area \f$a_f\f$ and edges \f$ e_0 \f$, \f$ e_1 \f$ and \f$ e_2 \f$,
 * compute three dot products \f$ \langle e_{i+2}, e_{i+1} \rangle / a_f \f$, for \f$ i = 0,1,2 \f$.
 */
template<typename ConfiguratorType>
void getDotProductsAndSquaredFaceAreas( const MeshTopologySaver& Topology, 
                                        const typename ConfiguratorType::VectorType& ReferenceGeometry, 
                                        std::vector<typename ConfiguratorType::VecType>& DotProducts,
                                        typename ConfiguratorType::VectorType& FaceAreasSqr ) {      
#ifdef DEBUGMODE
      if( hasNanEntries(ReferenceGeometry) )
          throw BasicException("getDotProductsAndSquaredFaceAreas(): NaN error in reference geometry!");
#endif   
      
    if( ReferenceGeometry.size() != 3 * Topology.getNumVertices() )  
        throw BasicException("getDotProductsAndSquaredFaceAreas(): geometry has wrong size!");
      
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
      
    DotProducts.resize( Topology.getNumFaces() ); 
    FaceAreasSqr.resize( Topology.getNumFaces() );
    
    // run over all faces
    for ( int faceIdx = 0; faceIdx < Topology.getNumFaces(); ++faceIdx ){

          int pi( Topology.getNodeOfTriangle(faceIdx,0) ),
              pj( Topology.getNodeOfTriangle(faceIdx,1) ),
              pk( Topology.getNodeOfTriangle(faceIdx,2) );

          // set up deformed vertices and edges
          VecType Ei, Ej, Ek, temp;
          getXYZCoord<VectorType, VecType>( ReferenceGeometry, temp, pi);
          getXYZCoord<VectorType, VecType>( ReferenceGeometry, Ej, pj);
          getXYZCoord<VectorType, VecType>( ReferenceGeometry, Ek, pk);
          // e_i = e_0 = x_2 - x_1
          Ei = Ek - Ej;
          // e_j = e_1 = x_0 - x_2
          Ej = temp - Ek;
          // e_k = -e_2 = (x_2 - x_1) + (x_0 - x_2) = - (x_1 - x_0)
          //CAUTION mind the signs! (Ek is actually -Ek here!)
          Ek = Ei + Ej;
      
          // compute volume and dot products
          temp.makeCrossProduct( Ei, Ej );
          FaceAreasSqr[faceIdx] = temp.normSqr() / 4.;
          DotProducts[faceIdx][0] = dotProduct(Ej,Ek);
          DotProducts[faceIdx][1] = dotProduct(Ek,Ei);
          DotProducts[faceIdx][2] = dotProduct(Ei,Ej);
    }
    
}


#endif