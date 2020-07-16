// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Fast version of simple shell model using a fixed reference domain.
 * \author Heeren
 *
 */
 #ifndef OPTIMIZEDSHELLDEFORMATIONS_HH
#define OPTIMIZEDSHELLDEFORMATIONS_HH


//== INCLUDES =================================================================
#include <goast/Core/Auxiliary.h>
#include <goast/Core/LocalMeshGeometry.h>
#include <goast/Core/Topology.h>
#include "AuxiliaryFunctions.h"

//==========================================================================================================
// FAST VERSION OF SUM OF MEMBRANE AND BENDING ENERGY IN A SINGLE CLASS (WITH FIXED REFERENCE DOMAIN)
//==========================================================================================================

/**
 * \brief Compute sum of nonlinear membrane and simple bending energy for fixed reference domain at once.
 * \author Heeren
 *
 * Literally the summation of the classes SimpleBendingFunctional and NonlinearMembraneFunctional.
 *
 * Here we benefit from the fact that we do not have to compute certain quantities twice, which are needed for both parts of the energy.
 */
template<typename ConfiguratorType>
class FastShellFunctional : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;

    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;

    const MeshTopologySaver& _topology;
    RealType _bendWeight, _mu, _lambdaQuarter, _muHalfPlusLambdaQuarter;
    VectorType _refDihedralAngles, _refFaceAreas, _refSqrEdgeLengths;

public:
    FastShellFunctional( const MeshTopologySaver& Topology, const VectorType& refShape, RealType BendingWeight ) 
     : _topology(Topology), 
     _bendWeight(BendingWeight), 
     _mu(1.), 
     _lambdaQuarter(0.25), 
     _muHalfPlusLambdaQuarter(0.5*_mu + _lambdaQuarter){
        
        // precompute reference quantities
          VectorType dummy;
          computeReferenceQuantitiesBending<ConfiguratorType>( _topology, refShape, _refDihedralAngles, dummy, _refSqrEdgeLengths);
          getFaceAreas<ConfiguratorType>( _topology, refShape, _refFaceAreas );
    }

    //
    void apply(const VectorType &Argument, RealType &Dest) const override {
        
      if (Argument.size() != 3 * _topology.getNumVertices() )
          throw BasicException("FastElasticMeanFunctional::apply: wrong size of dofs!");
      Dest.setZero();  
      
      // precompute quantities of argument
      VectorType dihedralAngles, sqrFaceAreas;
      getDihedralAngles<ConfiguratorType>( _topology, Argument, dihedralAngles );
      std::vector<VecType> dotProducts;
      getDotProductsAndSquaredFaceAreas<ConfiguratorType>( _topology, Argument, dotProducts, sqrFaceAreas );
    
        // bending contribution
        for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){
      
          if( !(_topology.isEdgeValid(edgeIdx)) )
	    continue;

          int pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ), pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );

          // no bending at boundary edges
          if( std::min( pl, pk) < 0 )
            continue;
          
          int fl(_topology.getAdjacentTriangleOfEdge(edgeIdx,0) ), fr(_topology.getAdjacentTriangleOfEdge(edgeIdx,1) );
          RealType delTheta = dihedralAngles[edgeIdx] - _refDihedralAngles[edgeIdx];             
          Dest[0] += _bendWeight * delTheta * delTheta * _refSqrEdgeLengths[edgeIdx] / (_refFaceAreas[fl] + _refFaceAreas[fr]);
        }
               
        // membrane contribution
        int numFaces = _topology.getNumFaces();
        for ( int faceIdx = 0; faceIdx < numFaces; ++faceIdx ){
          RealType volRef( _refFaceAreas[faceIdx] );
          // trace term = -1. *  \sum_{i =0,1,2} <e_{i+1}, e_{i+2}> |\bar e_i|^2, where \bar e_i are reference edges
          // note the signs! This is since we have dotProducts[faceIdx] = { <e_1, -e_2>, <-e_2, e_0>, <e_0, e_1> }
          RealType traceTerm( dotProducts[faceIdx][0] * _refSqrEdgeLengths[ _topology.getEdgeOfTriangle(faceIdx,0)] + dotProducts[faceIdx][1] * _refSqrEdgeLengths[_topology.getEdgeOfTriangle(faceIdx,1)] - dotProducts[faceIdx][2] * _refSqrEdgeLengths[_topology.getEdgeOfTriangle(faceIdx,2)] );      
          Dest[0] += (_mu/8. *  traceTerm + _lambdaQuarter * sqrFaceAreas[faceIdx]) / volRef -  ( _muHalfPlusLambdaQuarter * std::log( sqrFaceAreas[faceIdx] / (volRef*volRef) ) + _mu + _lambdaQuarter) * volRef;
        }         
        
    }

};

//! \brief Gradient of FastShellFunctional (see above)
//! \author Heeren
template<typename ConfiguratorType>
class FastShellGradient
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;

    const MeshTopologySaver& _topology;
    RealType _bendWeight, _mu, _lambdaQuarter, _muHalfPlusLambdaQuarter;
    VectorType _refDihedralAngles, _refFaceAreas, _refSqrEdgeLengths, _refFactors;

public:
    FastShellGradient( const MeshTopologySaver& Topology, const VectorType& refShape, RealType BendingWeight ) 
    : _topology(Topology), 
    _bendWeight(BendingWeight), 
    _mu(1.), 
    _lambdaQuarter(0.25), 
    _muHalfPlusLambdaQuarter(0.5*_mu + _lambdaQuarter){        
        // precompute reference quantities
        VectorType dummy;
        computeReferenceQuantitiesBending<ConfiguratorType>( _topology, refShape, _refDihedralAngles, dummy, _refSqrEdgeLengths);
        getFaceAreasAndFactors<ConfiguratorType>( _topology, refShape, _refFaceAreas, _refFactors );
    }


    //
    void apply(const VectorType &Argument, VectorType &Dest) const override {
        
      int numV = _topology.getNumVertices();  
      if (Argument.size() != 3 * numV )
          throw BasicException("FastShellGradient::apply: wrong size of dofs!");
      if (Dest.size() != 3 * numV )
          Dest.resize(  3 * numV  );
      Dest.setZero();
        
      // precompute quantities of argument
      std::vector<VecType> gradPi, gradPj, gradPk, gradPl;
      VectorType dihedralAngles;
      getDihedralAnglesGradients<ConfiguratorType>( _topology, Argument, dihedralAngles, gradPi, gradPj, gradPk, gradPl );
      
      // bending contribution         
      for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){
      
          if( !(_topology.isEdgeValid(edgeIdx)) )
	    continue;

	  int pi(_topology.getAdjacentNodeOfEdge(edgeIdx, 0)),
              pj(_topology.getAdjacentNodeOfEdge(edgeIdx, 1)),
	      pk(_topology.getOppositeNodeOfEdge(edgeIdx, 0)),
	      pl(_topology.getOppositeNodeOfEdge(edgeIdx, 1));

          // no bending at boundary edges
          if( std::min( pl, pk) < 0 )
            continue;
          
          // factor
          int fl(_topology.getAdjacentTriangleOfEdge(edgeIdx,0) ), fr(_topology.getAdjacentTriangleOfEdge(edgeIdx,1) );         
          RealType delTheta( -2. *  _refSqrEdgeLengths[edgeIdx] * (_refDihedralAngles[edgeIdx] - dihedralAngles[edgeIdx] ) / (_refFaceAreas[fl] + _refFaceAreas[fr]) );
          
          // assemble in global vector
	  for (int i = 0; i < 3; i++){
            Dest[i*numV + pi] += _bendWeight * delTheta * gradPi[edgeIdx][i];
            Dest[i*numV + pj] += _bendWeight * delTheta * gradPj[edgeIdx][i];
            Dest[i*numV + pk] += _bendWeight * delTheta * gradPk[edgeIdx][i];
            Dest[i*numV + pl] += _bendWeight * delTheta * gradPl[edgeIdx][i];
	  }
     }

     // get face areas and area gradients
     VectorType faceAreas;
     getAreaGradients<ConfiguratorType>( _topology, Argument, faceAreas, gradPi, gradPj, gradPk, gradPl );
    
     // membrane contribution     
     for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); ++faceIdx ){
       std::vector<int> nodesIdx(3);
       for (int j = 0; j < 3; j++)
         nodesIdx[j] = _topology.getNodeOfTriangle(faceIdx, j);

       RealType factor = 2. * (_lambdaQuarter * faceAreas[faceIdx] / _refFaceAreas[faceIdx] - _muHalfPlusLambdaQuarter * _refFaceAreas[faceIdx] / faceAreas[faceIdx]);
       for (int j = 0; j < 3; j++){
              
           // det part
           Dest[j*numV + nodesIdx[0]] += factor * gradPi[faceIdx][j];
           Dest[j*numV + nodesIdx[1]] += factor * gradPj[faceIdx][j];
           Dest[j*numV + nodesIdx[2]] += factor * gradPk[faceIdx][j];
              
           // trace part
           for (int i = 0; i < 3; i++){
             int nextIdx = (i + 1) % 3;
             int prevIdx = (i + 2) % 3;
             Dest[j*numV + nodesIdx[i]] += 0.25 * _mu * (_refFactors[3*faceIdx + prevIdx] * gradPl[3*faceIdx + prevIdx][j] - _refFactors[3*faceIdx + nextIdx] * gradPl[3*faceIdx + nextIdx][j]);
           }
       }                 
     }  
        
    }

};

#endif