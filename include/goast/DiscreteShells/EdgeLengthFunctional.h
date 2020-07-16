// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Edge length energy (for fixed reference domain) as a very simple membrane model.
 * \author Heeren
 *
 */
#ifndef __EDGELENGTHFUNCTIONAL_H
#define __EDGELENGTHFUNCTIONAL_H

#include <goast/Core/Auxiliary.h>
#include <goast/Core/LocalMeshGeometry.h>
#include <goast/Core/Topology.h>

//==========================================================================================================
// DISCRETE SHELLS EDGE LENGTH ENERGY (FOR FIXED REFERENCE DOMAIN)
//==========================================================================================================

/**
 * \brief Edge length functional (for fixed reference geometry).
 * \author Heeren
 *
 * Run over all edges \f$ e = (i,j) \f$ and compute \f$ w_e (\|x_i - x_j\| - d_e )^2 \f$ per edge.
 * where \f$ w_e \f$ is an integration weight and \f$ d_e \f$ some target edge length per edge, i.e.
 * \f[ E[x] = \sum_e w_e (l_e[x] - d_e )^2\, .\f]
 */
template< typename ConfiguratorType>
class EdgeLengthFunctional : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;

  const MeshTopologySaver& _topology;
  VectorType _targetLengths, _weights;

public:
  EdgeLengthFunctional( const MeshTopologySaver& topology, const VectorType& targetLengths, const VectorType& Weights ) 
  : _topology( topology), 
    _targetLengths( targetLengths ), 
    _weights( Weights ) {}
    
  EdgeLengthFunctional( const MeshTopologySaver& topology, const VectorType& referenceGeometry ) 
  : _topology( topology){
        computeIntegrationWeightsLengthFunctional<ConfiguratorType>( _topology, referenceGeometry, _weights );
        getEdgeLengths<ConfiguratorType>( _topology, referenceGeometry, _targetLengths );
    }

  // energy evaluation
  void apply( const VectorType& defGeometry, RealType & Dest ) const {
    if (3 * _topology.getNumVertices() != defGeometry.size() )
      throw BasicException("EdgeLengthFunctional::apply(): sizes dont match!");
    Dest = 0.;
    
    for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){
      int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) );

      VecType Pi, Pj;
      getXYZCoord<VectorType, VecType>( defGeometry, Pi, pi); 
      getXYZCoord<VectorType, VecType>( defGeometry, Pj, pj );
      RealType defEdgeLength = std::sqrt( dotProduct(Pj-Pi,Pj-Pi) ); 
      Dest += _weights[edgeIdx] * (defEdgeLength - _targetLengths[edgeIdx]) * (defEdgeLength - _targetLengths[edgeIdx]);
    }
  }

};

//==========================================================================================================
//! \brief First derivative of EdgeLengthFunctional (for fixed reference geometry).
//! \author Heeren
template< typename ConfiguratorType>
class EdgeLengthGradient : public BaseOp<typename ConfiguratorType::VectorType> {

  typedef typename ConfiguratorType::RealType   RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;

  const MeshTopologySaver& _topology;
  VectorType  _targetLengths, _weights;

public:
EdgeLengthGradient( const MeshTopologySaver& topology, const VectorType& targetLengths, const VectorType& Weights ) 
  : _topology( topology), 
    _targetLengths( targetLengths ), 
    _weights( Weights ) {}
    
EdgeLengthGradient( const MeshTopologySaver& topology, const VectorType& referenceGeometry ) 
  : _topology( topology){
        computeIntegrationWeightsLengthFunctional<ConfiguratorType>( _topology, referenceGeometry, _weights );
        getEdgeLengths<ConfiguratorType>( _topology, referenceGeometry, _targetLengths );
    }

void apply( const VectorType& defGeometry, VectorType& Dest ) const {

  if (3 * _topology.getNumVertices() != defGeometry.size() )
    throw BasicException("EdgeLengthGradient::apply(): sizes dont match!");

  if (Dest.size() != defGeometry.size())
    Dest.resize(defGeometry.size());
  Dest.setZero();
  
  for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){

        int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
            pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) );	    

      // set up vertices and edge        
      VecType Pi, Pj;      
      getXYZCoord<VectorType, VecType>( defGeometry, Pi, pi);
      getXYZCoord<VectorType, VecType>( defGeometry, Pj, pj);
      VecType edge = Pj-Pi;
      RealType defEdgeLength = std::sqrt( dotProduct(edge, edge) );       
      RealType factor = 2. * _weights[edgeIdx] * (defEdgeLength - _targetLengths[edgeIdx]);

      // assemble in global matrix
      for( int i = 0; i < 3; i++ ){
        Dest[i * _topology.getNumVertices() + pi] -= factor * edge[i] / defEdgeLength;
        Dest[i * _topology.getNumVertices() + pj] += factor * edge[i] / defEdgeLength;
      }

    }
  }

};

//==========================================================================================================
//! \brief Second derivative of EdgeLengthFunctional (for fixed reference domain).
//! \author Heeren
template <typename ConfiguratorType>
class EdgeLengthHessian : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::VecType     VecType;
  typedef typename ConfiguratorType::MatType     MatType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver& _topology;
  VectorType _targetLengths, _weights;
  RealType _factor;
  mutable int _rowOffset, _colOffset;

public:
EdgeLengthHessian( const MeshTopologySaver& topology,
                   const VectorType& targetLengths, 
                   const VectorType& Weights,
		   const RealType Factor = 1.,
                   int rowOffset = 0,
                   int colOffset = 0 ) 
  : _topology( topology), 
    _targetLengths( targetLengths ), 
    _weights( Weights ), 
    _factor( Factor ),
    _rowOffset(rowOffset), 
    _colOffset(colOffset){}
    
EdgeLengthHessian( const MeshTopologySaver& topology,
                   const VectorType& referenceGeometry,
                   const RealType Factor = 1.,
                   int rowOffset = 0,
                   int colOffset = 0 ) 
  : _topology( topology),
    _factor( Factor ),
    _rowOffset(rowOffset), 
    _colOffset(colOffset){
        computeIntegrationWeightsLengthFunctional<ConfiguratorType>( _topology, referenceGeometry, _weights );
        getEdgeLengths<ConfiguratorType>( _topology, referenceGeometry, _targetLengths );
    }
    
    void setRowOffset( int rowOffset ) const {
        _rowOffset = rowOffset;
    }
    
    void setColOffset( int colOffset ) const {
        _colOffset = colOffset;
    }

  //
  void apply( const VectorType& defGeometry, MatrixType& Dest ) const {    
    assembleHessian( defGeometry, Dest );
  }

  // assmeble Hessian matrix via triplet list
  void assembleHessian( const VectorType& defGeometry, MatrixType& Hessian ) const {
    int dofs = 3*_topology.getNumVertices();
    if( (Hessian.rows() != dofs) || (Hessian.cols() != dofs) )
        Hessian.resize( dofs, dofs );
    Hessian.setZero();
    
    // set up triplet list
    TripletListType tripletList;
    // per edge we have 2 active vertices, i.e. 4 combinations each producing a 3x3-matrix
    tripletList.reserve( 4 * 9 * _topology.getNumEdges() );
    
    pushTriplets( defGeometry, tripletList );
    
    // fill matrix from triplets
    Hessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }

  //
  void pushTriplets( const VectorType& defGeometry, TripletListType& tripletList ) const  {
      
  if (3 * _topology.getNumVertices() != defGeometry.size() )
    throw BasicException("EdgeLengthHessian::apply(): sizes dont match!");

    for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){

      int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) );

      // set up vertices and edges
      VecType Pi, Pj, P, edge;
      getXYZCoord<VectorType, VecType>( defGeometry, Pi, pi); 
      getXYZCoord<VectorType, VecType>( defGeometry, Pj, pj);
      edge = Pj-Pi;
      RealType defEdgeLengthSqr = dotProduct(edge, edge);    
      RealType defEdgeLength = std::sqrt( defEdgeLengthSqr );
      edge /= defEdgeLength;      

      // now compute second derivatives of dihedral angle
      MatType tensorProduct;
      tensorProduct.makeTensorProduct( edge, edge );
      tensorProduct *= 2. * _weights[edgeIdx] * ( 1. - (defEdgeLength - _targetLengths[edgeIdx]) / defEdgeLength );
      
      //ii  
      tensorProduct.addToDiagonal( 2. *  _weights[edgeIdx] * (defEdgeLength - _targetLengths[edgeIdx]) / defEdgeLength );
      localToGlobal( tripletList, pi, pi, tensorProduct );     
      localToGlobal( tripletList, pj, pj, tensorProduct ); 

      //ij & ji (Hij = Hji)
      tensorProduct *= -1.;
      localToGlobal( tripletList, pi, pj, tensorProduct );   
    }
  }

protected:
  void localToGlobal( TripletListType& tripletList, int k, int l, const MatType& localMatrix ) const {
    int numV = _topology.getNumVertices();
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          tripletList.push_back( TripletType( _rowOffset + i*numV + k, _colOffset + j*numV + l, _factor * localMatrix(i,j) ) );	
	
    if( k != l){
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          tripletList.push_back( TripletType( _rowOffset + i*numV + l, _colOffset + j*numV + k, _factor * localMatrix(j,i) ) );
    }
  }
  
};

#endif