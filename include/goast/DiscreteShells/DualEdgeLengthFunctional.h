// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Dual edge length energy (for fixed reference domain) as a simple and non-consistent bending approximation, which, however, is quadratic in vertex positions.
 * \author Heeren
 *
 */
 #ifndef __DUALEDGELENGTHFUNCTIONAL_H
#define __DUALEDGELENGTHFUNCTIONAL_H

#include <goast/Core/Auxiliary.h>
#include <goast/Core/LocalMeshGeometry.h>
#include <goast/Core/Topology.h>

//==========================================================================================================
// DUAL EDGE LENGTH ENERGY (FOR FIXED REFERENCE DOMAIN)
//==========================================================================================================

/**
 * \brief Dual edge length energy for fixed reference domain.
 * \author Heeren
 *
 * If \f$ e = (i,j) \f$ is an edge with adjacent faces \f$ (i,j,k) \f$ and \f$ (j,i,l) \f$, respectively,
 * then we denote \f$ (k,l) \f$ the dual edge of \f$ e = (i,j) \f$.
 *
 * For a closed mesh \f$ x = (x_1, \ldots, x_n) \f$, this class implements \f[ F[x] = \sum_{e \in x} w_e (\|x_k - x_l\| - d_e )^2 \, ,\f]
 * where \f$ w_e \f$ is an integration weight and \f$ d_e \f$ some dual target edge length for \f$ e \f$.
 */
template< typename ConfiguratorType>
class DualEdgeLengthFunctional : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;

  const MeshTopologySaver& _topology;
  VectorType  _targetLengths, _weights;

public:
  DualEdgeLengthFunctional( const MeshTopologySaver& topology, const VectorType& targetLengths, const VectorType& Weights ) 
  : _topology( topology), 
    _targetLengths( targetLengths ), 
    _weights( Weights ) {}
    
  DualEdgeLengthFunctional( const MeshTopologySaver& topology, const VectorType& referenceGeometry ) 
  : _topology( topology){
        computeIntegrationWeightsLengthFunctional<ConfiguratorType>( _topology, referenceGeometry, _weights );
        getDualEdgeLengths<ConfiguratorType>( _topology, referenceGeometry, _targetLengths );
    }

  // energy evaluation
  void apply( const VectorType& defGeometry, RealType & Dest ) const {
    if (3 * _topology.getNumVertices() != defGeometry.size() )
      throw BasicException("DualEdgeLengthFunctional::apply(): sizes dont match!");
    Dest = 0.;
    
    for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){
      int pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ),
          pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );
      if( std::min( pl, pk) < 0 )
        continue;

      VecType Pk, Pl;
      getXYZCoord<VectorType, VecType>( defGeometry, Pk, pk); 
      getXYZCoord<VectorType, VecType>( defGeometry, Pl, pl );
      RealType dualEdgeLength = std::sqrt( dotProduct(Pk-Pl,Pk-Pl) ); 
      Dest += _weights[edgeIdx] * (dualEdgeLength - _targetLengths[edgeIdx]) * (dualEdgeLength - _targetLengths[edgeIdx]);
    }
  }

};


//==========================================================================================================
//! \brief First derivative of DualEdgeLengthFunctional (for fixed reference geometry).
//! \author Heeren
template< typename ConfiguratorType>
class DualEdgeLengthGradient : public BaseOp<typename ConfiguratorType::VectorType> {

  typedef typename ConfiguratorType::RealType   RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;

  const MeshTopologySaver& _topology;
  VectorType _targetLengths, _weights;

public:
DualEdgeLengthGradient( const MeshTopologySaver& topology, const VectorType& targetLengths, const VectorType& Weights ) 
  : _topology( topology), 
    _targetLengths( targetLengths ), 
    _weights( Weights ) {}
    
DualEdgeLengthGradient( const MeshTopologySaver& topology, const VectorType& referenceGeometry ) 
  : _topology( topology){
        computeIntegrationWeightsLengthFunctional<ConfiguratorType>( _topology, referenceGeometry, _weights );
        getDualEdgeLengths<ConfiguratorType>( _topology, referenceGeometry, _targetLengths );
    }

void apply( const VectorType& defGeometry, VectorType& Dest ) const {

  if (3 * _topology.getNumVertices() != defGeometry.size() )
    throw BasicException("DualEdgeLengthGradient::apply(): sizes dont match!");

  if (Dest.size() != defGeometry.size())
    Dest.resize(defGeometry.size());
  Dest.setZero();
  
  for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){

        int pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ),
            pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );
        if( std::min( pl, pk) < 0 )
        continue;

      // set up vertices and edge        
      VecType Pk, Pl;      
      getXYZCoord<VectorType, VecType>( defGeometry, Pk, pk);
      getXYZCoord<VectorType, VecType>( defGeometry, Pl, pl);
      VecType edge = Pl-Pk;
      RealType dualEdgeLength = std::sqrt( dotProduct(edge, edge) );       
      RealType factor = 2. * _weights[edgeIdx] * (dualEdgeLength - _targetLengths[edgeIdx]);

      // assemble in global matrix
      for( int i = 0; i < 3; i++ ){
        Dest[i * _topology.getNumVertices() + pk] -= factor * edge[i] / dualEdgeLength;
        Dest[i * _topology.getNumVertices() + pl] += factor * edge[i] / dualEdgeLength;
      }

    }
  }

};

//==========================================================================================================
//! \brief Second derivative of DualEdgeLengthFunctional (for fixed reference domain).
//! \author Heeren
template <typename ConfiguratorType>
class DualEdgeLengthHessian : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

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
DualEdgeLengthHessian( const MeshTopologySaver& topology,
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
    
DualEdgeLengthHessian( const MeshTopologySaver& topology, 
                       const VectorType& referenceGeometry,
                       const RealType Factor = 1.,
                       int rowOffset = 0,
                       int colOffset = 0 ) 
  : _topology( topology),
    _factor( Factor ),
    _rowOffset(rowOffset), 
    _colOffset(colOffset){
        computeIntegrationWeightsLengthFunctional<ConfiguratorType>( _topology, referenceGeometry, _weights );
        getDualEdgeLengths<ConfiguratorType>( _topology, referenceGeometry, _targetLengths );
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

      int pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ),
          pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );
      if( std::min( pl, pk) < 0 )
        continue;

      // set up vertices and edges
      VecType Pk, Pl, edge;
      getXYZCoord<VectorType, VecType>( defGeometry, Pk, pk); 
      getXYZCoord<VectorType, VecType>( defGeometry, Pl, pl);
      edge = Pl-Pk;   
      RealType dualEdgeLength = std::sqrt( dotProduct(edge, edge) );
      edge /= dualEdgeLength;      

      // now compute second derivatives of dihedral angle
      MatType tensorProduct;
      tensorProduct.makeTensorProduct( edge, edge );
      tensorProduct *= 2. * _weights[edgeIdx] * ( 1. - (dualEdgeLength - _targetLengths[edgeIdx]) / dualEdgeLength );
      
      //ii  
      tensorProduct.addToDiagonal( 2. *  _weights[edgeIdx] * (dualEdgeLength - _targetLengths[edgeIdx]) / dualEdgeLength );
      localToGlobal( tripletList, pk, pk, tensorProduct );     
      localToGlobal( tripletList, pl, pl, tensorProduct ); 

      //ij & ji (Hij = Hji)
      tensorProduct *= -1.;
      localToGlobal( tripletList, pk, pl, tensorProduct );   
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
