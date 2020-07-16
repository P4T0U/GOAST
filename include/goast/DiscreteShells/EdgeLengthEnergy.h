// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Edge length energy as a very simple membrane model.
 * \author Heeren
 *
 */
 #ifndef __EDGELENGTHENERGY_H
#define __EDGELENGTHENERGY_H

#include <goast/Core/Auxiliary.h>
#include <goast/Core/LocalMeshGeometry.h>
#include <goast/Core/Topology.h>

//==========================================================================================================
// DISCRETE SHELLS EDGE LENGTH ENERGY
//==========================================================================================================

/**
 * \brief Edge length energy between two discrete shells given as triangular meshes.
 * \author Heeren
 *
 * This energy is a modified version of the membrane model introduced in Grinspun et al. \cite GrHiDeSc03
 *
 * Let \f$ x, \tilde x \f$ two meshes (resp. their geometries) in dense correspondence.
 * Then this class represents the membrane energy \f[ E[x, \tilde x] = \sum_e \frac{ (l_e - \tilde l_e)^2 }{ l_e^2 } d_e \, ,\f]
 * where the sum is over all edges \f$ e \in x \f$, where \f$ l_e \f$ is the edge length of \f$ e \f$
 * and \f$ d_e = |T_1| + |T_2| \f$, if \f$ T_1, T_2 \f$ are the triangles sharing edge \f$ e \f$.
 *
 * Note that the energy might either be thought of as \f$ x \mapsto E[x, \tilde x] \f$ (active shell is undeformed shell)
 * or \f$ \tilde x \mapsto E[x, \tilde x] \f$ (active shell is deformed shell).
 * The active shell is considered the argument whereas the inactive shell is given in the constructor.
 */
template< typename ConfiguratorType>
class EdgeLengthEnergy : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;

  const MeshTopologySaver& _topology;
  const VectorType&  _inactiveGeometry;
  const bool _activeShellIsDeformed;

public:

  EdgeLengthEnergy( const MeshTopologySaver& topology,
                    const VectorType& InactiveGeometry,
		    const bool ActiveShellIsDeformed ) 
  : _topology( topology), 
    _inactiveGeometry(InactiveGeometry), 
    _activeShellIsDeformed( ActiveShellIsDeformed ) {}

  // energy evaluation
  void apply( const VectorType& ActiveGeometry, RealType & Dest ) const {

    const VectorType* defShellP   = _activeShellIsDeformed ? &ActiveGeometry : &_inactiveGeometry;
    const VectorType* undefShellP = _activeShellIsDeformed ? &_inactiveGeometry : &ActiveGeometry;
    Dest.setZero();
    
    for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){

      int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
          pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ),
          pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );

      // set up vertices and edges
      VecType Pi, Pj, P, temp;

      getXYZCoord<VectorType, VecType>( *defShellP, Pi, pi); 
      getXYZCoord<VectorType, VecType>( *defShellP, Pj, pj );
      RealType defEdgeLength = std::sqrt( dotProduct(Pj-Pi,Pj-Pi) ); 
      
      getXYZCoord<VectorType, VecType>( *undefShellP, Pi, pi); 
      getXYZCoord<VectorType, VecType>( *undefShellP, Pj, pj );
      RealType undefEdgeLength = std::sqrt( dotProduct(Pj-Pi,Pj-Pi) );
      
      // compute volume, length of edge and theta difference
      RealType vol = 0.;
      if( pk != -1 ){
	getXYZCoord<VectorType, VecType>( *undefShellP, P, pk);
        temp.makeCrossProduct( P-Pj, Pi-P );
        vol += temp.norm() / 6.;
      }
      if( pl != -1 ){
	getXYZCoord<VectorType, VecType>( *undefShellP, P, pl);
        temp.makeCrossProduct( P-Pi, Pj-P );
        vol += temp.norm() / 6.;
      }
      
      RealType diff = ( undefEdgeLength - defEdgeLength ) / undefEdgeLength;
      Dest[0] += vol * diff * diff;

    }
  }

 };

 
//==========================================================================================================
//! \brief First derivative of EdgeLengthEnergy w.r.t. the deformed configuration
//! \author Heeren
template< typename ConfiguratorType>
class EdgeLengthGradientDef : public BaseOp<typename ConfiguratorType::VectorType> {

  typedef typename ConfiguratorType::RealType   RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;

  const MeshTopologySaver& _topology;
  const VectorType&  _undefShell;

public:
EdgeLengthGradientDef( const MeshTopologySaver& topology, const VectorType& undefShell ) : _topology( topology), _undefShell(undefShell) {}

void apply( const VectorType& defShell, VectorType& Dest ) const {

  if (_undefShell.size() != defShell.size()){
			std::cerr << "size of undef = " << _undefShell.size() << " vs. size of def = " << defShell.size() << std::endl;
			throw BasicException("EdgeLengthGradientDef::apply(): sizes dont match!");
  }

  if (Dest.size() != defShell.size())
    Dest.resize(defShell.size());
  Dest.setZero();
  
  for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){

        int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
            pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
            pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ),
            pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );	    

      // set up vertices and edges
      VecType Pi, Pj, P, edge;      
      getXYZCoord<VectorType, VecType>( _undefShell, Pi, pi); 
      getXYZCoord<VectorType, VecType>( _undefShell, Pj, pj); 
      RealType undefEdgeLength = std::sqrt( dotProduct(Pj-Pi, Pj-Pi) );
      
      // compute volume, length of edge and theta difference
      RealType vol = 0.;
      if( pk != -1 ){
	getXYZCoord<VectorType, VecType>( _undefShell, P, pk); 
        edge.makeCrossProduct( P-Pj, Pi-P );
        vol += edge.norm() / 6.;
      }
      if( pl != -1 ){
	getXYZCoord<VectorType, VecType>( _undefShell, P, pl);
        edge.makeCrossProduct( P-Pi, Pj-P );
        vol += edge.norm() / 6.;
      }      
            
      getXYZCoord<VectorType, VecType>( defShell, Pi, pi);
      getXYZCoord<VectorType, VecType>( defShell, Pj, pj);
      edge = Pj-Pi;
      RealType defEdgeLength = std::sqrt( dotProduct(edge, edge) );       
      
      vol *= 2. * (defEdgeLength - undefEdgeLength) / ( undefEdgeLength * undefEdgeLength );

      // assemble in global matrix
      for( int i = 0; i < 3; i++ ){
        Dest[i * _topology.getNumVertices() + pi] -= vol * edge[i] / defEdgeLength;
        Dest[i * _topology.getNumVertices() + pj] += vol * edge[i] / defEdgeLength;
      }

    }
  }

};


//==========================================================================================================
//! \brief Second derivative of EdgeLengthEnergy w.r.t. the deformed configuration
//! \author Heeren
template <typename ConfiguratorType>
class EdgeLengthHessianDef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::VecType     VecType;
  typedef typename ConfiguratorType::MatType     MatType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver& _topology;
  const VectorType& _undefShell;
  RealType _factor;
  mutable int _rowOffset, _colOffset;

public:
EdgeLengthHessianDef( const MeshTopologySaver& topology,
                             const VectorType& undefShell,
			     const RealType Factor = 1.,
                             int rowOffset = 0,
                             int colOffset = 0 ) 
  : _topology( topology), 
    _undefShell(undefShell), 
    _factor( Factor ),
    _rowOffset(rowOffset), 
    _colOffset(colOffset){}
    
    void setRowOffset( int rowOffset ) const {
        _rowOffset = rowOffset;
    }
    
    void setColOffset( int colOffset ) const {
        _colOffset = colOffset;
    }

  //
  void apply( const VectorType& defShell, MatrixType& Dest ) const {    
    assembleHessian( defShell, Dest );
  }

  // assmeble Hessian matrix via triplet list
  void assembleHessian( const VectorType& defShell, MatrixType& Hessian ) const {
    int dofs = 3*_topology.getNumVertices();
    if( (Hessian.rows() != dofs) || (Hessian.cols() != dofs) )
        Hessian.resize( dofs, dofs );
    Hessian.setZero();
    
    // set up triplet list
    TripletListType tripletList;
    // per edge we have 2 active vertices, i.e. 4 combinations each producing a 3x3-matrix
    tripletList.reserve( 4 * 9 * _topology.getNumEdges() );
    
    pushTriplets( defShell, tripletList );
    
    // fill matrix from triplets
    Hessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }

  //
  void pushTriplets( const VectorType& defShell, TripletListType& tripletList ) const  {
      
    if( _undefShell.size() != defShell.size() ){
      std::cerr << "size of undef = " << _undefShell.size() << " vs. size of def = " << defShell.size() << std::endl;
      throw BasicException( "EdgeLengthHessianDef::pushTriplets(): sizes dont match!");
    }

    for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){

      int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
          pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ),
          pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );


      // set up vertices and edges
      VecType Pi, Pj, P, edge;
      getXYZCoord<VectorType, VecType>( _undefShell, Pi, pi); 
      getXYZCoord<VectorType, VecType>( _undefShell, Pj, pj);
      RealType undefEdgeLengthSqr = dotProduct(Pj-Pi, Pj-Pi);
      RealType undefEdgeLength = std::sqrt( undefEdgeLengthSqr );
      
      // compute volume, length of edge and theta difference
      RealType vol = 0.;
      if( pk != -1 ){
	getXYZCoord<VectorType, VecType>( _undefShell, P, pk);
        edge.makeCrossProduct( P-Pj, Pi-P );
        vol += edge.norm() / 6.;
      }
      if( pl != -1 ){
	getXYZCoord<VectorType, VecType>( _undefShell, P, pl);
        edge.makeCrossProduct( P-Pi, Pj-P );
        vol += edge.norm() / 6.;
      }      
            
      // deformed quantities
      getXYZCoord<VectorType, VecType>( defShell, Pi, pi); 
      getXYZCoord<VectorType, VecType>( defShell, Pj, pj);
      edge = Pj-Pi;
      RealType defEdgeLengthSqr = dotProduct(edge, edge);    
      RealType defEdgeLength = std::sqrt( defEdgeLengthSqr );
      

      // now compute second derivatives of dihedral angle
      MatType tensorProduct;
      tensorProduct.makeTensorProduct( edge, edge );
      tensorProduct *= 2. * vol / ( defEdgeLengthSqr * defEdgeLength * undefEdgeLength );
      
      //ii  
      tensorProduct.addToDiagonal( 2. * vol * (defEdgeLength - undefEdgeLength) / ( undefEdgeLengthSqr * defEdgeLength ) );
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


//==========================================================================================================
//! \brief First derivative of EdgeLengthEnergy w.r.t. the undeformed configuration
//! \author Heeren
template< typename ConfiguratorType>
class EdgeLengthGradientUndef :  public BaseOp<typename ConfiguratorType::VectorType> {

  typedef typename ConfiguratorType::RealType   RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;

  const MeshTopologySaver& _topology;
  const VectorType&  _defShell;

public:
EdgeLengthGradientUndef( const MeshTopologySaver& topology, const VectorType& defShell ) : _topology( topology), _defShell(defShell) {}

void apply( const VectorType& undefShell, VectorType& Dest ) const {

  if (_defShell.size() != undefShell.size()){
			std::cerr << "size of undef = " << undefShell.size() << " vs. size of def = " << _defShell.size() << std::endl;
			throw BasicException("EdgeLengthGradientUndef::apply(): sizes dont match!");
  }

  if (Dest.size() != undefShell.size())
    Dest.resize(undefShell.size());
  Dest.setZero();

  for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){

        int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
            pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
            pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ),
            pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );

      // first get deformed quantities
      VecType Pi, Pj, Pk, Pl, temp;
      getXYZCoord<VectorType, VecType>( _defShell, Pi, pi); 
      getXYZCoord<VectorType, VecType>( _defShell, Pj, pj); 
      RealType defEdgeLengthSqr( dotProduct(Pj-Pi,Pj-Pi) );
      RealType defEdgeLength = std::sqrt( defEdgeLengthSqr );
      
      //!now get the undeformed values
       getXYZCoord<VectorType, VecType>( undefShell, Pi, pi); 
      getXYZCoord<VectorType, VecType>( undefShell, Pj, pj); 
  
      // compute Ak + Al, |e|^2 and diff. o dihedral angles
     RealType vol = 0.;
      if( pk != -1 ){
	getXYZCoord<VectorType, VecType>( undefShell, Pk, pk);
        temp.makeCrossProduct( Pk-Pj, Pi-Pk );
        vol += temp.norm() / 6.;
      }
      if( pl != -1 ){
	getXYZCoord<VectorType, VecType>( undefShell, Pl, pl);
        temp.makeCrossProduct( Pl-Pi, Pj-Pl );
        vol += temp.norm() / 6.;
      } 
      
      VecType edge( Pj-Pi );
      RealType undefEdgeLengthSqr( dotProduct(edge,edge) ); 
      RealType undefEdgeLength = std::sqrt( undefEdgeLengthSqr );
      //edge /= undefEdgeLength;
      
      // derivatives    
      VecType gradk, gradl, gradi, gradj, gradTheta, gradArea;  
      RealType diff = undefEdgeLength - defEdgeLength;
      RealType factorGradArea = diff * diff / (3. * undefEdgeLengthSqr);
      RealType factorGradEdgeLengthSqr =  2. * vol * diff * defEdgeLength / ( undefEdgeLengthSqr * undefEdgeLengthSqr );   
         
      // d_k
      if( pk != -1 ){
        getAreaGradK( Pi, Pj, Pk, gradk );
        gradk *= factorGradArea;
      }
      
      // d_l
      if( pl != -1 ){
        getAreaGradK( Pj, Pi, Pl, gradl );
        gradl *= factorGradArea;
      }
      
      // d_i
      if( pk != -1 ){
        getAreaGradK( Pj, Pk, Pi, gradArea );
        gradi.addMultiple( gradArea, factorGradArea );
      }
      if( pl != -1 ){
        getAreaGradK( Pl, Pj, Pi, gradArea );
        gradi.addMultiple( gradArea, factorGradArea );
      }
      gradi.addMultiple( edge, -1. * factorGradEdgeLengthSqr );
      
      // d_j
      if( pk != -1 ){
        getAreaGradK( Pk, Pi, Pj, gradArea );
        gradj.addMultiple( gradArea, factorGradArea );
      }
      if( pl != -1 ){
        getAreaGradK( Pi, Pl, Pj, gradArea );
        gradj.addMultiple( gradArea, factorGradArea );
      }
      gradj.addMultiple( edge, factorGradEdgeLengthSqr );      
      
      // assemble in global matrix
      for( int i = 0; i < 3; i++ ){
        Dest[i*_topology.getNumVertices() + pi] += gradi[i];
        Dest[i*_topology.getNumVertices() + pj] += gradj[i];
      }
      
      if( pl != -1 )
	for( int i = 0; i < 3; i++ )
          Dest[i*_topology.getNumVertices() + pl] += gradl[i];
	
      if( pk != -1 )
	for( int i = 0; i < 3; i++ )
          Dest[i*_topology.getNumVertices() + pk] += gradk[i];
  }
}

};


//==========================================================================================================
//! \brief Second derivative of EdgeLengthEnergy w.r.t. the undeformed configuration
//! \author Heeren
template <typename ConfiguratorType >
class EdgeLengthHessianUndef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::VecType     VecType;
  typedef typename ConfiguratorType::MatType     MatType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver& _topology;
  const VectorType& _defShell;
  RealType _factor;
  mutable int _rowOffset, _colOffset;

public:
EdgeLengthHessianUndef( const MeshTopologySaver& topology,
                             const VectorType& defShell,
			     const RealType Factor = 1.,
                             int rowOffset = 0,
                             int colOffset = 0 ) 
  : _topology( topology), 
    _defShell(defShell), 
    _factor( Factor ),
    _rowOffset(rowOffset), 
    _colOffset(colOffset){}
    
    void setRowOffset( int rowOffset ) const {
        _rowOffset = rowOffset;
    }
    
    void setColOffset( int colOffset ) const {
        _colOffset = colOffset;
    }

  //
  void apply( const VectorType& undefShell, MatrixType& Dest ) const {    
    assembleHessian( undefShell, Dest );
  }

  // assmeble Hessian matrix via triplet list
  void assembleHessian( const VectorType& undefShell, MatrixType& Hessian ) const {
    int dofs = 3*_topology.getNumVertices();
    if( (Hessian.rows() != dofs) || (Hessian.cols() != dofs) )
        Hessian.resize( dofs, dofs );
    Hessian.setZero();
    
    // set up triplet list
    TripletListType tripletList;
    // per edge we have 2 active vertices, i.e. 4 combinations each producing a 3x3-matrix
    tripletList.reserve( 4 * 9 * _topology.getNumEdges() );
    
    pushTriplets( undefShell, tripletList );
    
    // fill matrix from triplets
    Hessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }

  //
  void pushTriplets( const VectorType& undefShell, TripletListType& tripletList ) const  {
      
    if( _defShell.size() != undefShell.size() ){
      std::cerr << "size of undef = " << _defShell.size() << " vs. size of def = " << undefShell.size() << std::endl;
      throw BasicException( "EdgeLengthHessianUndef::pushTriplets(): sizes dont match!");
    }

    for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){

      int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
          pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ),
          pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );

      // set up vertices and edges
      VecType Pi, Pj, Pk, Pl, temp;

      // first get defomed quantities
      getXYZCoord<VectorType, VecType>( _defShell, Pi, pi); 
      getXYZCoord<VectorType, VecType>( _defShell, Pj, pj);
      RealType defEdgeLengthSqr( dotProduct(Pj-Pi,Pj-Pi) );
      RealType defEdgeLength = std::sqrt( defEdgeLengthSqr );
      
      // get undeformed vertex positions
      getXYZCoord<VectorType, VecType>( undefShell, Pi, pi); 
      getXYZCoord<VectorType, VecType>( undefShell, Pj, pj);  
            
      // compute Ak + Al, |e|^2 and dihedral angles
      RealType vol = 0.;
      if( pk != -1 ){
	getXYZCoord<VectorType, VecType>( undefShell, Pk, pk); 
        temp.makeCrossProduct( Pk-Pj, Pi-Pk );
        vol += temp.norm() / 6.;
      }
      if( pl != -1 ){
	getXYZCoord<VectorType, VecType>( undefShell, Pl, pl); 
        temp.makeCrossProduct( Pl-Pi, Pj-Pl );
        vol += temp.norm() / 6.;
      } 
     
      // compute first derivatives of area
      VecType areak, areal, areai, areaj;
      if( pk != -1 ){
        getAreaGradK( Pi, Pj, Pk, areak );
        getAreaGradK( Pj, Pk, Pi, areai );
	getAreaGradK( Pk, Pi, Pj, areaj );
      }
      if( pl != -1 ){
        getAreaGradK( Pj, Pi, Pl, areal );
        getAreaGradK( Pl, Pj, Pi, temp );
        areai += temp;
        getAreaGradK( Pi, Pl, Pj, temp );
        areaj += temp;
      }

      VecType edge( Pj-Pi );
      RealType undefEdgeLengthSqr( dotProduct(edge,edge) ); 
      RealType undefEdgeLength = std::sqrt( undefEdgeLengthSqr );
      RealType undefEdgeLengthSqrSqr = undefEdgeLengthSqr * undefEdgeLengthSqr;
      //edge /= undefEdgeLength;
      
      RealType diff = undefEdgeLength - defEdgeLength;
      RealType areaFactor = diff * diff / (3. * undefEdgeLengthSqr );
      RealType mixedFactor = 2.* diff * defEdgeLength / (3. * undefEdgeLengthSqrSqr );
      RealType diagonalFactor = 2. * vol * diff * defEdgeLength / undefEdgeLengthSqrSqr;
      RealType eCrossEFactor = 2. * vol * defEdgeLength * (4. * defEdgeLength - 3. * undefEdgeLength) / ( undefEdgeLengthSqrSqr * undefEdgeLengthSqr );
      
      // now compute second derivatives
      MatType H, auxMat, EcrossE;
      EcrossE.makeTensorProduct( edge, edge );
      
      //*k      
      if( pk != -1 ){
	//kk      
        H.setZero();
        getHessAreaKK( Pi, Pj, Pk, auxMat );
        H.addMultiple( auxMat, areaFactor ); 
        localToGlobal( tripletList, pk, pk, H ); 
        //ik
        H.setZero();     
        getHessAreaIK( Pi, Pj, Pk, auxMat );
        H.addMultiple( auxMat, areaFactor );     
        auxMat.makeTensorProduct( edge, areak );
        H.addMultiple( auxMat, -1. * mixedFactor );      
        localToGlobal( tripletList, pi, pk, H ); 

        //jk
        H.setZero();     
        getHessAreaIK( Pj, Pi, Pk, auxMat );
        H.addMultiple( auxMat, areaFactor );    
        auxMat.makeTensorProduct( edge, areak );
        H.addMultiple( auxMat, mixedFactor );         
        localToGlobal( tripletList, pj, pk, H ); 
      }
      
      //*l
      if( pl != -1 ){
        //ll      
        H.setZero(); 
        getHessAreaKK( Pj, Pi, Pl, auxMat );
        H.addMultiple( auxMat, areaFactor  );           
        localToGlobal( tripletList, pl, pl, H ); 
      
        //il
        H.setZero(); 
        getHessAreaIK( Pi, Pj, Pl, auxMat );
        H.addMultiple( auxMat, areaFactor );     
        auxMat.makeTensorProduct( edge, areal );
        H.addMultiple( auxMat, -1. * mixedFactor );  
        localToGlobal( tripletList, pi, pl, H ); 
      
        //jl
        H.setZero();    
        getHessAreaIK( Pj, Pi, Pl, auxMat );
        H.addMultiple( auxMat, areaFactor );   
        auxMat.makeTensorProduct( edge, areal );
        H.addMultiple( auxMat, mixedFactor );  
        localToGlobal( tripletList, pj, pl, H ); 
      }
      
      //*j
     
      //jj     
      H.setZero();
      auxMat.makeTensorProduct( areaj, edge );
      H.addMultiple( auxMat, mixedFactor ); 
      auxMat.makeTensorProduct( edge, areaj );
      H.addMultiple( auxMat, mixedFactor ); 
      
      if( pk != -1 ){
        getHessAreaKK( Pk, Pi, Pj, auxMat );
        H.addMultiple( auxMat,  areaFactor );  
      }
      if( pl != -1 ){
        getHessAreaKK( Pi, Pl, Pj, auxMat );
        H.addMultiple( auxMat, areaFactor ); 
      }
      
      H.addMultiple( EcrossE, eCrossEFactor );
      H.addToDiagonal( diagonalFactor );
      localToGlobal( tripletList, pj, pj, H );
      
      //ij     
      H.setZero();
      auxMat.makeTensorProduct( areai, edge );
      H.addMultiple( auxMat, mixedFactor );       
      auxMat.makeTensorProduct( edge, areaj );
      H.addMultiple( auxMat, -1. * mixedFactor ); 
      
      if( pk != -1 ){
        getHessAreaIK( Pi, Pk, Pj, auxMat );
        H.addMultiple( auxMat, areaFactor );  
      }
      if( pl != -1 ){
        getHessAreaIK( Pi, Pl, Pj, auxMat );
        H.addMultiple( auxMat, areaFactor );      
      }
      
      H.addMultiple( EcrossE, -1. * eCrossEFactor );
      H.addToDiagonal( -1. * diagonalFactor  );
      localToGlobal( tripletList, pi, pj, H );
      
      //*i
      
      //ii     
      H.setZero();
      auxMat.makeTensorProduct( areai, edge );
      H.addMultiple( auxMat, -1. * mixedFactor );
      auxMat.makeTensorProduct( edge, areai );
      H.addMultiple( auxMat, -1. * mixedFactor );   
      
      if( pl != -1 ){
        getHessAreaKK( Pl, Pj, Pi, auxMat );
        H.addMultiple( auxMat, areaFactor );  
      }
      if( pk != -1 ){
        getHessAreaKK( Pj, Pk, Pi, auxMat );
        H.addMultiple( auxMat, areaFactor );  
      }
      
      H.addMultiple( EcrossE, eCrossEFactor );
      H.addToDiagonal( diagonalFactor  );
      localToGlobal( tripletList, pi, pi, H );
     
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


//==========================================================================================================
//! \brief Mixed second derivative of EdgeLengthEnergy w.r.t. the undeformed and deformed configuration
//! \author Heeren
template <typename ConfiguratorType>
class EdgeLengthHessianMixed : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::VecType     VecType;
  typedef typename ConfiguratorType::MatType     MatType;
  
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver& _topology;
  const VectorType&  _inactiveGeometry;
  const bool _activeShellIsDeformed, _firstDerivWRTDef;
  const RealType _factor;
  mutable int _rowOffset, _colOffset;

public:
EdgeLengthHessianMixed( const MeshTopologySaver& topology,
                        const VectorType& InactiveGeometry,
		        const bool ActiveShellIsDeformed,
			const bool FirstDerivWRTDef,
			   const RealType Factor = 1.,
                           int rowOffset = 0, 
                           int colOffset = 0,
                           RealType Mu = 1., 
                           RealType Lambda = 1. ) 
: _topology( topology), 
  _inactiveGeometry(InactiveGeometry), 
  _activeShellIsDeformed(ActiveShellIsDeformed), 
  _firstDerivWRTDef( FirstDerivWRTDef ), 
  _factor(Factor),
  _rowOffset(rowOffset), 
  _colOffset(colOffset){}
    
  void setRowOffset( int rowOffset ) const {
    _rowOffset = rowOffset;
  }
    
  void setColOffset( int colOffset ) const {
    _colOffset = colOffset;
  }

  //
  void apply( const VectorType& ActiveGeometry, MatrixType& Dest ) const {    
    int dofs = 3*_topology.getNumVertices();
    if( dofs != ActiveGeometry.size() )
        throw BasicException("EdgeLengthHessianMixed::apply: sizes dont match!");        
    if( (Dest.rows() != dofs) || (Dest.cols() != dofs) )
        Dest.resize( dofs, dofs );
    assembleHessian( ActiveGeometry, Dest );
  }

  // assmeble Hessian matrix via triplet list
  void assembleHessian( const VectorType& ActiveGeometry, MatrixType& Hessian ) const {   
      
    // set up triplet list
    TripletListType tripletList;
    // per face we have 3 active vertices, i.e. 9 combinations each producing a 3x3-matrix
    tripletList.reserve( 9 * 9 * _topology.getNumFaces() );

    pushTriplets( ActiveGeometry, tripletList );
    
    // fill matrix from triplets
    Hessian.setZero();
    Hessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }

  //
  void pushTriplets( const VectorType& ActiveGeometry, TripletListType& tripletList ) const {
      
    const VectorType* undefShellP = _activeShellIsDeformed ? &_inactiveGeometry : &ActiveGeometry;
    const VectorType* defShellP   = _activeShellIsDeformed ? &ActiveGeometry    : &_inactiveGeometry;
      
    for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){

      int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
          pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ),
          pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );

      // set up vertices and edges
      VecType Pi, Pj, Pk, Pl, temp;
      MatType matrix;

      // first get defomed quantities
      getXYZCoord<VectorType, VecType>( *defShellP, Pi, pi); 
      getXYZCoord<VectorType, VecType>( *defShellP, Pj, pj);     
      VecType defGradi(Pi-Pj), defGradj(Pj-Pi);      
      RealType defEdgeLength = std::sqrt( dotProduct( defGradi, defGradi ) );
      defGradi /= defEdgeLength;
      defGradj /= defEdgeLength;      
      
      // get undeformed vertex positions
      getXYZCoord<VectorType, VecType>( *undefShellP, Pi, pi); 
      getXYZCoord<VectorType, VecType>( *undefShellP, Pj, pj);            
      VecType undefEdge( Pj-Pi );
      RealType undefEdgeLengthSqr =  dotProduct( undefEdge, undefEdge );
      RealType undefEdgeLength = std::sqrt( undefEdgeLengthSqr );
      undefEdge /= undefEdgeLength;
      
      // compute Ak + Al, |e|^2 and diff. o dihedral angles
      RealType vol = 0.;
      if( pk != -1 ){
	getXYZCoord<VectorType, VecType>( *undefShellP, Pk, pk);
        temp.makeCrossProduct( Pk-Pj, Pi-Pk );
        vol += temp.norm() / 6.;
      }
      if( pl != -1 ){
	getXYZCoord<VectorType, VecType>( *undefShellP, Pl, pl);
        temp.makeCrossProduct( Pl-Pi, Pj-Pl );
        vol += temp.norm() / 6.;
      } 
      
      // compute first derivatives of area
      VecType areak, areal, areai, areaj;
      if( pk != -1 ){
	getAreaGradK( Pi, Pj, Pk, areak );
	getAreaGradK( Pj, Pk, Pi, areai );
	getAreaGradK( Pk, Pi, Pj, areaj );
      }
      if( pl != -1 ){
	getAreaGradK( Pj, Pi, Pl, areal );
	getAreaGradK( Pl, Pj, Pi, temp );
        areai += temp;
	getAreaGradK( Pi, Pl, Pj, temp );
        areaj += temp;
      }               
      
      // derivatives                 
      RealType areaFactor = 2. * (defEdgeLength - undefEdgeLength) / (3. * undefEdgeLengthSqr);
      RealType edgeFactor = 2. * vol * ( 2.*defEdgeLength - undefEdgeLength) / (undefEdgeLengthSqr * undefEdgeLength);

      areak *= areaFactor;
      areal *= areaFactor;
      areai *= areaFactor;
      areai.addMultiple( undefEdge, edgeFactor );
      areaj *= areaFactor;
      areaj.addMultiple( undefEdge, -1. * edgeFactor );      
      
      
      // k*
      if( pk != -1 ){        
        matrix.makeTensorProduct( areak, defGradi );
        localToGlobal( tripletList, pk, pi, matrix );
        matrix.makeTensorProduct( areak, defGradj );
        localToGlobal( tripletList, pk, pj, matrix );
      }
      
      // l*
      if( pl != -1 ){ 
        matrix.makeTensorProduct( areal, defGradi );
        localToGlobal( tripletList, pl, pi, matrix );
        matrix.makeTensorProduct( areal, defGradj );
        localToGlobal( tripletList, pl, pj, matrix ); 
      }
      
      // i*
      matrix.makeTensorProduct( areai, defGradi );
      localToGlobal( tripletList, pi, pi, matrix );
      matrix.makeTensorProduct( areai, defGradj );
      localToGlobal( tripletList, pi, pj, matrix );
      
      // j*
      matrix.makeTensorProduct( areaj, defGradi );
      localToGlobal( tripletList, pj, pi, matrix );
      matrix.makeTensorProduct( areaj, defGradj );
      localToGlobal( tripletList, pj, pj, matrix );

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