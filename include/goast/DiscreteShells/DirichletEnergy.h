// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Classical Dirichlet's energy but implemented in a Computer Graphic's fashion.
 * \author Heeren
 *
 */
#ifndef HARMONICENERGY_HH
#define HARMONICENERGY_HH

//== INCLUDES =================================================================
#include <goast/Core/Auxiliary.h>
#include <goast/Core/LocalMeshGeometry.h>
#include <goast/Core/Topology.h>
#include <goast/Core/DeformationInterface.h>
#include <goast/Core/Auxiliary.h>
#include <goast/Core/LocalMeshGeometry.h>
#include <goast/Core/Topology.h>

//==========================================================================================================
/**
 * \brief Classical Dirichlet energy implemented as sum over edges.
 * \author Heeren
 *
 * Let \f$ x, \tilde x \f$ two meshes (resp. their geometries) in dense correspondence.
 * Then this class represents the classical Dirichlet energy,
 * however, here we make use of a representation frequently used in computer graphics, i.e.:
 * \f[ E[x, \tilde x] =  \sum_{e \in x} ( \cot \alpha_e + \cot \beta_e ) l_e^2 \, , \f]
 * where \f$ \alpha_e \f$ and \f$ \beta_e \f$ are the opposite inner angles in the two adjacent triangles of \f$ e \f$
 * and \f$ l_e \f$ is the length of the edge.
 *
 * Note that we might rewrite \f$ E \f$ as follows:
 *  \f[ E[x, \tilde x] = =  \sum_{f \in x} \sum_{i = 0}^2 \cot \alpha_i  \| e_i \|^2 \, .\f]
 * Here \f$ e_0, e_1, e_2 \f$ are the three edges of a particular face \f$ f \f$
 * and \f$ \alpha_i \f$ the inner angle in \f$ f \f$ opposite edge \f$ e_i \f$.
 * Note that we have the relationship
 * \f[ \cot \alpha_i = \frac{ \langle e_{i-1}, e_{i+1} \rangle }{a_f}\, , \f]
 * where the index is to be understood modulo 3.
 *
 * Note that the energy might either be thought of as \f$ x \mapsto E[x, \tilde x] \f$ (active shell is undeformed shell)
 * or \f$ \tilde x \mapsto E[x, \tilde x] \f$ (active shell is deformed shell).
 * The active shell is considered the argument whereas the inactive shell is given in the constructor.
 */
template<typename ConfiguratorType>
class DirichletEnergy : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::MatType    MatType;

  const MeshTopologySaver& _topology;
  const VectorType&  _inactiveGeometry;
  const bool _activeShellIsDeformed;
  const RealType _weight;

public:

  DirichletEnergy( const MeshTopologySaver& topology,
                           const VectorType& InactiveGeometry,
		           const bool ActiveShellIsDeformed,
                           RealType Weight = 1.) 
  : _topology( topology), 
    _inactiveGeometry(InactiveGeometry), 
    _activeShellIsDeformed( ActiveShellIsDeformed),
    _weight(Weight){}
    

  // energy evaluation
  void apply( const VectorType& ActiveGeometry, RealType & Dest ) const {

    const VectorType* defShellP   = _activeShellIsDeformed ? &ActiveGeometry : &_inactiveGeometry;
    const VectorType* undefShellP = _activeShellIsDeformed ? &_inactiveGeometry : &ActiveGeometry;
    
    Dest = 0.;
    
    for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); ++faceIdx ){

      int pi( _topology.getNodeOfTriangle(faceIdx,0) ),
          pj( _topology.getNodeOfTriangle(faceIdx,1) ),
          pk( _topology.getNodeOfTriangle(faceIdx,2) );

      // set up deformed vertices and edges
      VecType Ei, Ej, Ek, temp;
      getXYZCoord<VectorType, VecType>( *defShellP, temp, pi);
      getXYZCoord<VectorType, VecType>( *defShellP, Ej, pj);
      getXYZCoord<VectorType, VecType>( *defShellP, Ek, pk);
      Ei = Ek - Ej;
      Ej = temp - Ek;
      Ek = Ei + Ej;
      
      // compute edge lengths
      RealType liSqr = Ei.normSqr();
      RealType ljSqr = Ej.normSqr();
      RealType lkSqr = Ek.normSqr();   
      
      // set up undeformed vertices and edges
      getXYZCoord<VectorType, VecType>( *undefShellP, temp, pi);
      getXYZCoord<VectorType, VecType>( *undefShellP, Ej, pj);
      getXYZCoord<VectorType, VecType>( *undefShellP, Ek, pk);
      Ei = Ek - Ej;
      Ej = temp - Ek;
      Ek = Ei+Ej;

      // compute volume
      temp.makeCrossProduct( Ei, Ej );
      RealType volUndef = std::sqrt( temp.normSqr() / 4. );
      
      //CAUTION mind the signs! (Ek is actually -Ek here!)
      RealType traceTerm = ( dotProduct(Ej,Ek) * liSqr + dotProduct(Ek,Ei) * ljSqr - dotProduct(Ei,Ej) * lkSqr );
      
      // volume of triangle * evaluation of energy density  
      Dest +=  0.125 * _weight *  traceTerm  / volUndef;
    }
  }

  void setWeight(RealType Weight)
  {
		_weight = Weight;
  }
 
  
 };

//! \brief First derivative of DirichletEnergy w.r.t. the deformed configuration.
//! \author Heeren
template<typename ConfiguratorType>
class DirichletGradientDef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {

protected:    
  typedef typename ConfiguratorType::RealType   RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::MatType    MatType;

  const MeshTopologySaver& _topology;
  const VectorType&  _undefShell;
  const RealType _weight;
  
public:
DirichletGradientDef( const MeshTopologySaver& topology,
                              const VectorType& undefShell,
                              RealType Weight = 1. ) 
  : _topology( topology), 
    _undefShell(undefShell),
    _weight(Weight){}

//
void apply( const VectorType& defShell, VectorType& Dest ) const {

  if( _undefShell.size() != defShell.size() )
      throw BasicException( "DirichletGradientDef::apply(): sizes dont match!");
  
  if( Dest.size() != defShell.size() )
    Dest.resize( defShell.size() );
  
  Dest.setZero();

    for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); ++faceIdx ){

      std::vector<int> nodesIdx(3);
      std::vector<VecType> nodes(3), undefEdges(3), defEdges(3);
      VecType temp;
      for( int j = 0; j < 3; j++ )
        nodesIdx[j] = _topology.getNodeOfTriangle(faceIdx,j);        
      
      //! get undeformed quantities
      for( int j = 0; j < 3; j++ )
        getXYZCoord<VectorType, VecType>( _undefShell, nodes[j], nodesIdx[j]);
      for( int j = 0; j < 3; j++ )
          undefEdges[j] = nodes[(j+2)%3] - nodes[(j+1)%3];
      // compute volume
      temp.makeCrossProduct( nodes[2] - nodes[1], nodes[0] - nodes[2] );
      RealType volUndefSqr = temp.normSqr() / 4.;
      RealType volUndef = std::sqrt( volUndefSqr );      
      
      //! get deformed quantities
      for( int j = 0; j < 3; j++ )
        getXYZCoord<VectorType, VecType>( defShell, nodes[j], nodesIdx[j]); 
      for( int j = 0; j < 3; j++ )
          defEdges[j] = nodes[(j+2)%3] - nodes[(j+1)%3];
      
      //! trace part of gradient, energy_tr = volUndef * _mu/2. *  trace(DistTensor) 
      VecType factors;
      for( int i = 0; i < 3; i++ )
        factors[i] = -0.25 * _weight * dotProduct( undefEdges[(i+2)%3], undefEdges[(i+1)%3] ) / volUndef;  

      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          Dest[j*_topology.getNumVertices() + nodesIdx[i]] += factors[(i+1)%3] * defEdges[(i+1)%3][j] - factors[(i+2)%3] * defEdges[(i+2)%3][j];

  }
}

void setWeight(RealType Weight)
  {
		_weight = Weight;
  }
 
  
}; 

/**
 * \brief Second derivative of DirichletEnergy w.r.t. the deformed configuration.
 * \author Heeren
 *
 * \note The harmonic energy is quadratic wrt. its second (i.e. the deformed) argument, hence we pre-compute and store the constant Hessian matrix
 */
template<typename ConfiguratorType>
class DirichletHessianDef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

protected:
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::VecType     VecType;
  typedef typename ConfiguratorType::MatType     MatType;
  
  typedef std::vector<TripletType> TripletListType;
  
  const MeshTopologySaver& _topology;
  const VectorType& _undefShell;
  RealType _factor;
  int _rowOffset, _colOffset;
  
  TripletListType _constHessianTriplets;
  MatrixType _constHessian;
  int _nnz; // number of nonzero entries in Hessian (per face we have 3 active vertices, i.e. 9 combinations each producing a 3x3-matrix)

public:
DirichletHessianDef( const MeshTopologySaver& topology,
                             const VectorType& undefShell,
			     const RealType Factor = 1.,
                             int rowOffset = 0,
                             int colOffset = 0 ) 
  : _topology( topology), 
    _undefShell(undefShell), 
    _factor( Factor ),
    _rowOffset(rowOffset), 
    _colOffset(colOffset),
    _constHessian(3*_topology.getNumVertices(), 3*_topology.getNumVertices()),
    _nnz(9 * 9 * _topology.getNumFaces()){ // per face we have 3 active vertices, i.e. 9 combinations each producing a 3x3-matrix
        assembleHessian();
    }
    
    void setRowOffset( int rowOffset ) {
        _rowOffset = rowOffset;
        assembleHessian( );
    }
    
    void setColOffset( int colOffset ) {
        _colOffset = colOffset;
        assembleHessian( );
    }

  //
  void apply( const VectorType& /*defShell*/, MatrixType& Dest ) const {    
      int dofs = 3*_topology.getNumVertices();
      if( (Dest.rows() != dofs) || (Dest.cols() != dofs) )
        Dest.resize( dofs, dofs );
      Dest = _constHessian;
  }

  //
  void pushTriplets( const VectorType& /*defShell*/, TripletListType& tripletList ) const  {
    tripletList.clear();
    tripletList.resize( _constHessianTriplets.size() );
    for( uint i = 0; i < _constHessianTriplets.size(); i++ )
        tripletList[i] = _constHessianTriplets[i];
  }

  void setWeight(RealType Weight)
  {
		_factor = Weight;
		assembleHessian();
  }
 
  
protected:      
  // assemble Hessian matrix
  void assembleHessian( ) {    
    // per face we have 3 active vertices, i.e. 9 combinations each producing a 3x3-matrix
    _constHessianTriplets.clear();
    _constHessianTriplets.reserve( _nnz );
    
    // run over all faces and push triplets
    for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); ++faceIdx ){

      std::vector<int> nodesIdx(3);
      std::vector<VecType> nodes(3), undefEdges(3), defEdges(3);
      VecType temp;
      for( int j = 0; j < 3; j++ )
        nodesIdx[j] = _topology.getNodeOfTriangle(faceIdx,j);        
      
      //! get undeformed quantities
      for( int j = 0; j < 3; j++ )
        getXYZCoord<VectorType, VecType>( _undefShell, nodes[j], nodesIdx[j]);
      for( int j = 0; j < 3; j++ )
          undefEdges[j] = nodes[(j+2)%3] - nodes[(j+1)%3];
      // compute volume
      temp.makeCrossProduct( undefEdges[0], undefEdges[1] );
      RealType volUndefSqr = temp.normSqr() / 4.;
      RealType volUndef = std::sqrt( volUndefSqr );    
      
      VecType traceFactors;
      for( int i = 0; i < 3; i++ )
        traceFactors[i] = -0.25 * dotProduct( undefEdges[(i+2)%3], undefEdges[(i+1)%3] ) / volUndef;

      // compute local matrices
      MatType H;
      
      // i==j
      for( int i = 0; i < 3; i++ ){
        H.setZero();
        H.addToDiagonal( traceFactors[(i+1)%3] + traceFactors[(i+2)%3] );
        localToGlobal( nodesIdx[i], nodesIdx[i], H );
      }      
      
      // i!=j
      for( int i = 0; i < 3; i++ ){
        H.setZero();
        H.addToDiagonal( -traceFactors[(i+1)%3] );
        localToGlobal( nodesIdx[i], nodesIdx[(i+2)%3], H );
      }

    }
    
    // fill matrix from triplets
    _constHessian.setZero();
    _constHessian.setFromTriplets( _constHessianTriplets.cbegin(), _constHessianTriplets.cend() );
  }
  
  // push triplets
  void localToGlobal( int k, int l, const MatType& localMatrix ) {
    int numV = _topology.getNumVertices();
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          _constHessianTriplets.push_back( TripletType( _rowOffset + i*numV + k, _colOffset + j*numV + l, _factor * localMatrix(i,j) ) );	
	
    if( k != l){
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          _constHessianTriplets.push_back( TripletType( _rowOffset + i*numV + l, _colOffset + j*numV + k, _factor * localMatrix(j,i) ) );
    }
  }
  
};


//! \brief First derivative of DirichletEnergy w.r.t. the undeformed configuration.
//! \author Heeren
template< typename ConfiguratorType>
class DirichletGradientUndef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {

protected:    
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::VecType     VecType;

  const MeshTopologySaver& _topology;
  const VectorType& _defShell;
  const RealType _weight;
  
public:
DirichletGradientUndef( const MeshTopologySaver& topology, 
                                const VectorType& defShell, 
                                RealType Weight = 1. )
  : _topology( topology),
    _defShell(defShell),
    _weight(Weight){}

//
void apply( const VectorType& undefShell, VectorType& Dest ) const {
  
  if( Dest.size() != undefShell.size() )
    Dest.resize( undefShell.size() );
  
  Dest.setZero();

    for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); ++faceIdx ){

      std::vector<int> nodesIdx(3);
      std::vector<VecType> nodes(3), undefEdges(3), fixedEdges(3);
      VecType temp;
      for( int j = 0; j < 3; j++ )
        nodesIdx[j] = _topology.getNodeOfTriangle(faceIdx,j);     
      
      //! get fixed edgess      
      for( int j = 0; j < 3; j++ )
        getXYZCoord<VectorType, VecType>( _defShell, nodes[j], nodesIdx[j]);
      for( int j = 0; j < 3; j++ )
          fixedEdges[j] = nodes[(j+2)%3] - nodes[(j+1)%3];
      
      VecType defLengthSqr;
      for( int i = 0; i < 3; i++ )
          defLengthSqr[i] = fixedEdges[i].normSqr();
      
      //! get undeformed quantities
      for( int j = 0; j < 3; j++ )
        getXYZCoord<VectorType, VecType>( undefShell, nodes[j], nodesIdx[j]); 
      for( int j = 0; j < 3; j++ )
          undefEdges[j] = nodes[(j+2)%3] - nodes[(j+1)%3];
      
      // compute volume
      temp.makeCrossProduct(  nodes[1] - nodes[0], nodes[2] - nodes[1] );
      RealType volUndefSqr = temp.normSqr() / 4.;
      RealType volUndef = std::sqrt( volUndefSqr );
      
      RealType traceTerm = 0.;
      for( int i = 0; i < 3; i++ )
        traceTerm -= dotProduct( undefEdges[(i+1)%3], undefEdges[(i+2)%3] ) * defLengthSqr[i];     

      RealType factorAreaGrad  = 0.125 * _weight *  traceTerm / volUndefSqr;
      RealType factorTraceGrad = 0.125 * _weight / volUndef;
      
      std::vector<VecType> gradTrace(3);
      for( int i = 0; i < 3; i++ )
          for( int j = 0; j < 3; j++ )   
            gradTrace[i][j] = defLengthSqr[i] * (undefEdges[(i+1)%3][j] - undefEdges[(i+2)%3][j]) +  undefEdges[i][j] * (defLengthSqr[(i+1)%3] - defLengthSqr[(i+2)%3] );      
      
      // E = (0.125 * _mu * traceTerm ) / volUndef;
      // grad E[j] = 0.125 * _mu * grad traceTerm[j] / volUndef  - (0.125 * _mu * traceTerm )  grad volUndef[j] / volUndef^2
      for( int i = 0; i < 3; i++ ){
        getAreaGradient(  nodes[(i+1)%3], nodes[(i+2)%3], nodes[i], temp );
        for( int j = 0; j < 3; j++ )
          Dest[j*_topology.getNumVertices() + nodesIdx[i]] += factorTraceGrad * gradTrace[i][j] - factorAreaGrad * temp[j];
      }

  }
}

void setWeight(RealType Weight)
  {
		_weight = Weight;
  }
 
  
}; 


//! \brief Second derivative of DirichletEnergy w.r.t. the undeformed configuration.
//! \author Heeren
template <typename ConfiguratorType>
class DirichletHessianUndef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

protected:    
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::VecType     VecType;
  typedef typename ConfiguratorType::MatType     MatType;
  
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver& _topology;
  const VectorType& _defShell;
  const RealType _factor;
  int _rowOffset, _colOffset;
  
public:
DirichletHessianUndef( const MeshTopologySaver& topology, const VectorType& defShell, RealType Factor = 1., int rowOffset = 0, int colOffset = 0 ) 
  : _topology( topology),
    _defShell(defShell),
    _factor(Factor),
    _rowOffset(rowOffset), 
    _colOffset(colOffset){}
    
  void setRowOffset( int rowOffset ) {
    _rowOffset = rowOffset;
  }
    
  void setColOffset( int colOffset ) {
    _colOffset = colOffset;
  }

  //
  void apply( const VectorType& undefShell, MatrixType& Dest ) const {    
    int dofs = 3*_topology.getNumVertices();
    if( (Dest.rows() != dofs) || (Dest.cols() != dofs) )
        Dest.resize( dofs, dofs );
    Dest.setZero();
    assembleHessian( undefShell, Dest );
  }


  // assmeble Hessian matrix via triplet list
  void assembleHessian( const VectorType& undefShell, MatrixType& Hessian ) const {

    // set up triplet list
    TripletListType tripletList;
    // per face we have 3 active vertices, i.e. 9 combinations each producing a 3x3-matrix
    tripletList.reserve( 9 * 9 * _topology.getNumFaces() );

    pushTriplets( undefShell, tripletList );
    // fill matrix from triplets
    Hessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
  
  //
  void pushTriplets( const VectorType& undefShell, TripletListType& tripletList ) const  {
    // run over all faces
    for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); ++faceIdx ){

      std::vector<int> nodesIdx(3);
      std::vector<VecType> nodes(3), undefEdges(3), fixedEdges(3);
      VecType temp;
      for( int j = 0; j < 3; j++ )
        nodesIdx[j] = _topology.getNodeOfTriangle(faceIdx,j);        

      //! get fixed edgess      
      for( int j = 0; j < 3; j++ )
        getXYZCoord<VectorType, VecType>( _defShell, nodes[j], nodesIdx[j]);
      for( int j = 0; j < 3; j++ )
          fixedEdges[j] = nodes[(j+2)%3] - nodes[(j+1)%3];
      
      // compute volume
      temp.makeCrossProduct(  nodes[1] - nodes[0], nodes[2] - nodes[1] );
      RealType volDefSqr = temp.normSqr() / 4.;
      RealType volDef = std::sqrt( volDefSqr ); 
      
      VecType defLengthSqr;
      for( int i = 0; i < 3; i++ )
          defLengthSqr[i] = fixedEdges[i].normSqr();
      
      //! get deformed quantities
      for( int j = 0; j < 3; j++ )
        getXYZCoord<VectorType, VecType>( undefShell, nodes[j], nodesIdx[j]); 
      for( int j = 0; j < 3; j++ )
          undefEdges[j] = nodes[(j+2)%3] - nodes[(j+1)%3];
      // compute volume
      temp.makeCrossProduct( undefEdges[0], undefEdges[1] );
      RealType volUndefSqr = temp.normSqr() / 4.;
      RealType volUndef = std::sqrt( volUndefSqr );
      
      RealType traceTerm = 0.;
      for( int i = 0; i < 3; i++ )
        traceTerm -= dotProduct( undefEdges[(i+1)%3], undefEdges[(i+2)%3] ) * defLengthSqr[i];     
      
      std::vector<VecType> gradTrace(3);
      for( int i = 0; i < 3; i++ )
          for( int j = 0; j < 3; j++ )   
            gradTrace[i][j] = defLengthSqr[i] * (undefEdges[(i+1)%3][j] - undefEdges[(i+2)%3][j]) +  undefEdges[i][j] * (defLengthSqr[(i+1)%3] - defLengthSqr[(i+2)%3]);
      
      // precompute area gradients
      std::vector<VecType> gradArea(3);
      for( int i = 0; i < 3; i++ )
        getAreaGradient(  nodes[(i+1)%3], nodes[(i+2)%3], nodes[i], gradArea[i] );

      RealType areaFactor        = 0.125 * traceTerm;
      RealType negHessAreaFactor =  areaFactor / volUndefSqr;
      RealType mixedAreaFactor   = 2 * areaFactor / (volUndefSqr * volUndef);
      RealType mixedFactor       = -0.125 / volUndefSqr;
      RealType hessTraceFactor   =  0.125 / volUndef;
      
      
      // compute local matrices
      MatType tensorProduct, H, auxMat;
      
      // i==j
      for( int i = 0; i < 3; i++ ){
        getHessAreaKK( nodes[(i+1)%3], nodes[(i+2)%3], nodes[i], H );
        H *= -1. * negHessAreaFactor;
        
        tensorProduct.makeTensorProduct( gradArea[i], gradTrace[i] );
        H.addMultiple( tensorProduct, mixedFactor );
        tensorProduct.makeTensorProduct( gradTrace[i], gradArea[i] );
        H.addMultiple( tensorProduct, mixedFactor );
        
        tensorProduct.makeTensorProduct( gradArea[i], gradArea[i] );
        H.addMultiple( tensorProduct, mixedAreaFactor );
        
        H.addToDiagonal( hessTraceFactor * 2 * defLengthSqr[i] );
        localToGlobal( tripletList, nodesIdx[i], nodesIdx[i], H );
      }    
      
      // i!=j
      for( int i = 0; i < 3; i++ ){
        getHessAreaIK( nodes[i], nodes[(i+1)%3], nodes[(i+2)%3], H );
        H *= -1. * negHessAreaFactor;
        
        tensorProduct.makeTensorProduct( gradArea[i], gradTrace[(i+2)%3] );
        H.addMultiple( tensorProduct, mixedFactor );
        tensorProduct.makeTensorProduct( gradTrace[i], gradArea[(i+2)%3] );
        H.addMultiple( tensorProduct, mixedFactor );
        
        tensorProduct.makeTensorProduct( gradArea[i], gradArea[(i+2)%3] );
        H.addMultiple( tensorProduct, mixedAreaFactor );
        
        H.addToDiagonal( hessTraceFactor * (defLengthSqr[(i+1)%3] - defLengthSqr[i] - defLengthSqr[(i+2)%3]) );
        localToGlobal( tripletList, nodesIdx[i], nodesIdx[(i+2)%3], H );
      }

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


//! \brief Second (mixed) derivative of DirichletEnergy
//! \author Heeren
template <typename ConfiguratorType >
class DirichletHessianMixed : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

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
  int _rowOffset, _colOffset;

public:
DirichletHessianMixed( const MeshTopologySaver& topology,
                           const VectorType& InactiveGeometry,
		           const bool ActiveShellIsDeformed,
			   const bool FirstDerivWRTDef,
			   const RealType Factor = 1.,
                           int rowOffset = 0, 
                           int colOffset = 0 ) 
: _topology( topology), 
  _inactiveGeometry(InactiveGeometry), 
  _activeShellIsDeformed(ActiveShellIsDeformed), 
  _firstDerivWRTDef( FirstDerivWRTDef ), 
  _factor(Factor),
  _rowOffset(rowOffset), 
  _colOffset(colOffset){}
    
  void setRowOffset( int rowOffset ) {
    _rowOffset = rowOffset;
  }
    
  void setColOffset( int colOffset ) {
    _colOffset = colOffset;
  }

  //
  void apply( const VectorType& ActiveGeometry, MatrixType& Dest ) const {    
    int dofs = 3*_topology.getNumVertices();
    if( dofs != ActiveGeometry.size() )
        throw BasicException("DirichletHessianMixed::apply: sizes dont match!");        
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
    
    // run over all faces
    for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); ++faceIdx ){

      std::vector<int> nodesIdx(3);
      std::vector<VecType> nodes(3), undefEdges(3), defEdges(3);
      VecType temp;
      for( int j = 0; j < 3; j++ )
        nodesIdx[j] = _topology.getNodeOfTriangle(faceIdx,j);        
           
      
      //! get undeformed quantities
      for( int j = 0; j < 3; j++ )
        getXYZCoord<VectorType, VecType>( *undefShellP, nodes[j], nodesIdx[j]);
      for( int j = 0; j < 3; j++ )
          undefEdges[j] = nodes[(j+2)%3] - nodes[(j+1)%3];
      
      // compute volume
      temp.makeCrossProduct( undefEdges[0], undefEdges[1] );
      RealType volUndefSqr = temp.normSqr() / 4.;
      RealType volUndef = std::sqrt( volUndefSqr );      
            
      // precompute undeformed area gradients
      std::vector<VecType> gradUndefArea(3);
      for( int i = 0; i < 3; i++ )
        getAreaGradient(  nodes[(i+1)%3], nodes[(i+2)%3], nodes[i], gradUndefArea[i] );
      
      //
      VecType factors;
      for( int i = 0; i < 3; i++ )
        factors[i] = dotProduct( undefEdges[(i+1)%3], undefEdges[(i+2)%3]);      
      
      //! get deformed quantities     
      for( int j = 0; j < 3; j++ )
        getXYZCoord<VectorType, VecType>( *defShellP, nodes[j], nodesIdx[j]);
      for( int j = 0; j < 3; j++ )
          defEdges[j] = nodes[(j+2)%3] - nodes[(j+1)%3];
      
      // compute volume
      temp.makeCrossProduct(  nodes[1] - nodes[0], nodes[2] - nodes[1] );
      RealType volDefSqr = temp.normSqr() / 4.;
      RealType volDef = std::sqrt( volDefSqr ); 
        
      // precomputed deformed trace gradients
      std::vector<VecType> gradDefTrace(3);
      for( int i = 0; i < 3; i++ )
          getWeightedVectorSum<RealType>( -2 * factors[(i+1)%3], defEdges[(i+1)%3], 2 * factors[(i+2)%3], defEdges[(i+2)%3], gradDefTrace[i] );
   
      
      // compute local matrices
      MatType tensorProduct, H, auxMat;           
      RealType mixedTraceHessFactor = 0.125 / volUndef;
      RealType MixedFactor          = -0.125  / volUndefSqr;       

      // i!=j
      for( int i = 0; i < 3; i++ ){
        for( int j = 0; j < 3; j++ ){
          // Hess trace term
          if( i == j ){
            H.makeTensorProduct( undefEdges[i], defEdges[i] );
          }
          else{
            int k = (2*i+2*j)%3;        
            H.makeTensorProduct( undefEdges[j] - undefEdges[k], defEdges[i] );
            auxMat.makeTensorProduct( undefEdges[i], defEdges[k]  );
            H += auxMat;
          }
          H *= -2 * mixedTraceHessFactor;
        
          // mixed term
          tensorProduct.makeTensorProduct( gradUndefArea[i], gradDefTrace[j] );
          H.addMultiple( tensorProduct, MixedFactor );
 
          localToGlobal( tripletList, nodesIdx[i], nodesIdx[j], H );
        }
      }

    }
  }
  
protected:
  void localToGlobal( TripletListType& tripletList, int k, int l, const MatType& localMatrix ) const {
    int numV = _topology.getNumVertices();
    if( !_firstDerivWRTDef ){
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          tripletList.push_back( TripletType( _rowOffset + i*numV + k, _colOffset + j*numV + l, _factor * localMatrix(i,j) ) );	
    }
    else{
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          tripletList.push_back( TripletType( _rowOffset + i*numV + l, _colOffset + j*numV + k, _factor * localMatrix(j,i) ) );
    }
  }
  
};


//==========================================================================================================
//! \brief Deformation class for DirichletEnergy
//! \author Heeren
template <typename ConfiguratorType>
class DirichletDeformation : public DeformationBase<ConfiguratorType>{
    
protected:    
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;

  typedef typename ConfiguratorType::TripletType TripletType;  
  typedef std::vector<TripletType> TripletListType;
  
  const MeshTopologySaver& _topology;
  RealType _memWeight;

public:
  DirichletDeformation( const MeshTopologySaver& Topology, RealType memWeight ) : _topology( Topology ), _memWeight( memWeight ) {}
  
  void applyEnergy ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, RealType & Dest ) const {
    DirichletEnergy<ConfiguratorType>( _topology, UndeformedGeom, true ).apply( DeformedGeom, Dest );
    Dest *= _memWeight;
  }
  
  void applyUndefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest ) const {
    DirichletGradientUndef<ConfiguratorType>( _topology, DeformedGeom ).apply( UndeformedGeom, Dest );
    Dest *= _memWeight;
  }
  
  void applyDefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest ) const {
    DirichletGradientDef<ConfiguratorType>( _topology, UndeformedGeom ).apply( DeformedGeom, Dest );
    Dest *= _memWeight;
  }
  
  void pushTripletsDefHessian ( const VectorType& UndefGeom, const VectorType& DefGeom, TripletListType& triplets, int rowOffset, int colOffset, RealType factor = 1.0 ) const {
    DirichletHessianDef<ConfiguratorType>( _topology, UndefGeom, factor * _memWeight, rowOffset, colOffset ).pushTriplets( DefGeom, triplets );
  }
   
  void pushTripletsUndefHessian ( const VectorType& UndefGeom, const VectorType& DefGeom, TripletListType& triplets, int rowOffset, int colOffset, RealType factor = 1.0 ) const {
      DirichletHessianUndef<ConfiguratorType>( _topology, DefGeom, factor * _memWeight, rowOffset, colOffset ).pushTriplets( UndefGeom, triplets );
  }
  
  // mixed second derivative of deformation energy E[S_1, S_2], i.e. if "FirstDerivWRTDef" we have D_1 D_2 E[.,.], otherwise D_2 D_1 E[.,.]
  void pushTripletsMixedHessian ( const VectorType& UndefGeom, const VectorType& DefGeom, TripletListType& triplets, int rowOffset, int colOffset, const bool FirstDerivWRTDef, RealType factor = 1.0 ) const {
      DirichletHessianMixed<ConfiguratorType>( _topology, UndefGeom, true, FirstDerivWRTDef, factor * _memWeight, rowOffset, colOffset ).pushTriplets( DefGeom, triplets );
  }
  
  int numOfNonZeroHessianEntries () const {
    // per face we have 3 active vertices, i.e. 9 combinations each producing a 3x3-matrix
    return 9 * 9 * _topology.getNumFaces();    
  }
    
};

#endif
