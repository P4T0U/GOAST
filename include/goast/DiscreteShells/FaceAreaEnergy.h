// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Face area energy as a very simple membrane model (to be used together with edge length energy).
 * \author Heeren
 * \todo TO BE CONVERTED FROM QUOCMESH CODE!!!!
 */

#ifdef QUOCMESH_DISCRETE_SHELLS_FACE_ENERGY 

//!==========================================================================================================
//! DISCRETE SHELLS FACE AREA ENERGY
//!==========================================================================================================

//! \brief Face area energy between two thin shells given as triangular meshes.
//! \author Heeren
//!
//! Energy is taken from "Discrete shells" paper by Grinspun et al., 2003.
//! Here \f$ E[S, \tilde S] = \sum_t \frac{ (A_t - A_{\tilde t})^2 }{ A_t } \f$, 
//! where the sum is over all triangles \f$ t \in S \f$, where \f$ A_t \f$ and \f$ A_{\tilde t} \f$ are the volumes of \f$ t \f$ and \f$ \tilde t \f$, respectively.
//!
//! Note that the energy might either be thought of as \f$ S \mapsto E[S, \tilde S] \f$ (active shell is undeformed shell) or \f$ \tilde S \mapsto E[S, \tilde S] \f$ (active shell is deformed shell).
//! The active shell is considered the argument whereas the inactive shell is given in the constructor.
template< typename MeshType >
class FaceAreaEnergy : public aol::Op< aol::MultiVector<typename MeshType::RealType>, aol::Scalar<typename MeshType::RealType> > {

protected:
  typedef typename MeshType::RealType  RealType;
  typedef aol::MultiVector<RealType>   VectorType;

  const MeshTopologySaver<MeshType>& _topology;
  const aol::MultiVector<RealType>&  _inactiveGeometry;
  const bool _activeShellIsDeformed;

public:

  FaceAreaEnergy( const MeshTopologySaver<MeshType>& topology,
                       const VectorType& InactiveGeometry,
		       const bool ActiveShellIsDeformed ) 
  : _topology( topology), 
    _inactiveGeometry(InactiveGeometry), 
    _activeShellIsDeformed( ActiveShellIsDeformed) {}

  // energy evaluation
  void applyAdd( const VectorType& ActiveGeometry, aol::Scalar<RealType>& Dest ) const {

    const VectorType* defShellP   = _activeShellIsDeformed ? &ActiveGeometry : &_inactiveGeometry;
    const VectorType* undefShellP = _activeShellIsDeformed ? &_inactiveGeometry : &ActiveGeometry;
    
    for ( int faceIdx = 0; faceIdx < this->_topology.getNumFaces(); ++faceIdx ){

      int pi( this->_topology.getNodeOfTriangle(faceIdx,0) ),
          pj( this->_topology.getNodeOfTriangle(faceIdx,1) ),
          pk( this->_topology.getNodeOfTriangle(faceIdx,2) );

      // set up vertices and edges
      aol::Vec3<RealType> Pi, Pj, Pk, normal;

      defShellP->getTo( pi, Pi ); 
      defShellP->getTo( pj, Pj );
      defShellP->getTo( pk, Pk );
      normal.makeCrossProduct( Pk-Pj, Pi-Pk );      
      RealType defFaceArea = normal.norm() / 2.; 
      
      undefShellP->getTo( pi, Pi ); 
      undefShellP->getTo( pj, Pj );
      undefShellP->getTo( pk, Pk );
      normal.makeCrossProduct( Pk-Pj, Pi-Pk );      
      RealType undefFaceArea = normal.norm() / 2.; 
      
      Dest[0] += aol::Sqr( undefFaceArea - defFaceArea) / undefFaceArea;

    }
  }

 };
 
//==========================================================================================================
//! \brief First derivative of face area energy w.r.t. the deformed configuration (cf. class FaceAreaEnergy<> above).
//! \author Heeren
template< typename MeshType >
class FaceAreaGradientDef : public aol::Op< aol::MultiVector<typename MeshType::RealType> > {

  typedef typename MeshType::RealType RealType;
  typedef aol::MultiVector<RealType>  VectorType;
  typedef aol::Vec3<RealType>         Vec;

  const MeshTopologySaver<MeshType>& _topology;
  const aol::MultiVector<RealType>&  _undefShell;

public:
FaceAreaGradientDef( const MeshTopologySaver<MeshType>& topology,
                          const VectorType& undefShell ) : _topology( topology), _undefShell(undefShell) {}

FaceAreaGradientDef( const MeshType& /*dummy*/,
			  const MeshTopologySaver<MeshType>& topology,
                          const VectorType& undefShell ) : _topology( topology), _undefShell(undefShell) {}

void applyAdd( const VectorType& defShell, VectorType& Dest ) const {

  assert( _undefShell.getTotalSize() == defShell.getTotalSize() );
  if( (Dest.numComponents() != 3) || ( Dest[0].size() != defShell[0].size() ) )
    Dest.reallocate( 3, defShell[0].size() );
  
    for ( int faceIdx = 0; faceIdx < this->_topology.getNumFaces(); ++faceIdx ){

      int pi( this->_topology.getNodeOfTriangle(faceIdx,0) ),
          pj( this->_topology.getNodeOfTriangle(faceIdx,1) ),
          pk( this->_topology.getNodeOfTriangle(faceIdx,2) );

      // set up vertices and edges
      aol::Vec3<RealType> Pi, Pj, Pk, grad;
      
      _undefShell.getTo( pi, Pi ); 
      _undefShell.getTo( pj, Pj );
      _undefShell.getTo( pk, Pk );
      grad.makeCrossProduct( Pk-Pj, Pi-Pk );      
      RealType undefFaceArea = grad.norm() / 2.; 
      
      defShell.getTo( pi, Pi ); 
      defShell.getTo( pj, Pj );
      defShell.getTo( pk, Pk );
      grad.makeCrossProduct( Pk-Pj, Pi-Pk );      
      RealType defFaceArea = grad.norm() / 2.; 
      
      // pk
      getAreaGradK( Pi, Pj, Pk, grad );
      for( int i = 0; i < 3; i++ )
        Dest[i][pk] += 2. * ( defFaceArea - undefFaceArea ) * grad[i] / undefFaceArea;
      
      // d_i
      getAreaGradK( Pj, Pk, Pi, grad );
      for( int i = 0; i < 3; i++ )
        Dest[i][pi] += 2. * ( defFaceArea - undefFaceArea ) * grad[i] / undefFaceArea;
      
      // d_j
      getAreaGradK( Pk, Pi, Pj, grad );
      for( int i = 0; i < 3; i++ )
        Dest[i][pj] += 2. * ( defFaceArea - undefFaceArea ) * grad[i] / undefFaceArea;

    }
  }

};

//==========================================================================================================
//! \brief Second derivative of face area energy w.r.t. the deformed configuration (cf. class FaceAreaEnergy<> above).
//! \author Heeren
template <typename MeshType, typename BlockMatrixType = aol::SparseBlockMatrix< aol::SparseMatrix<typename MeshType::RealType> > >
class FaceAreaHessianDef : public aol::Op<aol::MultiVector<typename MeshType::RealType>, BlockMatrixType > {

protected:
  typedef typename MeshType::RealType RealType;
  typedef aol::Vec3<RealType>                Vec;
  typedef aol::MultiVector<RealType>         VectorType;

  const MeshTopologySaver<MeshType>& _topology;
  const VectorType& _undefShell;
  const RealType _factor;

public:
FaceAreaHessianDef( const MeshTopologySaver<MeshType>& topology,
                         const VectorType& undefShell,
			 const RealType Factor = 1. ) : _topology( topology), _undefShell(undefShell), _factor( Factor ) {} 



void applyAdd( const VectorType& defShell, BlockMatrixType& Dest ) const {

  for ( int faceIdx = 0; faceIdx < this->_topology.getNumFaces(); ++faceIdx ){

      int pi( this->_topology.getNodeOfTriangle(faceIdx,0) ),
          pj( this->_topology.getNodeOfTriangle(faceIdx,1) ),
          pk( this->_topology.getNodeOfTriangle(faceIdx,2) );

      // set up vertices and edges
      aol::Vec3<RealType> Pi, Pj, Pk, gradi, gradj, gradk;
      
      _undefShell.getTo( pi, Pi ); 
      _undefShell.getTo( pj, Pj );
      _undefShell.getTo( pk, Pk );
      gradi.makeCrossProduct( Pk-Pj, Pi-Pk );      
      RealType undefFaceArea = gradi.norm() / 2.; 
      
      defShell.getTo( pi, Pi ); 
      defShell.getTo( pj, Pj );
      defShell.getTo( pk, Pk );
      gradi.makeCrossProduct( Pk-Pj, Pi-Pk );      
      RealType defFaceArea = gradi.norm() / 2.; 
      
      RealType hessFactor = 2. * (defFaceArea - undefFaceArea) / undefFaceArea;
      RealType mixedFactor = 2. / undefFaceArea;
      
      // gradients
      getAreaGradK( Pi, Pj, Pk, gradk );
      getAreaGradK( Pj, Pk, Pi, gradi );
      getAreaGradK( Pk, Pi, Pj, gradj );
      
      //
      aol::Matrix33<RealType> H, auxMat;
	
      //*k      
      //kk      
      H.setZero();
      getHessAreaKK( Pi, Pj, Pk, auxMat );
      H.addMultiple( auxMat, hessFactor ); 
      auxMat.makeTensorProduct( gradk, gradk );
      H.addMultiple( auxMat, mixedFactor );
      localToGlobal( Dest, pk, pk, H ); 
      //ik
      H.setZero();     
      getHessAreaIK( Pi, Pj, Pk, auxMat );
      H.addMultiple( auxMat, hessFactor );     
      auxMat.makeTensorProduct( gradi, gradk );
      H.addMultiple( auxMat, mixedFactor );   
      localToGlobal( Dest, pi, pk, H ); 
      //jk
      H.setZero();     
      getHessAreaIK( Pj, Pi, Pk, auxMat );
      H.addMultiple( auxMat, hessFactor );    
      auxMat.makeTensorProduct( gradj, gradk );
      H.addMultiple( auxMat, mixedFactor );          
      localToGlobal( Dest, pj, pk, H ); 
      
      //*j      
      //jj      
      H.setZero();
      getHessAreaKK( Pk, Pi, Pj, auxMat );
      H.addMultiple( auxMat, hessFactor ); 
      auxMat.makeTensorProduct( gradj, gradj );
      H.addMultiple( auxMat, mixedFactor );
      localToGlobal( Dest, pj, pj, H ); 
      //ij
      H.setZero();     
      getHessAreaIK( Pi, Pk, Pj, auxMat );
      H.addMultiple( auxMat, hessFactor );     
      auxMat.makeTensorProduct( gradi, gradj );
      H.addMultiple( auxMat, mixedFactor );   
      localToGlobal( Dest, pi, pj, H ); 
    
      //ii  
      H.setZero();
      getHessAreaKK( Pj, Pk, Pi, auxMat );
      H.addMultiple( auxMat, hessFactor ); 
      auxMat.makeTensorProduct( gradi, gradi );
      H.addMultiple( auxMat, mixedFactor );
      localToGlobal( Dest, pi, pi, H );      

    }
  }

protected:
  void localToGlobal( BlockMatrixType& BlockOp, int k, int l, const aol::Matrix33<RealType>& localMatrix ) const {
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          BlockOp.getReference(i,j).add( k, l, _factor * localMatrix.get(i,j) );	
	
    if( k != l)
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          BlockOp.getReference(i,j).add( l, k, _factor * localMatrix.get(j,i) );	
  }
  
};

//==========================================================================================================
//! \brief First derivative of face area energy w.r.t. the undeformed configuration (cf. class FaceAreaEnergy<> above).
//! \author Heeren
template< typename MeshType >
class FaceAreaGradientUndef : public aol::Op< aol::MultiVector<typename MeshType::RealType> > {

  typedef typename MeshType::RealType   RealType;
  typedef aol::MultiVector<RealType>           VectorType;
  typedef aol::Vec3<RealType>                  Vec;

  const MeshTopologySaver<MeshType>& _topology;
  const aol::MultiVector<RealType>&  _defShell;

public:
FaceAreaGradientUndef( const MeshTopologySaver<MeshType>& topology,
                            const VectorType& defShell ) : _topology( topology), _defShell(defShell) {}

  void applyAdd( const VectorType& undefShell, VectorType& Dest ) const {

    assert( _defShell.getTotalSize() == undefShell.getTotalSize() );
    if( (Dest.numComponents() != 3) || ( Dest[0].size() != undefShell[0].size() ) )
      Dest.reallocate( 3, undefShell[0].size() );

    for ( int faceIdx = 0; faceIdx < this->_topology.getNumFaces(); ++faceIdx ){

      int pi( this->_topology.getNodeOfTriangle(faceIdx,0) ),
          pj( this->_topology.getNodeOfTriangle(faceIdx,1) ),
          pk( this->_topology.getNodeOfTriangle(faceIdx,2) );

      // set up vertices and edges
      aol::Vec3<RealType> Pi, Pj, Pk, grad;
      
      _defShell.getTo( pi, Pi ); 
      _defShell.getTo( pj, Pj );
      _defShell.getTo( pk, Pk );
      grad.makeCrossProduct( Pk-Pj, Pi-Pk );      
      RealType defFaceAreaSqr = grad.normSqr() / 4.;      
      
      undefShell.getTo( pi, Pi ); 
      undefShell.getTo( pj, Pj );
      undefShell.getTo( pk, Pk );
      grad.makeCrossProduct( Pk-Pj, Pi-Pk );      
      RealType undefFaceAreaSqr = grad.normSqr() / 4.; 
      
      // pk
      getAreaGradK( Pi, Pj, Pk, grad );
      for( int i = 0; i < 3; i++ )
        Dest[i][pk] +=  ( undefFaceAreaSqr - defFaceAreaSqr ) * grad[i] / undefFaceAreaSqr;
      
      // d_i
      getAreaGradK( Pj, Pk, Pi, grad );
      for( int i = 0; i < 3; i++ )
        Dest[i][pi] +=  ( undefFaceAreaSqr - defFaceAreaSqr ) * grad[i] / undefFaceAreaSqr;
      
      // d_j
      getAreaGradK( Pk, Pi, Pj, grad );
      for( int i = 0; i < 3; i++ )
        Dest[i][pj] +=  ( undefFaceAreaSqr - defFaceAreaSqr ) * grad[i] / undefFaceAreaSqr;

    }
  }

};

//==========================================================================================================
//! \brief Second derivative of face area energy w.r.t. the undeformed configuration (cf. class FaceAreaEnergy<> above).
//! \author Heeren
template <typename MeshType, typename BlockMatrixType = aol::SparseBlockMatrix< aol::SparseMatrix<typename MeshType::RealType> > >
class FaceAreaHessianUndef : public aol::Op<aol::MultiVector<typename MeshType::RealType>, BlockMatrixType > {

protected:
  typedef typename MeshType::RealType RealType;
  typedef aol::Vec3<RealType>                Vec;
  typedef aol::MultiVector<RealType>         VectorType;

  const MeshTopologySaver<MeshType>& _topology;
  const VectorType& _defShell;
  const RealType _factor;

public:
FaceAreaHessianUndef( const MeshTopologySaver<MeshType>& topology,
                           const VectorType& defShell,
			   const RealType Factor = 1. ) : _topology( topology), _defShell(defShell), _factor( Factor ) {}


void applyAdd( const VectorType& undefShell, BlockMatrixType& Dest ) const {
  
  for ( int faceIdx = 0; faceIdx < this->_topology.getNumFaces(); ++faceIdx ){

      int pi( this->_topology.getNodeOfTriangle(faceIdx,0) ),
          pj( this->_topology.getNodeOfTriangle(faceIdx,1) ),
          pk( this->_topology.getNodeOfTriangle(faceIdx,2) );

      // set up vertices and edges
      aol::Vec3<RealType> Pi, Pj, Pk, gradi, gradj, gradk;
      
      _defShell.getTo( pi, Pi ); 
      _defShell.getTo( pj, Pj );
      _defShell.getTo( pk, Pk );
      gradi.makeCrossProduct( Pk-Pj, Pi-Pk );      
      RealType defFaceAreaSqr = gradi.normSqr() / 4.;   
      
      undefShell.getTo( pi, Pi ); 
      undefShell.getTo( pj, Pj );
      undefShell.getTo( pk, Pk );
      gradi.makeCrossProduct( Pk-Pj, Pi-Pk );      
      RealType undefFaceAreaSqr = gradi.normSqr() / 4.; 
      
      RealType hessFactor = (undefFaceAreaSqr - defFaceAreaSqr) / undefFaceAreaSqr;
      RealType mixedFactor = 2. * defFaceAreaSqr / ( undefFaceAreaSqr * std::sqrt(undefFaceAreaSqr) );
      
      // gradients
      getAreaGradK( Pi, Pj, Pk, gradk );
      getAreaGradK( Pj, Pk, Pi, gradi );
      getAreaGradK( Pk, Pi, Pj, gradj );
      
      //
      aol::Matrix33<RealType> H, auxMat;
	
      //*k      
      //kk      
      H.setZero();
      getHessAreaKK( Pi, Pj, Pk, auxMat );
      H.addMultiple( auxMat, hessFactor ); 
      auxMat.makeTensorProduct( gradk, gradk );
      H.addMultiple( auxMat, mixedFactor );
      localToGlobal( Dest, pk, pk, H ); 
      //ik
      H.setZero();     
      getHessAreaIK( Pi, Pj, Pk, auxMat );
      H.addMultiple( auxMat, hessFactor );     
      auxMat.makeTensorProduct( gradi, gradk );
      H.addMultiple( auxMat, mixedFactor );   
      localToGlobal( Dest, pi, pk, H ); 
      //jk
      H.setZero();     
      getHessAreaIK( Pj, Pi, Pk, auxMat );
      H.addMultiple( auxMat, hessFactor );    
      auxMat.makeTensorProduct( gradj, gradk );
      H.addMultiple( auxMat, mixedFactor );          
      localToGlobal( Dest, pj, pk, H ); 
      
      //*j      
      //jj      
      H.setZero();
      getHessAreaKK( Pk, Pi, Pj, auxMat );
      H.addMultiple( auxMat, hessFactor ); 
      auxMat.makeTensorProduct( gradj, gradj );
      H.addMultiple( auxMat, mixedFactor );
      localToGlobal( Dest, pj, pj, H ); 
      //ij
      H.setZero();     
      getHessAreaIK( Pi, Pk, Pj, auxMat );
      H.addMultiple( auxMat, hessFactor );     
      auxMat.makeTensorProduct( gradi, gradj );
      H.addMultiple( auxMat, mixedFactor );   
      localToGlobal( Dest, pi, pj, H ); 
    
      //ii  
      H.setZero();
      getHessAreaKK( Pj, Pk, Pi, auxMat );
      H.addMultiple( auxMat, hessFactor ); 
      auxMat.makeTensorProduct( gradi, gradi );
      H.addMultiple( auxMat, mixedFactor );
      localToGlobal( Dest, pi, pi, H );  

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
//! \brief Mixed second derivative of face area energy w.r.t. the undeformed and deformed configuration (cf. class FaceAreaEnergy above).
//! \author Heeren
template <typename MeshType, typename BlockMatrixType = aol::SparseBlockMatrix< aol::SparseMatrix<typename MeshType::RealType> > >
class FaceAreaHessianMixed : public aol::Op<aol::MultiVector<typename MeshType::RealType>, BlockMatrixType > {

protected:
  typedef typename MeshType::RealType RealType;
  typedef aol::Vec3<RealType>                Vec;
  typedef aol::MultiVector<RealType>         VectorType;

  const MeshTopologySaver<MeshType>& _topology;
  const aol::MultiVector<RealType>&  _inactiveGeometry;
  const bool _activeShellIsDeformed, _firstDerivWRTDef;
  const RealType _factor;

public:
FaceAreaHessianMixed( const MeshTopologySaver<MeshType>& topology,
                           const VectorType& InactiveGeometry,
		           const bool ActiveShellIsDeformed,
			   const bool FirstDerivWRTDef,
			   const RealType Factor = 1. ) : _topology( topology), _inactiveGeometry(InactiveGeometry), _activeShellIsDeformed(ActiveShellIsDeformed), _firstDerivWRTDef( FirstDerivWRTDef ), _factor( Factor ) {}

//
void applyAdd( const VectorType& ActiveGeometry, BlockMatrixType& Dest ) const {
  
    const VectorType* undefShellP = _activeShellIsDeformed ? &_inactiveGeometry : &ActiveGeometry;
    const VectorType* defShellP   = _activeShellIsDeformed ? &ActiveGeometry    : &_inactiveGeometry;
      
  for ( int faceIdx = 0; faceIdx < this->_topology.getNumFaces(); ++faceIdx ){

      int pi( this->_topology.getNodeOfTriangle(faceIdx,0) ),
          pj( this->_topology.getNodeOfTriangle(faceIdx,1) ),
          pk( this->_topology.getNodeOfTriangle(faceIdx,2) );
      
      aol::Vec3<RealType> Pi, Pj, Pk;
      aol::Vec3<RealType> defGradk, defGradi, defGradj;
      aol::Vec3<RealType> undefGradk, undefGradi, undefGradj;
      
      // set up vertices and gradients
      defShellP->getTo( pi, Pi ); 
      defShellP->getTo( pj, Pj );
      defShellP->getTo( pk, Pk );
      defGradi.makeCrossProduct( Pk-Pj, Pi-Pk );      
      RealType factor = defGradi.norm();   
      
      // gradients      
      getAreaGradK( Pi, Pj, Pk, defGradk );
      getAreaGradK( Pj, Pk, Pi, defGradi );
      getAreaGradK( Pk, Pi, Pj, defGradj );
      
      undefShellP->getTo( pi, Pi ); 
      undefShellP->getTo( pj, Pj );
      undefShellP->getTo( pk, Pk );
      undefGradi.makeCrossProduct( Pk-Pj, Pi-Pk );      
      factor *= -4. / undefGradi.normSqr(); 
     
      // gradients      
      getAreaGradK( Pi, Pj, Pk, undefGradk );
      getAreaGradK( Pj, Pk, Pi, undefGradi );
      getAreaGradK( Pk, Pi, Pj, undefGradj );
      
      // scaling
      undefGradk *= factor;
      undefGradi *= factor;
      undefGradj *= factor;      
      
      // k*
      aol::Matrix33<RealType> matrix;
      matrix.makeTensorProduct( undefGradk, defGradk );
      localToGlobal( Dest, pk, pk, matrix );
      matrix.makeTensorProduct( undefGradk, defGradi );
      localToGlobal( Dest, pk, pi, matrix );
      matrix.makeTensorProduct( undefGradk, defGradj );
      localToGlobal( Dest, pk, pj, matrix );
      
      // i*
      matrix.makeTensorProduct( undefGradi, defGradi );
      localToGlobal( Dest, pi, pi, matrix );
      matrix.makeTensorProduct( undefGradi, defGradj );
      localToGlobal( Dest, pi, pj, matrix );
      matrix.makeTensorProduct( undefGradi, defGradk );
      localToGlobal( Dest, pi, pk, matrix );
      
      // j*
      matrix.makeTensorProduct( undefGradj, defGradj );
      localToGlobal( Dest, pj, pj, matrix );
      matrix.makeTensorProduct( undefGradj, defGradi );
      localToGlobal( Dest, pj, pi, matrix );
      matrix.makeTensorProduct( undefGradj, defGradk );
      localToGlobal( Dest, pj, pk, matrix );

    }
  }

protected:
  void localToGlobal( BlockMatrixType& BlockOp, int k, int l, const aol::Matrix33<RealType>& localMatrix ) const {    

    if( _firstDerivWRTDef )   
    {
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          BlockOp.getReference(i,j).add( k, l, _factor * localMatrix.get(i,j) );
    }
    else
    {
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          BlockOp.getReference(i,j).add( l, k, _factor * localMatrix.get(j,i) );
    }
  }
  
};

//! \brief wrapper in order to test mixed second derivatives of FaceAreaEnergy
//! \f$ S_2 \mapsto \partial_1 W[S_1, S_2] \f$
template <typename MeshType>
class FaceAreaGradientUndefWrapper : public aol::Op<aol::MultiVector<typename MeshType::RealType> > {

protected:
  typedef typename MeshType::RealType RealType;
  typedef aol::Vec3<RealType>                Vec;
  typedef aol::MultiVector<RealType>         VectorType;

  const MeshTopologySaver<MeshType>& _topology;
  const VectorType& _undefShell;

public:
FaceAreaGradientUndefWrapper( const MeshTopologySaver<MeshType>& topology,
                                   const VectorType& undefShell ) : _topology( topology), _undefShell(undefShell) {}



void applyAdd( const VectorType& defShell, VectorType& Dest ) const {
  FaceAreaGradientUndef<MeshType>( _topology, defShell ).applyAdd( _undefShell, Dest );
}

};

//! \brief wrapper in order to test mixed second derivatives of FaceAreaEnergy
//! \f$ S_1 \mapsto \partial_2 W[S_1, S_2] \f$
template <typename MeshType>
class FaceAreaGradientDefWrapper : public aol::Op<aol::MultiVector<typename MeshType::RealType> > {

protected:
  typedef typename MeshType::RealType RealType;
  typedef aol::Vec3<RealType>                Vec;
  typedef aol::MultiVector<RealType>         VectorType;

  const MeshTopologySaver<MeshType>& _topology;
  const VectorType& _defShell;

public:
FaceAreaGradientDefWrapper( const MeshTopologySaver<MeshType>& topology,
                              const VectorType& defShell ) : _topology( topology), _defShell(defShell) {}


void applyAdd( const VectorType& UndefShell, VectorType& Dest ) const {
  FaceAreaGradientDef<MeshType>( _topology, UndefShell ).applyAdd( _defShell, Dest );
}

};

//!===============================================================================================================================
//! \brief Deformation energy representing the discrete membrane energy given by the class FaceAreaEnergy<>
//! \author Heeren
template <typename MeshType>
class FaceAreaDeformation : public DeformationBase<MeshType> {

protected:
  typedef typename MeshType::RealType   RealType;
  typedef aol::MultiVector<RealType>           VectorType;
  typedef typename DeformationBase<MeshType>::BlockMatrixType BlockMatrixType;

public:
  FaceAreaDeformation( const MeshTopologySaver<MeshType>& Topology ) : DeformationBase<MeshType>( Topology ) {}
		            
  FaceAreaDeformation( const MeshTopologySaver<MeshType>& Topology, const aol::ParameterParser& /*Pparser*/ ) : DeformationBase<MeshType>( Topology ) {}
  
  void applyEnergy ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, aol::Scalar<RealType>& Dest  ) const {
    FaceAreaEnergy<MeshType>( this->_topology, UndeformedGeom, true ).apply( DeformedGeom, Dest );
  }
  
  void applyUndefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest ) const {
    FaceAreaGradientUndef<MeshType>( this->_topology, DeformedGeom ).apply( UndeformedGeom, Dest );
  }
  
  void applyDefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest ) const {
    FaceAreaGradientDef<MeshType>( this->_topology, UndeformedGeom ).apply( DeformedGeom, Dest );
  }
  
  void assembleAddDefHessian  ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, BlockMatrixType& Dest, RealType factor = 1.0 ) const {
    this->checkIfMatrixIsAllocated( Dest );
    FaceAreaHessianDef<MeshType>( this->_topology, UndeformedGeom, factor ).applyAdd( DeformedGeom, Dest ); 
  }
  
  void assembleAddUndefHessian  ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, BlockMatrixType& Dest, RealType factor = 1.0 ) const {
    this->checkIfMatrixIsAllocated( Dest );
    FaceAreaHessianUndef<MeshType>( this->_topology, DeformedGeom, factor ).applyAdd( UndeformedGeom, Dest ); 
  }
  
  void assembleAddMixedHessian ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, BlockMatrixType& Dest, const bool FirstDerivWRTDef, RealType factor = 1.0 ) const {
    this->checkIfMatrixIsAllocated( Dest );
    bool ActiveShellIsDeformed = true;
    FaceAreaHessianMixed<MeshType>( this->_topology, UndeformedGeom, ActiveShellIsDeformed, FirstDerivWRTDef, factor ).applyAdd( DeformedGeom, Dest );
  }

};

//!==========================================================================================================
//! DISCRETE SHELLS MEMBRANE DEFORMATION ENERGY
//!==========================================================================================================

//! \brief Deformation energy representing the discrete membrane energy in the "Discrete shells" paper by Grinspun et al., 2003.
//! \author Heeren
//! The energy is given as a (weighted) sum of the EdgeLengthDeformation and the FaceAreaDeformation, see documentation thereof above.
//! \f$ E[S, \tilde S] = \mu \sum_t \frac{ (A_t - A_{\tilde t})^2 }{ A_t } + \lambda \sum_e \frac{ (l_e - l_{\tilde e})^2 }{ l_e^2 } A_e \f$
template <typename MeshType>
class DiscreteShellsMembraneDeformation : public DeformationBase<MeshType> {

protected:
  typedef typename MeshType::RealType   RealType;
  typedef aol::MultiVector<RealType>           VectorType;
  typedef typename DeformationBase<MeshType>::BlockMatrixType BlockMatrixType;
  
  FaceAreaDeformation<MeshType> _faceAreaDeformation;
  EdgeLengthDeformation<MeshType> _edgeLengthDeformation;
  
  RealType _lengthWeight, _volWeight;

public:
  DiscreteShellsMembraneDeformation( const MeshTopologySaver<MeshType>& Topology, RealType lengthWeight = 1., RealType volWeight = 1. ) 
    : DeformationBase<MeshType>( Topology ), 
      _faceAreaDeformation( Topology ),
      _edgeLengthDeformation( Topology ),
      _lengthWeight( lengthWeight ),
      _volWeight( volWeight ){}
		            
  DiscreteShellsMembraneDeformation( const MeshTopologySaver<MeshType>& Topology, const aol::ParameterParser& Pparser ) 
    : DeformationBase<MeshType>( Topology ),
      _faceAreaDeformation( Topology ),
      _edgeLengthDeformation( Topology ),
      _lengthWeight( Pparser.getDoubleOrDefault("lengthWeight", 1.0) ),
      _volWeight( Pparser.getDoubleOrDefault("volWeight", 1.0) ){}
  
  void applyEnergy ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, aol::Scalar<RealType>& Dest  ) const {
    Dest.setZero();
    _faceAreaDeformation.applyAddEnergy( UndeformedGeom, DeformedGeom, Dest, _volWeight );
    _edgeLengthDeformation.applyAddEnergy( UndeformedGeom, DeformedGeom, Dest, _lengthWeight );
  }
  
  void applyUndefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest ) const {
    Dest.setZero();
    _faceAreaDeformation.applyAddUndefGradient( UndeformedGeom, DeformedGeom, Dest, _volWeight );
    _edgeLengthDeformation.applyAddUndefGradient( UndeformedGeom, DeformedGeom, Dest, _lengthWeight );
  }
  
  void applyDefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest ) const {
    Dest.setZero();
    _faceAreaDeformation.applyAddDefGradient( UndeformedGeom, DeformedGeom, Dest, _volWeight );
    _edgeLengthDeformation.applyAddDefGradient( UndeformedGeom, DeformedGeom, Dest, _lengthWeight );
  }
  
  void assembleAddDefHessian  ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, BlockMatrixType& Dest, RealType factor = 1.0 ) const {
    this->checkIfMatrixIsAllocated( Dest );
    _faceAreaDeformation.assembleAddDefHessian( UndeformedGeom, DeformedGeom, Dest, _volWeight * factor );
    _edgeLengthDeformation.assembleAddDefHessian( UndeformedGeom, DeformedGeom, Dest, _lengthWeight * factor );
  }
  
  void assembleAddUndefHessian  ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, BlockMatrixType& Dest, RealType factor = 1.0 ) const {
    this->checkIfMatrixIsAllocated( Dest );
    _faceAreaDeformation.assembleAddUndefHessian( UndeformedGeom, DeformedGeom, Dest, _volWeight * factor );
    _edgeLengthDeformation.assembleAddUndefHessian( UndeformedGeom, DeformedGeom, Dest, _lengthWeight * factor );
  }
  
  void assembleAddMixedHessian ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, BlockMatrixType& Dest, const bool FirstDerivWRTDef, RealType factor = 1.0 ) const {
    this->checkIfMatrixIsAllocated( Dest );
    _faceAreaDeformation.assembleAddMixedHessian( UndeformedGeom, DeformedGeom, Dest, FirstDerivWRTDef, _volWeight * factor );
    _edgeLengthDeformation.assembleAddMixedHessian( UndeformedGeom, DeformedGeom, Dest, FirstDerivWRTDef, _lengthWeight * factor );
  }

};
#endif
 
