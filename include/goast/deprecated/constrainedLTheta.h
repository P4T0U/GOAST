// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef CONSTRAINEDLTHETA_H
#define CONSTRAINEDLTHETA_H

#include <iostream>
#include <chrono>
#include <ctime>
#include <string>

#include <goast/Core.h>

//==========================================================================================================
//! \brief For given weights (w_e)_e, prescribed angles (theta_e^*)_e and and augmented Lagrange weight \mu, 
//!        compute E[X] = \sum_e w_e (theta_e(X) - theta_e^*) + \mu \sum_e (theta_e(X) - theta_e^*)^2
//! \author Heeren
template<typename ConfiguratorType>
class DihedralSumEnergy : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::MatType    MatType;

  const MeshTopologySaver& _topology;
  const VectorType&  _prescribedAngles, _weights;
  RealType _augmWeight;

public:

  DihedralSumEnergy( const MeshTopologySaver& topology,
		     const VectorType& Weights,
                     const VectorType& prescribedAngles,
                     RealType AugmentedWeight = 0. ) 
  : _topology( topology), 
    _prescribedAngles(prescribedAngles),
    _weights(Weights),
    _augmWeight(AugmentedWeight){
      if( _weights.size() != _topology.getNumEdges() )
	throw BasicException( "DihedralSumEnergy: weights have wrong size!");
      if( _prescribedAngles.size() != _topology.getNumEdges() )
	throw BasicException( "DihedralSumEnergy: prescribed thetas have wrong size!");
    }

  // energy evaluation
  void apply( const VectorType& Geometry, RealType & Dest ) const {
    VectorType temp;
    evaluateIndivualTerms( Geometry, temp );
    Dest[0] = temp.dot( _weights );
    if( _augmWeight > 1.e-15 )
        Dest[0] += _augmWeight * temp.squaredNorm();
  }
  
  //
  void setWeight( RealType AugmentedWeight ) {
    _augmWeight = AugmentedWeight;    
  }
  
  //
  void evaluateIndivualTerms( const VectorType& Geometry, VectorType& Dest ) const {
     
    if( Geometry.size() != 3 * _topology.getNumVertices() )
      throw BasicException( "DihedralSumEnergy::apply(): geometry has wrong size!");
    
    Dest.resize( _topology.getNumEdges() );    
    Dest.setZero();
    
    for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){
      
      if( !(_topology.isEdgeValid(edgeIdx)) )
	continue;

      int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
          pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ),
          pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );

      // no bending at boundary edges
      if( std::min( pl, pk) < 0 )
        continue;

       // get nodal positions
      VecType Pi, Pj, Pk, Pl;     
      getXYZCoord<VectorType, VecType>( Geometry, Pi, pi);
      getXYZCoord<VectorType, VecType>( Geometry, Pj, pj);
      getXYZCoord<VectorType, VecType>( Geometry, Pk, pk);
      getXYZCoord<VectorType, VecType>( Geometry, Pl, pl);

      RealType Theta = getDihedralAngle( Pi, Pj, Pk, Pl );
      Dest[edgeIdx] =  Theta - _prescribedAngles[edgeIdx];      
    } 
  }

 };

//==========================================================================================================
//! \brief Gradient of DihedralSumEnergy (see above)
//! \author Heeren
template<typename ConfiguratorType>
class DihedralSumGradient : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {

  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::VectorType VectorType;  
  typedef typename ConfiguratorType::VecType     VecType;
  typedef typename ConfiguratorType::MatType     MatType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver& _topology;
  const VectorType&  _prescribedAngles, _weights;
  RealType _augmWeight;
  BitVector _bdryMask;

public:

  DihedralSumGradient( const MeshTopologySaver& topology,
		       const VectorType& Weights,
                       const VectorType& prescribedAngles,
                       RealType AugmentedWeight = 0. ) 
  : _topology( topology), 
    _prescribedAngles(prescribedAngles),
    _weights(Weights),
    _augmWeight(AugmentedWeight),
    _bdryMask( _topology.getNumVertices() ){
      if( _weights.size() != _topology.getNumEdges() )
	throw BasicException( "DihedralSumGradient: weights have wrong size!");
    }
    
  //
  void setWeight( RealType AugmentedWeight ) {
    _augmWeight = AugmentedWeight;    
  }
  
  void setBoundaryMask( const std::vector<int>& mask ) {
    _bdryMask.setAll( false );
    for( int i = 0; i < mask.size(); i++ )
        if( mask[i] < _topology.getNumVertices() )
            _bdryMask.set( mask[i], true );
  }

  void apply(const VectorType& Geometry, VectorType& Dest) const {

    if( Geometry.size() != 3 * _topology.getNumVertices() )
      throw BasicException( "DihedralSumGradient::apply(): geometry has wrong size!");

    if (Dest.size() != Geometry.size())
      Dest.resize(Geometry.size());
    Dest.setZero();

	for (int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx){

	    if (!(_topology.isEdgeValid(edgeIdx)))
				continue;

			int pi(_topology.getAdjacentNodeOfEdge(edgeIdx, 0)),
				pj(_topology.getAdjacentNodeOfEdge(edgeIdx, 1)),
				pk(_topology.getOppositeNodeOfEdge(edgeIdx, 0)),
				pl(_topology.getOppositeNodeOfEdge(edgeIdx, 1));

			// no bending at boundary edges
			if (std::min(pl, pk) < 0)
				continue;

			//! get nodal positions
			VecType Pi, Pj, Pk, Pl, temp;
			getXYZCoord<VectorType, VecType>(Geometry, Pi, pi);
			getXYZCoord<VectorType, VecType>(Geometry, Pj, pj);
			getXYZCoord<VectorType, VecType>(Geometry, Pk, pk);
			getXYZCoord<VectorType, VecType>(Geometry, Pl, pl);

			// compute first derivatives of dihedral angle
			VecType thetak, thetal, thetai, thetaj;
			getThetaGradK(Pi, Pj, Pk, thetak);
			getThetaGradK(Pj, Pi, Pl, thetal);
			getThetaGradI(Pi, Pj, Pk, Pl, thetai);
			getThetaGradJ(Pi, Pj, Pk, Pl, thetaj);
                        
                        // scalar factor
                        RealType factor = _weights[edgeIdx];
                        if( _augmWeight > 1.e-15 ){
                            RealType Theta = getDihedralAngle( Pi, Pj, Pk, Pl );
                            factor += 2. * _augmWeight * (Theta - _prescribedAngles[edgeIdx]);
                        }

			// assemble in global vector
			for (int i = 0; i < 3; i++){
				Dest[i*_topology.getNumVertices() + pi] += factor * thetai[i];
				Dest[i*_topology.getNumVertices() + pj] += factor * thetaj[i];
				Dest[i*_topology.getNumVertices() + pk] += factor * thetak[i];
				Dest[i*_topology.getNumVertices() + pl] += factor * thetal[i];
			}

		}
  }
	
  // Push m x N block, where m = #edges and N = 3 * #vertices.
  // In the ith row the gradient wrt. the ith edge is pushed, i.e. the ith row will have 4 * 3 = 12 nonzero column entries
  // One can specify to offset the block as well as to add the transposed block (subject to a transposed offset).
  // TODO make efficient, reduce code duplication of apply function above!
  void pushTriplets( const VectorType& Geometry, int rowOffset, int colOffset, TripletListType& tripletList, bool pushTransposed = false ) const {
    if( Geometry.size() != 3 * _topology.getNumVertices() )
      throw BasicException( "DihedralSumGradient::pushTriplets(): geometry has wrong size!");

    for (int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx){

	    if (!(_topology.isEdgeValid(edgeIdx)))
              continue;

			int pi(_topology.getAdjacentNodeOfEdge(edgeIdx, 0)),
				pj(_topology.getAdjacentNodeOfEdge(edgeIdx, 1)),
				pk(_topology.getOppositeNodeOfEdge(edgeIdx, 0)),
				pl(_topology.getOppositeNodeOfEdge(edgeIdx, 1));

			// no bending at boundary edges
			if (std::min(pl, pk) < 0)
				continue;

			//! get nodal positions
			VecType Pi, Pj, Pk, Pl, temp;
			getXYZCoord<VectorType, VecType>(Geometry, Pi, pi);
			getXYZCoord<VectorType, VecType>(Geometry, Pj, pj);
			getXYZCoord<VectorType, VecType>(Geometry, Pk, pk);
			getXYZCoord<VectorType, VecType>(Geometry, Pl, pl);

			// compute first derivatives of dihedral angle
			VecType thetak, thetal, thetai, thetaj;
			getThetaGradK(Pi, Pj, Pk, thetak);
			getThetaGradK(Pj, Pi, Pl, thetal);
			getThetaGradI(Pi, Pj, Pk, Pl, thetai);
			getThetaGradJ(Pi, Pj, Pk, Pl, thetaj);
                        
                        // scalar factor
                        RealType factor = _weights[edgeIdx];
                        if( _augmWeight > 1.e-15 ){
                            RealType Theta = getDihedralAngle( Pi, Pj, Pk, Pl );
                            factor += 2. * _augmWeight * (Theta - _prescribedAngles[edgeIdx]);
                        }

			// push triplets
			for (int i = 0; i < 3; i++){
                            if( !_bdryMask[pi] )
				tripletList.push_back( TripletType( rowOffset + edgeIdx, colOffset + i*_topology.getNumVertices() + pi, factor * thetai[i] ) );
                            if( !_bdryMask[pj] )
				tripletList.push_back( TripletType( rowOffset + edgeIdx, colOffset + i*_topology.getNumVertices() + pj, factor * thetaj[i] ) );
                            if( !_bdryMask[pk] )
				tripletList.push_back( TripletType( rowOffset + edgeIdx, colOffset + i*_topology.getNumVertices() + pk, factor * thetak[i] ) );
                            if( !_bdryMask[pl] )
				tripletList.push_back( TripletType( rowOffset + edgeIdx, colOffset + i*_topology.getNumVertices() + pl, factor * thetal[i] ) );
			}
			
			if( !pushTransposed )
                            continue;
                        
                        // push transposed triplets
			for (int i = 0; i < 3; i++){
                            if( !_bdryMask[pi] )
				tripletList.push_back( TripletType( colOffset + i*_topology.getNumVertices() + pi, rowOffset + edgeIdx, factor * thetai[i] ) );
                            if( !_bdryMask[pj] )
				tripletList.push_back( TripletType( colOffset + i*_topology.getNumVertices() + pj, rowOffset + edgeIdx, factor * thetaj[i] ) );
                            if( !_bdryMask[pk] )
				tripletList.push_back( TripletType( colOffset + i*_topology.getNumVertices() + pk, rowOffset + edgeIdx, factor * thetak[i] ) );
                            if( !_bdryMask[pl] )
				tripletList.push_back( TripletType( colOffset + i*_topology.getNumVertices() + pl, rowOffset + edgeIdx, factor * thetal[i] ) );
			}

    }
  }	
};

//==========================================================================================================
//! \brief Hessian of DihedralSumEnergy (see above)
//! \author Heeren
template<typename ConfiguratorType>
class DihedralSumHessian : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

protected:
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::VecType     VecType;
  typedef typename ConfiguratorType::MatType     MatType;
  
  typedef std::vector<TripletType> TripletListType;
  
  const MeshTopologySaver& _topology;
  const VectorType&  _prescribedAngles, _weights;
  RealType _augmWeight, _factor;
  mutable int _rowOffset, _colOffset;
  
public:
  DihedralSumHessian( const MeshTopologySaver& topology,
		       const VectorType& Weights,
                       const VectorType& prescribedAngles,
                       const RealType AugmentedWeight = 0.,
                       const RealType Factor = 1.,
                       int rowOffset = 0, 
                       int colOffset = 0 ) : 
    _topology( topology),     
    _prescribedAngles(prescribedAngles),
    _weights(Weights),
    _augmWeight(AugmentedWeight),
    _factor( Factor ), 
    _rowOffset(rowOffset), 
    _colOffset(colOffset){}
    
  //
  void setWeight( RealType AugmentedWeight ) {
    _augmWeight = AugmentedWeight;    
  }  
    
  void setRowOffset( int rowOffset ) const {
        _rowOffset = rowOffset;
  }
    
  void setColOffset( int colOffset ) const {
        _colOffset = colOffset;
  }

  //
  void apply( const VectorType& Geometry, MatrixType& Dest ) const {    
    assembleHessian( Geometry, Dest );
  }

  // assmeble Hessian matrix via triplet list
  void assembleHessian( const VectorType& Geometry, MatrixType& Hessian ) const {
    int dofs = 3*_topology.getNumVertices();
    if( (Hessian.rows() != dofs) || (Hessian.cols() != dofs) )
        Hessian.resize( dofs, dofs );
    Hessian.setZero();
    
    // set up triplet list
    TripletListType tripletList;
    // per edge we have 4 active vertices, i.e. 16 combinations each producing a 3x3-matrix
    tripletList.reserve( 16 * 9 * _topology.getNumEdges() );   
    // fill matrix from triplets
    pushTriplets( Geometry, tripletList );
    Hessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
  
  // fill triplets
  void pushTriplets( const VectorType& Geometry, TripletListType& tripletList ) const {

    if( Geometry.size() != 3 * _topology.getNumVertices() )
      throw BasicException( "DihedralSumHessian::apply(): geometry has wrong size!");
  
    // run over all edges and fill triplets
    for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){
      
      if( !(_topology.isEdgeValid(edgeIdx)) )
	continue;

      int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
          pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ),
          pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );

      // no bending at boundary edges
      if( std::min( pl, pk) < 0 )
        continue;

      //! get nodal positions
      VecType Pi, Pj, Pk, Pl;
      getXYZCoord<VectorType, VecType>( Geometry, Pi, pi);
      getXYZCoord<VectorType, VecType>( Geometry, Pj, pj);
      getXYZCoord<VectorType, VecType>( Geometry, Pk, pk);
      getXYZCoord<VectorType, VecType>( Geometry, Pl, pl);

      // now compute second derivatives of dihedral angle
      MatType H;
      
      // scalar factor
      RealType HessFactor = _weights[edgeIdx];
      if( _augmWeight > 1.e-15 ){
        RealType Theta = getDihedralAngle( Pi, Pj, Pk, Pl );
        HessFactor += 2. * _augmWeight * (Theta - _prescribedAngles[edgeIdx]);
      }
      
      //kk
      getHessThetaKK( Pi, Pj, Pk, H );
      localToGlobal( tripletList, pk, pk, H, HessFactor );
            
      //ik & ki (Hki = Hik)
      getHessThetaIK( Pi, Pj, Pk, H);
      localToGlobal( tripletList, pi, pk, H, HessFactor );
      
      //jk & kj (Hkj = Hjk)
      getHessThetaJK( Pi, Pj, Pk, H );
      localToGlobal( tripletList, pj, pk, H, HessFactor );      
      
      //ll
      getHessThetaKK( Pj, Pi, Pl, H );
      localToGlobal( tripletList, pl, pl, H, HessFactor );     
      
      //il & li (Hli = Hil)
      getHessThetaJK( Pj, Pi, Pl, H );
      localToGlobal( tripletList, pi, pl, H, HessFactor );
      
      //jl & lj (Hlj = Hjl)
      getHessThetaIK( Pj, Pi, Pl, H );
      localToGlobal( tripletList, pj, pl, H, HessFactor );           
      
      //kl/lk: Hkl = 0 and Hlk = 0
        
      //ii  
      getHessThetaII( Pi, Pj, Pk, Pl, H );
      localToGlobal( tripletList, pi, pi, H, HessFactor );             

      //jj
      getHessThetaII( Pj, Pi, Pl, Pk, H );       
      localToGlobal( tripletList, pj, pj, H, HessFactor );

      //ij & ji (Hij = Hji)
      getHessThetaJI( Pi, Pj, Pk, Pl, H );     
      localToGlobal( tripletList, pi, pj, H, HessFactor ); 

      
      // add mixed terms
      if( _augmWeight < 1.e-15 )
          continue;
      
      // compute gradients
      VecType thetak, thetal, thetai, thetaj;
      getThetaGradK(Pi, Pj, Pk, thetak);
      getThetaGradK(Pj, Pi, Pl, thetal);
      getThetaGradI(Pi, Pj, Pk, Pl, thetai);
      getThetaGradJ(Pi, Pj, Pk, Pl, thetaj);
            
      //kk
      H.makeTensorProduct( thetak, thetak );
      localToGlobal( tripletList, pk, pk, H, 2. * _augmWeight );
            
      //ik & ki (Hki = Hik)
      H.makeTensorProduct( thetai, thetak );
      localToGlobal( tripletList, pi, pk, H, 2. * _augmWeight );
      
      //jk & kj (Hkj = Hjk)
      H.makeTensorProduct( thetaj, thetak );
      localToGlobal( tripletList, pj, pk, H, 2. * _augmWeight );      
            
      //kl & lk
      H.makeTensorProduct( thetal, thetak );
      localToGlobal( tripletList, pl, pk, H, 2. * _augmWeight );
      
      //ll
      H.makeTensorProduct( thetal, thetal );
      localToGlobal( tripletList, pl, pl, H, 2. * _augmWeight );     
      
      //il & li (Hli = Hil)
      H.makeTensorProduct( thetai, thetal );
      localToGlobal( tripletList, pi, pl, H, 2. * _augmWeight );
      
      //jl & lj (Hlj = Hjl)
      H.makeTensorProduct( thetaj, thetal );
      localToGlobal( tripletList, pj, pl, H, 2. * _augmWeight );        
        
      //ii  
      H.makeTensorProduct( thetai, thetai );
      localToGlobal( tripletList, pi, pi, H, 2. * _augmWeight );             

      //jj
      H.makeTensorProduct( thetaj, thetaj );       
      localToGlobal( tripletList, pj, pj, H, 2. * _augmWeight );

      //ij & ji (Hij = Hji)
      H.makeTensorProduct( thetai, thetaj );     
      localToGlobal( tripletList, pi, pj, H, 2. * _augmWeight ); 
      
    }
  }


protected:
  void localToGlobal( TripletListType& tripletList, int k, int l, const MatType& localMatrix, RealType Weight ) const {
    int numV = _topology.getNumVertices();
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          tripletList.push_back( TripletType( _rowOffset + i*numV + k, _colOffset + j*numV + l, _factor * Weight * localMatrix(i,j) ) );	
	
    if( k != l){
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          tripletList.push_back( TripletType( _rowOffset + i*numV + l, _colOffset + j*numV + k, _factor * Weight * localMatrix(j,i) ) );
    }
  } 
  
};
 
//==========================================================================================================
//==========================================================================================================
//==========================================================================================================

//==========================================================================================================
//! \brief For given weights (w_e)_e, prescribed lengths (l_e^*)_e and and augmented Lagrange weight \mu,
//!        compute E[X] = \sum_e w_e (l_e(X) - l_e^*) + \mu \sum_e (l_e(X) - l_e^*)^2
//! \author Heeren
template< typename ConfiguratorType>
class EdgeLengthSumEnergy : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::MatType    MatType;

  const MeshTopologySaver& _topology;
  const VectorType&  _prescribedLengths, _weights;
  RealType _augmWeight;

public:

  EdgeLengthSumEnergy( const MeshTopologySaver& topology,
		       const VectorType& Weights,
                       const VectorType& prescribedLengths,
                       RealType AugmentedWeight = 0. ) 
  : _topology( topology), 
    _prescribedLengths(prescribedLengths),
    _weights(Weights),
    _augmWeight(AugmentedWeight){
      if( _weights.size() != _topology.getNumEdges() )
	throw BasicException( "EdgeLengthSumEnergy: weights have wrong size!");
      if( _prescribedLengths.size() != _topology.getNumEdges() )
	throw BasicException( "EdgeLengthSumEnergy: prescribed thetas have wrong size!");
    }

  //
  void setWeight( RealType AugmentedWeight ) {
    _augmWeight = AugmentedWeight;    
  }  
    
  // energy evaluation
  void apply( const VectorType& Geometry, RealType & Dest ) const {
    VectorType temp;
    evaluateIndivualTerms( Geometry, temp );
    Dest[0] = temp.dot( _weights );
    if( _augmWeight > 1.e-15 )
        Dest[0] += _augmWeight * temp.squaredNorm();
  }
  
  // returns vector (l_e(X) - l_e^*)_e, where (l_e^*)_e is the vector of prescribed edge lengths
  void evaluateIndivualTerms( const VectorType& Geometry, VectorType& Dest ) const {

    if( Geometry.size() != 3 * _topology.getNumVertices() )
      throw BasicException( "EdgeLengthSumEnergy::apply(): geometry has wrong size!");
    
    Dest.resize( _topology.getNumEdges() );
    Dest.setZero();
    
    for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){

      int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) );

      VecType Pi, Pj;
      getXYZCoord<VectorType, VecType>( Geometry, Pi, pi);
      getXYZCoord<VectorType, VecType>( Geometry, Pj, pj);
      RealType EdgeLength = std::sqrt( dotProduct(Pj-Pi,Pj-Pi) );      
      Dest[edgeIdx] =  EdgeLength - _prescribedLengths[edgeIdx]; 
    }
  }

 };
 

//==========================================================================================================
//! \brief Gradient of EdgeLengthSumEnergy (see above)
//! \author Heeren
template< typename ConfiguratorType>
class EdgeLengthSumGradient : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {

  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::MatType    MatType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver& _topology;
  const VectorType&  _prescribedLengths, _weights;
  RealType _augmWeight;
  BitVector _bdryMask;

public:

  EdgeLengthSumGradient( const MeshTopologySaver& topology,
		         const VectorType& Weights,
                         const VectorType& prescribedLengths,
                         RealType AugmentedWeight = 0. ) 
  : _topology( topology), 
    _prescribedLengths(prescribedLengths),
    _weights(Weights),
    _augmWeight(AugmentedWeight),
    _bdryMask( _topology.getNumVertices() ){
      if( _weights.size() != _topology.getNumEdges() )
	throw BasicException( "EdgeLengthSumGradient: weights have wrong size!");
    }
    
      //
  void setWeight( RealType AugmentedWeight ) {
    _augmWeight = AugmentedWeight;    
  }  
  
  void setBoundaryMask( const std::vector<int>& mask ) {
    _bdryMask.setAll( false );
    for( int i = 0; i < mask.size(); i++ )
        if( mask[i] < _topology.getNumVertices() )
            _bdryMask.set( mask[i], true );
  }

  void apply(const VectorType& Geometry, VectorType& Dest) const {

    if( Geometry.size() != 3 * _topology.getNumVertices() )
      throw BasicException( "EdgeLengthSumGradient::apply(): geometry has wrong size!");

    if (Dest.size() != Geometry.size())
      Dest.resize(Geometry.size());
    Dest.setZero();

	for (int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx){

      int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) );

      // set up vertices and edges
      VecType Pi, edge;
      getXYZCoord<VectorType, VecType>( Geometry, Pi, pi);
      getXYZCoord<VectorType, VecType>( Geometry, edge, pj);   
      edge -= Pi;
      RealType EdgeLength = std::sqrt( dotProduct(edge, edge) );
      RealType factor = ( _weights[edgeIdx] + 2. * _augmWeight * (EdgeLength - _prescribedLengths[edgeIdx]) ) / EdgeLength;  
      
      // assemble in global matrix
      for( int i = 0; i < 3; i++ ){
        Dest[i * _topology.getNumVertices() + pi] -= factor * edge[i];
        Dest[i * _topology.getNumVertices() + pj] += factor * edge[i];
      }

    }
  }
  
  // Push m x N block, where m = #edges and N = 3 * #vertices.
  // In the ith row the gradient wrt. the ith edge is pushed, i.e. the ith row will have 2 * 3 = 6 nonzero column entries
  // One can specify to offset the block as well as to add the transposed block (subject to a transposed offset).
  // TODO make efficient, reduce code duplication of apply function above!
  void pushTriplets( const VectorType& Geometry, int rowOffset, int colOffset, TripletListType& tripletList, bool pushTransposed = false ) const {
    if( Geometry.size() != 3 * _topology.getNumVertices() )
      throw BasicException( "EdgeLengthSumGradient::pushTriplets(): geometry has wrong size!");

    for (int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx){

      int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) );

      // set up vertices and edges
      VecType Pi, edge;
      getXYZCoord<VectorType, VecType>( Geometry, Pi, pi);
      getXYZCoord<VectorType, VecType>( Geometry, edge, pj);   
      edge -= Pi;
      RealType EdgeLength = std::sqrt( dotProduct(edge, edge) );
      RealType factor = ( _weights[edgeIdx] + 2. * _augmWeight * (EdgeLength - _prescribedLengths[edgeIdx]) ) / EdgeLength;  
      
      // push triplets
      for( int i = 0; i < 3; i++ ){
        if( !_bdryMask[pi] )
            tripletList.push_back( TripletType( rowOffset + edgeIdx, colOffset + i * _topology.getNumVertices() + pi, -1. * factor * edge[i] ) );
        if( !_bdryMask[pj] )
            tripletList.push_back( TripletType( rowOffset + edgeIdx, colOffset + i * _topology.getNumVertices() + pj,  1. * factor * edge[i] ) );
      }
      
      if(!pushTransposed)
          continue;
      
      // push transposed triplets
      for( int i = 0; i < 3; i++ ){
          if( !_bdryMask[pi] )
              tripletList.push_back( TripletType( colOffset + i * _topology.getNumVertices() + pi, rowOffset + edgeIdx, -1. * factor * edge[i] ) );
          if( !_bdryMask[pj] )
              tripletList.push_back( TripletType( colOffset + i * _topology.getNumVertices() + pj, rowOffset + edgeIdx,  1. * factor * edge[i] ) );
      }

    }
    
  }

};


//==========================================================================================================
//! \brief Hessian of EdgeLengthSumEnergy (see above)
//! \author Heeren
template <typename ConfiguratorType>
class EdgeLengthSumHessian :  public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

protected:
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::VecType     VecType;
  typedef typename ConfiguratorType::MatType     MatType;
  
  typedef std::vector<TripletType> TripletListType;
  
  const MeshTopologySaver& _topology;
  const VectorType&  _prescribedLengths, _weights;
  RealType _augmWeight, _factor;
  mutable int _rowOffset, _colOffset;
  
public:
  EdgeLengthSumHessian( const MeshTopologySaver& topology,
		       const VectorType& Weights,
                       const VectorType& prescribedLengths,
                       const RealType AugmentedWeight = 0.,
                       const RealType Factor = 1.,
                       int rowOffset = 0, 
                       int colOffset = 0 ) : 
    _topology( topology),     
    _prescribedLengths(prescribedLengths),
    _weights(Weights),
    _augmWeight(AugmentedWeight),
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
  void setWeight( RealType AugmentedWeight ) {
    _augmWeight = AugmentedWeight;    
  }  

  //
  void apply( const VectorType& Geometry, MatrixType& Dest ) const {    
    assembleHessian( Geometry, Dest );
  }

  // assmeble Hessian matrix via triplet list
  void assembleHessian( const VectorType& Geometry, MatrixType& Hessian ) const {
    int dofs = 3*_topology.getNumVertices();
    if( (Hessian.rows() != dofs) || (Hessian.cols() != dofs) )
        Hessian.resize( dofs, dofs );
    Hessian.setZero();
    
    // set up triplet list
    TripletListType tripletList;
    // per edge we have 4 active vertices, i.e. 16 combinations each producing a 3x3-matrix
    tripletList.reserve( 16 * 9 * _topology.getNumEdges() );   
    // fill matrix from triplets
    pushTriplets( Geometry, tripletList );
    Hessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
  
  // fill triplets
  void pushTriplets( const VectorType& Geometry, TripletListType& tripletList ) const {

    for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){

      int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) );

      // set up vertices and edges
      VecType Pi, edge;
      getXYZCoord<VectorType, VecType>( Geometry, Pi, pi);
      getXYZCoord<VectorType, VecType>( Geometry, edge, pj);   
      edge -= Pi;
      RealType EdgeLengthSqr = dotProduct(edge, edge); 
      RealType EdgeLength    = std::sqrt( EdgeLengthSqr );  
      edge /= EdgeLength;

      // now compute second derivatives of dihedral angle
      MatType tensorProduct;
      tensorProduct.makeTensorProduct( edge, edge );      
      RealType factor = (_weights[edgeIdx] + 2. * _augmWeight * (EdgeLength - _prescribedLengths[edgeIdx]) ) / EdgeLength;
      tensorProduct *= (2. * _augmWeight - factor);
      tensorProduct.addToDiagonal( factor );
      
      //ii        
      localToGlobal( tripletList, pi, pi, tensorProduct, 1. );     
      localToGlobal( tripletList, pj, pj, tensorProduct, 1. ); 

      //ij & ji (Hij = Hji)
      localToGlobal( tripletList, pi, pj, tensorProduct, -1. );  
    }
  }

protected:
  void localToGlobal( TripletListType& tripletList, int k, int l, const MatType& localMatrix, RealType Weight ) const {
    int numV = _topology.getNumVertices();
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          tripletList.push_back( TripletType( _rowOffset + i*numV + k, _colOffset + j*numV + l, _factor * Weight * localMatrix(i,j) ) );	
	
    if( k != l){
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          tripletList.push_back( TripletType( _rowOffset + i*numV + l, _colOffset + j*numV + k, _factor * Weight * localMatrix(j,i) ) );
    }
  } 
  
};


//==========================================================================================================
//==========================================================================================================
//==========================================================================================================


//==========================================================================================================
//! \brief For given weights (w_e)_e and prescribed quantities (z_e^*)_e compute E[Z] = \sum_e w_e (z_e - z_e^*)^2
//! \author Heeren
template<typename ConfiguratorType>
class QuadraticMismatchFunctional : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::MatType    MatType;

  const int _numDOFs;
  const VectorType&  _prescribedQuantities, _weights;

public:

  QuadraticMismatchFunctional( const VectorType& Weights, const VectorType& prescribedQuantities ) 
  : _numDOFs(  prescribedQuantities.size() ), 
    _prescribedQuantities(prescribedQuantities),
    _weights(Weights){
      if( _weights.size() != _numDOFs )
	throw BasicException( "QuadraticMismatchFunctional: weights have wrong size!");
    }

  // energy evaluation
  void apply( const VectorType& Quantities, RealType & Dest ) const {

    if( Quantities.size() != _numDOFs )
      throw BasicException( "QuadraticMismatchFunctional::apply(): argument has wrong size!");
    
    Dest.setZero();
    
    for ( int k = 0; k < _numDOFs; ++k )      
      Dest[0] +=  0.5 * _weights[k] * (Quantities[k] - _prescribedQuantities[k]) * (Quantities[k] - _prescribedQuantities[k]); 
 
  }

 };

//==========================================================================================================
//! \brief Gradient of QuadraticMismatchFunctional (see above)
//! \author Heeren
template<typename ConfiguratorType>
class QuadraticMismatchGradient : public BaseOp<typename ConfiguratorType::VectorType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::MatType    MatType;

  const int _numDOFs;
  const VectorType&  _prescribedQuantities, _weights;

public:

  QuadraticMismatchGradient( const VectorType& Weights, const VectorType& prescribedQuantities ) 
  : _numDOFs(  prescribedQuantities.size() ), 
    _prescribedQuantities(prescribedQuantities),
    _weights(Weights){
      if( _weights.size() != _numDOFs )
	throw BasicException( "QuadraticMismatchGradient: weights have wrong size!");
    }

  // energy evaluation
  void apply( const VectorType& Quantities, VectorType& Gradient ) const {

    if( Quantities.size() != _numDOFs )
      throw BasicException( "QuadraticMismatchGradient::apply(): argument has wrong size!");
    
    if( Gradient.size() != _numDOFs )
        Gradient.resize( _numDOFs );
    
     for ( int k = 0; k < _numDOFs; ++k )      
       Gradient[k] = _weights[k] * (Quantities[k] - _prescribedQuantities[k]); 
   
  }

 };
 
//==========================================================================================================
//! \brief Hessian of QuadraticMismatchFunctional (see above)
//! \author Heeren
template<typename ConfiguratorType>
class QuadraticMismatchHessian : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::MatType    MatType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  const int _numDOFs;
  const VectorType&   _weights;
  RealType _factor;
  mutable int _rowOffset, _colOffset;
  mutable MatrixType _constHessian;
  
public:

  QuadraticMismatchHessian( const VectorType& Weights,
		            const RealType Factor = 1.,
                            int rowOffset = 0, 
                            int colOffset = 0 ) : _numDOFs( Weights.size() ), _weights(Weights), _factor( Factor ), _rowOffset(rowOffset), _colOffset(colOffset){}
    
  void setRowOffset( int rowOffset ) const {
        _rowOffset = rowOffset;
  }
    
  void setColOffset( int colOffset ) const {
        _colOffset = colOffset;
  }

  //
  void apply( const VectorType& Quantities, MatrixType& Dest ) const {        
    if( Quantities.size() != _numDOFs )
      throw BasicException( "QuadraticMismatchHessian::apply(): argument has wrong size!");
    if( _constHessian.rows() != _numDOFs )  
      assembleHessian( );
    Dest.resize( _numDOFs, _numDOFs );
    Dest = _constHessian;
  }
  
  void pushTriplets( TripletListType &tripletList ) const {
      for( int i = 0; i < _numDOFs; i++ )
        tripletList.push_back( TripletType( _rowOffset + i, _colOffset + i, _factor * _weights[i]) );
  }
  
protected:    
  // assmeble constant Hessian 
  void assembleHessian() const { 
    _constHessian.resize( _numDOFs, _numDOFs );    
    // set up triplet list
    TripletListType tripletList;
    // Hessian is diagonal matrix
    tripletList.reserve( _numDOFs );   
    for( int i = 0; i < _numDOFs; i++ )
        tripletList.push_back( TripletType(i,i, _factor * _weights[i]) );
    _constHessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }

 };
 
//==========================================================================================================
//==========================================================================================================
//==========================================================================================================

//==========================================================================================================
//! \brief L[X,Z,Lambda] = E[Z] + Lambda * F[X,Z], where E is (sum of different) quadratice edge based matching energies and F is sum of DihedralSumEnergy and EdgeLengthSumEnergy
//! \author Heeren
template<typename ConfiguratorType>
class ConstrainedLThetaLagrangeFunctional : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::MatType    MatType;

  const MeshTopologySaver& _topology;
  int _numGeomDofs, _dimLTheta, _numDOFs;
  const VectorType& _prescribedLengths, _prescribedAngles;
  const VectorType& _lengthsWeights, _angleWeights; 
  
public:
  ConstrainedLThetaLagrangeFunctional( const MeshTopologySaver& topology,
                                       const VectorType& prescribedLengths,
                                       const VectorType& prescribedAngles,
                                       const VectorType& lengthsWeights,
                                       const VectorType& angleWeights )  
  : _topology(topology), 
    _numGeomDofs( 3 * _topology.getNumVertices() ), 
    _dimLTheta( 2 * _topology.getNumEdges() ),
    _numDOFs( _numGeomDofs + 2 * _dimLTheta ),
    _prescribedLengths(prescribedLengths),
    _prescribedAngles(prescribedAngles),
    _lengthsWeights(lengthsWeights),
    _angleWeights(angleWeights){  }

  // argument is (X,Z,Lambda), where X is geometry, Z is LTheta variable and Lambda are Lagrange multipliers
  void apply( const VectorType& Argument, RealType & Dest ) const {
    if( Argument.size() != _numDOFs )
      throw BasicException( "ConstrainedLThetaLagrangeFunctional::apply(): argument has wrong size!");
    
    // get segments (TODO make efficient!)
    int numEdges = _topology.getNumEdges();
    VectorType geometry = Argument.segment(0, _numGeomDofs );
    VectorType lengths  = Argument.segment( _numGeomDofs, numEdges );
    VectorType angles   = Argument.segment( _numGeomDofs + numEdges, numEdges );
    VectorType lambda1  = Argument.segment( _numGeomDofs + _dimLTheta, numEdges );
    VectorType lambda2  = Argument.segment( _numGeomDofs + _dimLTheta + numEdges, numEdges );
    
    Dest.setZero();    
    
    // add LTheta mismatch
    QuadraticMismatchFunctional<ConfiguratorType>( _lengthsWeights, _prescribedLengths ).applyAdd( lengths, Dest );
    QuadraticMismatchFunctional<ConfiguratorType>( _angleWeights, _prescribedAngles ).applyAdd( angles, Dest );
    
    // add constraints
    EdgeLengthSumEnergy<ConfiguratorType>( _topology, lambda1, lengths).applyAdd( geometry, Dest );
    DihedralSumEnergy<ConfiguratorType>( _topology, lambda2, angles).applyAdd( geometry, Dest );
  }

};
 

//==========================================================================================================
//! \brief Gradient of ConstrainedLThetaLagrangeFunctional (see above)
//! \author Heeren
template<typename ConfiguratorType>
class ConstrainedLThetaLagrangeGradient : public BaseOp<typename ConfiguratorType::VectorType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::MatType    MatType;

  const MeshTopologySaver& _topology;
  int _numGeomDofs, _dimLTheta, _numDOFs;
  const VectorType& _prescribedLengths, _prescribedAngles;
  const VectorType& _lengthsWeights, _angleWeights;
  
  
public:
  ConstrainedLThetaLagrangeGradient( const MeshTopologySaver& topology,
                                     const VectorType& prescribedLengths,
                                     const VectorType& prescribedAngles,
                                     const VectorType& lengthsWeights,
                                     const VectorType& angleWeights )  
  : _topology(topology), 
    _numGeomDofs( 3 * _topology.getNumVertices() ), 
    _dimLTheta( 2 * _topology.getNumEdges() ),
    _numDOFs( _numGeomDofs + 2 * _dimLTheta ),
    _prescribedLengths(prescribedLengths),
    _prescribedAngles(prescribedAngles),
    _lengthsWeights(lengthsWeights),
    _angleWeights(angleWeights){  }

  // argument is (X,Z,Lambda), where X is geometry, Z is LTheta variable and Lambda are Lagrange multipliers
  void apply( const VectorType& Argument, VectorType& Gradient ) const {
    if( Argument.size() != _numDOFs )
      throw BasicException( "ConstrainedLThetaLagrangeGradient::apply(): argument has wrong size!");
    if( Gradient.size() != _numDOFs )
        Gradient.resize( _numDOFs );
    Gradient.setZero();
    
    // get segments (TODO make efficient!)
    int numEdges = _topology.getNumEdges();
    VectorType geometry = Argument.segment(0, _numGeomDofs );
    VectorType lengths  = Argument.segment( _numGeomDofs, numEdges );
    VectorType angles   = Argument.segment( _numGeomDofs + numEdges, numEdges );
    VectorType lambda1  = Argument.segment( _numGeomDofs + _dimLTheta, numEdges );
    VectorType lambda2  = Argument.segment( _numGeomDofs + _dimLTheta + numEdges, numEdges );
    
    VectorType partialGrad;
    
    // gradient wrt. geometry
    EdgeLengthSumGradient<ConfiguratorType>( _topology, lambda1, lengths ).apply( geometry, partialGrad );
    DihedralSumGradient<ConfiguratorType>( _topology, lambda2, angles ).applyAdd( geometry, partialGrad );
    for( int i = 0; i < _numGeomDofs; i++ )
        Gradient[i] = partialGrad[i];
    
    // gradient wrt. ltheta
    QuadraticMismatchGradient<ConfiguratorType>( _lengthsWeights, _prescribedLengths ).apply( lengths, partialGrad);
    for( int i = 0; i < numEdges; i++ )
        Gradient[_numGeomDofs + i] = partialGrad[i] - lambda1[i];
    QuadraticMismatchGradient<ConfiguratorType>( _angleWeights, _prescribedAngles ).apply( angles, partialGrad );
    for( int i = 0; i < numEdges; i++ )
        Gradient[_numGeomDofs + numEdges + i] = partialGrad[i] - lambda2[i];
    
    // gradient wrt. multipliers
    EdgeLengthSumEnergy<ConfiguratorType>( _topology, lambda1, lengths).evaluateIndivualTerms( geometry, partialGrad );
    for( int i = 0; i < numEdges; i++ )
        Gradient[_numGeomDofs + _dimLTheta + i] = partialGrad[i];
    DihedralSumEnergy<ConfiguratorType>( _topology, lambda2, angles).evaluateIndivualTerms( geometry, partialGrad );
    for( int i = 0; i < numEdges; i++ )
        Gradient[_numGeomDofs + _dimLTheta + numEdges + i] = partialGrad[i];
  }

};

//==========================================================================================================
//! \brief Gradient of ConstrainedLThetaLagrangeFunctional (see above)
//! \author Heeren
template<typename ConfiguratorType>
class ConstrainedLThetaLagrangeHessian : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::MatType    MatType;

  const MeshTopologySaver& _topology;
  int _numGeomDofs, _dimLTheta, _numDOFs;
  const VectorType& _prescribedLengths, _prescribedAngles;
  const VectorType& _lengthsWeights, _angleWeights;
  std::vector<int> _bdryMask;
  
public:
  ConstrainedLThetaLagrangeHessian( const MeshTopologySaver& topology,
                                     const VectorType& prescribedLengths,
                                     const VectorType& prescribedAngles,
                                     const VectorType& lengthsWeights,
                                     const VectorType& angleWeights )  
  : _topology(topology), 
    _numGeomDofs( 3 * _topology.getNumVertices() ), 
    _dimLTheta( 2 * _topology.getNumEdges() ),
    _numDOFs( _numGeomDofs + 2 * _dimLTheta ),
    _prescribedLengths(prescribedLengths),
    _prescribedAngles(prescribedAngles),
    _lengthsWeights(lengthsWeights),
    _angleWeights(angleWeights){  }
    
    
  void apply( const VectorType& Argument, MatrixType& Hessian ) const {
    if( Argument.size() != _numDOFs )
      throw BasicException( "ConstrainedLThetaLagrangeHessian::apply(): argument has wrong size!");
    Hessian.resize( _numDOFs, _numDOFs );
    Hessian.setZero();
    
    // set up triplet list
    TripletListType tripletList;
    // fill matrix from triplets
    pushTriplets( Argument, tripletList );
    Hessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
  
  //
  void setBoundaryMask( const std::vector<int>& mask ) {
    _bdryMask.resize( mask.size() );
    _bdryMask = mask;    
  }
  
  // fill triplets
  void pushTriplets( const VectorType& Argument, TripletListType& tripletList ) const {
    if( Argument.size() != _numDOFs )
      throw BasicException( "ConstrainedLThetaLagrangeHessian::pushTriplets(): argument has wrong size!");
    
    // get segments (TODO make efficient!)
    int numEdges = _topology.getNumEdges();
    VectorType geometry = Argument.segment(0, _numGeomDofs );
    VectorType lengths  = Argument.segment( _numGeomDofs, numEdges );
    VectorType angles   = Argument.segment( _numGeomDofs + numEdges, numEdges );
    VectorType lambda1  = Argument.segment( _numGeomDofs + _dimLTheta, numEdges );
    VectorType lambda2  = Argument.segment( _numGeomDofs + _dimLTheta + numEdges, numEdges );
    
    
    // push diagonal blocks
    EdgeLengthSumHessian<ConfiguratorType>( _topology, lambda1, lengths).pushTriplets( geometry, tripletList );
    DihedralSumHessian<ConfiguratorType>( _topology, lambda2, angles).pushTriplets( geometry, tripletList );
     
    QuadraticMismatchHessian<ConfiguratorType>( _lengthsWeights, 1., _numGeomDofs, _numGeomDofs ).pushTriplets( tripletList );
    QuadraticMismatchHessian<ConfiguratorType>( _angleWeights, 1., _numGeomDofs + numEdges, _numGeomDofs + numEdges ).pushTriplets( tripletList );
    
    // push Lagrange multiplier blocks
    VectorType ones = VectorType::Ones( numEdges );
    EdgeLengthSumGradient<ConfiguratorType> DE( _topology, ones, lengths );
    if( _bdryMask.size() > 0 )
        DE.setBoundaryMask( _bdryMask );
    DE.pushTriplets( geometry, _numGeomDofs + _dimLTheta, 0, tripletList, true );    
    
    DihedralSumGradient<ConfiguratorType> DF( _topology, ones, angles );
    if( _bdryMask.size() > 0 )
        DF.setBoundaryMask( _bdryMask );
    DF.pushTriplets( geometry, _numGeomDofs + _dimLTheta + numEdges, 0, tripletList, true );  
    
    // push identity blocks
    for( int i = 0; i < _dimLTheta; i++ ){
        tripletList.push_back( TripletType( _numGeomDofs + i, _numGeomDofs + _dimLTheta + i, -1.) );
        tripletList.push_back( TripletType( _numGeomDofs + _dimLTheta + i, _numGeomDofs + i, -1.) );
    }

  }    
    
};

//==========================================================================================================
//==========================================================================================================
//==========================================================================================================

//==========================================================================================================
//! \brief L[X,Z] = Q[Z] + <Lambda, C[X,Z]> + \mu  * <C[X,Z], C[X,Z]>, for given scaler \mu and vector Lambda
//!  Here Q is (sum of different) quadratic edge based matching energies, and C is a vector-valued constraint energy 
//! \author Heeren
template<typename ConfiguratorType>
class ConstrainedLThetaAugmentedLagrangeFunctional : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::MatType    MatType;

  const MeshTopologySaver& _topology;
  int _numGeomDofs, _dimLTheta, _numDOFs;
  const VectorType& _prescribedLengths, _prescribedAngles;
  const VectorType& _lengthsWeights, _angleWeights; 
  const VectorType& _lengthsMultipliers, _angleMultipliers;
  RealType _augmWeight;
  
public:
  ConstrainedLThetaAugmentedLagrangeFunctional( const MeshTopologySaver& topology,
                                                const VectorType& prescribedLengths,
                                                const VectorType& prescribedAngles,
                                                const VectorType& lengthsWeights,
                                                const VectorType& angleWeights,
                                                const VectorType& lengthsMultipliers,
                                                const VectorType& angleMultipliers,
                                                RealType AugmentedWeight )  
  : _topology(topology), 
    _numGeomDofs( 3 * _topology.getNumVertices() ), 
    _dimLTheta( 2 * _topology.getNumEdges() ),
    _numDOFs( _numGeomDofs + _dimLTheta ),
    _prescribedLengths(prescribedLengths),
    _prescribedAngles(prescribedAngles),
    _lengthsWeights(lengthsWeights),
    _angleWeights(angleWeights),
    _lengthsMultipliers(lengthsMultipliers), 
    _angleMultipliers(angleMultipliers),
    _augmWeight(AugmentedWeight){  }
    
  void setAugmentedWeight( RealType AugmentedWeight ) {
    _augmWeight = AugmentedWeight;    
  }  

  // argument is (X,Z), where X is geometry and Z is LTheta variable
  void apply( const VectorType& Argument, RealType & Dest ) const {
    if( Argument.size() != _numDOFs )
      throw BasicException( "ConstrainedLThetaAugmentedLagrangeFunctional::apply(): argument has wrong size!");
    
    // get segments (TODO make efficient!)
    int numEdges = _topology.getNumEdges();
    VectorType geometry = Argument.segment(0, _numGeomDofs );
    VectorType lengths  = Argument.segment( _numGeomDofs, numEdges );
    VectorType angles   = Argument.segment( _numGeomDofs + numEdges, numEdges );
    
    Dest.setZero();    
    
    // add LTheta mismatch
    QuadraticMismatchFunctional<ConfiguratorType>( _lengthsWeights, _prescribedLengths ).applyAdd( lengths, Dest );
    QuadraticMismatchFunctional<ConfiguratorType>( _angleWeights, _prescribedAngles ).applyAdd( angles, Dest );
    
    // add constraints
    EdgeLengthSumEnergy<ConfiguratorType>( _topology, _lengthsMultipliers, lengths, _augmWeight ).applyAdd( geometry, Dest );
    DihedralSumEnergy<ConfiguratorType>( _topology, _angleMultipliers, angles, _augmWeight ).applyAdd( geometry, Dest );
  }
  
  //
  void getLengthsConstraintVector( const VectorType& Argument, VectorType& Dest ) const {
    if( Argument.size() != _numDOFs )
      throw BasicException( "ConstrainedLThetaAugmentedLagrangeFunctional::apply(): argument has wrong size!");
    
    // get segments (TODO make efficient!)
    int numEdges = _topology.getNumEdges();
    VectorType geometry = Argument.segment(0, _numGeomDofs );
    VectorType lengths  = Argument.segment( _numGeomDofs, numEdges );
    VectorType dummyWeights = VectorType::Zero( numEdges );

    EdgeLengthSumEnergy<ConfiguratorType>( _topology, dummyWeights, lengths, 0. ).evaluateIndivualTerms( geometry, Dest );
  }
  
  //
  void getAnglesConstraintVector( const VectorType& Argument, VectorType& Dest ) const {
    if( Argument.size() != _numDOFs )
      throw BasicException( "ConstrainedLThetaAugmentedLagrangeFunctional::apply(): argument has wrong size!");
    
    // get segments (TODO make efficient!)
    int numEdges = _topology.getNumEdges();
    VectorType geometry = Argument.segment(0, _numGeomDofs );
    VectorType angles   = Argument.segment( _numGeomDofs + numEdges, numEdges );
    VectorType dummyWeights = VectorType::Zero( numEdges );

    DihedralSumEnergy<ConfiguratorType>( _topology, dummyWeights, angles, 0. ).evaluateIndivualTerms( geometry, Dest );
  }

};


//==========================================================================================================
//! \brief Gradient of ConstrainedLThetaAugmentedLagrangeFunctional (see above)
//! \author Heeren
template<typename ConfiguratorType>
class ConstrainedLThetaAugmentedLagrangeGradient : public BaseOp<typename ConfiguratorType::VectorType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::MatType    MatType;

  const MeshTopologySaver& _topology;
  int _numGeomDofs, _dimLTheta, _numDOFs;
  const VectorType& _prescribedLengths, _prescribedAngles;
  const VectorType& _lengthsWeights, _angleWeights;
  const VectorType& _lengthsMultipliers, _angleMultipliers;
  RealType _augmWeight;
  
  
public:
  ConstrainedLThetaAugmentedLagrangeGradient( const MeshTopologySaver& topology,
                                     const VectorType& prescribedLengths,
                                     const VectorType& prescribedAngles,
                                     const VectorType& lengthsWeights,
                                     const VectorType& angleWeights,
                                     const VectorType& lengthsMultipliers,
                                     const VectorType& angleMultipliers,
                                     RealType AugmentedWeight )  
  : _topology(topology), 
    _numGeomDofs( 3 * _topology.getNumVertices() ), 
    _dimLTheta( 2 * _topology.getNumEdges() ),
    _numDOFs( _numGeomDofs + _dimLTheta ),
    _prescribedLengths(prescribedLengths),
    _prescribedAngles(prescribedAngles),
    _lengthsWeights(lengthsWeights),
    _angleWeights(angleWeights),
    _lengthsMultipliers(lengthsMultipliers), 
    _angleMultipliers(angleMultipliers),
    _augmWeight(AugmentedWeight){  }

  //  
  void setAugmentedWeight( RealType AugmentedWeight ) {
    _augmWeight = AugmentedWeight;    
  }    
    
  // argument is (X,Z), where X is geometry and Z is LTheta variable
  void apply( const VectorType& Argument, VectorType& Gradient ) const {
    if( Argument.size() != _numDOFs )
      throw BasicException( "ConstrainedLThetaAugmentedLagrangeGradient::apply(): argument has wrong size!");
    if( Gradient.size() != _numDOFs )
        Gradient.resize( _numDOFs );
    Gradient.setZero();
    
    // get segments (TODO make efficient!)
    int numEdges = _topology.getNumEdges();
    VectorType geometry = Argument.segment(0, _numGeomDofs );
    VectorType lengths  = Argument.segment( _numGeomDofs, numEdges );
    VectorType angles   = Argument.segment( _numGeomDofs + numEdges, numEdges );
    
    VectorType partialGrad;
    
    // gradient wrt. geometry
    EdgeLengthSumGradient<ConfiguratorType>( _topology, _lengthsMultipliers, lengths, _augmWeight ).apply( geometry, partialGrad );
    DihedralSumGradient<ConfiguratorType>( _topology, _angleMultipliers, angles, _augmWeight ).applyAdd( geometry, partialGrad );
    for( int i = 0; i < _numGeomDofs; i++ )
        Gradient[i] = partialGrad[i];
    
    // gradient wrt. ltheta
    QuadraticMismatchGradient<ConfiguratorType>( _lengthsWeights, _prescribedLengths ).apply( lengths, partialGrad);
    for( int i = 0; i < numEdges; i++ )
        Gradient[_numGeomDofs + i] = partialGrad[i] - _lengthsMultipliers[i];
    QuadraticMismatchGradient<ConfiguratorType>( _angleWeights, _prescribedAngles ).apply( angles, partialGrad );
    for( int i = 0; i < numEdges; i++ )
        Gradient[_numGeomDofs + numEdges + i] = partialGrad[i] - _angleMultipliers[i];
  }

};
 

//==========================================================================================================
//! \brief Gradient of ConstrainedLThetaAugmentedLagrangeFunctional (see above)
//! \author Heeren
template<typename ConfiguratorType>
class ConstrainedLThetaAugmentedLagrangeHessian : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::MatType    MatType;

  const MeshTopologySaver& _topology;
  int _numGeomDofs, _dimLTheta, _numDOFs;
  const VectorType& _prescribedLengths, _prescribedAngles;
  const VectorType& _lengthsWeights, _angleWeights;
  const VectorType& _lengthsMultipliers, _angleMultipliers;
  RealType _augmWeight;
  std::vector<int> _bdryMask;
  
public:
  ConstrainedLThetaAugmentedLagrangeHessian( const MeshTopologySaver& topology,
                                     const VectorType& prescribedLengths,
                                     const VectorType& prescribedAngles,
                                     const VectorType& lengthsWeights,
                                     const VectorType& angleWeights,
                                     const VectorType& lengthsMultipliers,
                                     const VectorType& angleMultipliers,
                                     RealType AugmentedWeight )  
  : _topology(topology), 
    _numGeomDofs( 3 * _topology.getNumVertices() ), 
    _dimLTheta( 2 * _topology.getNumEdges() ),
    _numDOFs( _numGeomDofs + _dimLTheta ),
    _prescribedLengths(prescribedLengths),
    _prescribedAngles(prescribedAngles),
    _lengthsWeights(lengthsWeights),
    _angleWeights(angleWeights),
    _lengthsMultipliers(lengthsMultipliers), 
    _angleMultipliers(angleMultipliers),
    _augmWeight(AugmentedWeight){  }
    
  //  
  void setAugmentedWeight( RealType AugmentedWeight ) {
    _augmWeight = AugmentedWeight;    
  }    
  
  // argument is (X,Z), where X is geometry and Z is LTheta variable
  void apply( const VectorType& Argument, MatrixType& Hessian ) const {
    if( Argument.size() != _numDOFs )
      throw BasicException( "ConstrainedLThetaAugmentedLagrangeHessian::apply(): argument has wrong size!");
    Hessian.resize( _numDOFs, _numDOFs );
    Hessian.setZero();
    
    // set up triplet list
    TripletListType tripletList;
    // fill matrix from triplets
    pushTriplets( Argument, tripletList );
    Hessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
  
  //
  void setBoundaryMask( const std::vector<int>& mask ) {
    _bdryMask.resize( mask.size() );
    _bdryMask = mask;    
  }
  
  // argument is (X,Z), where X is geometry and Z is LTheta variable
  void pushTriplets( const VectorType& Argument, TripletListType& tripletList ) const {
    if( Argument.size() != _numDOFs )
      throw BasicException( "ConstrainedLThetaAugmentedLagrangeHessian::pushTriplets(): argument has wrong size!");
    
    // get segments (TODO make efficient!)
    int numEdges = _topology.getNumEdges();
    VectorType geometry = Argument.segment(0, _numGeomDofs );
    VectorType lengths  = Argument.segment( _numGeomDofs, numEdges );
    VectorType angles   = Argument.segment( _numGeomDofs + numEdges, numEdges );    
    
    // push diagonal blocks
    EdgeLengthSumHessian<ConfiguratorType>( _topology, _lengthsMultipliers, lengths, _augmWeight).pushTriplets( geometry, tripletList );
    DihedralSumHessian<ConfiguratorType>( _topology, _angleMultipliers, angles, _augmWeight).pushTriplets( geometry, tripletList );
     
    QuadraticMismatchHessian<ConfiguratorType>( _lengthsWeights, 1., _numGeomDofs, _numGeomDofs ).pushTriplets( tripletList );
    QuadraticMismatchHessian<ConfiguratorType>( _angleWeights, 1., _numGeomDofs + numEdges, _numGeomDofs + numEdges ).pushTriplets( tripletList );
  }    
    
};

#endif //

