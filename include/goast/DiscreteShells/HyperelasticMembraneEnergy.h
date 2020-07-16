// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Hyperelastic membrane energy and derivatives.
 * \author Heeren
 */

#ifndef HYPERELASTICENERGY_HH
#define HYPERELASTICENERGY_HH

#include <goast/Core/Auxiliary.h>
#include <goast/Core/LocalMeshGeometry.h>
#include <goast/Core/Topology.h>
#include <goast/Core/DeformationInterface.h>

/**
 * \brief Hyperelastic membrane energy (see eq. (8) in \cite HeRuWa12 )
 * \author Heeren
 *
 * Let \f$ x, \tilde x \f$ be the geometries of two meshes that are in dense correspondence, let \f $ F \f$ be the set of faces.
 * For a face \f$ f \in F \f$ we consider discrete first fundamental forms \f$ g_f \f$ and \f$ \tilde g_f \f$, respectively.
 *
 * In detail, if \f$ x_i, x_j x_k \f$ are vertex positions of three nodes of \f$ f \f$ in \f$ x \f$,
 * we define the \f$ 2 \times 2 \f$ - matrix \f$ g_f = [x_j - x_i | x_k - x_i]^T [x_j - x_i | x_k - x_i] \f$.
 *
 * Furthermore, for \f$ f \in F \f$ we define the geometric distortion tensor as \f$ G_f = g_f^{-1}\tilde g_f \f$.
 * Then the hyperelastic energy is given as \f[ E[x, \tilde x] = \sum_{f \in F} a_f W( \tr G_f, \det G_f) \f]
 * where \f$ a_f \f$ is the face area and the hyperelastic energy density
 * \f[ W(a,d) = \frac\mu2 a + \frac\lambda4 d - (\frac\mu2 + \frac\lambda4) \log d - \mu -\frac\lambda4 \f]
 * for physical parameters \f$ \mu, \lambda > 0 \f$.
 *
 * Note that the energy might either be thought of as \f$ x \mapsto E[x, \tilde x] \f$ (active shell is undeformed shell)
 * or \f$ \tilde x \mapsto E[x, \tilde x] \f$ (active shell is deformed shell).
 * The active shell is considered the argument whereas the inactive shell is given in the constructor.
 */
template<typename ConfiguratorType>
class NonlinearMembraneEnergy : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::MatType    MatType;

  const MeshTopologySaver& _topology;
  const VectorType&  _inactiveGeometry;
  const bool _activeShellIsDeformed;
  RealType _mu, _lambdaQuarter, _muHalfPlusLambdaQuarter;

public:

  NonlinearMembraneEnergy( const MeshTopologySaver& topology,
                           const VectorType& InactiveGeometry,
		           const bool ActiveShellIsDeformed,
                           RealType Mu = 1., 
                           RealType Lambda = 1. ) 
  : _topology( topology), 
    _inactiveGeometry(InactiveGeometry), 
    _activeShellIsDeformed( ActiveShellIsDeformed),
    _mu(Mu),
    _lambdaQuarter(Lambda/4.),
    _muHalfPlusLambdaQuarter(_mu/2. + _lambdaQuarter){}
    

  // energy evaluation
  void apply( const VectorType& ActiveGeometry, RealType & Dest ) const {
      
    if( ActiveGeometry.size() != _inactiveGeometry.size() ){
      std::cerr << "size of active = " << ActiveGeometry.size() << " vs. size of inactive = " << _inactiveGeometry.size() << std::endl;
      throw BasicException( "NonlinearMembraneEnergy::apply(): sizes dont match!");
    }

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
      // compute volume
      temp.makeCrossProduct( Ei, Ej );
      RealType volDefSqr = temp.normSqr() / 4.;
      
      // check whether area is finite
      if( std::sqrt( volDefSqr ) < 1e-15 ){
          std::cerr << "WARNING thrown in NonlinearMembraneEnergy! Probably face " << faceIdx << " in deformed mesh is degenerated! " << std::endl;
          Dest = std::numeric_limits<RealType>::infinity();
          return;
      }
      
      // set up undeformed vertices and edges
      getXYZCoord<VectorType, VecType>( *undefShellP, temp, pi);
      getXYZCoord<VectorType, VecType>( *undefShellP, Ej, pj);
      getXYZCoord<VectorType, VecType>( *undefShellP, Ek, pk);
      Ei = Ek - Ej;
      Ej = temp - Ek;
      Ek = Ei+Ej;

      // compute volume
      temp.makeCrossProduct( Ei, Ej );
      RealType volUndefSqr = temp.normSqr() / 4.;
      RealType volUndef = std::sqrt( volUndefSqr );
      
      // check whether area is finite
      if( volUndef < 1e-15 ){
          std::cerr << "WARNING thrown in NonlinearMembraneEnergy! Probably face " << faceIdx << " in undeformed mesh is degenerated! " << std::endl;
          Dest = std::numeric_limits<RealType>::infinity();
          return;
      }      
      //CAUTION mind the signs! (Ek is actually -Ek here!)
      RealType traceTerm = ( dotProduct(Ej,Ek) * liSqr + dotProduct(Ek,Ei) * ljSqr - dotProduct(Ei,Ej) * lkSqr );
      
      // volume of triangle * evaluation of energy density  
      Dest +=  (_mu/8. *  traceTerm + _lambdaQuarter * volDefSqr) / volUndef -  ( _muHalfPlusLambdaQuarter * std::log( volDefSqr / volUndefSqr ) + _mu + _lambdaQuarter) * volUndef;
      
#ifdef DEBUGMODE
      if( std::isnan( Dest ) ){
          std::cerr << "NaN in membrane energy in face " << faceIdx << "! " << std::endl;
          if( hasNanEntries(ActiveGeometry) )
            std::cerr << "Argument has NaN entries! " << std::endl;
          else{
            std::cerr << "traceTerm = " << traceTerm << std::endl;
            std::cerr << "volUndefSqr = " << volUndefSqr << std::endl;
            std::cerr << "volDefSqr = " << volDefSqr << std::endl;
          }
          throw BasicException("NonlinearMembraneEnergy::apply(): NaN Error!");
      }
#endif
    }
  }

 };
   
 
//! \brief First derivative of NonlinearMembraneEnergy w.r.t. the deformed configuration.
//! \author Heeren
template<typename ConfiguratorType>
class NonlinearMembraneGradientDef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {

protected:
	typedef typename ConfiguratorType::RealType   RealType;

	typedef typename ConfiguratorType::VectorType VectorType;
	typedef typename ConfiguratorType::VecType    VecType;
	typedef typename ConfiguratorType::MatType    MatType;

	const MeshTopologySaver& _topology;
	const VectorType&  _undefShell;
	RealType _mu, _lambdaQuarter, _muHalfPlusLambdaQuarter;

public:
	NonlinearMembraneGradientDef(const MeshTopologySaver& topology,
		const VectorType& undefShell,
		RealType Mu = 1.,
		RealType Lambda = 1.)
		: _topology(topology),
		_undefShell(undefShell),
		_mu(Mu),
		_lambdaQuarter(Lambda / 4.),
		_muHalfPlusLambdaQuarter(_mu / 2. + _lambdaQuarter){}

	//
	void apply(const VectorType& defShell, VectorType& Dest) const {

		if (_undefShell.size() != defShell.size()){
			std::cerr << "size of undef = " << _undefShell.size() << " vs. size of def = " << defShell.size() << std::endl;
			throw BasicException("NonlinearMembraneGradientDef::apply(): sizes dont match!");
		}

		if (Dest.size() != defShell.size())
			Dest.resize(defShell.size());

		Dest.setZero();

		for (int faceIdx = 0; faceIdx < _topology.getNumFaces(); ++faceIdx){

			std::vector<int> nodesIdx(3);
			std::vector<VecType> nodes(3), undefEdges(3), defEdges(3);
			VecType temp;
			for (int j = 0; j < 3; j++)
				nodesIdx[j] = _topology.getNodeOfTriangle(faceIdx, j);

			//! get undeformed quantities
			for (int j = 0; j < 3; j++)
				getXYZCoord<VectorType, VecType>(_undefShell, nodes[j], nodesIdx[j]);
			for (int j = 0; j < 3; j++)
				undefEdges[j] = nodes[(j + 2) % 3] - nodes[(j + 1) % 3];
			// compute volume
			temp.makeCrossProduct(nodes[2] - nodes[1], nodes[0] - nodes[2]);
			RealType volUndefSqr = temp.normSqr() / 4.;
			RealType volUndef = std::sqrt(volUndefSqr);

			//! get deformed quantities
			for (int j = 0; j < 3; j++)
				getXYZCoord<VectorType, VecType>(defShell, nodes[j], nodesIdx[j]);
			for (int j = 0; j < 3; j++)
				defEdges[j] = nodes[(j + 2) % 3] - nodes[(j + 1) % 3];
			// compute volume
			temp.makeCrossProduct(nodes[2] - nodes[1], nodes[0] - nodes[2]);
			RealType volDefSqr = temp.normSqr() / 4.;
			RealType volDef = std::sqrt(volDefSqr);

			//! trace part of gradient, e_tr = volUndef * _mu/2. *  traceDistTensor   
			VecType factors;
			for (int i = 0; i < 3; i++)
				factors[i] = -0.25 * _mu * dotProduct(undefEdges[(i + 2) % 3], undefEdges[(i + 1) % 3]) / volUndef;
			RealType factor = 2. * (_lambdaQuarter * volDef / volUndef - _muHalfPlusLambdaQuarter * volUndef / volDef);

			for (int i = 0; i < 3; i++){
				getAreaGradient(nodes[(i + 1) % 3], nodes[(i + 2) % 3], nodes[i], temp);
				for (int j = 0; j < 3; j++)
					Dest[j*_topology.getNumVertices() + nodesIdx[i]] += factor * temp[j] + factors[(i + 1) % 3] * defEdges[(i + 1) % 3][j] - factors[(i + 2) % 3] * defEdges[(i + 2) % 3][j];
			}

		}
	}

};


//! \brief Second derivative of NonlinearMembraneEnergy w.r.t. the deformed configuration.
//! \author Heeren
template<typename ConfiguratorType>
class NonlinearMembraneHessianDef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

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
  const RealType _mu, _lambda, _muHalfPlusLambdaQuarter;
  mutable int _rowOffset, _colOffset;

public:
NonlinearMembraneHessianDef( const MeshTopologySaver& topology,
                             const VectorType& undefShell,
			     const RealType Factor = 1.,
                             int rowOffset = 0,
                             int colOffset = 0,
                             RealType Mu = 1., 
                             RealType Lambda = 1. ) 
  : _topology( topology), 
    _undefShell(undefShell), 
    _factor( Factor ),
    _mu(Mu),
    _lambda(Lambda),
    _muHalfPlusLambdaQuarter(_mu/2. + _lambda/4.),
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
    // per face we have 3 active vertices, i.e. 9 combinations each producing a 3x3-matrix
    tripletList.reserve( 9 * 9 * _topology.getNumFaces() );
    
    pushTriplets( defShell, tripletList );
    
    // fill matrix from triplets
    Hessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }

  //
  void pushTriplets( const VectorType& defShell, TripletListType& tripletList ) const  {
      
    if( _undefShell.size() != defShell.size() ){
      std::cerr << "size of undef = " << _undefShell.size() << " vs. size of def = " << defShell.size() << std::endl;
      throw BasicException( "NonlinearMembraneHessianDef::pushTriplets(): sizes dont match!");
    }
      
    // run over all faces
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

      //! get deformed quantities
      for( int j = 0; j < 3; j++ )
        getXYZCoord<VectorType, VecType>( defShell, nodes[j], nodesIdx[j]); 
      for( int j = 0; j < 3; j++ )
          defEdges[j] = nodes[(j+2)%3] - nodes[(j+1)%3];
      // compute volume
      temp.makeCrossProduct( defEdges[0], defEdges[1] );
      RealType volDefSqr = temp.normSqr() / 4.;
      RealType volDef = std::sqrt( volDefSqr );
      
      VecType traceFactors;
      for( int i = 0; i < 3; i++ )
        traceFactors[i] = -0.25 * _mu * dotProduct( undefEdges[(i+2)%3], undefEdges[(i+1)%3] ) / volUndef;
      
      RealType mixedFactor = 0.5 * _lambda / volUndef + 2. * _muHalfPlusLambdaQuarter * volUndef / volDefSqr;
      RealType areaFactor  = 0.5 * _lambda * volDef / volUndef - 2. * _muHalfPlusLambdaQuarter * volUndef / volDef;
      
      // precompute area gradients
      std::vector<VecType> gradArea(3);
      for( int i = 0; i < 3; i++ )
        getAreaGradient(  nodes[(i+1)%3], nodes[(i+2)%3], nodes[i], gradArea[i] );
        
      // compute local matrices
      MatType tensorProduct, H, auxMat;
      
      // i==j
      for( int i = 0; i < 3; i++ ){
        getHessAreaKK( nodes[(i+1)%3], nodes[(i+2)%3], nodes[i], auxMat );
        tensorProduct.makeTensorProduct( gradArea[i], gradArea[i] );
        getWeightedMatrixSum( areaFactor, auxMat, mixedFactor, tensorProduct, H );
        H.addToDiagonal( traceFactors[(i+1)%3] + traceFactors[(i+2)%3] );
        localToGlobal( tripletList, nodesIdx[i], nodesIdx[i], H );
      }      
      
      // i!=j
      for( int i = 0; i < 3; i++ ){
        getHessAreaIK( nodes[i], nodes[(i+1)%3], nodes[(i+2)%3], auxMat );
        tensorProduct.makeTensorProduct( gradArea[i], gradArea[(i+2)%3] );
        getWeightedMatrixSum( areaFactor, auxMat, mixedFactor, tensorProduct, H );
        H.addToDiagonal( -traceFactors[(i+1)%3] );
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


//! \brief First derivative of NonlinearMembraneEnergy w.r.t. the undeformed configuration.
//! \author Heeren
template< typename ConfiguratorType>
class NonlinearMembraneGradientUndef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {

protected:    
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::VecType     VecType;

  const MeshTopologySaver& _topology;
  const VectorType& _defShell;
  const RealType _mu, _lambdaQuarter, _muHalfPlusLambdaQuarter, _const;
  
public:
NonlinearMembraneGradientUndef( const MeshTopologySaver& topology, 
                                const VectorType& defShell, 
                                RealType Mu = 1., 
                                RealType Lambda = 1. )
  : _topology( topology),
    _defShell(defShell),
    _mu(Mu),
    _lambdaQuarter(Lambda/4.),
    _muHalfPlusLambdaQuarter(_mu/2. + _lambdaQuarter),
    _const(_mu + _lambdaQuarter){}

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
      
      // compute volume
      temp.makeCrossProduct(  nodes[1] - nodes[0], nodes[2] - nodes[1] );
      RealType volDefSqr = temp.normSqr() / 4.;
      RealType volDef = std::sqrt( volDefSqr ); 
      
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
      
      RealType factor1 = (0.125 * _mu *  traceTerm + _lambdaQuarter * volDefSqr) / volUndefSqr;
      RealType factor2 = _muHalfPlusLambdaQuarter * std::log( volDefSqr / volUndefSqr ) + _const - 2 * _muHalfPlusLambdaQuarter;
      RealType factorAreaGrad  = factor1  + factor2;
      RealType factorTraceGrad = 0.125 * _mu / volUndef;
      
      std::vector<VecType> gradTrace(3);
      for( int i = 0; i < 3; i++ )
          for( int j = 0; j < 3; j++ )   
            gradTrace[i][j] = defLengthSqr[i] * (undefEdges[(i+1)%3][j] - undefEdges[(i+2)%3][j]) +  undefEdges[i][j] * (defLengthSqr[(i+1)%3] - defLengthSqr[(i+2)%3] );      
      
      // E = (0.125 * _mu * traceTerm + _lambdaQuarter * volDefSqr ) / volUndef;
      // grad E[j] = 0.125 * _mu * grad traceTerm[j] / volUndef  - (0.125 * _mu * traceTerm + _lambdaQuarter * volDefSqr )  grad volUndef[j] / volUndef^2
      for( int i = 0; i < 3; i++ ){
        getAreaGradient(  nodes[(i+1)%3], nodes[(i+2)%3], nodes[i], temp );
        for( int j = 0; j < 3; j++ )
          Dest[j*_topology.getNumVertices() + nodesIdx[i]] += factorTraceGrad * gradTrace[i][j] - factorAreaGrad * temp[j];
      }

  }
}

}; 


//! \brief Second derivative of NonlinearMembraneEnergy w.r.t. the undeformed configuration.
//! \author Heeren
template <typename ConfiguratorType>
class NonlinearMembraneHessianUndef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

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
  const RealType _mu, _lambdaQuarter, _muHalfPlusLambdaQuarter, _const;
  const RealType _factor;
  mutable int _rowOffset, _colOffset;
  
public:
NonlinearMembraneHessianUndef( const MeshTopologySaver& topology, const VectorType& defShell, RealType Factor = 1., int rowOffset = 0, int colOffset = 0, RealType Mu = 1., RealType Lambda = 1. ) 
  : _topology( topology),
    _defShell(defShell),
    _mu(Mu),
    _lambdaQuarter(Lambda/4.),
    _muHalfPlusLambdaQuarter(_mu/2. + _lambdaQuarter),
    _const(_mu + _lambdaQuarter),
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

      RealType areaFactor        = 0.125 * _mu * traceTerm + _lambdaQuarter * volDefSqr;
      RealType negHessAreaFactor =  _muHalfPlusLambdaQuarter * (std::log(volDefSqr / volUndefSqr) - 2.) + _const + areaFactor / volUndefSqr;
      RealType mixedAreaFactor   = 2 * (areaFactor / volUndefSqr + _muHalfPlusLambdaQuarter) / volUndef;
      RealType mixedFactor       = -0.125 * _mu / volUndefSqr;
      RealType hessTraceFactor   =  0.125 * _mu / volUndef;
      
      
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


//! \brief Second (mixed) derivative of NonlinearMembraneEnergy.
//! \author Heeren
template <typename ConfiguratorType >
class NonlinearMembraneHessianMixed : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

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
  const RealType _mu, _lambdaQuarter, _muHalfPlusLambdaQuarter, _const;
  const RealType _factor;
  mutable int _rowOffset, _colOffset;

public:
NonlinearMembraneHessianMixed( const MeshTopologySaver& topology,
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
  _mu(Mu),
  _lambdaQuarter(Lambda/4.),
  _muHalfPlusLambdaQuarter(_mu/2. + _lambdaQuarter),
  _const(_mu + _lambdaQuarter),
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
        throw BasicException("NonlinearMembraneHessianMixed::apply: sizes dont match!");        
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
      
      // precompute deformed area gradients
      std::vector<VecType> gradDefArea(3);
      for( int i = 0; i < 3; i++ )
        getAreaGradient(  nodes[(i+1)%3], nodes[(i+2)%3], nodes[i], gradDefArea[i] );      
      
      // compute local matrices
      MatType tensorProduct, H, auxMat;           
      RealType mixedTraceHessFactor = 0.125 * _mu / volUndef;
      RealType MixedAreaFactor      = -2. * (_lambdaQuarter * volDef / volUndef + _muHalfPlusLambdaQuarter * volUndef/ volDef ) / volUndef;
      RealType MixedFactor          = -0.125 * _mu / volUndefSqr;       

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
          
          // mixed area term
          tensorProduct.makeTensorProduct( gradUndefArea[i], gradDefArea[j] );
          H.addMultiple( tensorProduct, MixedAreaFactor );
        
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



//! \brief Deformation class for NonlinearMembraneEnergy
//! \author Heeren
template <typename ConfiguratorType>
class NonlinearMembraneDeformation : public DeformationBase<ConfiguratorType>{
    
protected:    
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;

  typedef typename ConfiguratorType::TripletType TripletType;  
  typedef std::vector<TripletType> TripletListType;
  
  const MeshTopologySaver& _topology;
  RealType _memWeight;

public:
  NonlinearMembraneDeformation( const MeshTopologySaver& Topology, RealType memWeight ) : _topology( Topology ), _memWeight( memWeight ) {}
  
  void applyEnergy ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, RealType & Dest ) const {
    NonlinearMembraneEnergy<ConfiguratorType>( _topology, UndeformedGeom, true ).apply( DeformedGeom, Dest );
    Dest *= _memWeight;
  }
  
  void applyUndefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest ) const {
    NonlinearMembraneGradientUndef<ConfiguratorType>( _topology, DeformedGeom ).apply( UndeformedGeom, Dest );
    Dest *= _memWeight;
  }
  
  void applyDefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest ) const {
    NonlinearMembraneGradientDef<ConfiguratorType>( _topology, UndeformedGeom ).apply( DeformedGeom, Dest );
    Dest *= _memWeight;
  }
  
  void pushTripletsDefHessian ( const VectorType& UndefGeom, const VectorType& DefGeom, TripletListType& triplets, int rowOffset, int colOffset, RealType factor = 1.0 ) const {
    NonlinearMembraneHessianDef<ConfiguratorType>( _topology, UndefGeom, factor * _memWeight, rowOffset, colOffset ).pushTriplets( DefGeom, triplets );
  }
   
  void pushTripletsUndefHessian ( const VectorType& UndefGeom, const VectorType& DefGeom, TripletListType& triplets, int rowOffset, int colOffset, RealType factor = 1.0 ) const {
      NonlinearMembraneHessianUndef<ConfiguratorType>( _topology, DefGeom, factor * _memWeight, rowOffset, colOffset ).pushTriplets( UndefGeom, triplets );
  }
  
  // mixed second derivative of deformation energy E[S_1, S_2], i.e. if "FirstDerivWRTDef" we have D_1 D_2 E[.,.], otherwise D_2 D_1 E[.,.]
  void pushTripletsMixedHessian ( const VectorType& UndefGeom, const VectorType& DefGeom, TripletListType& triplets, int rowOffset, int colOffset, const bool FirstDerivWRTDef, RealType factor = 1.0 ) const {
      NonlinearMembraneHessianMixed<ConfiguratorType>( _topology, UndefGeom, true, FirstDerivWRTDef, factor * _memWeight, rowOffset, colOffset ).pushTriplets( DefGeom, triplets );
  }
  
  int numOfNonZeroHessianEntries () const {
    // per face we have 3 active vertices, i.e. 9 combinations each producing a 3x3-matrix
    return 9 * 9 * _topology.getNumFaces();    
  }
    
};


#endif
