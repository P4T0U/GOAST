// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Original Discrete Shell's bending energy and derivatives.
 * \author Heeren
 * \cite GrHiDeSc03
 */
 #ifndef SIMPLEBENDINGENERGY_HH
#define SIMPLEBENDINGENERGY_HH


//== INCLUDES =================================================================
#include <goast/Core/Auxiliary.h>
#include <goast/Core/LocalMeshGeometry.h>
#include <goast/Core/Topology.h>
#include <goast/Core/DeformationInterface.h>

//#define DEBUGMODE

//==========================================================================================================
// DISCRETE SHELLS BENDING ENERGY (BASED ON GRINSPUN ET AL. 2003)
//==========================================================================================================

/**
 * \brief Discrete Shell's bending energy taken from  Grinspun et al. \cite GrHiDeSc03
 * \author Heeren
 *
 * Let \f$ x, \tilde x \f$ be the geometries of two meshes that are in dense correspondence, let \f $ E \f$ be the set of edges.
 * Then this class realizes the energy \f[ E[x, \tilde x] = |sum_{e \in E}  \frac{ (\theta_e[x] - \theta_e[\tilde x])^2 }{d_e[x]} l_e^2[x]\, , \f]
 * where \f$ \theta_e \f$ is the dihedral angle at the edge, \f$ l_e \f$ the length of the edge
 * and \f$ d_e = \frac13 (a_1 + a_2) \f$, if \f$ a_1 \f$ and \f$ a_2 \f$ denote the face area of the two adjacent triangles, respectively.
 *
 * Note that the energy might either be thought of as \f$ x \mapsto E[x, \tilde x] \f$ (active shell is undeformed shell)
 * or \f$ \tilde x \mapsto E[x, \tilde x] \f$ (active shell is deformed shell).
 * The active shell is considered the argument whereas the inactive shell is given in the constructor.
 */
template<typename ConfiguratorType>
class SimpleBendingEnergy : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::MatType    MatType;

  const MeshTopologySaver& _topology;
  const VectorType&  _inactiveGeometry;
  const bool _activeShellIsDeformed;
  RealType _weight;

public:

  SimpleBendingEnergy( const MeshTopologySaver& topology,
                       const VectorType& InactiveGeometry,
		       const bool ActiveShellIsDeformed,
                       RealType Weight = 1. ) 
  : _topology( topology), 
    _inactiveGeometry(InactiveGeometry), 
    _activeShellIsDeformed( ActiveShellIsDeformed),
    _weight(Weight){}

  // energy evaluation
  void apply( const VectorType& ActiveGeometry, RealType & Dest ) const {

    if( ActiveGeometry.size() != _inactiveGeometry.size() ){
      std::cerr << "size of active = " << ActiveGeometry.size() << " vs. size of inactive = " << _inactiveGeometry.size() << std::endl;
      throw BasicException( "SimpleBendingEnergy::apply(): sizes dont match!");
    }
    
    const VectorType* defShellP   = _activeShellIsDeformed ? &ActiveGeometry : &_inactiveGeometry;
    const VectorType* undefShellP = _activeShellIsDeformed ? &_inactiveGeometry : &ActiveGeometry;
    
    Dest = 0.;
    
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

      // set up vertices and edges
      VecType Pi, Pj, Pk, Pl, temp;

      // get deformed geometry
      getXYZCoord<VectorType, VecType>( *defShellP, Pi, pi);
      getXYZCoord<VectorType, VecType>( *defShellP, Pj, pj);
      getXYZCoord<VectorType, VecType>( *defShellP, Pk, pk);
      getXYZCoord<VectorType, VecType>( *defShellP, Pl, pl);

      // compute deformed dihedral angle
      RealType delTheta = getDihedralAngle( Pi, Pj, Pk, Pl );
      
      // get undeformed geometry
      getXYZCoord<VectorType, VecType>( *undefShellP, Pi, pi);
      getXYZCoord<VectorType, VecType>( *undefShellP, Pj, pj);
      getXYZCoord<VectorType, VecType>( *undefShellP, Pk, pk);
      getXYZCoord<VectorType, VecType>( *undefShellP, Pl, pl);

      // compute volume, length of edge and theta difference
      // A_e = h_e * l_e = (|T_1| + |T_2|)/3 if T_1 and T_2 share edge e
      // Furthermore, |T| = 0.5 * |e_1 x e_2|, if e_1 and e_2 are edges of T
      temp.makeCrossProduct( Pk-Pj, Pi-Pk );
      RealType vol = temp.norm() / 2.;
      temp.makeCrossProduct( Pl-Pi, Pj-Pl );
      vol += temp.norm() / 2.;
      RealType elengthSqr( dotProduct(Pj-Pi,Pj-Pi) );       
      delTheta -= getDihedralAngle( Pi, Pj, Pk, Pl );      
      // CAUTION We omitted a factor 3 here!
      Dest += _weight * delTheta * delTheta * elengthSqr / vol;
      
#ifdef DEBUGMODE
      if( std::isnan( Dest ) ){
          std::cerr << "NaN in simple bending energy in edge " << edgeIdx << "! " << std::endl;
          if( hasNanEntries(ActiveGeometry) ){
            std::cerr << "Argument has NaN entries! " << std::endl;
          }
          else{
            std::cerr << "delTheta = " << delTheta << std::endl;            
            std::cerr << "elengthSqr = " << elengthSqr << std::endl;
            std::cerr << "vol = " << vol << std::endl;  
          }
          throw BasicException("SimpleBendingEnergy::apply(): NaN Error!");
      }
#endif
    }
  }

  void setWeight(double Eta)
  {
	  _weight = Eta;
  }
 };

//==========================================================================================================
//! \brief First derivative of SimpleBendingEnergy w.r.t. the deformed configuration
//! \author Heeren
template<typename ConfiguratorType>
class SimpleBendingGradientDef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {

	typedef typename ConfiguratorType::RealType   RealType;
	typedef typename ConfiguratorType::VectorType VectorType;
	typedef typename ConfiguratorType::VecType    VecType;
	typedef typename ConfiguratorType::MatType    MatType;

	const MeshTopologySaver& _topology;
	const VectorType&  _undefShell;
	RealType _weight;

public:
	SimpleBendingGradientDef(const MeshTopologySaver& topology,
		const VectorType& undefShell,
		RealType Weight = 1.) : _topology(topology), _undefShell(undefShell), _weight(Weight) {}

	void apply(const VectorType& defShell, VectorType& Dest) const {

		if (_undefShell.size() != defShell.size()){
			std::cerr << "size of undef = " << _undefShell.size() << " vs. size of def = " << defShell.size() << std::endl;
			throw BasicException("SimpleBendingGradientDef::apply(): sizes dont match!");
		}

		if (Dest.size() != defShell.size())
			Dest.resize(defShell.size());

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

			//! first get undefomed quantities
			VecType Pi, Pj, Pk, Pl, temp;
			getXYZCoord<VectorType, VecType>(_undefShell, Pi, pi);
			getXYZCoord<VectorType, VecType>(_undefShell, Pj, pj);
			getXYZCoord<VectorType, VecType>(_undefShell, Pk, pk);
			getXYZCoord<VectorType, VecType>(_undefShell, Pl, pl);
			RealType delTheta = getDihedralAngle(Pi, Pj, Pk, Pl);


			// compute Ak + Al
			temp.makeCrossProduct(Pk - Pj, Pi - Pk);
			RealType vol = temp.norm() / 2.;
			temp.makeCrossProduct(Pl - Pi, Pj - Pl);
			vol += temp.norm() / 2.;
			RealType elengthSqr(dotProduct(Pj - Pi, Pj - Pi));

			//! now get the deformed values
			getXYZCoord<VectorType, VecType>(defShell, Pi, pi);
			getXYZCoord<VectorType, VecType>(defShell, Pj, pj);
			getXYZCoord<VectorType, VecType>(defShell, Pk, pk);
			getXYZCoord<VectorType, VecType>(defShell, Pl, pl);

			// compute weighted differnce of dihedral angles
			delTheta -= getDihedralAngle(Pi, Pj, Pk, Pl);
			delTheta *= -2. * _weight * elengthSqr / vol;

			// compute first derivatives of dihedral angle
			VecType thetak, thetal, thetai, thetaj;
			getThetaGradK(Pi, Pj, Pk, thetak);
			getThetaGradK(Pj, Pi, Pl, thetal);
			getThetaGradI(Pi, Pj, Pk, Pl, thetai);
			getThetaGradJ(Pi, Pj, Pk, Pl, thetaj);

			// assemble in global vector
			for (int i = 0; i < 3; i++){
				Dest[i*_topology.getNumVertices() + pi] += delTheta * thetai[i];
				Dest[i*_topology.getNumVertices() + pj] += delTheta * thetaj[i];
				Dest[i*_topology.getNumVertices() + pk] += delTheta * thetak[i];
				Dest[i*_topology.getNumVertices() + pl] += delTheta * thetal[i];
			}


#ifdef DEBUGMODE
			if( hasNanEntries( Dest ) ){
				std::cerr << "NaN in simple bending gradient deformed in edge " << edgeIdx << "! " << std::endl;
				if( hasNanEntries(defShell) ){
					std::cerr << "Argument has NaN entries! " << std::endl;
				}
				else{
					std::cerr << "delTheta = " << delTheta << std::endl;            
					std::cerr << "elengthSqr = " << elengthSqr << std::endl;
					std::cerr << "vol = " << vol << std::endl;  
				}
				throw BasicException("SimpleBendingGradientDef::apply(): NaN Error!");
			}
#endif

		}
	}

	void setWeight(double Eta)
	{
		_weight = Eta;
	}
};


//! \brief Second derivative of SimpleBendingEnergy w.r.t. the deformed configuration
//! \author Heeren
template<typename ConfiguratorType>
class SimpleBendingHessianDef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

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
  mutable int _rowOffset, _colOffset;
  
public:
  SimpleBendingHessianDef( const MeshTopologySaver& topology,
                           const VectorType& undefShell,
			   const RealType Weight = 1.,
                           int rowOffset = 0, 
                           int colOffset = 0 ) : 
                    _topology( topology),
                    _undefShell(undefShell), 
                    _factor( Weight ), 
                    _rowOffset(rowOffset), 
                    _colOffset(colOffset)  {}
    
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
    // per edge we have 4 active vertices, i.e. 16 combinations each producing a 3x3-matrix
    tripletList.reserve( 16 * 9 * _topology.getNumEdges() );   
    // fill matrix from triplets
    pushTriplets( defShell, tripletList );
    Hessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
  
  // fill triplets
  void pushTriplets( const VectorType& defShell, TripletListType& tripletList ) const {
      
    if( _undefShell.size() != defShell.size() ){
      std::cerr << "size of undef = " << _undefShell.size() << " vs. size of def = " << defShell.size() << std::endl;
      throw BasicException( "SimpleBendingHessianDef::pushTriplets(): sizes dont match!");
    }
  
  
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

      //! first get undefomed quantities
      VecType Pi, Pj, Pk, Pl, temp;
      getXYZCoord<VectorType, VecType>( _undefShell, Pi, pi);
      getXYZCoord<VectorType, VecType>( _undefShell, Pj, pj);
      getXYZCoord<VectorType, VecType>( _undefShell, Pk, pk);
      getXYZCoord<VectorType, VecType>( _undefShell, Pl, pl);
      RealType delThetaDouble = getDihedralAngle( Pi, Pj, Pk, Pl );
      
      // compute Ak + Al
      temp.makeCrossProduct( Pk-Pj, Pi-Pk );
      RealType vol = temp.norm() / 2.;
      temp.makeCrossProduct( Pl-Pi, Pj-Pl );
      vol += temp.norm() / 2.;
      RealType elengthSqr = VecType(Pj-Pi).normSqr();  
      
      // get deformed vertex positions
      getXYZCoord<VectorType, VecType>( defShell, Pi, pi);
      getXYZCoord<VectorType, VecType>( defShell, Pj, pj);
      getXYZCoord<VectorType, VecType>( defShell, Pk, pk);
      getXYZCoord<VectorType, VecType>( defShell, Pl, pl);
      
      // compute difference in dihedral angles
      delThetaDouble -= getDihedralAngle( Pi, Pj, Pk, Pl );
      delThetaDouble *= -2. * elengthSqr / vol;
      RealType factor = 2. * elengthSqr / vol;

      // compute first derivatives of dihedral angle
      VecType thetak, thetal, thetai, thetaj;
      getThetaGradK( Pi, Pj, Pk, thetak );
      getThetaGradK( Pj, Pi, Pl, thetal );
      getThetaGradI( Pi, Pj, Pk, Pl, thetai );
      getThetaGradJ( Pi, Pj, Pk, Pl, thetaj );

      // now compute second derivatives of dihedral angle
      MatType tensorProduct, H, aux;
      
      //kk
      getHessThetaKK( Pi, Pj, Pk, aux );
      tensorProduct.makeTensorProduct( thetak, thetak );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H);
      localToGlobal( tripletList, pk, pk, H ); 
            
      //ik & ki (Hki = Hik)
      getHessThetaIK( Pi, Pj, Pk, aux);
      tensorProduct.makeTensorProduct( thetai, thetak );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H);
      localToGlobal( tripletList, pi, pk, H );
      
      //jk & kj (Hkj = Hjk)
      getHessThetaJK( Pi, Pj, Pk, aux );
      tensorProduct.makeTensorProduct( thetaj, thetak );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H);
      localToGlobal( tripletList, pj, pk, H );       
      
      //ll
      getHessThetaKK( Pj, Pi, Pl, aux );
      tensorProduct.makeTensorProduct( thetal, thetal );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H);
      localToGlobal( tripletList, pl, pl, H );      
      
      //il & li (Hli = Hil)
      getHessThetaJK( Pj, Pi, Pl, aux);
      tensorProduct.makeTensorProduct( thetai, thetal );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H);
      localToGlobal( tripletList, pi, pl, H );
      
      //jl & lj (Hlj = Hjl)
      getHessThetaIK( Pj, Pi, Pl, aux);
      tensorProduct.makeTensorProduct( thetaj, thetal );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H);
      localToGlobal( tripletList, pj, pl, H );            
      
      //kl/lk: Hkl = 0 and Hlk = 0
      tensorProduct.makeTensorProduct( thetak, thetal );
      tensorProduct *= factor;
      localToGlobal( tripletList, pk, pl, tensorProduct );
        
      //ii  
      getHessThetaII( Pi, Pj, Pk, Pl, aux );
      tensorProduct.makeTensorProduct( thetai, thetai );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H);
      localToGlobal( tripletList, pi, pi, H );              

      //jj
      getHessThetaII( Pj, Pi, Pl, Pk, aux );       
      tensorProduct.makeTensorProduct( thetaj, thetaj );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H);
      localToGlobal( tripletList, pj, pj, H );

      //ij & ji (Hij = Hji)
      getHessThetaJI( Pi, Pj, Pk, Pl, H );     
      H *= delThetaDouble;
      tensorProduct.makeTensorProduct( thetai, thetaj );
      H.addMultiple( tensorProduct, factor );
      localToGlobal( tripletList, pi, pj, H );  

    }
  }

  void setWeight(double Eta)
  {
	  _factor = Eta;
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
//! \brief First derivative of SimpleBendingEnergy w.r.t. the undeformed configuration
//! \author Heeren
template<typename ConfiguratorType>
class SimpleBendingGradientUndef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {

  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::MatType    MatType;

  const MeshTopologySaver& _topology;
  const VectorType&  _defShell;
  const RealType _weight;

public:
SimpleBendingGradientUndef( const MeshTopologySaver& topology,
                            const VectorType& defShell,
                            const RealType Weight = 1. ) : _topology( topology), _defShell(defShell), _weight(Weight) {}

void apply( const VectorType& undefShell, VectorType& Dest ) const {

  if( undefShell.size() != _defShell.size() ){
      std::cerr << "size of undef = " << undefShell.size() << " vs. size of def = " << _defShell.size() << std::endl;
      throw BasicException( "SimpleBendingGradientUndef::apply(): sizes dont match!");
  }
      
  if( Dest.size() != undefShell.size() )
    Dest.resize( undefShell.size() );

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

      // first get defomed quantities
      VecType Pi, Pj, Pk, Pl, temp;
      getXYZCoord<VectorType, VecType>( _defShell, Pi, pi);
      getXYZCoord<VectorType, VecType>( _defShell, Pj, pj);
      getXYZCoord<VectorType, VecType>( _defShell, Pk, pk);
      getXYZCoord<VectorType, VecType>( _defShell, Pl, pl);
      RealType delTheta = getDihedralAngle( Pi, Pj, Pk, Pl );

      //!now get the undeformed values
      getXYZCoord<VectorType, VecType>( undefShell, Pi, pi);
      getXYZCoord<VectorType, VecType>( undefShell, Pj, pj);
      getXYZCoord<VectorType, VecType>( undefShell, Pk, pk);
      getXYZCoord<VectorType, VecType>( undefShell, Pl, pl);

      // compute Ak + Al, |e|^2 and diff. o dihedral angles
      temp.makeCrossProduct( Pk-Pj, Pi-Pk );
      RealType vol = temp.norm() / 2.;
      temp.makeCrossProduct( Pl-Pi, Pj-Pl );
      vol += temp.norm() / 2.;
      RealType elengthSqr( dotProduct(Pj-Pi,Pj-Pi) ); 
      // note signe here!
      delTheta -= getDihedralAngle( Pi, Pj, Pk, Pl );      

      // derivatives    
      VecType gradk, gradl, gradi, gradj, gradTheta, gradArea;
      RealType factorGradTheta = -2. * delTheta * elengthSqr / vol;    
      RealType factorGradArea = -1. * delTheta * delTheta * elengthSqr / (vol*vol);
      RealType factorGradEdgeLengthSqr = 2. * delTheta * delTheta / vol;   
         
      // d_k
      getThetaGradK( Pi, Pj, Pk, gradTheta );
      getAreaGradK( Pi, Pj, Pk, gradArea );
      getWeightedVectorSum<RealType>( factorGradTheta, gradTheta, factorGradArea, gradArea, gradk );
      
      // d_l
      getThetaGradK( Pj, Pi, Pl, gradTheta );
      getAreaGradK( Pj, Pi, Pl, gradArea );
      getWeightedVectorSum<RealType>( factorGradTheta, gradTheta, factorGradArea, gradArea, gradl );
      
      // d_i
      getThetaGradI( Pi, Pj, Pk, Pl, gradTheta );
      getAreaGradK( Pj, Pk, Pi, gradArea );
      getWeightedVectorSum<RealType>( factorGradTheta, gradTheta, factorGradArea, gradArea, gradi );
      getAreaGradK( Pl, Pj, Pi, gradArea );
      gradi.addMultiple( gradArea, factorGradArea );
      gradi.addMultiple( Pi-Pj, factorGradEdgeLengthSqr );
      
      // d_j
      getThetaGradJ( Pi, Pj, Pk, Pl, gradTheta );
      getAreaGradK( Pk, Pi, Pj, gradArea );
      getWeightedVectorSum<RealType>( factorGradTheta, gradTheta, factorGradArea, gradArea, gradj );
      getAreaGradK( Pi, Pl, Pj, gradArea );
      gradj.addMultiple( gradArea, factorGradArea );
      gradj.addMultiple( Pj-Pi, factorGradEdgeLengthSqr );      
      
      // assemble in global vector
      for( int i = 0; i < 3; i++ ){
        Dest[i*_topology.getNumVertices() + pi] += _weight * gradi[i];
        Dest[i*_topology.getNumVertices() + pj] += _weight * gradj[i];
        Dest[i*_topology.getNumVertices() + pk] += _weight * gradk[i];
        Dest[i*_topology.getNumVertices() + pl] += _weight * gradl[i];
      }

  }
}

};

//! \brief Second derivative of SimpleBendingEnergy w.r.t. the undeformed configuration
//! \author Heeren
template<typename ConfiguratorType>
class SimpleBendingHessianUndef : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

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
  mutable int _rowOffset, _colOffset;
  
public:
  SimpleBendingHessianUndef( const MeshTopologySaver& topology,
                             const VectorType& defShell,
			     const RealType Factor = 1., 
                             int rowOffset = 0, 
                             int colOffset = 0 ) : _topology( topology), _defShell(defShell), _factor( Factor ), _rowOffset(rowOffset),  _colOffset(colOffset){}
    
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
    // per edge we have 4 active vertices, i.e. 16 combinations each producing a 3x3-matrix
    tripletList.reserve( 16 * 9 * _topology.getNumEdges() );
    
    // fill matrix from triplets
    pushTriplets( undefShell, tripletList );    
    Hessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
  
  //
  void pushTriplets( const VectorType& undefShell, TripletListType& tripletList ) const {
          
    if( undefShell.size() != _defShell.size() ){
      std::cerr << "size of undef = " << undefShell.size() << " vs. size of def = " << _defShell.size() << std::endl;
      throw BasicException( "SimpleBendingHessianUndef::pushTriplets(): sizes dont match!");
    }
    
    // fill triplets
    for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){
      
      if( !(this->_topology.isEdgeValid(edgeIdx)) )
	continue;

      int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
          pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ),
          pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );

      // no bending at boundary edges
      if( std::min( pl, pk) < 0 )
        continue;

      // set up vertices and edges
      VecType Pi, Pj, Pk, Pl, temp;

      // first get defomed quantities
      getXYZCoord<VectorType, VecType>( _defShell, Pi, pi);
      getXYZCoord<VectorType, VecType>( _defShell, Pj, pj);
      getXYZCoord<VectorType, VecType>( _defShell, Pk, pk);
      getXYZCoord<VectorType, VecType>( _defShell, Pl, pl);
      RealType delTheta = getDihedralAngle( Pi, Pj, Pk, Pl );

      
      // get undeformed vertex positions
      getXYZCoord<VectorType, VecType>( undefShell, Pi, pi);
      getXYZCoord<VectorType, VecType>( undefShell, Pj, pj);
      getXYZCoord<VectorType, VecType>( undefShell, Pk, pk);
      getXYZCoord<VectorType, VecType>( undefShell, Pl, pl);    
            
      // compute Ak + Al, |e|^2 and dihedral angles
      temp.makeCrossProduct( Pk-Pj, Pi-Pk );
      RealType vol = temp.norm() / 2.;
      temp.makeCrossProduct( Pl-Pi, Pj-Pl );
      vol += temp.norm() / 2.;
      RealType elengthSqr = VecType(Pj-Pi).normSqr();  
      // note the sign!
      delTheta -= getDihedralAngle( Pi, Pj, Pk, Pl );
      delTheta *= -1.;
      
      // compute first derivatives of dihedral angle
      VecType thetak, thetal, thetai, thetaj;
      getThetaGradK( Pi, Pj, Pk, thetak );
      getThetaGradK( Pj, Pi, Pl, thetal );
      getThetaGradI( Pi, Pj, Pk, Pl, thetai );
      getThetaGradJ( Pi, Pj, Pk, Pl, thetaj );
      
      // compute first derivatives of area
      VecType areak, areal, areai, areaj;
      getAreaGradK( Pi, Pj, Pk, areak );
      getAreaGradK( Pj, Pi, Pl, areal );
      getAreaGradK( Pj, Pk, Pi, areai );
      getAreaGradK( Pl, Pj, Pi, temp );
      areai += temp;
      getAreaGradK( Pk, Pi, Pj, areaj );
      getAreaGradK( Pi, Pl, Pj, temp );
      areaj += temp;

      // now compute second derivatives
      MatType H, auxMat;
      VecType auxVec, e(Pj-Pi);
      
      //*k
      getWeightedVectorSum<RealType>( elengthSqr, thetak, -1. * delTheta * elengthSqr / vol, areak,  auxVec );
      getWeightedVectorSum<RealType>( 2., thetak, -1. * delTheta / vol, areak,  temp );
      
      //kk      
      H.makeTensorProduct( thetak, auxVec );            
      getHessThetaKK( Pi, Pj, Pk, auxMat );
      H.addMultiple( auxMat,  delTheta * elengthSqr );     
      auxMat.makeTensorProduct( areak, auxVec );
      H.addMultiple( auxMat, -1. * delTheta / vol );  
      getHessAreaKK( Pi, Pj, Pk, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );           
      H *= 2./vol;
      localToGlobal( tripletList, pk, pk, H ); 
      
      //lk      
      H.makeTensorProduct( thetal, auxVec );            
      auxMat.makeTensorProduct( areal, auxVec );
      H.addMultiple( auxMat, -1. * delTheta / vol );               
      H *= 2./vol;
      localToGlobal( tripletList, pl, pk, H ); 

      //ik
      H.makeTensorProduct( thetai, auxVec );            
      getHessThetaIK( Pi, Pj, Pk, auxMat );
      H.addMultiple( auxMat,  delTheta * elengthSqr );      
      auxMat.makeTensorProduct( areai, auxVec );
      H.addMultiple( auxMat, -1. * delTheta / vol );    
      getHessAreaIK( Pi, Pj, Pk, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );     
      auxMat.makeTensorProduct( e, temp );
      H.addMultiple( auxMat, -1.*delTheta );      
      H *= 2./vol;
      localToGlobal( tripletList, pi, pk, H ); 

      //jk
      H.makeTensorProduct( thetaj, auxVec );            
      getHessThetaJK( Pi, Pj, Pk, auxMat );
      H.addMultiple( auxMat,  delTheta * elengthSqr );      
      auxMat.makeTensorProduct( areaj, auxVec );
      H.addMultiple( auxMat, -1. * delTheta / vol );   
      getHessAreaIK( Pj, Pi, Pk, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );    
      auxMat.makeTensorProduct( e, temp );
      H.addMultiple( auxMat, delTheta );         
      H *= 2./vol;
      localToGlobal( tripletList, pj, pk, H ); 
      
      //*l
      getWeightedVectorSum<RealType>( elengthSqr, thetal, -1. * delTheta * elengthSqr / vol, areal,  auxVec );
      getWeightedVectorSum<RealType>( 2., thetal, -1. * delTheta / vol, areal,  temp );
      
      //ll      
      H.makeTensorProduct( thetal, auxVec );            
      getHessThetaKK( Pj, Pi, Pl, auxMat );
      H.addMultiple( auxMat,  delTheta * elengthSqr );     
      auxMat.makeTensorProduct( areal, auxVec );
      H.addMultiple( auxMat, -1. * delTheta / vol ); 
      getHessAreaKK( Pj, Pi, Pl, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );           
      H *= 2./vol;
      localToGlobal( tripletList, pl, pl, H ); 
      
      //il
      H.makeTensorProduct( thetai, auxVec );            
      getHessThetaJK( Pj, Pi, Pl, auxMat );
      H.addMultiple( auxMat,  delTheta * elengthSqr );      
      auxMat.makeTensorProduct( areai, auxVec );
      H.addMultiple( auxMat, -1. * delTheta / vol );   
      getHessAreaIK( Pi, Pj, Pl, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );     
      auxMat.makeTensorProduct( e, temp );
      H.addMultiple( auxMat, -1.*delTheta );  
      H *= 2./vol;
      localToGlobal( tripletList, pi, pl, H ); 
      
      //jl
      H.makeTensorProduct( thetaj, auxVec );            
      getHessThetaIK( Pj, Pi, Pl, auxMat );
      H.addMultiple( auxMat,  delTheta * elengthSqr );      
      auxMat.makeTensorProduct( areaj, auxVec );
      H.addMultiple( auxMat, -1. * delTheta / vol );   
      getHessAreaIK( Pj, Pi, Pl, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );   
      auxMat.makeTensorProduct( e, temp );
      H.addMultiple( auxMat, delTheta );  
      H *= 2./vol;
      localToGlobal( tripletList, pj, pl, H ); 
      
      //*j
      getWeightedVectorSum<RealType>( elengthSqr, thetaj, -1. * delTheta * elengthSqr / vol, areaj,  auxVec );
      auxVec.addMultiple( e, 2.*delTheta );// caution with factor 2!!!!!!!
      getWeightedVectorSum<RealType>( 2., thetaj, -1. * delTheta / vol, areaj,  temp );
     
      //jj     
      H.makeTensorProduct( thetaj, auxVec );   
      getHessThetaII( Pj, Pi, Pl, Pk, auxMat );
      H.addMultiple( auxMat,  delTheta * elengthSqr );   
      auxVec.addMultiple( e, -1.*delTheta );
      auxMat.makeTensorProduct( areaj, auxVec );
      H.addMultiple( auxMat, -1. * delTheta / vol );       
      getHessAreaKK( Pk, Pi, Pj, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );  
      getHessAreaKK( Pi, Pl, Pj, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol ); 
      auxMat.makeTensorProduct( e, temp );
      H.addMultiple( auxMat, delTheta ); 
      H.addToDiagonal( delTheta * delTheta );
      H *= 2./vol;
      localToGlobal( tripletList, pj, pj, H );
      
      //ij     
      auxVec.addMultiple( e, delTheta );
      H.makeTensorProduct( thetai, auxVec );   
      getHessThetaJI( Pi, Pj, Pk, Pl, auxMat );
      H.addMultiple( auxMat,  delTheta * elengthSqr );
      auxVec.addMultiple( e, -1.*delTheta );
      auxMat.makeTensorProduct( areai, auxVec );
      H.addMultiple( auxMat, -1. * delTheta / vol );       
      getHessAreaIK( Pi, Pk, Pj, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );  
      getHessAreaIK( Pi, Pl, Pj, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );      
      auxMat.makeTensorProduct( e, temp );
      H.addMultiple( auxMat, -1.*delTheta ); 
      H.addToDiagonal( -1. * delTheta * delTheta );
      H *= 2./vol;
      localToGlobal( tripletList, pi, pj, H );
      
      //*i
      getWeightedVectorSum<RealType>( elengthSqr, thetai, -1. * delTheta * elengthSqr / vol, areai,  auxVec );
      auxVec.addMultiple( e, -2.*delTheta ); // caution with factor 2!!!!!!!
      getWeightedVectorSum<RealType>( 2., thetai, -1. * delTheta / vol, areai,  temp );
      
      //ii     
      H.makeTensorProduct( thetai, auxVec );   
      getHessThetaII( Pi, Pj, Pk, Pl, auxMat );
      H.addMultiple( auxMat,  delTheta * elengthSqr ); 
      auxVec.addMultiple( e, 1.*delTheta );
      auxMat.makeTensorProduct( areai, auxVec );
      H.addMultiple( auxMat, -1. * delTheta / vol );       
      getHessAreaKK( Pl, Pj, Pi, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );  
      getHessAreaKK( Pj, Pk, Pi, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );  
      auxMat.makeTensorProduct( e, temp );
      H.addMultiple( auxMat, -1.*delTheta ); 
      H.addToDiagonal( delTheta * delTheta );
      H *= 2./vol;
      localToGlobal( tripletList, pi, pi, H );     
    }
    
  }

protected:
  void localToGlobal( TripletListType& tripletList, int k, int l, const MatType& localMatrix ) const {
    int numV = _topology.getNumVertices();
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          tripletList.push_back( TripletType( _rowOffset + i*numV + k, _colOffset + j*numV + l, _factor * localMatrix(i,j) ) );	
	
    if( k != l)
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          tripletList.push_back( TripletType( _rowOffset + i*numV + l, _colOffset + j*numV + k, _factor * localMatrix(j,i) ) );
  }
  
};


//! \brief Mixed second derivative of SimpleBendingEnergy
//! \author Heeren
template<typename ConfiguratorType>
class SimpleBendingHessianMixed : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

protected:
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
SimpleBendingHessianMixed( const MeshTopologySaver& topology,
                           const VectorType& InactiveGeometry,
		           const bool ActiveShellIsDeformed,
			   const bool FirstDerivWRTDef,
			   const RealType Factor = 1., 
                           int rowOffset = 0, 
                           int colOffset = 0 ) : _topology( topology), _inactiveGeometry(InactiveGeometry), _activeShellIsDeformed(ActiveShellIsDeformed), _firstDerivWRTDef( FirstDerivWRTDef ), _factor( Factor ),_rowOffset(rowOffset), _colOffset(colOffset){}
    
  void setRowOffset( int rowOffset ) const {
        _rowOffset = rowOffset;
  }
    
  void setColOffset( int colOffset ) const {
        _colOffset = colOffset;
  }

  //
  void apply( const VectorType& ActiveGeometry, MatrixType& Dest ) const {    
    int dofs = 3*_topology.getNumVertices();
    if( (Dest.rows() != dofs) || (Dest.cols() != dofs) )
        Dest.resize( dofs, dofs );
    Dest.setZero();
    assembleHessian( ActiveGeometry, Dest );
  }
  
  // assmeble Hessian matrix via triplet list
  void assembleHessian( const VectorType& ActiveGeometry, MatrixType& Hessian ) const {
      
    // set up triplet list
    TripletListType tripletList;
    // per edge we have 4 active vertices, i.e. 16 combinations each producing a 3x3-matrix
    tripletList.reserve( 16 * 9 * _topology.getNumEdges() );
    
    pushTriplets( ActiveGeometry, tripletList );

    // fill matrix from triplets
    Hessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
  
  //
  void pushTriplets( const VectorType& ActiveGeometry, TripletListType& tripletList ) const {
      
    const VectorType* undefShellP = _activeShellIsDeformed ? &_inactiveGeometry : &ActiveGeometry;
    const VectorType* defShellP   = _activeShellIsDeformed ? &ActiveGeometry    : &_inactiveGeometry;
    
        // fill triplets
    for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){
      
      if( !(this->_topology.isEdgeValid(edgeIdx)) )
	continue;

      int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
          pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ),
          pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );
          
      std::vector<int> idx; 
      idx.push_back( pi );
      idx.push_back( pj );
      idx.push_back( pk );
      idx.push_back( pl );

      // no bending at boundary edges
      if( std::min( pl, pk) < 0 )
        continue;

      // set up vertices and edges
      VecType Pi, Pj, Pk, Pl, gradArea, gradTheta, temp;
      std::vector<VecType> undefGrads(4), defGrads(4); 

      // first get deformed quantities
      getXYZCoord<VectorType, VecType>( *defShellP, Pi, pi);
      getXYZCoord<VectorType, VecType>( *defShellP, Pj, pj);
      getXYZCoord<VectorType, VecType>( *defShellP, Pk, pk);
      getXYZCoord<VectorType, VecType>( *defShellP, Pl, pl);
      RealType delTheta = getDihedralAngle( Pi, Pj, Pk, Pl );
      
      getThetaGradI( Pi, Pj, Pk, Pl, defGrads[0] );
      getThetaGradJ( Pi, Pj, Pk, Pl, defGrads[1] );
      getThetaGradK( Pi, Pj, Pk, defGrads[2] );
      getThetaGradK( Pj, Pi, Pl, defGrads[3] );

      
      // get undeformed vertex positions
      getXYZCoord<VectorType, VecType>( *undefShellP, Pi, pi);
      getXYZCoord<VectorType, VecType>( *undefShellP, Pj, pj);
      getXYZCoord<VectorType, VecType>( *undefShellP, Pk, pk);
      getXYZCoord<VectorType, VecType>( *undefShellP, Pl, pl);   
            
      // compute Ak + Al, |e|^2 and diff. o dihedral angles
      temp.makeCrossProduct( Pk-Pj, Pi-Pk );
      RealType vol = temp.norm() / 2.;
      temp.makeCrossProduct( Pl-Pi, Pj-Pl );
      vol += temp.norm() / 2.;
      RealType elengthSqr( dotProduct(Pj-Pi,Pj-Pi) ); 
      // note signe here!
      delTheta -= getDihedralAngle( Pi, Pj, Pk, Pl );      

      // factors   
      RealType factorGradTheta = -2. * elengthSqr / vol;    
      RealType factorGradArea = -2. * delTheta * elengthSqr / (vol*vol);
      RealType factorGradEdgeLengthSqr = 4. * delTheta / vol;       
            
      // d_i
      getThetaGradI( Pi, Pj, Pk, Pl, gradTheta );
      getAreaGradK( Pj, Pk, Pi, gradArea );
      getWeightedVectorSum<RealType>( factorGradTheta, gradTheta, factorGradArea, gradArea, undefGrads[0] );
      getAreaGradK( Pl, Pj, Pi, gradArea );
      undefGrads[0].addMultiple( gradArea, factorGradArea );
      undefGrads[0].addMultiple( Pi-Pj, factorGradEdgeLengthSqr );
      
      // d_j
      getThetaGradJ( Pi, Pj, Pk, Pl, gradTheta );
      getAreaGradK( Pk, Pi, Pj, gradArea );
      getWeightedVectorSum<RealType>( factorGradTheta, gradTheta, factorGradArea, gradArea, undefGrads[1] );
      getAreaGradK( Pi, Pl, Pj, gradArea );
      undefGrads[1].addMultiple( gradArea, factorGradArea );
      undefGrads[1].addMultiple( Pj-Pi, factorGradEdgeLengthSqr );      
      
      // d_k
      getThetaGradK( Pi, Pj, Pk, gradTheta );
      getAreaGradK( Pi, Pj, Pk, gradArea );
      getWeightedVectorSum<RealType>( factorGradTheta, gradTheta, factorGradArea, gradArea, undefGrads[2] );
      
      // d_l
      getThetaGradK( Pj, Pi, Pl, gradTheta );
      getAreaGradK( Pj, Pi, Pl, gradArea );
      getWeightedVectorSum<RealType>( factorGradTheta, gradTheta, factorGradArea, gradArea, undefGrads[3] );
 
      // local to global
      for( int m = 0; m < 4; m++ ){
        for( int n = 0; n < 4; n++ ){
          MatType matrix;
          matrix.makeTensorProduct( defGrads[m], undefGrads[n] );
          localToGlobal( tripletList, idx[m], idx[n], matrix );
        }
      }

      
    }
    
  }

protected:
  void localToGlobal( TripletListType& tripletList, int k, int l, const MatType& localMatrix ) const {
    int numV = _topology.getNumVertices();
    if( _firstDerivWRTDef )   
    {
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          tripletList.push_back( TripletType( _rowOffset + i*numV + k, _colOffset + j*numV + l, _factor * localMatrix(i,j) ) );
    }
    else
    {
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          tripletList.push_back( TripletType( _rowOffset + i*numV + l, _colOffset + j*numV + k, _factor * localMatrix(j,i) ) );
    }
  }
  
};


//==========================================================================================================

//! \brief Deformation class for SimpleBendingEnergy
//! \author Heeren
template <typename ConfiguratorType>
class SimpleBendingDeformation : public DeformationBase<ConfiguratorType>{
    
protected:   
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;

  typedef typename ConfiguratorType::TripletType TripletType;  
  typedef std::vector<TripletType> TripletListType;
  
  const MeshTopologySaver& _topology;
  RealType _bendWeight;

public:
  SimpleBendingDeformation( const MeshTopologySaver& Topology, RealType bendWeight ) : _topology( Topology ), _bendWeight( bendWeight ) {}
  
  void applyEnergy ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, RealType & Dest ) const {
    SimpleBendingEnergy<ConfiguratorType>( _topology, UndeformedGeom, true ).apply( DeformedGeom, Dest );
    Dest *= _bendWeight;
  }
  
  void applyUndefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest ) const {
    SimpleBendingGradientUndef<ConfiguratorType>( _topology, DeformedGeom ).apply( UndeformedGeom, Dest );
    Dest *= _bendWeight;
  }
  
  void applyDefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest ) const {
    SimpleBendingGradientDef<ConfiguratorType>( _topology, UndeformedGeom ).apply( DeformedGeom, Dest );
    Dest *= _bendWeight;
  }
  
  void pushTripletsDefHessian ( const VectorType& UndefGeom, const VectorType& DefGeom, TripletListType& triplets, int rowOffset, int colOffset, RealType factor = 1.0 ) const {
    SimpleBendingHessianDef<ConfiguratorType>( _topology, UndefGeom, factor * _bendWeight, rowOffset, colOffset ).pushTriplets( DefGeom, triplets );
  }
   
  void pushTripletsUndefHessian ( const VectorType& UndefGeom, const VectorType& DefGeom, TripletListType& triplets, int rowOffset, int colOffset, RealType factor = 1.0 ) const {
    SimpleBendingHessianUndef<ConfiguratorType>( _topology, DefGeom, factor * _bendWeight, rowOffset, colOffset ).pushTriplets( UndefGeom, triplets );
  }
  
  // mixed second derivative of deformation energy E[S_1, S_2], i.e. if "FirstDerivWRTDef" we have D_1 D_2 E[.,.], otherwise D_2 D_1 E[.,.]
  void pushTripletsMixedHessian ( const VectorType& UndefGeom, const VectorType& DefGeom, TripletListType& triplets, int rowOffset, int colOffset, const bool FirstDerivWRTDef, RealType factor = 1.0 ) const {
    SimpleBendingHessianMixed<ConfiguratorType>( _topology, UndefGeom, true, FirstDerivWRTDef, factor * _bendWeight, rowOffset, colOffset ).pushTriplets( DefGeom, triplets );
  }
  
  int numOfNonZeroHessianEntries () const {
    // per edge we have 4 active vertices, i.e. 16 combinations each producing a 3x3-matrix
    return 16 * 9 * _topology.getNumEdges();
  }
    
};

#endif 
