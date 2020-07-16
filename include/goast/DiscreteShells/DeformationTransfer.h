// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Deformation transfer based on deformation gradients.
 * \author Heeren
 * \cite BoSuPa06
 */
#ifndef __DEFORMATIONTRANSFEROP_H
#define __DEFORMATIONTRANSFEROP_H

#include <goast/Optimization.h>

/**
 * \brief Gradient based deformation transfer based on \cite BoSuPa06 (see all \cite SuPo04 )
 * \author Heeren
 *
 * Let \f$ S \f$ be a source mesh and \f$ S' = \Phi(S) \f$ a deformed source mesh, let \f$ T \f$ be the target mesh.
 * Output is a deformed target mesh \f$ T' = \Phi(T) \f$.
 * All meshes share the same topology with n vertices and m faces.
 * Compute nodal positions \f$ q = (q_x, q_y, q_z) \f$ of \f$ T' \f$ by solving a least squares problem, i.e.
 * \f[ A[T]^T D[T] A[T] q_i = A[T]^T D[T] G_i(\Phi) \f] for  \f$ i \in \{x,y,z\} \f$ , with \f$ q_i \in \R^n \f$.
 * Here \f$ A[T] \in \R^{3m,n} \f$ and \f$ D[T] \in \R^{3m,3m} \f$ is a diagonal matrix containing face areas on the diagonal (repeated three times).
 * The gradients \f$ G_i \in \R^{3m} \f$ on the right depend on \f$ \Phi \f$, i.e. the geometry of the source \f$ S \f$ AND the deformed source \f$ S' \f$.
 */
template<typename ConfiguratorType>
class DeformationTransferOperator{
  
   typedef typename ConfiguratorType::RealType RealType;
   typedef typename ConfiguratorType::VectorType VectorType;    
   typedef typename ConfiguratorType::VecType VecType;    
   typedef typename ConfiguratorType::MatType MatType; 
   
   typedef typename ConfiguratorType::SparseMatrixType MatrixType;
   typedef typename ConfiguratorType::TripletType TripletType;  
   typedef std::vector<TripletType> TripletListType;
   
   
   const MeshTopologySaver& _topology;
   std::vector<int> _fixedIndices; 
   bool _quiet;   

public:  
  DeformationTransferOperator( const MeshTopologySaver& Topology, bool quiet = false )  : _topology( Topology ),  _quiet( quiet ){ }
  
  //! constructor with fixed handle indices and positions
  DeformationTransferOperator( const MeshTopologySaver& Topology, const std::vector<int>& fixedIndices, bool quiet = false ) 
  : _topology( Topology ), _fixedIndices(fixedIndices), _quiet( quiet ){ }
  
  //! perform deformation transfer
  void execute( const VectorType& source, const VectorType& defSource, const VectorType& target, VectorType& defTarget ) const {
    
    // matrices: gradient matrix G and diagonal area matrix D
    MatrixType gradMat, diagMat, systemMat;
    if( !_quiet ) std::cerr << "Assemble matrices..." << std::endl;
    // temporarily, the area matrix is stored in systemMat
    computeGradientSystemMatrix( target, gradMat, diagMat );
    MatrixType tempMat = gradMat.transpose() * diagMat;
    
    // compute source gradients
    if( !_quiet ) std::cerr << "Prepare rhs..." << std::endl;
    std::vector<VectorType> sourceGradients(3), rhs(3); 
    computeDeformationGradients( source, defSource, sourceGradients );
    // rhs = G^T D S, where S are the source gradients
    for( int i = 0; i < 3; i++ )
      rhs[i] = tempMat * sourceGradients[i];
/* 
    // TEST
    MatrixType gradMatSource, diagMatSource;    
    computeGradientSystemMatrix( source, gradMatSource, diagMatSource );
    for( int i = 0; i < 3; i++ )
      rhs[i] = tempMat * gradMatSource * defSource.segment( i * _topology.getNumVertices(), _topology.getNumVertices() );
*/
    // system mat A = G^T D G
    if( !_quiet ) std::cerr << "Assemble system matrix..." << std::endl;
    systemMat = tempMat * gradMat;
 
    // bc
    if( _fixedIndices.size() > 0 ){
      if( !_quiet ) std::cerr << "Manipulate system matrix to account for bc..." << std::endl;
      // manipulate matrix
      for( int i = 0; i < _fixedIndices.size(); i++ ){
        for (typename MatrixType::InnerIterator it(systemMat,_fixedIndices[i]); it; ++it)
          it.valueRef() = (it.row() == it.col()) ? 1. : 0.;
        systemMat.coeffRef(_fixedIndices[i], _fixedIndices[i]) =  1.;
      }       
      // prescribe fixed positions (obtained from target) on rhs   
      for( int k = 0; k < _fixedIndices.size(); k++ ){
          VecType targetNode;
          getXYZCoord<VectorType, VecType>(target, targetNode, _fixedIndices[k]);
          for( int j = 0; j < 3; j++ )
            rhs[j][ _fixedIndices[k] ] = targetNode[j];
      }
    }

    // solving
    int dofs = _topology.getNumVertices();
    defTarget.resize( 3 * dofs );
    if( !_quiet ) std::cerr << "Factorize system matrix..."; 
    auto t_start = std::chrono::high_resolution_clock::now();
    LinearSolver<ConfiguratorType> directSolver; 
    directSolver.prepareSolver( systemMat );
    auto t_end = std::chrono::high_resolution_clock::now();
    if( !_quiet ) std::cerr<<"done in " << std::chrono::duration<double, std::ratio<1> >(t_end-t_start).count() << " seconds." << std::endl;
    VectorType temp(dofs);
    //TODO avoid copying by means of references!
    for( int j = 0; j < 3; j++ ){
      directSolver.backSubstitute( rhs[j], temp );
      for( int k = 0; k < dofs; k++ )
          defTarget[j*dofs + k] = temp[k];
    }
    if( !_quiet ) std::cerr << "Deformation transfer done." << std::endl;
  }
 
  //! perform inverse deformation transfer
  void executeInverse( const VectorType& target, const VectorType& defTarget, const VectorType& source, VectorType& defSource ) const {
    
    int numVerts = _topology.getNumVertices();  
      
    // matrices: gradient matrix G and diagonal area matrix D
    if( !_quiet ) std::cerr << "Assemble target and source matrices..." << std::endl;  
    MatrixType gradMatTarget, diagMatTarget;    
    computeGradientSystemMatrix( target, gradMatTarget, diagMatTarget );
    MatrixType gradMatSource, diagMatSource;    
    computeGradientSystemMatrix( source, gradMatSource, diagMatSource );
    
    // A[T]^T * D[T]
    MatrixType tempMat = gradMatTarget.transpose() * diagMatTarget;
    
    // compute rhs = A[T]^T * D[T] * A[T] * T'
    if( !_quiet ) std::cerr << "Prepare rhs..." << std::endl;
    // get x-, y-, z-components of target
    std::vector<VectorType> rhs(3); 
    for( int i = 0; i < 3; i++ )
        rhs[i] = tempMat * gradMatTarget * defTarget.segment( i * numVerts, numVerts );

    // system mat A = A[T]^T * D[T] * A[S]
    if( !_quiet ) std::cerr << "Assemble system matrix..." << std::endl;
    MatrixType systemMat = tempMat * gradMatSource;
 
    // bc
    if( _fixedIndices.size() > 0 ){
      if( !_quiet ) std::cerr << "Manipulate system matrix to account for bc..." << std::endl;
      // manipulate matrix
      for( int i = 0; i < _fixedIndices.size(); i++ ){
        for (typename MatrixType::InnerIterator it(systemMat,_fixedIndices[i]); it; ++it)
          it.valueRef() = (it.row() == it.col()) ? 1. : 0.;
        systemMat.coeffRef(_fixedIndices[i], _fixedIndices[i]) =  1.;
      }       
      // prescribe fixed positions (obtained from target) on rhs   
      for( int k = 0; k < _fixedIndices.size(); k++ ){
          VecType sourceNode;
          getXYZCoord<VectorType, VecType>(source, sourceNode, _fixedIndices[k]);
          for( int j = 0; j < 3; j++ )
            rhs[j][ _fixedIndices[k] ] = sourceNode[j];
      }
    }

    // solving
    defSource.resize( 3 * numVerts );
    if( !_quiet ) std::cerr << "Factorize system matrix..."; 
    auto t_start = std::chrono::high_resolution_clock::now();
    LinearSolver<ConfiguratorType> directSolver; 
    directSolver.prepareSolver( systemMat );
    auto t_end = std::chrono::high_resolution_clock::now();
    if( !_quiet ) std::cerr<<"done in " << std::chrono::duration<double, std::ratio<1> >(t_end-t_start).count() << " seconds." << std::endl;
    VectorType solutionComp(numVerts);
    //TODO avoid copying by means of references!
    for( int j = 0; j < 3; j++ ){
      directSolver.backSubstitute( rhs[j], solutionComp );
      for( int k = 0; k < numVerts; k++ )
          defSource[j*numVerts + k] = solutionComp[k];
    }
    if( !_quiet ) std::cerr << "Inverse deformation transfer done." << std::endl;
  }
 
/* 
  //! test two optons to compute deformation gradients
  void test( const VectorType& source, const VectorType& defSource ) const {

      int numVerts = _topology.getNumVertices();
    if( !_quiet ) std::cerr << "Start test..." << std::endl;  
    MatrixType gradMatSource, diagMatSource;    
    computeGradientSystemMatrix( source, gradMatSource, diagMatSource );

    std::vector<VectorType> sourceGradients(3); 
    computeDeformationGradients( source, defSource, sourceGradients );
    for( int i = 0; i < 3; i++ ){
        sourceGradients[i] -= gradMatSource * defSource.segment( i * numVerts, numVerts );
        std::cerr << i << "th comp error = " << sourceGradients[i].norm() << std::endl;
    }
  }
*/

  //! compute deformation gradients G_i[\Phi], i \in \{x,y,z\} of source mesh and deformed source mesh
  void computeDeformationGradients( const VectorType& source, const VectorType& defSource, std::vector<VectorType>& defGradientsTrans ) const {
    
    if( defGradientsTrans.size() != 3 )
      defGradientsTrans.resize(3);
    for( int i = 0; i < 3 ; i++ )
      defGradientsTrans[i].resize( 3 * _topology.getNumFaces() );
    
    // run over all faces
    for( int faceIdx = 0; faceIdx < _topology.getNumFaces(); faceIdx++ ){
      
      // get nodal positions of triangle
      std::vector< VecType > p(3), q(3);
      for( int j = 0; j < 3; j++ ){
        getXYZCoord<VectorType, VecType>( source,    p[j], _topology.getNodeOfTriangle(faceIdx,j) );
	getXYZCoord<VectorType, VecType>( defSource, q[j], _topology.getNodeOfTriangle(faceIdx,j) );
      }
      
      // compute edges ei = pi - p0, fi = qi - q0, i = 1,2
      p[1] -= p[0]; p[2] -= p[0];
      q[1] -= q[0]; q[2] -= q[0];     
            
      // compute normals
      p[0].makeCrossProduct( p[1], p[2] );
      p[0] /= p[0].norm();
      q[0].makeCrossProduct( q[1], q[2] );
      q[0] /= q[0].norm();
      
      // compute inverse matrix B = (e1, e2, pn)^{-1}
      MatType pMat, qMat;
      for( int j = 0; j < 3; j++ )
        pMat.setCol( j, p[(j+1)%3] );
      
      // compute matrix A = (f1, f2, qn)
      for( int j = 0; j < 3; j++ )
        qMat.setCol( j, q[(j+1)%3] );
      // C = A * B
      //std::cerr<< faceIdx << " " << pMat.det() << std::endl;
      qMat *= pMat.inverse();   
      
      // local to global, note transposition here!
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
	  defGradientsTrans[j][3*faceIdx+i] = qMat.get(j,i);
      
    }
  }
    
  // compute system matrix A[T] and diagonal Matrix D[T] corresponding to geometry of the target mesh T 
  void computeGradientSystemMatrix( const VectorType& targetGeometry, MatrixType& systemMat, MatrixType& diagMat ) const {
      
      // allocation
      if( (diagMat.rows() != 3 * _topology.getNumFaces()) || (diagMat.cols() != 3 * _topology.getNumFaces() ) )
        diagMat.resize( 3 * _topology.getNumFaces(), 3 * _topology.getNumFaces() );
      if( (systemMat.rows() != 3 * _topology.getNumFaces()) || (systemMat.cols() != _topology.getNumVertices() ) )
        systemMat.resize( 3 * _topology.getNumFaces(), _topology.getNumVertices() );
      
      // set up triplet lists
      TripletListType diagTriplets, systemTriplets;
      diagTriplets.reserve( 3 * _topology.getNumFaces() );
      systemTriplets.reserve( 9 * _topology.getNumFaces() );
	  
      //     |-1 1 0 |
      // C = |-1 0 1 |
      //     | 0 0 0 |
      MatType constMat(-1., 1., 0., -1., 0., 1., 0., 0., 0.);
      
      // run over all faces
      for( int faceIdx = 0; faceIdx < _topology.getNumFaces(); faceIdx++ ){
	
	// compute base function gradient
	MatType baseFncGradient, temp;
	
	// get nodal positions pi, edges and unit normal n
	std::vector< VecType > nodalPos( 3 );
	for( int i = 0; i < 3; i++ ){
	  getXYZCoord<VectorType, VecType>( targetGeometry, nodalPos[i],  _topology.getNodeOfTriangle( faceIdx, i ) );
	  if( i > 0 )
	    nodalPos[i] -= nodalPos[0];
	}
	// write unit normal in first entyr
	nodalPos[0].makeCrossProduct( nodalPos[1], nodalPos[2] );
	RealType vol = nodalPos[0].norm();
	nodalPos[0] /= vol;
	vol /= 2.;
	
	// write area volume of triangle in diagonal matrix
	for( int j = 0; j < 3; j++ )
          diagTriplets.push_back( TripletType( 3 * faceIdx + j, 3 * faceIdx + j, vol ) );
	
	// assemble transposed (!!) local matrix
	for( int i = 0; i < 3; i++ )
	  temp.setCol( i, nodalPos[(i+1)%3] );
	
	// B = [ p1-p0 | p2-p0 | n ]^{-T} * C
	baseFncGradient.makeProductAtransposedB( temp.inverse(), constMat );	
	
	// local to global
	for( int i = 0; i < 3; i++ )
	  for( int j = 0; j < 3; j++ )
	      systemTriplets.push_back( TripletType( 3 * faceIdx + i, _topology.getNodeOfTriangle( faceIdx, j ), baseFncGradient.get(i,j) ) );
	
      }
      
      diagMat.setFromTriplets( diagTriplets.cbegin(), diagTriplets.cend() );
      systemMat.setFromTriplets( systemTriplets.cbegin(), systemTriplets.cend() );
  }
 
};

#endif