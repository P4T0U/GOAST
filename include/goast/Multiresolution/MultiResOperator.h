// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef __MULTIRESOPERATOR_H
#define __MULTIRESOPERATOR_H

#include "MeshDecimation.h"
#include <goast/Optimization.h>

//! \brief Multi-resolution operator based on deformation transfer as described in \cite FrBo11 (see also appendix of Behrend's thesis)
//! \author Heeren
//!\todo add documentation
template<typename ConfiguratorType, typename QuadricType = SimultaneousQuadricErrorMetric<ConfiguratorType> >
class MultiResolutionOperator{
  
   typedef typename ConfiguratorType::RealType RealType;
   typedef typename ConfiguratorType::VectorType VectorType;    
   typedef typename ConfiguratorType::VecType VecType;    
   typedef typename ConfiguratorType::MatType MatType; 
   
   typedef typename ConfiguratorType::SparseMatrixType MatrixType;
   typedef typename ConfiguratorType::TripletType TripletType;  
   typedef std::vector<TripletType> TripletListType;

   const MeshTopologySaver& _fineTopology;
   
   int _numOfVerticesInCoarseMesh;   
   VectorType _coarseRefGeom;
   mutable VectorType _fineBaseGeom; 
   std::vector<int> _anchors;   
   mutable bool _quiet, _factorized;
   mutable MatrixType _systemMat;
   mutable LinearSolver<ConfiguratorType> _directSolverBiLaplace;
   std::vector<int> _fixedVertexIndices;

public:  
  MultiResolutionOperator( const MeshTopologySaver& Topology, bool quiet = false )
  : _fineTopology( Topology ),
    _numOfVerticesInCoarseMesh( -1 ),
    _quiet( quiet ),
    _factorized( false ){ }
    
  MultiResolutionOperator( const MeshTopologySaver& Topology, const std::vector<int>& FixedVertexIndices, bool quiet = false ) 
  : _fineTopology( Topology ),
    _numOfVerticesInCoarseMesh( -1 ),
    _quiet( quiet ),
    _factorized( false ),
    _fixedVertexIndices(FixedVertexIndices){ }
  
    //! decimate meshes, compute prolongation and anchors
  void computeDecimation( const std::vector<TriMesh>& FineMeshes, 
			  RealType theta,
			  std::vector<TriMesh>& CoarseMeshes ) { 
      DecimationInfo<ConfiguratorType> decInfo;
      computeDecimation( FineMeshes, theta, CoarseMeshes, decInfo );
  }
  
  //! decimate meshes, compute prolongation and anchors
  void computeDecimation( const std::vector<TriMesh>& FineMeshes, 
			  RealType theta,
			  std::vector<TriMesh>& CoarseMeshes,
                          DecimationInfo<ConfiguratorType>& decInfo ) {    
    
    if( !_quiet ) std::cerr << "Decimation..." << std::endl;
    decInfo.resize( _fineTopology.getNumVertices() );
    SimultaneousDecimater<ConfiguratorType, QuadricType> simDec( FineMeshes, _fixedVertexIndices, _quiet );
    auto t_start = std::chrono::high_resolution_clock::now();
    simDec.calcDecimation( theta, decInfo ); 
    auto t_end = std::chrono::high_resolution_clock::now();
    if( !_quiet ) std::cerr<<"Done in " << std::chrono::duration<double, std::ratio<1> >(t_end-t_start).count() << " seconds." << std::endl;
    
    if( !_quiet ) std::cerr << "Flip edges with bad angles: " << simDec.flipEdgesWithBadAngles( 3.0 ) << std::endl;
    
    // get coarse reference mesh
    getGeometry( simDec.getCoarseMesh( 0 ), _coarseRefGeom );
    _numOfVerticesInCoarseMesh = simDec.getCoarseMesh( 0 ).n_vertices();    
    
    // get coarse meshes (omit first one which is the reference mesh!)
    CoarseMeshes.clear();
    for( uint j = 0; j < FineMeshes.size(); j++ )
      CoarseMeshes.push_back( simDec.getCoarseMesh( j ) );     
    
    // set anchors
    if( !_quiet ) std::cerr << "Set anchors..." << std::endl;
    _anchors.resize( _numOfVerticesInCoarseMesh );
    for( int i = 0; i < _numOfVerticesInCoarseMesh; i++ )
      _anchors[i] = decInfo.vertexMap[i];
  }

  void setCoarseData(const VectorType &coarseRefGeom, const std::vector<int> &vertexMap) {
    _coarseRefGeom = coarseRefGeom;
    _numOfVerticesInCoarseMesh = coarseRefGeom.size() / 3;
    if( !_quiet ) std::cerr << "Set anchors..." << std::endl;
    _anchors.resize( _numOfVerticesInCoarseMesh );
    for( int i = 0; i < _numOfVerticesInCoarseMesh; i++ )
      _anchors[i] = vertexMap[i];
  }
    
  const std::vector<int>& getAnchors() const {
    return _anchors;
  }

  const VectorType& getCoarseRefGeom() const {
    return _coarseRefGeom;
  }
  
  void setCoarseRefGeom( const VectorType& coarseRefGeom ) {
      _coarseRefGeom.resize( coarseRefGeom.size() );
      _coarseRefGeom = coarseRefGeom;
      _numOfVerticesInCoarseMesh = _coarseRefGeom[0].size();
  }
  
  // compute coarse boundary mask from fine boundary mask
  void computeDecimatedBoundaryMask( const BitVector& fineMask, BitVector& coarseMask ) const {
    coarseMask.resize( _numOfVerticesInCoarseMesh );
    coarseMask.setAll( false );
    for( int i = 0; i < _numOfVerticesInCoarseMesh; i++ )
      if( fineMask[_anchors[i]] )
	coarseMask.set( i, true );
  }
  
  void saveBaseMesh( std::string filename ) const {
    
    if( _fineBaseGeom.size() != 3*_fineTopology.getNumVertices() )
      throw BasicException ( "MultiResolutionOperator::saveBaseMesh(): wrong size!!" );
     
    TriMesh tempMesh( _fineTopology.getGrid() );
    setGeometry( tempMesh,  _fineBaseGeom );
    if( !OpenMesh::IO::write_mesh(tempMesh, filename) )
      throw BasicException ( "MultiResolutionOperator::saveBaseMesh(): could not write base mesh!!" );
  }
 
  
  //! prolongation by detail transfer
  void prolongate( const std::vector<VectorType>& defCoarse, std::vector<VectorType>& defTarget, bool useParallelComp = true ) const {   
    
    if( defCoarse.size() != defTarget.size() )
      throw BasicException( "MultiResolutionOperator::performDetailTransfer(): sizes don't match!!" );
    
    //auto t_start_total = std::chrono::high_resolution_clock::now();
    VectorType target;
    getGeometry<VectorType>( _fineTopology.getGrid(), target );  
        
    if( _fineBaseGeom.size() == 0 )
      computeFineBaseMesh();
    
    if( useParallelComp ){
      bool quiet = _quiet;
      if( !_quiet ) std::cerr << "\n-------------------------------" << std::endl;
      if( !_quiet ) std::cerr << "Compute " << defCoarse.size() << " detail transfer parallely..." << std::endl; 
      _quiet = true;
      
#ifdef _OPENMP
#pragma omp parallel for
#endif    
      for( int i = 0; i < defCoarse.size(); i++ ){
        VectorType deformedFineBaseMesh;
        solveThinPlateEquation( defCoarse[i], deformedFineBaseMesh ); 
        performDetailTransfer( _fineBaseGeom, deformedFineBaseMesh, target, defCoarse[i], defTarget[i] );
      }
      _quiet = quiet;
    }
    else{
      for( int i = 0; i < defCoarse.size(); i++ ){
        VectorType deformedFineBaseMesh;
        auto t_start = std::chrono::high_resolution_clock::now();

        if( !_quiet ) std::cerr << "\n-------------------------------" << std::endl;
        if( !_quiet ) std::cerr << "Compute deformed base mesh (B')..." << std::endl; 
        solveThinPlateEquation( defCoarse[i], deformedFineBaseMesh ); 
    
        if( !_quiet ) std::cerr << "\n-------------------------------" << std::endl;
        if( !_quiet ) std::cerr << "Perform detail transfer (M')..." << std::endl; 
        performDetailTransfer( _fineBaseGeom, deformedFineBaseMesh, target, defCoarse[i], defTarget[i] );    
	
        auto t_end = std::chrono::high_resolution_clock::now();       
        if( !_quiet ) std::cerr << "Computed B' and M' in " << std::chrono::duration<double, std::ratio<1> >(t_end-t_start).count() << " seconds." << std::endl << std::endl;
      }
    }
    //auto t_end_total = std::chrono::high_resolution_clock::now();
    //if( !_quiet ) std::cerr << "Total detail transfer done in " << std::chrono::duration<double, std::ratio<1> >(t_end_total-t_start_total).count() << " seconds." << std::endl << std::endl; 
    
  }
     
protected:   
  //! solve thin plate problem on fine mesh with anchors
  void computeFineBaseMesh(  ) const {
    
    if( !_quiet ) std::cerr << "\n-------------------------------" << std::endl;
    if( !_quiet ) std::cerr << "Compute base mesh (i.e. solve thin plate eq.)..." << std::endl; 
    
    if( _coarseRefGeom.size() == 0 )
      throw BasicException ( "MultiResolutionOperator::computeFineBaseMesh(): coarse reference mesh not set!!" );
    
    solveThinPlateEquation( _coarseRefGeom, _fineBaseGeom );   
  }
  
  //
  void solveThinPlateEquation( const VectorType& AnchorPositions, VectorType& Solution ) const {
    
    if( AnchorPositions.size() != 3*_numOfVerticesInCoarseMesh )  
      throw BasicException ( "MultiResolutionOperator::solveThinPlateEquation(): Passed AnchorPositions has wrong size!!" );
    
    if( _anchors.size() != _numOfVerticesInCoarseMesh )
      throw BasicException ( "MultiResolutionOperator::solveThinPlateEquation(): anchors have wrong size!!" );     
    
    auto t_start_total = std::chrono::high_resolution_clock::now();
    // assemble system matrix and factorize
    if( !_factorized ){        
      if( !_quiet ) std::cerr << "Assemble system matrix..." << std::endl;

      if( !_quiet ) std::cerr << "Use cotan Laplace approximation..." << std::endl;
      approximateBiLaplaceMatrixByCotanFormula ( _systemMat );
      
      // account for bc
      if( !_quiet ) std::cerr << "Manipulate system matrix to account for bc..." << std::endl;
      applyMaskToRow( _anchors, _systemMat );

      // factorize system amatrix
      if( !_quiet ) std::cerr << "Factorize system matrix...";
      auto t_start = std::chrono::high_resolution_clock::now();
      _directSolverBiLaplace.prepareSolver( _systemMat );
      _factorized = true;
      auto t_end = std::chrono::high_resolution_clock::now();
      if( !_quiet ) std::cerr<<"done in " << std::chrono::duration<double, std::ratio<1> >(t_end-t_start).count() << " seconds." << std::endl;
    }

    // solve system       
    if( !_quiet ) std::cerr << "Back substitution..." << std::endl;
     int numFineDofs = _fineTopology.getNumVertices();
    Solution.resize( 3 * numFineDofs );    
    for( int j = 0; j < 3; j++ ){
      // assemble rhs, i.e. fix all anchor positions
      VectorType rhs( numFineDofs ); 
      rhs.setZero(); 
      for( int k = 0; k < _numOfVerticesInCoarseMesh; k++ )	
	  rhs[ _anchors[k] ] = AnchorPositions[ j*_numOfVerticesInCoarseMesh + k ];
      // back substitute
      VectorType temp( numFineDofs );
      _directSolverBiLaplace.backSubstitute( rhs, temp );
      for( int k = 0; k < numFineDofs; k++ )
          Solution[j*numFineDofs + k] = temp[k];
    }
    
/*
    for( int j = 0; j < 3; j++ ){
      // assemble rhs, i.e. fix all anchor positions
      VectorType rhs( numFineDofs + _anchors.size() ); 
      rhs.setZero(); 
      for( int k = 0; k < _numOfVerticesInCoarseMesh; k++ )	
	  rhs[ numFineDofs + k ] = AnchorPositions[ j*_numOfVerticesInCoarseMesh + k ];
      // back substitute
      VectorType temp( numFineDofs + _anchors.size() );
      _directSolverBiLaplace.backSubstitute( rhs, temp );
      for( int k = 0; k < numFineDofs; k++ )
          Solution[j*numFineDofs + k] = temp[k];
    }
*/      
    auto t_end_total = std::chrono::high_resolution_clock::now();
    if( !_quiet ) std::cerr<<"Thin plate eq. solved in " << std::chrono::duration<double, std::ratio<1> >(t_end_total-t_start_total).count() << " seconds." << std::endl;
  }
  
 
  //! compute deformation gradients of source mesh and deformed source mesh
  void computeDeformationGradients( const VectorType& source, const VectorType& defSource, std::vector<VectorType>& defGradientsTrans ) const {
    
    if( defGradientsTrans.size() != 3 )
      defGradientsTrans.resize(3);
    for( int i = 0; i < 3 ; i++ )
      defGradientsTrans[i].resize( 3 * _fineTopology.getNumFaces() );
    
    // run over all faces
    for( int faceIdx = 0; faceIdx < _fineTopology.getNumFaces(); faceIdx++ ){
      
      // get nodal positions of triangle
      std::vector< VecType > p(3), q(3);
      for( int j = 0; j < 3; j++ ){
        getXYZCoord<VectorType, VecType>( source,    p[j], _fineTopology.getNodeOfTriangle(faceIdx,j) );
	getXYZCoord<VectorType, VecType>( defSource, q[j], _fineTopology.getNodeOfTriangle(faceIdx,j) );
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
    
  //
  void computeGradientSystemMatrix( const VectorType& geom, MatrixType& systemMat, MatrixType& diagMat ) const {
      
      // allocation
      if( (diagMat.rows() != 3 * _fineTopology.getNumFaces()) || (diagMat.cols() != 3 * _fineTopology.getNumFaces() ) )
        diagMat.resize( 3 * _fineTopology.getNumFaces(), 3 * _fineTopology.getNumFaces() );
      if( (systemMat.rows() != 3 * _fineTopology.getNumFaces()) || (systemMat.cols() != _fineTopology.getNumVertices() ) )
        systemMat.resize( 3 * _fineTopology.getNumFaces(), _fineTopology.getNumVertices() );
      
      // set up triplet lists
      TripletListType diagTriplets, systemTriplets;
      diagTriplets.reserve( 3 * _fineTopology.getNumFaces() );
      systemTriplets.reserve( 9 * _fineTopology.getNumFaces() );
	  
      //     |-1 1 0 |
      // C = |-1 0 1 |
      //     | 0 0 0 |
      MatType constMat(-1., 1., 0., -1., 0., 1., 0., 0., 0.);
      
      // run over all faces
      for( int faceIdx = 0; faceIdx < _fineTopology.getNumFaces(); faceIdx++ ){
	
	// compute base function gradient
	MatType baseFncGradient, temp;
	
	// get nodal positions pi, edges and unit normal n
	std::vector< VecType > nodalPos( 3 );
	for( int i = 0; i < 3; i++ ){
	  getXYZCoord<VectorType, VecType>( geom, nodalPos[i],  _fineTopology.getNodeOfTriangle( faceIdx, i ) );
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
	      systemTriplets.push_back( TripletType( 3 * faceIdx + i, _fineTopology.getNodeOfTriangle( faceIdx, j ), baseFncGradient.get(i,j) ) );
	
      }
      
      diagMat.setFromTriplets( diagTriplets.cbegin(), diagTriplets.cend() );
      systemMat.setFromTriplets( systemTriplets.cbegin(), systemTriplets.cend() );
  }

  //! perform detail transfer
  void performDetailTransfer( const VectorType& source, const VectorType& defSource, const VectorType& target, const VectorType& AnchorPositions, VectorType& defTarget ) const {
    
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

    // system mat A = G^T D G
    if( !_quiet ) std::cerr << "Assemble system matrix..." << std::endl;
    systemMat = tempMat * gradMat;
 
    // bc
    if( !_quiet ) std::cerr << "Manipulate system matrix to account for bc..." << std::endl;
    applyMaskToRow( _anchors, systemMat );
    // prescribe anchor positions    
    for( int j = 0; j < 3; j++ )
      for( int k = 0; k < _numOfVerticesInCoarseMesh; k++ )
        rhs[j][ _anchors[k] ] = AnchorPositions[ j*_numOfVerticesInCoarseMesh + k];

    // solving
    int dofs = _fineTopology.getNumVertices();
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
    //if( !_quiet ) std::cerr << "Done." << std::endl;
  }

  // Laplace-Beltrami as mean urvature vector, as e.g. in Meyer et al 2003, eq. (8) or Botsch and Kobbelt 2004, eq. (4) 
  // Laplace( p_i ) = 1/A_i \sum_j (cot a_{ij} + cot b_{ij}) (p_j - p_i), where we sum over the vertex 1-ring of p_i
  void approximateBiLaplaceMatrixByCotanFormula ( MatrixType &BiLaplaceMat ) const {
    
    int dimOfMatrix = _fineTopology.getNumVertices(); // + _anchors.size();
    MatrixType LaplaceMat( dimOfMatrix, dimOfMatrix );
    
    TripletListType tripletList;
    // number of nonzero entries = number of vertices * maximal valence
    tripletList.reserve( 10 * _fineTopology.getNumVertices() );
    
    VectorType geometry;
    getGeometry<VectorType>( _fineTopology.getGrid(), geometry );
					
    // assemble Laplace-Beltrami matrix
    for( int i = 0; i < _fineTopology.getNumVertices(); i++ ){
      // get vertex p_i
      VecType vertex;
      getXYZCoord<VectorType, VecType>( geometry, vertex, i );
      
      // get one ring of vertices
      std::vector<int> indices;
      _fineTopology.getVertex1RingVertices( i, indices );
      int valence = indices.size();
      
      // vertex positions p_j
      std::vector< VecType > vertices( valence );
      for( int j = 0; j < valence; j++ )
	getXYZCoord<VectorType, VecType>( geometry, vertices[j], indices[j] );
      
       // compute edges P_j - p_i
      std::vector< VecType > edges( vertices );
      for( int j = 0; j < valence; j++ )
	edges[j] -= vertex;
      
      // cotan arrays  cot a_{ij} + cot b_{ij}
      VectorType alpha( valence ), beta( valence );
      alpha.setZero();
      beta.setZero();
      for( int j = 0; j < valence; j++ ){
	// outer edge = p_{j+1} - p_j
	VecType outerEdge( vertices[(j+1)%valence] - vertices[j] );	
	alpha[(j+1)%valence] += getCotan( -1. * getAngle( outerEdge, edges[j] ) ) / 2.;
	beta[j] += getCotan( getAngle( outerEdge, edges[(j+1)%valence] ) ) / 2.;
      }
      
      // compute (mixed) Voronoi area (Meyer et al 2003, formula (7) and Fig. 4 for pseudo code)
      // A_i = 1/8 \sum_j (cot a_{ij} + cot b_{ij}) (p_j - p_i)^2, where we sum over the vertex 1-ring of p_i
      RealType mixedArea = 0.;
      for( int j = 0; j < valence; j++ ){	
	//vertices[j].makeCrossProduct( edges[j], edges[(j+1)%valence] );
	//mixedArea += vertices[j].norm() / 6.;
	mixedArea += (alpha[j] + beta[j]) * edges[j].normSqr() / 8.;
      }
     
      // assemble matrix     
      //!TODO check for factor 2 (cf. eq. (4) in Botsch and Kobbelt paper)
      tripletList.push_back( TripletType(i, i, -1. *  (alpha.sum() + beta.sum()) / mixedArea) );
      for( int j = 0; j < valence; j++ )
        tripletList.push_back( TripletType( i, indices[j], (alpha[j] + beta[j]) / mixedArea ) );
    }
    
    LaplaceMat.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
    
    // BiLaplace matrix is squared Laplace matrix
    BiLaplaceMat.resize(  dimOfMatrix, dimOfMatrix );
    BiLaplaceMat = LaplaceMat * LaplaceMat; // TODO improve efficiency here?
/*    
    // set constraints
    int offset = _fineTopology.getNumVertices();
    for( int i = 0; i < _anchors.size(); i++ ){
        BiLaplaceMat.coeffRef( _anchors[i], offset + i ) = 1.;
        BiLaplaceMat.coeffRef( offset + i, _anchors[i] ) = 1.;
    }
*/    
  }
 
  // interior angle in triangle is smaller than pi!
  RealType getAngle( const VecType& p, const VecType& q ) const { 
    RealType angle = std::acos( p*q / ( p.norm() * q.norm() ) );
    if( angle < 0 )
      throw BasicException ( "MultiResolutionOperator::getAngle(): negative angle!!" );
    return( angle );     
  }
  
  //
  RealType getCotan( RealType x) const { 
    return( 1. / std::tan(x) );     
  }
 
};

#endif