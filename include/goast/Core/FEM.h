// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef FEM_HH
#define FEM_HH


//== INCLUDES =================================================================
#include "Auxiliary.h"
#include "LocalMeshGeometry.h"
#include "Topology.h"
#include "LinearSolver.h"

//=================================================================================
//! \brief Compute grid size statistics, i.e. minimal, maximal and average edge lengths
//! \author Heeren
template<typename ConfiguratorType>
void computeGridsize( const MeshTopologySaver& Topology, const typename ConfiguratorType::VectorType& Geometry, typename ConfiguratorType::VecType& gridSizeStats ){
    
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    
    // gridSizeStats = (min, max, average)
    gridSizeStats[0] = 1e+12;
    gridSizeStats[1] = 0.;
    gridSizeStats[2] = 0.;

    // compute minimal/maximal edge lengths
    for( int edgeIdx = 0; edgeIdx < Topology.getNumEdges(); edgeIdx++ ){
        
      int pi( Topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( Topology.getAdjacentNodeOfEdge(edgeIdx,1) );

      VecType Pi, Pj;
      getXYZCoord<VectorType, VecType>( Geometry, Pi, pi);
      getXYZCoord<VectorType, VecType>( Geometry, Pj, pj);
      
      RealType edgeLength = std::sqrt( dotProduct(Pi-Pj,Pi-Pj) );
      
      if( edgeLength < gridSizeStats[0] )
          gridSizeStats[0] = edgeLength;
      if( edgeLength > gridSizeStats[1] )
          gridSizeStats[1] = edgeLength;
      gridSizeStats[2] += edgeLength;
    }   
    
    gridSizeStats[2] *= 1. / Topology.getNumEdges();
}

//=================================================================================
/**
 * \brief Compute vertex weights
 * \author Heeren
 * For the $i$-th vertex the weight is defined as $w_i = 1/3 \sum_f a_f$,
 * where the sum is over all faces $f$ such that $i$ is part of $f$.
 */
template<typename ConfiguratorType>
void computeNodalAreas(const MeshTopologySaver& Topol, 
                       const typename ConfiguratorType::VectorType& Geometry, 
                       typename ConfiguratorType::VectorType& Areas ){
    
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::RealType RealType;
    
    Areas.resize( Topol.getNumVertices() );
    Areas.setZero();
    for ( int faceIdx = 0; faceIdx < Topol.getNumFaces(); ++faceIdx ){

      int pi( Topol.getNodeOfTriangle(faceIdx,0) ),
          pj( Topol.getNodeOfTriangle(faceIdx,1) ),
          pk( Topol.getNodeOfTriangle(faceIdx,2) );

      // set up deformed vertices and edges
      VecType Ei, Ej, Ek, temp;
      getXYZCoord<VectorType, VecType>( Geometry, temp, pi);
      getXYZCoord<VectorType, VecType>( Geometry, Ej, pj);
      getXYZCoord<VectorType, VecType>( Geometry, Ek, pk);
      Ei = Ek - Ej;
      Ej = temp - Ek;

      // compute volume
      temp.makeCrossProduct( Ei, Ej );
      RealType vol = std::sqrt( temp.normSqr() / 4. );       
      Areas[pi] += vol/3.;
      Areas[pj] += vol/3.;
      Areas[pk] += vol/3.;
    }
}


//=================================================================================
/**
 * \brief Compute vertex normals
 * \author Heeren
 * For the \f$i\f$-th vertex the normal is defined as $\f$n_i =  \sum_f n_f\f$,
 * where the sum is over all unit face normals \f$n_f\f$ of adjacent faces \f$f\f$ (such that \f$i\f$ is part of \f$f\f$).
 * If normalize = true, all vertex normals are normalized.
 */
template<typename ConfiguratorType>
void computeVertexNormals( const MeshTopologySaver& Topol, 
                           const typename ConfiguratorType::VectorType& Geometry, 
                           std::vector<typename ConfiguratorType::VecType>& vertexNormals,
                           bool normalize = true ){
    
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::RealType RealType;
    
    vertexNormals.clear();
    vertexNormals.resize( Topol.getNumVertices() );
    for( int i = 0; i < Topol.getNumVertices(); i++ )
        vertexNormals[i].setZero();

    for ( int faceIdx = 0; faceIdx < Topol.getNumFaces(); ++faceIdx ){
        
      // compute unit face normal and add to three vertices
      VecType Pi, Pj, Pk, faceNormal;  
      int pi( Topol.getNodeOfTriangle(faceIdx,0) ),
          pj( Topol.getNodeOfTriangle(faceIdx,1) ),
          pk( Topol.getNodeOfTriangle(faceIdx,2) );      
      getXYZCoord<VectorType, VecType>( Geometry, Pi, pi);
      getXYZCoord<VectorType, VecType>( Geometry, Pj, pj);
      getXYZCoord<VectorType, VecType>( Geometry, Pk, pk); 
      
      getNormal<RealType>( Pi, Pj, Pk, faceNormal );     
      vertexNormals[pi] += faceNormal; 
      vertexNormals[pj] += faceNormal; 
      vertexNormals[pk] += faceNormal;       
    }
    
    // normalize vertex normals
    if( normalize ){
      for( int i = 0; i < Topol.getNumVertices(); i++ )
        vertexNormals[i].normalize();
    }
}


//=================================================================================
/**
 * \brief Compute mass matrix for linear FEM
 * \author Heeren
 * Mass matrix \f$ M = (M_{ij})_{ij}\f$, with \f$ M_{ij} = \int_S \phi_i \phi_j \d a\f$,
 * where \f$\phi_i\f$ is the nodal basis function sitting at the i-th vertex, and the integration is over the mesh \f$ S \f$
 * Let \f$ n \f$ be the number of vertices in \f$ S \f$, then \f$ M \in \R^{n,n}\f$ resp. \f$ M \in \R^{3n,3n} \f$ if the matrix is supposed to be "extended".
 * In the latter case, the usual mass matrix is simply copied twice.
 */
template<typename ConfiguratorType>
void computeMassMatrix(const MeshTopologySaver& Topol, 
                       const typename ConfiguratorType::VectorType& Geometry, 
                       typename ConfiguratorType::SparseMatrixType& MassMatrix,
                       bool extended = true ){
    
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::TripletType TripletType;
    
    int extIdx = extended ? 3 : 1;
    MassMatrix.resize( extIdx * Topol.getNumVertices(), extIdx * Topol.getNumVertices() );  
    MassMatrix.setZero();
    std::vector<TripletType>  tripletList;
    tripletList.reserve( extIdx * 9 * Topol.getNumFaces() );   
        
    // run over all faces
    for ( int faceIdx = 0; faceIdx < Topol.getNumFaces(); ++faceIdx ){

      // get indices of vertices
      std::vector<int> verts(3);
      for( int j = 0; j < 3; j++ )
        verts[j] = Topol.getNodeOfTriangle(faceIdx,j);

      // get coordinates of vertices
      VecType Pi, Pj, Pk;
      getXYZCoord<VectorType, VecType>( Geometry, Pi, verts[0]);
      getXYZCoord<VectorType, VecType>( Geometry, Pj, verts[1]);
      getXYZCoord<VectorType, VecType>( Geometry, Pk, verts[2]);
      RealType area = getArea( Pi, Pj, Pk ) / 12.;
      
      // add local contribution
      for( int i = 0; i < 3; i++ )
          for( int j = 0; j < 3; j++ ){
              RealType factor = (i==j) ? 2. : 1.;
              for( int k = 0; k < extIdx; k++ )
                  tripletList.push_back( TripletType( k*Topol.getNumVertices() + verts[i], k*Topol.getNumVertices() + verts[j], factor * area) );
          }
    }
    
    MassMatrix.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
}


//=================================================================================
/**
 * \brief Compute lumped mass matrix for linear FEM
 * \author Heeren
 * Compute lumped mass matrix for mesh with n vertices.
 * The diagonal entries are given by the vertex weights.
 * If the matrix is supposed to be "extended", the output is a (3n)x(3n)-matrix.
 */
template<typename ConfiguratorType>
void computeLumpedMassMatrix( const MeshTopologySaver& Topol, 
                              const typename ConfiguratorType::VectorType& Geometry, 
                              typename ConfiguratorType::SparseMatrixType& LumpedMassMatrix,
                              bool extended = true ){
    
    int extIdx = extended ? 3 : 1;
    LumpedMassMatrix.resize( extIdx * Topol.getNumVertices(), extIdx * Topol.getNumVertices() );
    LumpedMassMatrix.setIdentity();
    
    typename ConfiguratorType::VectorType areas;
    computeNodalAreas<ConfiguratorType>( Topol, Geometry, areas );
    
    for( int i = 0; i < Topol.getNumVertices(); i++ )
      for( int j = 0; j < extIdx; j++ )
        LumpedMassMatrix.coeffRef( j* Topol.getNumVertices()+i, j* Topol.getNumVertices()+i ) = areas[i];
}

//=================================================================================
//! \brief Compute inverse lumped mass matrix for linear FEM
//! \author Heeren
template<typename ConfiguratorType>
void computeInverseLumpedMassMatrix( const MeshTopologySaver& Topol, 
                                     const typename ConfiguratorType::VectorType& Geometry, 
                                     typename ConfiguratorType::SparseMatrixType& InvLumpedMassMatrix,
                                     bool extended = true ){
    
    int extIdx = extended ? 3 : 1;
    InvLumpedMassMatrix.resize( extIdx * Topol.getNumVertices(), extIdx * Topol.getNumVertices() );
    InvLumpedMassMatrix.setIdentity();
    
    typename ConfiguratorType::VectorType areas;
    computeNodalAreas<ConfiguratorType>( Topol, Geometry, areas );
    
    for( int i = 0; i < Topol.getNumVertices(); i++ )
      for( int j = 0; j < extIdx; j++ )
        InvLumpedMassMatrix.coeffRef( j* Topol.getNumVertices()+i, j* Topol.getNumVertices()+i ) = 1. / areas[i];
}


//=================================================================================
/**
 * \brief Compute triplets for stiffness matrix for linear FEM
 * \author Heeren
 *
 * Stiffness matrix \f$ L = (L_{ij})_{ij}\f$, with \f$ L_{ij} = \int_S \nabla \phi_i \cdot \nabla \phi_j \d a\f$,
 * where \f$ \phi_i \f$ is the nodal basis function sitting at the i-th vertex, and the integration is over the mesh \f$ S \f$
 * Let \f$ n \f$ be the number of vertices in \f$ S \f$, then \f$ L \in \R^{n,n} \f$ resp. \f$ L \in \R^{3n,3n} \f$ if the matrix is supposed to be "extended".
 * In the latter case, the usual stiffness matrix is simply copied twice.
 *
 * Here, we make use of the formula provided e.g. in Corman et al. 2016, "Functional Characterization of Intrinsic and Extrinsic Geometry", cf. Fig. 1
 * This concides exactly with the cotan-formula (see class below), using the following simple geometric observations.
 * In a triangle \f$ t \f$ with edge lengths \f$(a,b,c)\f$, let \f$ \phi \f$ be the angle opposite edge \f$ a \f$ and \f$ vol(t) \f$ the area of \f$ t \f$.
 * Then \f$ vol(t) = 1/2 * \sin\phi * b * c\f$ and the cosine law \f$ a^2 = b^2 + c^2 - 2*b*c*\cos\phi \f$, which yield the equivalence.
 * HOWEVER, assembling the stiffness matrix this way is approximately twice as fast.
 */
template<typename ConfiguratorType>
void getStiffnessMatrixTriplets( const MeshTopologySaver& Topol, 
                                 const typename ConfiguratorType::VectorType& Geometry, 
                                 std::vector<typename ConfiguratorType::TripletType>& tripletListStiffnessMatrix,
                                 bool extended = true ){
    
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::TripletType TripletType;
    
    int extIdx = extended ? 3 : 1;
    tripletListStiffnessMatrix.clear();
    tripletListStiffnessMatrix.reserve( extIdx * 12 * Topol.getNumFaces() );   

    // run over all faces
    for ( int faceIdx = 0; faceIdx < Topol.getNumFaces(); ++faceIdx ){

      // get indices of vertices
      std::vector<int> verts(3);
      std::vector<VecType> nodes(3);
      for( int j = 0; j < 3; j++ ){
        verts[j] = Topol.getNodeOfTriangle(faceIdx,j);
        getXYZCoord<VectorType, VecType>( Geometry, nodes[j], verts[j]);
      }
      
      // get area and squared edge lengths
      RealType area = getArea( nodes[0], nodes[1], nodes[2] );      
      VecType lengthSqr;
      for( int j = 0; j < 3; j++ )
          lengthSqr[j] = dotProduct( nodes[(j+1)%3] -  nodes[(j+2)%3], nodes[(j+1)%3] -  nodes[(j+2)%3] );
      
      // add local contribution
      for( int j = 0; j < 3; j++ ){
          RealType entry =  0.125 * (lengthSqr[j] -  lengthSqr[(j+1)%3] -  lengthSqr[(j+2)%3]) / area;
          for( int k = 0; k < extIdx; k++ ){
              tripletListStiffnessMatrix.push_back( TripletType( k*Topol.getNumVertices() + verts[(j+1)%3], k*Topol.getNumVertices() + verts[(j+2)%3], entry) );
              tripletListStiffnessMatrix.push_back( TripletType( k*Topol.getNumVertices() + verts[(j+2)%3], k*Topol.getNumVertices() + verts[(j+1)%3], entry) );
              tripletListStiffnessMatrix.push_back( TripletType( k*Topol.getNumVertices() + verts[(j+1)%3], k*Topol.getNumVertices() + verts[(j+1)%3], -1. * entry) );
              tripletListStiffnessMatrix.push_back( TripletType( k*Topol.getNumVertices() + verts[(j+2)%3], k*Topol.getNumVertices() + verts[(j+2)%3], -1. * entry) );
          }
      }
    }
}

//=================================================================================
//! \brief Compute stiffness matrix for linear FEM (see documentation of triplet class above)
//! \author Heeren
template<typename ConfiguratorType>
void computeStiffnessMatrix( const MeshTopologySaver& Topol, 
                             const typename ConfiguratorType::VectorType& Geometry, 
                             typename ConfiguratorType::SparseMatrixType& StiffnessMatrix,
                             bool extended = true ){    
    int extIdx = extended ? 3 : 1;
    StiffnessMatrix.resize( extIdx * Topol.getNumVertices(), extIdx * Topol.getNumVertices() );    
    StiffnessMatrix.setZero();
    std::vector<typename ConfiguratorType::TripletType>  tripletList;
    getStiffnessMatrixTriplets<ConfiguratorType>( Topol, Geometry, tripletList, extended );   
    StiffnessMatrix.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
}

//=================================================================================
//! \brief Compute stiffness matrix for linear FEM (with fixed zeroth momentum to make it invertible)
//! \author Heeren
template<typename ConfiguratorType>
void computeStiffnessMatrixWithConstraint( const MeshTopologySaver& Topol, 
                                           const typename ConfiguratorType::VectorType& Geometry, 
                                           typename ConfiguratorType::SparseMatrixType& regularStiffnessMatrix ){    
    
    typedef typename ConfiguratorType::TripletType TripletType;
    typedef typename ConfiguratorType::VectorType VectorType;
    
    regularStiffnessMatrix.resize( Topol.getNumVertices() + 1, Topol.getNumVertices() + 1);    
    regularStiffnessMatrix.setZero();
    
    // first add triplets for stiffness matrix
    std::vector<TripletType>  tripletList;
    getStiffnessMatrixTriplets<ConfiguratorType>( Topol, Geometry, tripletList, false );   
    
    // now add constraint entries, i.e. L_{n,i} = L_{i,n} = m_i for i=0,...,n-1 where m_i is mass of ith vertex.
    VectorType areas;
    computeNodalAreas<ConfiguratorType>( Topol, Geometry, areas );
    for( int k = 0; k < areas.size(); k++ ){
        tripletList.push_back( TripletType( Topol.getNumVertices(), k, areas[k] ) );
        tripletList.push_back( TripletType( k, Topol.getNumVertices(), areas[k] ) );
    }
    
    regularStiffnessMatrix.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
} 
 
//=================================================================================
/**
 * \brief Compute stiffness matrix for linear FEM via cotan formula (cf. Pinkall and Polthier, 1993)
 * \author Heeren
 * Stiffness matrix \f$ L = (L_{ij})_{ij}\f$ using the cotan-formula - coincides with result of computeStiffnessMatrix() above.
 * Here \f$L_{ij} = -0.5 * (\cot \alpha_{ij} + \cot \beta_{ij})\f$, for \f$ i \neq j \f$, and \f$ L_{ii} = - \sum_{j} L_{ij} \f$
 * If \f$ e = (i,j) \f$ is an edge and \f$ t_k = (i,j,k)\f$ and \f$ t_l = (j,i,l)\f$ two adjacent triangles sharing edge \f$ e \f$,
 * Then \f$\alpha_{ij}\f$ is the interior triangle angle in \f$ t_k\f$ at vertex \f$ k \f$,
 * and \f$\beta_{ij}\f$ the interior triangle angle in \f$t_l\f$ at vertex \f$ l \f$.
 */
template<typename ConfiguratorType>
void computeStiffnessMatrixCotan( const MeshTopologySaver& Topology, 
                                  const typename ConfiguratorType::VectorType& Geometry, 
                                  typename ConfiguratorType::SparseMatrixType& StiffnessMatrix,
                                  bool extended = true ){
    
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::VecType VecType;
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::TripletType TripletType;
    
    int extIdx = extended ? 3 : 1;
    int numV = Topology.getNumVertices();
    StiffnessMatrix.resize( extIdx * numV, extIdx * numV );    
    StiffnessMatrix.setZero();
    std::vector<TripletType>  tripletList;
    tripletList.reserve( extIdx * 4 * Topology.getNumEdges() );  

    // assemble Laplace-Beltrami matrix
    for( int edgeIdx = 0; edgeIdx < Topology.getNumEdges(); edgeIdx++ ){
        
      int pi( Topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( Topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
          pk( Topology.getOppositeNodeOfEdge(edgeIdx,0) ),
          pl( Topology.getOppositeNodeOfEdge(edgeIdx,1) );

      // get vertices
      VecType Pi, Pj, Pk, Pl;
      getXYZCoord<VectorType, VecType>( Geometry, Pi, pi);
      getXYZCoord<VectorType, VecType>( Geometry, Pj, pj);
      getXYZCoord<VectorType, VecType>( Geometry, Pk, pk);
      getXYZCoord<VectorType, VecType>( Geometry, Pl, pl);

      // cotan(x) = 1 / tan(x)
      RealType cotAlphaij = 1. / std::tan( getInnerAngle<RealType>( Pk - Pi, Pk - Pj ) );
      RealType cotBetaij  = 1. / std::tan( getInnerAngle<RealType>( Pl - Pi, Pl - Pj ) );
     
      // assemble matrix     
      for( int k = 0; k < extIdx; k++ ){
        tripletList.push_back( TripletType( k*numV + pi,  k*numV + pj, -0.5 * (cotAlphaij + cotBetaij) ) );
        tripletList.push_back( TripletType( k*numV + pj,  k*numV + pi, -0.5 * (cotAlphaij + cotBetaij) ) );
        tripletList.push_back( TripletType( k*numV + pi,  k*numV + pi, 0.5 * (cotAlphaij + cotBetaij) ) );
        tripletList.push_back( TripletType( k*numV + pj,  k*numV + pj, 0.5 * (cotAlphaij + cotBetaij) ) );
      }
    }    
    StiffnessMatrix.setFromTriplets( tripletList.cbegin(), tripletList.cend() ); 
}

//=================================================================================
//! \brief Compute lowest eigenvalues and eigenfunctions of Laplace-Beltrami
//! \author Heeren
//! \todo Check if we need to do normalization wrt. mass matrix!
template<typename ConfiguratorType>
void computeLaplaceBeltramiEigenfunctions( const MeshTopologySaver& Topol, 
                                           const typename ConfiguratorType::VectorType& Geometry, 
                                           int numEigenvalues,
                                           std::vector<typename ConfiguratorType::RealType>& Eigenvalues,
                                           std::vector<typename ConfiguratorType::VectorType>& Eigenfunctions,
                                           int maxIterations = 1000,
                                           bool useL2Metric = true,
                                           bool quiet = true ){
    
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::RealType RealType;  
  
  // compute pairs (\lamba_i, u_i) with Lu_i = \lambda_i Mu_i
  typename ConfiguratorType::SparseMatrixType MassMatrix, extStiffnessMatrix;
  computeMassMatrix<ConfiguratorType>( Topol, Geometry, MassMatrix, false );
  computeStiffnessMatrixWithConstraint<ConfiguratorType>( Topol, Geometry, extStiffnessMatrix );  
  int dim = Topol.getNumVertices();
  
  // initialize solver   
  if(!quiet) std::cerr << "Factorize matrix..." << std::endl;
  LinearSolver<ConfiguratorType> linSolver;
  auto t_start = std::chrono::high_resolution_clock::now(); 
  linSolver.prepareSolver( extStiffnessMatrix );
  auto t_end = std::chrono::high_resolution_clock::now();
  if(!quiet) std::cerr << std::fixed << "done in " << std::chrono::duration<double, std::ratio<1> >(t_end - t_start).count() << " seconds." << std::endl;  
  
  // clear variables
  Eigenvalues.clear();
  Eigenvalues.resize( numEigenvalues );
  Eigenfunctions.clear();
  Eigenfunctions.resize( numEigenvalues );

  // compute eigenvalues via inverse power iteration
  for ( int i = 0; i < numEigenvalues; ++i ){
	  
	    if(!quiet)
	        std::cerr << "Start to compute " << i+1 << "th of " << numEigenvalues << " eigenvectors..." << std::endl;

	    if( Eigenfunctions[i].size() == 0 ){
          Eigenfunctions[i] = VectorType::Constant( dim, 1. );
          VectorType temp = MassMatrix * Eigenfunctions[i];
          RealType sqrNorm = useL2Metric ? Eigenfunctions[i].dot( temp ) : Eigenfunctions[i].squaredNorm();
	      Eigenfunctions[i] /= std::sqrt(sqrNorm);
	    }

        RealType err = 1.;
	    int iter = 0;
        Eigenvalues[i] = (i>0) ? Eigenvalues[i-1] : 0;
        VectorType rhs;
           
	    // start inverse vector iteration
        for (; iter < maxIterations && err > 1e-10; ++iter){
            		
                // rhs = [Mx^k; 0], where M is mass matrix
                VectorType rhs =  MassMatrix * Eigenfunctions[i];
                rhs.conservativeResize( dim + 1 );
                rhs[dim] = 0.;

                // solve L x^{k+1} = Mx^k
                VectorType newEV;
                linSolver.backSubstitute( rhs, newEV );
                newEV.conservativeResize( dim );
                rhs.conservativeResize( dim );

                // projection
                VectorType factors(i);
                for (int j = 0; j < i; ++j)
                    factors[j] = useL2Metric ? newEV.dot( MassMatrix * Eigenfunctions[j] ) : newEV.dot( Eigenfunctions[j] );
                for (int j = 0; j < i; ++j)
                  newEV -= factors[j] * Eigenfunctions[j];
                
                // approximation by Rayleigh coefficients         
                RealType newEigenValueInverse = useL2Metric ? newEV.dot( rhs ) : newEV.dot( Eigenfunctions[i] );
                RealType newEigenValue = 1. / newEigenValueInverse;                 
                                
                // normalize
                RealType sqrNorm = useL2Metric ? newEV.dot(MassMatrix * newEV) : newEV.squaredNorm();
		        newEV /= std::sqrt( sqrNorm );
		
		        // update error
                err = std::abs( (newEigenValue - Eigenvalues[i]) / Eigenvalues[i] );
                
                // update
                Eigenvalues[i] = newEigenValue;
                Eigenfunctions[i] = newEV;
        }
        
        // final console output
        if(!quiet) std::cerr << iter << " iterations: lambda[" << i << "] = " << Eigenvalues[i] << ", error = " << err << std::endl;
	    if(!quiet) std::cerr << "-----------------------------------------------------------" << std::endl;
  }
 
}


//=================================================================================
/**
 * \brief Solve heat equation with FEM and implicit Euler time-stepping
 * \author Heeren
 * For initial value \f$ u_0 = u \f$ solve iteratively \f$(Mu^k - Mu^{k-1})/ tau = -Lu^k \f$ for  \f$ u^k \f$ and \f$k > 0\f$,
 * where \f$ M \f$ and \f$ L \f$ denote the mass and stiffness matrices, respectively.
 * Dirichlet boundary data is encoded via an integer mask.
 */
template<typename ConfiguratorType>
void solveHeatFlowImplicit( const MeshTopologySaver& Topol, 
                            const typename ConfiguratorType::VectorType& Geometry, 
                            typename ConfiguratorType::RealType tau,
                            int steps,
                            std::vector<int>& DirichletMask,
                            typename ConfiguratorType::VectorType& u ){
    
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::RealType RealType;  
    
  typename ConfiguratorType::SparseMatrixType MassMatrix, StiffnessMatrix;
  computeMassMatrix<ConfiguratorType>( Topol, Geometry, MassMatrix, false );
  computeStiffnessMatrixCotan<ConfiguratorType>( Topol, Geometry, StiffnessMatrix, false );
  
  // set up solver to compute M^{-1}
  LinearSolver<ConfiguratorType> solver;
  typename ConfiguratorType::SparseMatrixType SystemMatrix = MassMatrix + tau * StiffnessMatrix;
  applyMaskToRow( DirichletMask, SystemMatrix );
  solver.prepareSolver( SystemMatrix );
  applyMaskToRow( DirichletMask, MassMatrix );
  
  // implicit time-stepping, i.e. solve (M + tau*L) u^k = Mu^{k-1} for u^k
  for( int i = 0; i < steps; i++ ){
    typename ConfiguratorType::VectorType rhs = MassMatrix * u;
    solver.backSubstitute( rhs, u ); 
  }
}

//=================================================================================
//! \brief Solve mean curvature flow with FEM and implicit Euler time-stepping
//! \author Heeren
template<typename ConfiguratorType>
void solveMeanCurvatureFlowImplicit( const MeshTopologySaver& Topol, 
                                     typename ConfiguratorType::VectorType& Geometry, 
                                     typename ConfiguratorType::RealType tau,
                                     int steps,
                                     std::vector<int>& DirichletMask,
                                     bool quiet = true ){    
  typedef typename ConfiguratorType::VectorType VectorType;    
  
  std::vector<bool> mask( Topol.getNumVertices() );
  for( int i = 0; i < Topol.getNumVertices(); i++ )
      mask[i] = false;
  for( int i = 0; i < DirichletMask.size(); i++ )
      if( DirichletMask[i] < Topol.getNumVertices() )
          mask[DirichletMask[i]] = true;
 
  // implicit time-stepping, i.e. solve (M + tau*L) u^k = Mu^{k-1} for u^k
  for( int i = 0; i < steps; i++ ){
      
      if(!quiet) std::cerr << "Start " << i+1 << "th MCM iteration of " << steps << std::endl;
      
      // stiffness matrix
      typename ConfiguratorType::SparseMatrixType MassMatrix, StiffnessMatrix;
      // full mass matrix
      //computeMassMatrix<ConfiguratorType>( Topol, Geometry, MassMatrix, false );
      // lumped mass matrix
      computeLumpedMassMatrix<ConfiguratorType>( Topol, Geometry, MassMatrix, false );
      computeStiffnessMatrixCotan<ConfiguratorType>( Topol, Geometry, StiffnessMatrix, false );
  
      // set up solver to compute (M + tau*L)^{-1}
      LinearSolver<ConfiguratorType> solver;
      typename ConfiguratorType::SparseMatrixType SystemMatrix = MassMatrix + tau * StiffnessMatrix;
      //applyMaskToSymmetricMatrix( DirichletMask, SystemMatrix );
      solver.prepareSolver( SystemMatrix );
      //applyMaskToSymmetricMatrix( DirichletMask, MassMatrix );
  
      // apply to (x,y,z)-components 
      VectorType oldGeom( Geometry );
      for(int j = 0; j < 3; j++){
        // get geometry component
        VectorType geomComp( Topol.getNumVertices() );
        for( int k = 0; k < Topol.getNumVertices(); k++ )
            geomComp[k] = Geometry[ j*Topol.getNumVertices() + k ];
        VectorType rhs = MassMatrix * geomComp;
        solver.backSubstitute( rhs, geomComp ); 
        
        // write back component
        for( int k = 0; k < Topol.getNumVertices(); k++ )
          if( !mask[k] )
              Geometry[j*Topol.getNumVertices() + k] = geomComp[k];
      }
  }
}

/*
//! Import from old QuocMesh code
template<typename ConfiguratorType>
class MeanCurvatureSmoother{
  
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  
  const MeshTopologySaver& _topology;

public:
  MeanCurvatureSmoother( const MeshTopologySaver& Topology ) : _topology( Topology ) {}
  
  void execute( VectorType& Geometry, RealType tau, int numSteps ) const {

    const int size = _topology.getNumVertices();
    std::vector<int> mask;
    _topology.fillFullBoundaryMask( mask );
    std::vector<bool> boolMask(size);
    for( int i = 0; i < size; i++ )
        boolMask[i] = false;
    for( int i = 0; i < mask.size(); i++ )
        boolMask[mask[i]] = true;

    for ( int i = 0; i < numSteps; ++i )
      mcmStep ( Geometry, tau, boolMask );
  }

protected:
void mcmStep ( VectorType& Geometry, RealType tau, const std::vector<bool>& mask ) const {

    std::cerr << "bla" << std::endl;
  const int size = _topology.getNumVertices();

  typename ConfiguratorType::SparseMatrixType MassMatrix, StiffnessMatrix;
  computeMassMatrix<ConfiguratorType>( _topology, Geometry, MassMatrix, false );
  computeStiffnessMatrixCotan<ConfiguratorType>( _topology, Geometry, StiffnessMatrix, false );
  
  // set up solver to compute (M + tau*L)^{-1}
  LinearSolver<ConfiguratorType> solver;
  typename ConfiguratorType::SparseMatrixType SystemMatrix = MassMatrix + tau * StiffnessMatrix;
  solver.prepareSolver( SystemMatrix );

  std::cerr << "bla" << std::endl;
  for( int j = 0; j < 3; j++ ){
    VectorType surf(size);
    for( int i = 0; i < size; i++ )
        surf[i] = Geometry[j*size+i];
    VectorType surf_old = surf;
    VectorType rhs = MassMatrix * surf;
    
    solver.backSubstitute ( rhs, surf_old );
    
    for( int i = 0; i < size; i++ )
      if ( !mask[i] )
        Geometry[j*size+i] = surf_old[i];
  }

}


};
*/

#endif