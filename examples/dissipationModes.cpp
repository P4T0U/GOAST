// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include <iostream>
#include <chrono>
#include <ctime>
#include <string>

#include <goast/Core.h>
#include <goast/DiscreteShells.h>
#include <goast/GeodesicCalculus.h>

//==============================================================================================================
typedef ShellDeformation<DefaultConfigurator, NonlinearMembraneDeformation<DefaultConfigurator>, SimpleBendingDeformation<DefaultConfigurator> > ShellDeformationType;


//==============================================================================================================
//! \brief Shoot eigenmodes by means of discrete exponential map
//! \author Heeren
template< typename ConfiguratorType>
void shootEigenmodes( const MeshTopologySaver& Topology,
                      const typename ConfiguratorType::VectorType& ReferenceShape, 
                      const std::vector<typename ConfiguratorType::VectorType>& Eigenvectors, 
                      const std::vector<int>& Mask,
                      const DeformationBase<ConfiguratorType>& W,
                      const OptimizationParameters<ConfiguratorType>& optPars,
                      int numShootingSteps, 
                      double scaling,
                      std::string saveNameStem, 
                      int saveInterval = 1 ) {
         	
        int numEigenvalues = Eigenvectors.size();
        std::cerr << "Start to shoot "  << numEigenvalues << " eigenmodes " << numShootingSteps << " steps." << std::endl;        
        int numSavingShapes = numShootingSteps / saveInterval;
        
        
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for( int i = 0; i < numEigenvalues; i++ ){
          // save reference mesh
          TriMesh tempMesh( Topology.getGrid() );  
          setGeometry( tempMesh, ReferenceShape );
          std::ostringstream refMeshName;
          refMeshName << saveNameStem << "_mode" << i << "_" << numSavingShapes << ".ply";
          OpenMesh::IO::write_mesh( tempMesh, refMeshName.str() );  
          
          std::vector<typename ConfiguratorType::VectorType> Path;
          int iter = 0;
          typename ConfiguratorType::VectorType Variation( ReferenceShape );
          Variation += scaling * Eigenvectors[i];
          
          // shoot in positive direction          
	  integerExtrapolation( Topology, ReferenceShape, Variation, Mask, W, optPars, numShootingSteps, Path, true, true );
          
          // save           
          iter = 0;
          for( int j = 1; j <= numSavingShapes; j++ ){
              iter += saveInterval;
              setGeometry( tempMesh, Path[iter] );
              std::ostringstream saveMeshName;
              saveMeshName << saveNameStem << "_mode" << i << "_" << numSavingShapes + j << ".ply";
              OpenMesh::IO::write_mesh( tempMesh, saveMeshName.str() );
          }
          
          // shoot in negative direction  
	  integerExtrapolation( Topology, ReferenceShape, Variation, Mask, W, optPars, numShootingSteps + 1, Path, false, true );
          
          // save           
          iter = 1;
          for( int j = 1; j <= numSavingShapes; j++ ){
              iter += saveInterval;
              setGeometry( tempMesh, Path[iter] );
              std::ostringstream saveMeshName;
              saveMeshName << saveNameStem << "_mode" << i << "_" << numSavingShapes - j << ".ply";
              OpenMesh::IO::write_mesh( tempMesh, saveMeshName.str() );
          }          
          
	}
      
}

/**
 * \brief Dissipation modes as eigenvectors of the Hessian of an elastic shell energy
 * \author Heeren
 * 
 * First, compute pairs \f$ (\lambda_i, v_i) \f$ of eigenmodes such that \f$ Hv_i = \lambda_i Mv_i \f$, \f$ i =1, \ldots, n \f$,
 * and the orthonormality property \f$ v_j^T M v_i = \delta_{ij} \f$.
 * Here, \f$ H = W_{_22}[\bar x, \bar x]\f$ denotes the Hessian of the elastic shell energy W of a reference shape \$ \bar x \f$ 
 * and \f$ M = M[\bar x] \f$ is a mass matrix encoding the underlying scalar product.
 * One can either chose the \f$ L^2 \f$ scalar product or the Euclidean one.
 * In the latter case, \f$ M \f$ is simply the identity matrix.
 * 
 * Second, compute discrete exponential shooting of eigenmodes \f$ \pm v_i \f$ starting at \$ \bar x \f$.
 * In detail, compute \f$ x_{k+1} = exp_{x_k}( \pm \alpha * v_i ) \f$ for \$ k = 0,1, \ldots \$, 
 * where \f$ x_0 = \bar x \f$ and \f$ \alpha > 0 \$ is a scaling factor.
 * 
 * For a reference, see e.g. Figure 8 and 9 in \cite HeRuSc14
 * 
 */
int main(int argc, char *argv[])
{

try{

  std::cerr << "=================================================================================" << std::endl;
  std::cerr << "EXAMPLE DISSIPATION MODES" << std::endl;
  std::cerr << "=================================================================================" << std::endl << std::endl;
  
  //! 
  
  TriMesh mesh;
  OpenMesh::IO::read_mesh(mesh, "../../data/cactus/cactus0.ply"); 
   
  //! load reference shape 
  typename DefaultConfigurator::VectorType sourceGeom, Eigenvalues;
  getGeometry( mesh, sourceGeom );
  
  
  MeshTopologySaver Topol( mesh );
  typename DefaultConfigurator::RealType bendingWeight = 0.01;
  ShellDeformationType W(Topol, bendingWeight );
  ShellHessianOperator<DefaultConfigurator> hessian(  W );
  
  //! specify metric  
  EuclideanMetric<DefaultConfigurator> metric;
  L2Metric<DefaultConfigurator> l2Metric( Topol, sourceGeom );
  EigenmodesOp<DefaultConfigurator> eigenOp( Topol, hessian, metric, false );
  
  
  //! One can compute dissispation modes with or without fixing a subset of the vertices!
  std::vector<int> Mask;
  bool fixDirichletNodes = true;
  if( fixDirichletNodes ){
     //! boundary mask (fix 5-ring around bottom vertex)
    Topol.getVertexNRingVertices(5, 619, Mask); 
    std::cerr << "#bdry nodes = " << Mask.size() << std::endl; 
    extendBoundaryMask(Topol.getNumVertices(), Mask );  
    eigenOp.setBoundaryMask( Mask );
  }
  
  std::cerr << "COMPUTE DISSIPATION MODES....." << std::endl; 
  int numEigenvalues = 10;
  std::vector<typename DefaultConfigurator::VectorType> Eigenvectors;
  eigenOp.execute( sourceGeom, numEigenvalues, Eigenvalues, Eigenvectors );
  for( int i = 0; i < numEigenvalues; i++ )
      std::cerr << i << "th eigenvalue is " << Eigenvalues[i] << std::endl;
  std::cerr << std::endl; 
  
  std::cerr << "SHOOT EIGENVACTORS....." << std::endl; 
  //! shoot ten steps in both directions, 
  //! i.e. x_{k+1} = exp_{x_k}( \pm \alpha * v_i ) for k = 0,1,...,9,
  //! where v_i is the ith mode, x_0 the reference shape and \alpha > 0 a scaling factor
  int numShootingSteps = 10;
  double alpha = 1.; 
  // save each shooting step
  int saveInterval = 1;   
  OptimizationParameters<DefaultConfigurator> optPars;
  optPars.setNewtonIterations( 25 );  
  auto shoot_start = std::chrono::high_resolution_clock::now();
  shootEigenmodes<DefaultConfigurator>( Topol, sourceGeom, Eigenvectors, Mask, W, optPars, numShootingSteps, alpha, "cactus_0p01mu_l2Metric_DirichletBC", saveInterval );
  auto shoot_end = std::chrono::high_resolution_clock::now();
  std::cout << std::fixed << "Shooting done in " << std::chrono::duration<double, std::ratio<1> >(shoot_end - shoot_start).count() << "seconds." << std::endl;
 
  
  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}