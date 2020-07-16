// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Example demonstrating how to compute elastic means in NRIC
 * \author Sassen
 *
 * Based on Sassen, J., Heeren, B., Hildebrandt, K., & Rumpf, M. (2020). Geometric optimization using nonlinear
 * rotation-invariant coordinates. Computer Aided Geometric Design, 77, 101829.
 */

#include <goast/Core.h>
#include <goast/NRIC.h>
#include <goast/Optimization.h>
#include <goast/GeodesicCalculus.h>


int main( int argc, char *argv[] ) {
  // Import short handles for default types
  using RealType = DefaultConfigurator::RealType;
  using VectorType = DefaultConfigurator::VectorType;
  using VecType = DefaultConfigurator::VecType;

  // Used deformation energies
  using NRICDeformationType = ShellDeformation<DefaultConfigurator, NRICMembraneDeformation<DefaultConfigurator>, NRICBendingDeformation<DefaultConfigurator>>;

  // Additional type for frame
  using FrameType = DefaultConfigurator::FrameType;

  std::cout << std::endl;
  std::cout << "                           Example: Elastic Mean in NRIC                           " << std::endl;
  std::cout << "-----------------------------------------------------------------------------------" << std::endl;
  std::cout << std::endl;

  // ========================================================================
  // ============================== Parameters ==============================
  // ========================================================================

  // Files
  std::vector<std::string> Files{ "../../data/cactus/cactus0.ply",
                                  "../../data/cactus/cactus1.ply",
                                  "../../data/cactus/cactus2.ply",
                                  "../../data/cactus/cactus3.ply",
                                  "../../data/cactus/cactus4.ply",
                                  "../../data/cactus/cactus5.ply",
                                  "../../data/cactus/cactus6.ply",
                                  "../../data/cactus/cactus7.ply",
                                  "../../data/cactus/cactus8.ply",
                                  "../../data/cactus/cactus9.ply" };

  std::string outputPrefix; // to specify potential output folder

  // Deformation model
  RealType bendingWeight = 0.001;

  // Reconstruction NRIC -> vertex position
  int gaussNewtonIterations = 100;
  int initVertexIndex = 0; // Fixed initial vertex for frame-based reconstruction
  int initFrameIndex = -1; // Fixed initial triangle for frame-based reconstruction, -1 = automatic choice

  // Geodesic
  int lengthOfGeodesic = 15;

  // NRIC Optimization
  int constrainedIterations = 100; // Augmented-Lagrange iterations

  // ========================================================================
  // ================================= Setup ================================
  // ========================================================================
  // Read first mesh which determines the topology etc.
  TriMesh baseMesh;
  if ( !OpenMesh::IO::read_mesh( baseMesh, Files[0] ))
    throw std::runtime_error( "Failed to read file: " + Files[0] );

  MeshTopologySaver Topology( baseMesh );
  const int numEdges = Topology.getNumEdges();

  std::cout << std::endl << " - Mesh size:" << std::endl;
  std::cout << " -- Number of vertices: " << Topology.getNumVertices() << std::endl;
  std::cout << " -- Number of edges: " << Topology.getNumEdges() << std::endl;
  std::cout << " -- Number of faces: " << Topology.getNumFaces() << std::endl;

  // NRIC Map
  NRICMap<DefaultConfigurator> Z( Topology );

  // Read all meshes
//  std::vector<TriMesh> Meshes;
  std::vector<VectorType> Geometries;
  std::vector<VectorType> GeometriesNRIC;
  VectorType stackedGeometriesNRIC( Files.size() * 2 * numEdges );


  for ( int i = 0; i < Files.size(); i++ ) {
    const auto &File = Files[i];
    TriMesh Mesh;
    if ( !OpenMesh::IO::read_mesh( Mesh, File ))
      throw std::runtime_error( "Failed to read file: " + File );

    // Extract geometry and NRIC
    VectorType Geom, GeomNRIC;
    getGeometry( Mesh, Geom );
    GeomNRIC = Z( Geom );

    Geometries.push_back( Geom );
    GeometriesNRIC.push_back( GeomNRIC );

    stackedGeometriesNRIC.segment( i * 2 * numEdges, 2 * numEdges ) = GeomNRIC;
  }

  // Extract geometry
  VectorType baseGeometry;
  getGeometry( baseMesh, baseGeometry );
  VectorType baseGeometryNRIC = Z( baseGeometry );

  // Deformation energies
  NRICDeformationType W( Topology, bendingWeight );


  // Prepare variable for disk output
  TriMesh output( Topology.getGrid());

  // Get weights for the quadratic energies
  VectorType squaredWeights, Weights;
  getReconstructionWeights<DefaultConfigurator>( Topology, baseGeometry, squaredWeights, true );
  getReconstructionWeights<DefaultConfigurator>( Topology, baseGeometry, Weights, false );

  QuaternionIntegrabilityOp<DefaultConfigurator> intOp( Topology );
  QuaternionIntegrabilityGradient<DefaultConfigurator> intGrad( Topology );
  QuaternionIntegrabilityHessian<DefaultConfigurator> intHess( Topology );

  // Set up frame-based reconstruction
  // Seed for reconstruction
  VecType initVertexPosition;
  getXYZCoord<VectorType, VecType>( baseGeometry, initVertexPosition, initVertexIndex );

  // Automatic choice of triangle if none is specified
  if ( initFrameIndex == -1 ) {
    for ( int faceIdx = 0; faceIdx < Topology.getNumFaces(); faceIdx++ ) {
      if ( Topology.getNodeOfTriangle( faceIdx, 0 ) == initVertexIndex ||
           Topology.getNodeOfTriangle( faceIdx, 1 ) == initVertexIndex ||
           Topology.getNodeOfTriangle( faceIdx, 2 ) == initVertexIndex ) {
        initFrameIndex = faceIdx;
        break;
      }
    }
  }

  FrameType initFrame = extractFrame<RealType, VectorType, FrameType>( baseGeometry, Topology, initFrameIndex );

  DirectReconstruction<DefaultConfigurator> dirRec( Topology, initVertexIndex, initVertexPosition,
                                                    initFrameIndex, initFrame );

  // Set up least-squares reconstruction

  // Adopt dirichlet boundary values from initial triangle of frame-based reconstruction
  std::vector<int> bndMask{ Topology.getNodeOfTriangle( initFrameIndex, 0 ),
                            Topology.getNodeOfTriangle( initFrameIndex, 1 ),
                            Topology.getNodeOfTriangle( initFrameIndex, 2 ) };
  extendBoundaryMask( Topology.getNumVertices(), bndMask );

  VectorType reconstructionTarget( baseGeometryNRIC );
  VectorType reconstructionResult( baseGeometry );
  LinearReconstructionFunctional<DefaultConfigurator> L( Topology, reconstructionTarget, Weights, 1., bendingWeight );
  LinearReconstructionDerivative<DefaultConfigurator> DL( Topology, reconstructionTarget, Weights, 1., bendingWeight );
  GaussNewtonAlgorithm<DefaultConfigurator> GNOp( 2 * Topology.getNumEdges(), L, DL, gaussNewtonIterations );
//  GNOp.setBoundaryMask( bndMask );
  GNOp.setQuiet();


  std::cout << std::endl;
  std::cout << "=====================================================" << std::endl;
  std::cout << "  Linear mean in NRIC" << std::endl;
  std::cout << "=====================================================" << std::endl;
  VectorType linearMeanNRIC = VectorType::Zero( 2 * numEdges );

  for ( const auto &Geometry : GeometriesNRIC )
    linearMeanNRIC += Geometry;
  linearMeanNRIC /= GeometriesNRIC.size();

  std::cout << " - Reconstruction (Gauss-Newton): " << std::flush;
  reconstructionTarget = linearMeanNRIC;

  auto t_start = std::chrono::high_resolution_clock::now();
  GNOp.solve( reconstructionResult, reconstructionResult );
  auto t_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done." << std::endl;

  std::cout << std::fixed << " -- Computation time: "
            << std::chrono::duration<double, std::ratio<1> >( t_end - t_start ).count() << "seconds." << std::endl;

  std::cout << std::endl;


  std::cout << " - Writing to disk... ";
  setGeometry( output, reconstructionResult );
  if ( !OpenMesh::IO::write_mesh( output, outputPrefix + "linear_mean.ply" ))
    throw std::runtime_error( "Failed to write file: " + outputPrefix + "linear_mean.ply" );

  std::cout << "Done." << std::endl;


  std::cout << std::endl;
  std::cout << "=====================================================" << std::endl;
  std::cout << "  NRIC Elastic Mean" << std::endl;
  std::cout << "=====================================================" << std::endl;
  // Variables for results
  VectorType constrainedMeanNRIC = linearMeanNRIC;

  // Define functionals
  VectorType Alphas( GeometriesNRIC.size());
  Alphas.setConstant( 1. / GeometriesNRIC.size());

  // Abuse elastic mean functionals to get energy to a fixed mesh as simple objective
  ElasticMeanFunctional<DefaultConfigurator> EMF( W, stackedGeometriesNRIC, Alphas, GeometriesNRIC.size());
  ElasticMeanFunctionalGradient<DefaultConfigurator> EMG( W, stackedGeometriesNRIC, Alphas, GeometriesNRIC.size());
  ElasticMeanFunctionalHessian<DefaultConfigurator> EMH( W, stackedGeometriesNRIC, Alphas, GeometriesNRIC.size());

  // Start and end point as starting point
  VectorType startingPoint = linearMeanNRIC;//GeometriesNRIC[0];


  // Initialize solver
  AugmentedLagrangeMethod<DefaultConfigurator> constrainedSolverAL( EMF, EMG, EMH, intOp, intGrad, intHess,
                                                                    constrainedIterations, 1e-8, false );
  constrainedSolverAL.setParameter( "penalty_factor", 100. );
  constrainedSolverAL.setParameter( "constraint_decrease_exponent", 0.9 );
  constrainedSolverAL.setParameter( "constraint_tolerance", 5e-10 );
  constrainedSolverAL.setParameter( "optimality_tolerance", 5e-8 );
  constrainedSolverAL.setParameter( "inner_method", "LSN" );
  constrainedSolverAL.setParameter( "inner__tau_increase", 10. );
  constrainedSolverAL.setParameter( "inner__direction_beta", 1.e-3 );
  constrainedSolverAL.setParameter( "inner__minimal_stepsize", 1.e-15 );


  std::cout << " - Computing geodesic using augmented lagrange... ";
  t_start = std::chrono::high_resolution_clock::now();
  constrainedSolverAL.solve( startingPoint, constrainedMeanNRIC );
  t_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done in " << std::fixed << std::chrono::duration<double, std::ratio<1> >( t_end - t_start ).count()
            << "seconds."
            << std::endl;
  for ( const auto &timing : constrainedSolverAL.status.additionalIterations )
    std::cout << " -- Iterations (" << timing.first << "): " << timing.second << std::endl;
  for ( const auto &timing : constrainedSolverAL.status.additionalTimings )
    std::cout << " -- Time (" << timing.first << "): " << timing.second << "ms" << std::endl;

  RealType optTime = std::chrono::duration<RealType, std::ratio<1> >( t_end - t_start ).count();


  std::cout << std::endl;

  std::cout << " - Reconstruction (Frame-based): " << std::endl;
  t_start = std::chrono::high_resolution_clock::now();

  dirRec.apply( constrainedMeanNRIC, reconstructionResult );
//    GNOp.solve( reconstructionResult, reconstructionResult );
  t_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done." << std::endl;

  std::cout << std::fixed << " -- Computation time: "
            << std::chrono::duration<double, std::ratio<1> >( t_end - t_start ).count() << "seconds." << std::endl;

  std::cout << std::endl;

  std::cout << " - Writing to disk... ";
  setGeometry( output, reconstructionResult );
  if ( !OpenMesh::IO::write_mesh( output, outputPrefix + "constrained_mean.ply" ))
    throw std::runtime_error( "Failed to write file: " + outputPrefix + "constrained_mean.ply" );
}