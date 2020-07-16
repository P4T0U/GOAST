// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Example demonstrating how to compute various kind of interpolating geodesics in NRIC
 * \author Sassen
 *
 * Based on Sassen, J., Heeren, B., Hildebrandt, K., & Rumpf, M. (2020). Geometric optimization using nonlinear
 * rotation-invariant coordinates. Computer Aided Geometric Design, 77, 101829.
 */

#include <goast/Core.h>
#include <goast/NRIC.h>
#include <goast/Optimization.h>
#include <goast/GeodesicCalculus/GeodesicInterpolation.h>


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
  std::cout << "                           Example: Geodesics in NRIC                           " << std::endl;
  std::cout << "--------------------------------------------------------------------------------" << std::endl;
  std::cout << std::endl;

  // ========================================================================
  // ============================== Parameters ==============================
  // ========================================================================

  // Files
  std::string startFile = "../../data/plate/plate121_rolledA.ply";
  std::string endFile = "../../data/plate/plate121_rolledB.ply";

  std::string outputPrefix; // to specify potential output folder

  // Deformation model
  RealType bendingWeight = 0.001;

  // Reconstruction NRIC -> vertex position
  int gaussNewtonIterations = 100;
  int initVertexIndex = 60; // Fixed initial vertex for frame-based reconstruction
  int initFrameIndex = -1; // Fixed initial triangle for frame-based reconstruction, -1 = automatic choice

  // Geodesic
  int lengthOfGeodesic = 15;

  // NRIC Optimization
  int constrainedIterations = 100; // Augmented-Lagrange iterations

  // ========================================================================
  // ================================= Setup ================================
  // ========================================================================
  TriMesh startMesh, endMesh;
  if ( !OpenMesh::IO::read_mesh( startMesh, startFile ))
    throw std::runtime_error( "Failed to read file: " + startFile );
  if ( !OpenMesh::IO::read_mesh( endMesh, endFile ))
    throw std::runtime_error( "Failed to read file: " + endFile );

  std::cout << " - Read meshes:" << std::endl;
  std::cout << " -- Start: " << startFile << std::endl;
  std::cout << " -- End: " << endFile << std::endl;

  MeshTopologySaver Topology( startMesh );
  const int numEdges = Topology.getNumEdges();
  const int numVertices = Topology.getNumVertices();
  const int numNodalDOFs = 3 * Topology.getNumVertices();

  std::cout << std::endl << " - Mesh size:" << std::endl;
  std::cout << " -- Number of vertices: " << Topology.getNumVertices() << std::endl;
  std::cout << " -- Number of edges: " << Topology.getNumEdges() << std::endl;
  std::cout << " -- Number of faces: " << Topology.getNumFaces() << std::endl;

  // Extract geometry
  VectorType startGeom, endGeom;
  getGeometry( startMesh, startGeom );
  getGeometry( endMesh, endGeom );

  // Go to NRIC
  NRICMap<DefaultConfigurator> Z( Topology );
  VectorType startGeomNRIC = Z( startGeom ), endGeomNRIC = Z( endGeom );

  // Deformation energies
  NRICDeformationType W( Topology, bendingWeight );


  // Prepare variable for disk output
  TriMesh output( Topology.getGrid());

  // Get weights for the quadratic energies
  VectorType squaredWeights, Weights;
  getReconstructionWeights<DefaultConfigurator>( Topology, startGeom, squaredWeights, true );
  getReconstructionWeights<DefaultConfigurator>( Topology, startGeom, Weights, false );

  QuaternionIntegrabilityOp<DefaultConfigurator> intOp( Topology );
  QuaternionIntegrabilityGradient<DefaultConfigurator> intGrad( Topology );
  QuaternionIntegrabilityHessian<DefaultConfigurator> intHess( Topology );

  // Set up frame-based reconstruction
  // Seed for reconstruction
  VecType initVertexPosition;
  getXYZCoord<VectorType, VecType>( endGeom, initVertexPosition, initVertexIndex );

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

  FrameType initFrame = extractFrame<RealType, VectorType, FrameType>( endGeom, Topology, initFrameIndex );

  DirectReconstruction<DefaultConfigurator> dirRec( Topology, initVertexIndex, initVertexPosition,
                                                    initFrameIndex, initFrame );

  // Set up least-squares reconstruction

  // Adopt dirichlet boundary values from initial triangle of frame-based reconstruction
  std::vector<int> bndMask{ Topology.getNodeOfTriangle( initFrameIndex, 0 ),
                            Topology.getNodeOfTriangle( initFrameIndex, 1 ),
                            Topology.getNodeOfTriangle( initFrameIndex, 2 ) };
  extendBoundaryMask( Topology.getNumVertices(), bndMask );

  VectorType reconstructionTarget( startGeomNRIC );
  VectorType reconstructionResult( startGeom );
  LinearReconstructionFunctional<DefaultConfigurator> L( Topology, reconstructionTarget, Weights, 1., bendingWeight );
  LinearReconstructionDerivative<DefaultConfigurator> DL( Topology, reconstructionTarget, Weights, 1., bendingWeight );
  GaussNewtonAlgorithm<DefaultConfigurator> GNOp( 2 * Topology.getNumEdges(), L, DL, gaussNewtonIterations );
//  GNOp.setBoundaryMask( bndMask );
  GNOp.setQuiet();


  std::cout << std::endl;
  std::cout << "=====================================================" << std::endl;
  std::cout << "  Linear interpolation in NRIC" << std::endl;
  std::cout << "=====================================================" << std::endl;
  VectorType linearGeodesicLT(( lengthOfGeodesic - 2 ) * 2 * numEdges );
  VectorType reconstructedLinearGeodesic(( lengthOfGeodesic - 2 ) * numNodalDOFs );
  VectorType reconstructedLinearGeodesicNRIC(( lengthOfGeodesic - 2 ) * 2 * numEdges );

  // Linear interpolation
  for ( int i = 0; i < lengthOfGeodesic - 2; i++ ) {
    linearGeodesicLT.segment( i * 2 * numEdges, 2 * numEdges ) =
            startGeomNRIC + ( i + 1 ) * ( endGeomNRIC - startGeomNRIC ) / ( lengthOfGeodesic - 1 );
  }

  std::cout << " - Reconstruction (Gauss-Newton): " << std::endl;
  double recTime = 0;
  for ( int i = 0; i < lengthOfGeodesic - 2; i++ ) {
    reconstructionTarget = linearGeodesicLT.segment( i * 2 * numEdges, 2 * numEdges );

    std::cout << " -- Reconstructing step " << i + 1 << "... ";
    auto t_start = std::chrono::high_resolution_clock::now();
    GNOp.solve( reconstructionResult, reconstructionResult );
    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout << "Done." << std::endl;

    reconstructedLinearGeodesic.segment( i * numNodalDOFs, numNodalDOFs ) = reconstructionResult;
    reconstructedLinearGeodesicNRIC.segment( i * 2 * numEdges, 2 * numEdges ) = Z( reconstructionResult );

    recTime += std::chrono::duration<double, std::ratio<1> >( t_end - t_start ).count();
  }

  std::cout << std::fixed << " -- Computation time: " << recTime << "seconds." << std::endl;

  std::cout << std::endl;

  std::cout << " - Quantitative measures: " << std::endl;
  std::cout << " -- Distribution of energy along reconstruction geodesic: ";
  std::cout << W( startGeomNRIC, reconstructedLinearGeodesicNRIC.segment( 0, 2 * numEdges )) << " ";
  for ( int i = 1; i < lengthOfGeodesic - 2; i++ ) {
    std::cout << W( reconstructedLinearGeodesicNRIC.segment(( i - 1 ) * 2 * numEdges, 2 * numEdges ),
                    reconstructedLinearGeodesicNRIC.segment( i * 2 * numEdges, 2 * numEdges )) << " ";
  }
  std::cout << W( reconstructedLinearGeodesicNRIC.segment(( lengthOfGeodesic - 3 ) * 2 * numEdges, 2 * numEdges ),
                  endGeomNRIC ) << " ";

  std::cout << std::endl;

  std::cout << " - Writing to disk... ";
  setGeometry( output, startGeom );
  if ( !OpenMesh::IO::write_mesh( output, outputPrefix + "linear_geodesic_0.ply" ))
    throw std::runtime_error( "Failed to write file: " + outputPrefix + "linear_geodesic_0.ply" );


  for ( int i = 0; i < lengthOfGeodesic - 2; i++ ) {
    setGeometry( output, reconstructedLinearGeodesic.segment( i * numNodalDOFs, numNodalDOFs ));
    if ( !OpenMesh::IO::write_mesh( output, outputPrefix + "linear_geodesic_" + std::to_string( i + 1 ) + ".ply" ))
      throw std::runtime_error( "Failed to write file: " + outputPrefix + "linear_geodesic_" +
                                std::to_string( i + 1 ) + ".ply" );
  }

  setGeometry( output, endGeom );
  if ( !OpenMesh::IO::write_mesh( output, outputPrefix + "linear_geodesic_" + std::to_string( lengthOfGeodesic - 1 ) +
                                          ".ply" ))
    throw std::runtime_error( "Failed to write file: " + outputPrefix + "linear_geodesic_" +
                              std::to_string( lengthOfGeodesic - 1 ) + ".ply" );

  std::cout << "Done." << std::endl;


  std::cout << std::endl;
  std::cout << "=====================================================" << std::endl;
  std::cout << "  NRIC geodesic" << std::endl;
  std::cout << "=====================================================" << std::endl;
  // Variables for results
  VectorType constrainedGeodesicLT( linearGeodesicLT );
  VectorType reconstructedConstrainedGeodesic( reconstructedLinearGeodesic );
  VectorType reconstructedConstrainedGeodesicNRIC( reconstructedLinearGeodesicNRIC );

  // Define functionals
  DiscretePathEnergy<DefaultConfigurator> DPE( W, lengthOfGeodesic - 1, startGeomNRIC, endGeomNRIC );
  DiscretePathEnergyGradient<DefaultConfigurator> DPG( W, lengthOfGeodesic - 1, startGeomNRIC, endGeomNRIC );
  DiscretePathEnergyHessian<DefaultConfigurator> DPH( W, lengthOfGeodesic - 1, startGeomNRIC, endGeomNRIC );

  // Apply integrability operator for all shapes along the path
  PathConstraintOp<DefaultConfigurator> multipleIntegrabilityOp( intOp, ( lengthOfGeodesic - 2 ));
  PathConstraintGradient<DefaultConfigurator> multipleIntegrabilityGrad( intGrad, ( lengthOfGeodesic - 2 ));
  PathConstraintHessian<DefaultConfigurator> multipleIntegrabilityHess( intHess, ( lengthOfGeodesic - 2 ));

  // Start and end point as starting point
  VectorType startingPoint = linearGeodesicLT;
  for ( int i = 0; i < lengthOfGeodesic - 2; i++ ) {
    if ( i < (lengthOfGeodesic - 2) / 2 )
      startingPoint.segment( i * 2 * numEdges, 2 * numEdges ) = startGeomNRIC;
    else
      startingPoint.segment( i * 2 * numEdges, 2 * numEdges ) = endGeomNRIC;
  }


  // Initialize solver
  AugmentedLagrangeMethod<DefaultConfigurator> constrainedSolverAL( DPE, DPG, DPH,
                                                                    multipleIntegrabilityOp,
                                                                    multipleIntegrabilityGrad,
                                                                    multipleIntegrabilityHess,
                                                                    constrainedIterations,
                                                                    1e-8,
                                                                    false );
  constrainedSolverAL.setParameter("penalty_factor", 100.);
  constrainedSolverAL.setParameter("constraint_decrease_exponent", 0.9);
  constrainedSolverAL.setParameter("constraint_tolerance", 5e-10);
  constrainedSolverAL.setParameter("optimality_tolerance", 5e-7);
  constrainedSolverAL.setParameter("inner_method", "LSN");
  constrainedSolverAL.setParameter("inner__tau_increase", 10.);
  constrainedSolverAL.setParameter("inner__direction_beta", 1.e-3);
  constrainedSolverAL.setParameter("inner__minimal_stepsize", 1.e-15);


  std::cout << " - Computing geodesic using augmented lagrange... ";
  auto t_start = std::chrono::high_resolution_clock::now();
  constrainedSolverAL.solve( startingPoint, constrainedGeodesicLT );
  auto t_end = std::chrono::high_resolution_clock::now();
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
  reconstructionResult = startGeom;
  recTime = 0;
  for ( int i = 0; i < lengthOfGeodesic - 2; i++ ) {
    reconstructionTarget = constrainedGeodesicLT.segment( i * 2 * numEdges, 2 * numEdges );

    std::cout << " -- Reconstructing step " << i + 1 << "... ";
    t_start = std::chrono::high_resolution_clock::now();

    dirRec.apply( reconstructionTarget, reconstructionResult );
//    GNOp.solve( reconstructionResult, reconstructionResult );
    t_end = std::chrono::high_resolution_clock::now();
    std::cout << "Done." << std::endl;

    reconstructedConstrainedGeodesic.segment( i * numNodalDOFs, numNodalDOFs ) = reconstructionResult;

    reconstructedConstrainedGeodesicNRIC.segment( i * 2 * numEdges, 2 * numEdges ) = Z( reconstructionResult );

    recTime += std::chrono::duration<double, std::ratio<1> >( t_end - t_start ).count();
  }
  std::cout << std::fixed << " -- Computation time: " << recTime << "seconds." << std::endl;

  std::cout << std::endl;


  std::cout << " - Quantitative measures: " << std::endl;
  std::cout << " -- Distribution of energy along reconstructed geodesic: " << std::scientific;
  std::cout << W( startGeomNRIC, reconstructedConstrainedGeodesicNRIC.segment( 0, 2 * numEdges )) << " ";
  for ( int i = 1; i < lengthOfGeodesic - 2; i++ ) {
    std::cout << W( reconstructedConstrainedGeodesicNRIC.segment(( i - 1 ) * 2 * numEdges, 2 * numEdges ),
                    reconstructedConstrainedGeodesicNRIC.segment( i * 2 * numEdges, 2 * numEdges )) << " ";
  }
  std::cout << W( reconstructedConstrainedGeodesicNRIC.segment(( lengthOfGeodesic - 3 ) * 2 * numEdges, 2 * numEdges ),
                  endGeomNRIC ) << std::endl;

  std::cout << " - Writing to disk... ";
  setGeometry( output, startGeom );
  if ( !OpenMesh::IO::write_mesh( output, outputPrefix + "constrained_geodesic_0.ply" ))
    throw std::runtime_error( "Failed to write file: " + outputPrefix + "constrained_geodesic_0.ply" );


  for ( int i = 0; i < lengthOfGeodesic - 2; i++ ) {
    setGeometry( output, reconstructedConstrainedGeodesic.segment( i * numNodalDOFs, numNodalDOFs ));
    if ( !OpenMesh::IO::write_mesh( output, outputPrefix + "constrained_geodesic_" + std::to_string( i + 1 ) + ".ply" ))
      throw std::runtime_error( "Failed to write file: " + outputPrefix + "constrained_geodesic_" +
                                std::to_string( i + 1 ) + ".ply" );
  }

  setGeometry( output, endGeom );
  if ( !OpenMesh::IO::write_mesh( output, outputPrefix + "constrained_geodesic_" +
                                          std::to_string( lengthOfGeodesic - 1 ) + ".ply" ))
    throw std::runtime_error( "Failed to write file: " + outputPrefix + "constrained_geodesic_" +
                              std::to_string( lengthOfGeodesic - 1 ) + ".ply" );

  std::cout << "Done." << std::endl;

  std::cout << std::fixed << " - Computation time: " << optTime + recTime << "seconds." << std::endl;


  // Check if lengths of start and end shape agree. Then also do isometric interpolation.
  VectorType lengthDifference = startGeomNRIC.head( numEdges ) - endGeomNRIC.head( numEdges );
  if ( lengthDifference.lpNorm<Eigen::Infinity>() < 1e-6 ) {
    std::cout << std::endl;
    std::cout << "=====================================================" << std::endl;
    std::cout << "  Isometric geodesic" << std::endl;
    std::cout << "=====================================================" << std::endl;

    std::vector<int> bndMaskAL;
    for ( int i = 0; i < lengthOfGeodesic - 2; i++ ) {
      for ( int j = 0; j < numEdges; j++ ) {
        bndMaskAL.push_back( i * 2 * numEdges + j );
        startingPoint[i * 2 * numEdges + j] = startGeomNRIC[j];
      }
    }

    constrainedSolverAL.setParameter("inner__reduced_direction", 1);
    constrainedSolverAL.setParameter("relative_stopping", 0);
    constrainedSolverAL.setParameter("optimality_tolerance", 5e-8);
    constrainedSolverAL.setBoundaryMask( bndMaskAL );

    std::cout << " - Computing geodesic using augmented lagrange... ";
    t_start = std::chrono::high_resolution_clock::now();
    constrainedSolverAL.solve( startingPoint, constrainedGeodesicLT );
    t_end = std::chrono::high_resolution_clock::now();
    std::cout << "Done in " << std::fixed << std::chrono::duration<double, std::ratio<1> >( t_end - t_start ).count()
              << "seconds."
              << std::endl;

    std::cout << std::endl;

    std::cout << " - Reconstruction (Frame-based): " << std::endl;
    reconstructionResult = startGeom;
    recTime = 0;
    for ( int i = 0; i < lengthOfGeodesic - 2; i++ ) {
      reconstructionTarget = constrainedGeodesicLT.segment( i * 2 * numEdges, 2 * numEdges );

      std::cout << " -- Reconstructing step " << i + 1 << "... ";
      t_start = std::chrono::high_resolution_clock::now();

      dirRec.apply( reconstructionTarget, reconstructionResult );
//    GNOp.solve( reconstructionResult, reconstructionResult );
      t_end = std::chrono::high_resolution_clock::now();
      std::cout << "Done." << std::endl;

      reconstructedConstrainedGeodesic.segment( i * numNodalDOFs, numNodalDOFs ) = reconstructionResult;

      reconstructedConstrainedGeodesicNRIC.segment( i * 2 * numEdges, 2 * numEdges ) = Z( reconstructionResult );

      recTime += std::chrono::duration<double, std::ratio<1> >( t_end - t_start ).count();
    }
    std::cout << std::fixed << " -- Computation time: " << recTime << "seconds." << std::endl;

    std::cout << std::endl;


    std::cout << " - Quantitative measures: " << std::endl;
    std::cout << " -- Distribution of energy along reconstructed geodesic: " << std::scientific;
    std::cout << W( startGeomNRIC, reconstructedConstrainedGeodesicNRIC.segment( 0, 2 * numEdges )) << " ";
    for ( int i = 1; i < lengthOfGeodesic - 2; i++ ) {
      std::cout << W( reconstructedConstrainedGeodesicNRIC.segment(( i - 1 ) * 2 * numEdges, 2 * numEdges ),
                      reconstructedConstrainedGeodesicNRIC.segment( i * 2 * numEdges, 2 * numEdges )) << " ";
    }
    std::cout << W( reconstructedConstrainedGeodesicNRIC.segment(( lengthOfGeodesic - 3 ) * 2 * numEdges, 2 * numEdges ),
                    endGeomNRIC ) << std::endl;

    std::cout << " - Writing to disk... ";
    setGeometry( output, startGeom );
    if ( !OpenMesh::IO::write_mesh( output, outputPrefix + "isometric_geodesic_0.ply" ))
      throw std::runtime_error( "Failed to write file: " + outputPrefix + "isometric_geodesic_0.ply" );


    for ( int i = 0; i < lengthOfGeodesic - 2; i++ ) {
      setGeometry( output, reconstructedConstrainedGeodesic.segment( i * numNodalDOFs, numNodalDOFs ));
      if ( !OpenMesh::IO::write_mesh( output, outputPrefix + "isometric_geodesic_" + std::to_string( i + 1 ) + ".ply" ))
        throw std::runtime_error( "Failed to write file: " + outputPrefix + "isometric_geodesic_" +
                                  std::to_string( i + 1 ) + ".ply" );
    }

    setGeometry( output, endGeom );
    if ( !OpenMesh::IO::write_mesh( output, outputPrefix + "isometric_geodesic_" +
                                            std::to_string( lengthOfGeodesic - 1 ) + ".ply" ))
      throw std::runtime_error( "Failed to write file: " + outputPrefix + "isometric_geodesic_" +
                                std::to_string( lengthOfGeodesic - 1 ) + ".ply" );

    std::cout << "Done." << std::endl;
  }
}