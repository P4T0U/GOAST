// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Example demonstrating how to compute constriction of meshes, i.e. reduction of edge lengths along a curve, in NRIC
 * \author Sassen
 *
 * Based on Sassen, J., Heeren, B., Hildebrandt, K., & Rumpf, M. (2020). Geometric optimization using nonlinear
 * rotation-invariant coordinates. Computer Aided Geometric Design, 77, 101829.
 */

#include <goast/Core.h>
#include <goast/NRIC.h>
#include <goast/Optimization.h>
#include <goast/GeodesicCalculus/GeodesicInterpolation.h>
#include <goast/GeodesicCalculus/ElasticMean.h>


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
  std::cout << "                           Example: Constrictions using NRIC                           " << std::endl;
  std::cout << "---------------------------------------------------------------------------------------" << std::endl;
  std::cout << std::endl;

  // ========================================================================
  // ============================== Parameters ==============================
  // ========================================================================

  // Files
  std::string File = "../../data/plate/paperCrissCross.ply";

  std::string outputPrefix; // to specify potential output folder

  // Deformation model
  RealType bendingWeight = 1.;

  // Reconstruction NRIC -> vertex position
  int initVertexIndex = 17; // Fixed initial vertex for frame-based reconstruction
  int initFrameIndex = -1; // Fixed initial triangle for frame-based reconstruction, -1 = automatic choice

  // Folding
  RealType foldingAngle = M_PI / 2.;
  // List of edges to constrict, each edge is given by the two indices of its vertices
  // This format ist for example used by OpenFlipper
  std::vector<int> foldingList{ 559, 558, 560, 559, 558, 557, 557, 556, 1006, 1039, 1039, 1072, 973, 1006, 940,
                                973 };

  bool isometricInterpolation = true;


  // NRIC Optimization
  int constrainedIterations = 100; // Augmented-Lagrange iterations

  // ========================================================================
  // ================================= Setup ================================
  // ========================================================================
  TriMesh Mesh;
  if ( !OpenMesh::IO::read_mesh( Mesh, File ))
    throw std::runtime_error( "Failed to read file: " + File );


  std::cout << " - Read mesh: " << File << std::endl;

  MeshTopologySaver Topology( Mesh );
  const int numEdges = Topology.getNumEdges();
  const int numVertices = Topology.getNumVertices();
  const int numNodalDOFs = 3 * Topology.getNumVertices();

  std::cout << std::endl << " - Mesh size:" << std::endl;
  std::cout << " -- Number of vertices: " << Topology.getNumVertices() << std::endl;
  std::cout << " -- Number of edges: " << Topology.getNumEdges() << std::endl;
  std::cout << " -- Number of faces: " << Topology.getNumFaces() << std::endl;

  // Extract geometry
  VectorType Geometry;
  getGeometry( Mesh, Geometry );

  // Go to NRIC
  NRICMap<DefaultConfigurator> Z( Topology );
  VectorType GeometryNRIC = Z( Geometry );

  // Deformation energies
  NRICDeformationType W( Topology, bendingWeight );


  // Prepare variable for disk output
  TriMesh output( Topology.getGrid());

  // Get weights for the quadratic energies
  VectorType squaredWeights, Weights;
  getReconstructionWeights<DefaultConfigurator>( Topology, Geometry, squaredWeights, true );
  getReconstructionWeights<DefaultConfigurator>( Topology, Geometry, Weights, false );

  QuaternionIntegrabilityOp<DefaultConfigurator> intOp( Topology );
  QuaternionIntegrabilityGradient<DefaultConfigurator> intGrad( Topology );
  QuaternionIntegrabilityHessian<DefaultConfigurator> intHess( Topology );

  // Set up frame-based reconstruction
  // Seed for reconstruction
  VecType initVertexPosition;
  getXYZCoord<VectorType, VecType>( Geometry, initVertexPosition, initVertexIndex );

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

  FrameType initFrame = extractFrame<RealType, VectorType, FrameType>( Geometry, Topology, initFrameIndex );

  DirectReconstruction<DefaultConfigurator> dirRec( Topology, initVertexIndex, initVertexPosition,
                                                    initFrameIndex, initFrame );

  std::cout << std::endl;
  std::cout << "=====================================================" << std::endl;
  std::cout << "  Constriction" << std::endl;
  std::cout << "=====================================================" << std::endl;

  // Extract actual edges from vertex indices (a bit inefficient, but meh)
  std::vector<int> foldedEdges;
  for ( int i = 0; i < foldingList.size() / 2; i++ ) {
    int v_a = foldingList[2 * i];
    int v_b = foldingList[2 * i + 1];

    for ( int edgeIdx = 0; edgeIdx < numEdges; edgeIdx++ ) {
      int e1 = Topology.getAdjacentNodeOfEdge( edgeIdx, 0 );
      int e2 = Topology.getAdjacentNodeOfEdge( edgeIdx, 1 );

      if (( e1 == v_a && e2 == v_b ) || ( e1 == v_b && e2 == v_a ))
        foldedEdges.push_back( edgeIdx );
    }
  }

  std::cout << " - Number of folded edges: " << foldedEdges.size() << std::endl;

  VectorType constrictedGeometryNRIC = GeometryNRIC;
  for ( const int edgeIdx : foldedEdges )
    constrictedGeometryNRIC[numEdges + edgeIdx] = foldingAngle;

  // Define functionals
  VectorType Alphas( 1 );
  Alphas.setConstant( 1. );

  // Abuse elastic mean functionals to get energy to a fixed mesh as simple objective
  ElasticMeanFunctional<DefaultConfigurator> EMF( W, GeometryNRIC, Alphas, 1 );
  ElasticMeanFunctionalGradient<DefaultConfigurator> EMG( W, GeometryNRIC, Alphas, 1 );
  ElasticMeanFunctionalHessian<DefaultConfigurator> EMH( W, GeometryNRIC, Alphas, 1 );

  std::vector<int> bndMask;
  for ( const int edgeIdx : foldedEdges )
    bndMask.emplace_back( numEdges + edgeIdx );
  if ( isometricInterpolation )
    for ( int edgeIdx = 0; edgeIdx < numEdges; edgeIdx++ )
      bndMask.emplace_back( edgeIdx );


  AugmentedLagrangeMethod<DefaultConfigurator> constrainedSolverAL( EMF, EMG, EMH, intOp, intGrad, intHess,
                                                                    constrainedIterations, 1e-8, false );
  constrainedSolverAL.setParameter( "penalty_factor", 5. );
  constrainedSolverAL.setParameter( "constraint_decrease_exponent", 0.5 );
  constrainedSolverAL.setParameter( "constraint_tolerance", 1e-8 );
  constrainedSolverAL.setParameter( "optimality_tolerance", 5e-7 );
  constrainedSolverAL.setParameter( "functional_scaling", 0.1 );
  constrainedSolverAL.setParameter( "inner_method", "LSN" );
  constrainedSolverAL.setParameter( "inner__tau_increase", 10. );
  constrainedSolverAL.setParameter( "inner__direction_beta", 1.e-4 );
  constrainedSolverAL.setParameter( "inner__minimal_stepsize", 1.e-15 );
  constrainedSolverAL.setParameter( "inner__reduced_direction", 1 );

  constrainedSolverAL.setBoundaryMask( bndMask );

  std::cout << " - Computing constriction using augmented lagrange... " << std::endl;
  auto t_start = std::chrono::high_resolution_clock::now();
  constrainedSolverAL.solve( constrictedGeometryNRIC, constrictedGeometryNRIC );
  auto t_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done in " << std::fixed << std::chrono::duration<double, std::ratio<1> >( t_end - t_start ).count()
            << "seconds."
            << std::endl;

  for ( const auto &timing : constrainedSolverAL.status.additionalIterations )
    std::cout << " -- Iterations (" << timing.first << "): " << timing.second << std::endl;
  for ( const auto &timing : constrainedSolverAL.status.additionalTimings )
    std::cout << " -- Time (" << timing.first << "): " << timing.second << "ms" << std::endl;

  std::cout << " - Reconstruction (Frame-based): " << std::flush;
  VectorType constrictedGeometry = dirRec( constrictedGeometryNRIC );
  std::cout << "Done." << std::endl;

  std::cout << " - Writing to disk... ";
  setGeometry( output, constrictedGeometry );
  if ( !OpenMesh::IO::write_mesh( output, outputPrefix + "folded_mesh.ply" ))
    throw std::runtime_error( "Failed to write file: " + outputPrefix + "folded_mesh.ply" );
  std::cout << "Done." << std::endl;

}