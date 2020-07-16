// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Example demonstrating how to check for isometric deformation modes using NRIC
 * \author Sassen
 *
 * In this example, we will check if, for a given mesh, there is an isometric mode in the tangent space of the NRIC
 * manifold, i.e. a tangent vector with zero length component. Then we will extrapolate it in a simple fashion.
 * Note that the criterion checked is only necessary not sufficient!
 * For more information, see Sassen, J., Heeren, B., Hildebrandt, K., & Rumpf, M. (2020). Geometric optimization using
 * nonlinear rotation-invariant coordinates. Computer Aided Geometric Design, 77, 101829.
 */

#include <goast/Core.h>
#include <goast/NRIC.h>
#include <goast/Optimization.h>
#include <goast/external/vtkIO.h>


int main( int argc, char *argv[] ) {
  // Import short handles for default types
  using RealType = DefaultConfigurator::RealType;
  using VectorType = DefaultConfigurator::VectorType;
  using MatrixType = DefaultConfigurator::SparseMatrixType;
  using FullMatrixType = DefaultConfigurator::FullMatrixType;

  // Used deformation energies
  using NRICDeformationType = ShellDeformation<DefaultConfigurator, NRICMembraneDeformation<DefaultConfigurator>, NRICBendingDeformation<DefaultConfigurator>>;

  // Additional type for frame
  using FrameType = DefaultConfigurator::FrameType;

  // Output format
  Eigen::IOFormat AsRowFmt( 8, 0, ", ", ",", "", "", "[ ", " ]" );

  std::cout << std::endl;
  std::cout << "                           Example: Isometric tangent modes                           " << std::endl;
  std::cout << "--------------------------------------------------------------------------------------" << std::endl;
  std::cout << std::endl;

  // ========================================================================
  // ============================== Parameters ==============================
  // ========================================================================

  // Files
  std::string File = "../../data/simple/SteffensPolyhedron.ply";

  std::string outputPrefix; // to specify potential output folder

  // Deformation model
  RealType bendingWeight = 0.001;

  // Reconstruction NRIC -> vertex position
  int gaussNewtonIterations = 100;

  // Tangent space processing
  RealType tangentThreshold = 1e-8; // Threshold for computing basis of tangent space using SVD
  RealType isometricThreshold = 1e-6; // Threshold for computing basis of isometric space using SVD


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
  int numBoundaryVertices = 0;
  for ( int vertexIdx = 0; vertexIdx < Topology.getNumVertices(); vertexIdx++ ) {
    for ( int edgeIdx = 0; edgeIdx < Topology.getNumEdges(); edgeIdx++ ) {
      if ( Topology.getAdjacentNodeOfEdge( edgeIdx, 0 ) == vertexIdx ||
           Topology.getAdjacentNodeOfEdge( edgeIdx, 1 ) == vertexIdx ) {
        int f0 = Topology.getOppositeNodeOfEdge( edgeIdx, 0 );
        int f1 = Topology.getOppositeNodeOfEdge( edgeIdx, 1 );

        if ( f0 == -1 || f1 == -1 )
          numBoundaryVertices++;
      }
    }
  }
  std::cout << " -- Number of boundary vertices: " << numBoundaryVertices << std::endl;

  int numBoundaryEdges = 0;
  std::vector<int> boundaryEdges;
  for ( int edgeIdx = 0; edgeIdx < numEdges; edgeIdx++ ) {
    int f0 = Topology.getOppositeNodeOfEdge( edgeIdx, 0 );
    int f1 = Topology.getOppositeNodeOfEdge( edgeIdx, 1 );

    if ( f0 == -1 || f1 == -1 ) {
      numBoundaryEdges++;
      boundaryEdges.push_back( edgeIdx );
    }
  }
  std::cout << " -- Number of boundary edges: " << numBoundaryEdges << std::endl;


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

  // Set up least-squares reconstruction
  VectorType reconstructionTarget( GeometryNRIC );
  VectorType reconstructionResult( Geometry );
  LinearReconstructionFunctional<DefaultConfigurator> L( Topology, reconstructionTarget, Weights, 1., bendingWeight );
  LinearReconstructionDerivative<DefaultConfigurator> DL( Topology, reconstructionTarget, Weights, 1., bendingWeight );
  GaussNewtonAlgorithm<DefaultConfigurator> GNOp( 2 * Topology.getNumEdges(), L, DL, gaussNewtonIterations );
//  GNOp.setBoundaryMask( bndMask );
  GNOp.setQuiet();

  std::cout << std::endl;
  std::cout << "=====================================================" << std::endl;
  std::cout << "  Tangent space" << std::endl;
  std::cout << "=====================================================" << std::endl;
  // Compute Jacobian
  MatrixType Jacobian = intGrad( GeometryNRIC );
  FullMatrixType denseJacobian = Jacobian;

  std::cout << " - Shape of Jacobian: " << denseJacobian.rows() << "x" << denseJacobian.cols() << std::endl;
  std::cout << " - Norm of Jacobian: " << denseJacobian.norm() << std::endl;

  FullMatrixType tangentBasis;
  int tangentDimension;

  std::cout << " - Computing SVD... " << std::flush;
  auto t_start = std::chrono::high_resolution_clock::now();
  Eigen::BDCSVD<FullMatrixType> svdSolver( denseJacobian, Eigen::ComputeFullV );
  auto t_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done in " << std::chrono::duration<double, std::ratio<1> >( t_end - t_start ).count()
            << "seconds." << std::endl;

  svdSolver.setThreshold( tangentThreshold );
  std::cout << " -- Rank of Jacobian: " << svdSolver.rank() << " with threshold " << svdSolver.threshold()
            << std::endl;
  std::cout << " -- Singular values: " << svdSolver.singularValues().format( AsRowFmt ) << std::endl;

  const FullMatrixType &V = svdSolver.matrixV();

  tangentDimension = V.cols() - svdSolver.rank();
  tangentBasis = svdSolver.matrixV().block( 0, svdSolver.rank(), V.rows(), tangentDimension ).sparseView();

  std::cout << std::endl;
  std::cout << "=====================================================" << std::endl;
  std::cout << "  Isometric modes" << std::endl;
  std::cout << "=====================================================" << std::endl;
  FullMatrixType dihedralBasis = FullMatrixType::Zero( 2 * numEdges, numEdges - numBoundaryEdges );

  // Basis of tangent space and isometric space in one matrix
  FullMatrixType combinedBases = FullMatrixType::Zero( 2 * numEdges,
                                                       tangentDimension + numEdges - numBoundaryEdges );
  combinedBases.block( 0, 0, 2 * numEdges, tangentDimension ) = tangentBasis;
  int interiorIdx = -1;
  for ( int edgeIdx = 0; edgeIdx < numEdges; edgeIdx++ ) {
    if ( std::find( boundaryEdges.begin(), boundaryEdges.end(), edgeIdx ) != boundaryEdges.end())
      continue;
    interiorIdx++;
    combinedBases( numEdges + edgeIdx, tangentDimension + interiorIdx ) = -1.;
    dihedralBasis( numEdges + edgeIdx, interiorIdx ) = -1.;
  }


  std::cout << " - Shape of combined bases: " << combinedBases.rows() << "x" << combinedBases.cols() << std::endl;
  std::cout << " - Norm of combined bases: " << combinedBases.norm() << std::endl;

  // Compute zero space of combined bases matrix
  std::cout << " - Computing SVD... " << std::flush;
  t_start = std::chrono::high_resolution_clock::now();
  Eigen::BDCSVD<FullMatrixType> svdSolver_isometric( combinedBases, Eigen::ComputeFullV );
  t_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done in " << std::chrono::duration<double, std::ratio<1> >( t_end - t_start ).count()
            << "seconds." << std::endl;
  svdSolver_isometric.setThreshold( isometricThreshold );
  std::cout << " -- Rank of combined bases: " << svdSolver_isometric.rank() << " with threshold "
            << svdSolver_isometric.threshold()
            << std::endl;

  std::cout << " -- Singular values: " << svdSolver_isometric.singularValues().format( AsRowFmt ) << std::endl;
  std::cout << " -- Smallest singular value: "
            << svdSolver_isometric.singularValues()[std::min( 2 * numEdges - numBoundaryEdges,
                                                              tangentDimension + numEdges - numBoundaryEdges ) - 1]
            << std::endl;


  // Compute basis of intersection
  const FullMatrixType &V_iso = svdSolver_isometric.matrixV();
  int intersectionDimension = V_iso.cols() - svdSolver_isometric.rank();

  FullMatrixType intersectionBasis =
          tangentBasis * V_iso.block( 0, svdSolver_isometric.rank(), tangentDimension, intersectionDimension );

  std::cout << " - Shape of intersection basis: " << intersectionBasis.rows() << "x" << intersectionBasis.cols()
            << std::endl;

  if ( intersectionDimension >= 1 ) {
    std::cout << std::endl;
    std::cout << "=====================================================" << std::endl;
    std::cout << "  Simple extrapolation" << std::endl;
    std::cout << "=====================================================" << std::endl;
    for ( int i = 0; i < intersectionDimension; i++ ) {
      std::cout << " - Mode " << i << std::endl;

      // Extract tangent mode and take care of possible boundary edges
      VectorType mode( 2 * numEdges );
      mode.setZero();
      interiorIdx = -1;
      for ( int edgeIdx = 0; edgeIdx < numEdges; edgeIdx++ ) {
        mode[edgeIdx] = intersectionBasis( edgeIdx, i );
        if ( std::find( boundaryEdges.begin(), boundaryEdges.end(), edgeIdx ) != boundaryEdges.end())
          continue;
        interiorIdx++;
        mode[numEdges + edgeIdx] = intersectionBasis( numEdges + interiorIdx, i );
      }

      for ( int t = 0; t <= 5; t++ ) {
        std::cout << " -- Timestep " << t << std::endl;
        // Reconstruct
        reconstructionTarget = GeometryNRIC + t / 5. * mode;
        GNOp.solve( reconstructionResult, reconstructionResult );

        // Create some analysis output
        std::map<std::string, VectorType> colorings;
        colorings["Isometry violation"] = VectorType::Zero( Topology.getNumFaces());

        for ( int faceIdx = 0; faceIdx < Topology.getNumFaces(); faceIdx++ ) {
          for ( int e : { 0, 1, 2 } ) {
            int edgeIdx = Topology.getEdgeOfTriangle( faceIdx, e );
            colorings["Isometry violation"][faceIdx] +=
                    std::abs( GeometryNRIC[edgeIdx] - Z( reconstructionResult )[edgeIdx] ) /
                    std::abs( GeometryNRIC[edgeIdx] );
          }
        }

#ifdef GOAST_WITH_VTK
        saveAsVTP( Topology, reconstructionResult,
                   outputPrefix + "mode_" + std::to_string( i ) + "_" + std::to_string( t ) + ".vtp",
                   colorings );
#else
        setGeometry( output, reconstructionResult );
        if ( !OpenMesh::IO::write_mesh( output, outputPrefix +
                                                "mode_" + std::to_string( i ) + "_" + std::to_string( t ) + ".ply" ))
          throw std::runtime_error( "Failed to write file: " + outputPrefix +
                                    "mode_" + std::to_string( i ) + "_" + std::to_string( t ) + ".ply" );
#endif
      }


    }

  }

}