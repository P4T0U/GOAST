// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Test of reconstruction operators for NRIC
 * \author Sassen
 */
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include "catch.hpp"

#include <goast/Core.h>
#include <goast/NRIC.h>

// ========================== Definitions ========================== //
typedef DefaultConfigurator ConfiguratorType;

// ========================== Tests ========================== //
TEST_CASE( "Frame-based reconstruction is inverse to NRIC map", "[NRIC]" ) {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::FrameType FrameType;

  const RealType tolerance = 1e-8;
  const std::string dataPath = "../../data/";

  std::vector<std::string> examples{ dataPath + "cactus/cactus0.ply", dataPath + "finger/finger0.ply",
                                     dataPath + "plate/plate121.ply" };

  for ( const auto &file : examples ) {
    SECTION( file ) {
      // Read mesh
      TriMesh mesh;
      if ( !OpenMesh::IO::read_mesh( mesh, file ))
        throw std::runtime_error( "Failed to read file: " + file );

      VectorType Geometry;
      getGeometry( mesh, Geometry );

      // Topology of the mesh
      MeshTopologySaver Topology( mesh );

      // Construct NRIC map and integrability operator
      NRICMap<ConfiguratorType> Z( Topology );

      VectorType GeometryNRIC = Z( Geometry );

      // Construct reconstruction operator
      int initVertexIndex = 0; // Just start with any vertex
      VecType initVertexPosition;
      getXYZCoord<VectorType, VecType>( Geometry, initVertexPosition, initVertexIndex );

      // Automatically determine a triangle matching the vertex
      int initFrameIndex = -1;
      for ( int faceIdx = 0; faceIdx < Topology.getNumFaces(); faceIdx++ ) {
        if ( Topology.getNodeOfTriangle( faceIdx, 0 ) == initVertexIndex ||
             Topology.getNodeOfTriangle( faceIdx, 1 ) == initVertexIndex ||
             Topology.getNodeOfTriangle( faceIdx, 2 ) == initVertexIndex ) {
          initFrameIndex = faceIdx;
          break;
        }
      }

      FrameType initFrame = extractFrame<RealType, VectorType, FrameType>( Geometry, Topology, initFrameIndex );

      DirectReconstruction<ConfiguratorType> R( Topology, initVertexIndex, initVertexPosition, initFrameIndex,
                                                initFrame );

      // Reconstructing yields the same
      SECTION ( "NRIC Map -> Reconstruction" ) {
        VectorType testGeometry = R( GeometryNRIC );

        REQUIRE(( testGeometry - Geometry ).norm() < tolerance );
      }

      SECTION ( "Reconstruction -> NRIC Map" ) {
        // Invert dihedral angles, preserves integrability while changing the geometry
        GeometryNRIC.tail( Topology.getNumEdges()) *= -1.;
        VectorType testGeometry = R( GeometryNRIC );
        VectorType testGeometryNRIC = Z( testGeometry );

        REQUIRE(( testGeometryNRIC - GeometryNRIC ).norm() < tolerance );
      }
    }
  }

}

TEST_CASE( "Least-squares reconstruction is inverse to NRIC map", "[NRIC]" ) {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::FrameType FrameType;

  const RealType bendingWeight = 0.001;
  const RealType tolerance = 1e-8;
  const std::string dataPath = "../../data/";

  std::vector<std::string> examples{ dataPath + "cactus/cactus0.ply", dataPath + "finger/finger0.ply",
                                     dataPath + "plate/plate121.ply" };

  for ( const auto &file : examples ) {
    SECTION( file ) {
      // Read mesh
      TriMesh mesh;
      if ( !OpenMesh::IO::read_mesh( mesh, file ))
        throw std::runtime_error( "Failed to read file: " + file );

      VectorType Geometry;
      getGeometry( mesh, Geometry );

      // Topology of the mesh
      MeshTopologySaver Topology( mesh );

      // Construct NRIC map and integrability operator
      NRICMap<ConfiguratorType> Z( Topology );
      VectorType GeometryNRIC = Z( Geometry );

      // Construct reconstruction operator
      // Fixed vertices
      int initFrameIndex = 0; // Just fix any triangle
      std::vector<int> bndMask{ Topology.getNodeOfTriangle( initFrameIndex, 0 ),
                                Topology.getNodeOfTriangle( initFrameIndex, 1 ),
                                Topology.getNodeOfTriangle( initFrameIndex, 2 ) };
      extendBoundaryMask( Topology.getNumVertices(), bndMask );

      // Get weights for the quadratic energies
      VectorType Weights;
      getReconstructionWeights<DefaultConfigurator>( Topology, Geometry, Weights, false );

      // Construct operator
      VectorType reconstructionTarget = GeometryNRIC;

      LinearReconstructionFunctional<DefaultConfigurator> L( Topology, reconstructionTarget, Weights, 1.,
                                                             bendingWeight );
      LinearReconstructionDerivative<DefaultConfigurator> DL( Topology, reconstructionTarget, Weights, 1.,
                                                              bendingWeight );
      GaussNewtonAlgorithm<DefaultConfigurator> GNOp( 2 * Topology.getNumEdges(), L, DL, 100 );
      GNOp.setBoundaryMask( bndMask );
      GNOp.setQuiet();


      // Reconstructing yields the same
      SECTION ( "NRIC Map -> Reconstruction" ) {
        // Initial point for Gauß-Newton, adding some random noise to the correct solution
        // Not too much, as otherwise local fold-overs etc might appear
        VectorType startGeometry = Geometry;
        startGeometry += VectorType::Random( 3 * Topology.getNumVertices()) *
                         startGeometry.array().abs().mean() * 0.005;
        for ( const int idx : bndMask )
          startGeometry[idx] = Geometry[idx];

        VectorType testGeometry;
        GNOp.solve( startGeometry, testGeometry );

        REQUIRE(( testGeometry - Geometry ).lpNorm<Eigen::Infinity>() < tolerance );
      }
    }
  }

}