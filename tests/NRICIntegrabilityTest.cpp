// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Test of integrability operators for NRIC
 * \author Sassen
 */
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include "catch.hpp"

#include <goast/Core.h>
#include <goast/NRIC.h>

// ========================== Definitions ========================== //
typedef DefaultConfigurator ConfiguratorType;

// ========================== Tests ========================== //
TEST_CASE( "NRIC obtained from meshes are integrable", "[NRIC]" ) {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;

  const RealType tolerance = 1e-8;
  const std::string dataPath = "../../data/";

  std::vector<std::string> examples{ dataPath + "cactus/cactus0.ply", dataPath + "finger/finger0.ply",
                                     dataPath + "plate/plate121.ply" };

  SECTION( "AngleIntegrability" ) {
    for ( const auto &file : examples ) {
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
      AngleIntegrabilityOp<ConfiguratorType> I( Topology );

      // Compute residual of integrability
      VectorType error = I( Z( Geometry )).array() - 3.;

      REQUIRE( error.lpNorm<Eigen::Infinity>() < tolerance );
    }
  }

  SECTION( "EulerIntegrability" ) {
    for ( const auto &file : examples ) {
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
      EulerIntegrabilityOp<ConfiguratorType> I( Topology );

      // Compute residual of integrability
      VectorType error = I( Z( Geometry ));

      REQUIRE( error.lpNorm<Eigen::Infinity>() < tolerance );
    }
  }

  SECTION( "QuaternionIntegrability" ) {
    for ( const auto &file : examples ) {
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
      QuaternionIntegrabilityOp<ConfiguratorType> I( Topology );

      // Compute residual of integrability
      VectorType error = I( Z( Geometry ));

      REQUIRE( error.lpNorm<Eigen::Infinity>() < tolerance );
    }
  }

}