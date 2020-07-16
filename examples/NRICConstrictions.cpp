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
  std::string File = "../../data/sphere/unitSphere3.ply";

  std::string outputPrefix; // to specify potential output folder

  // Deformation model
  RealType bendingWeight = 0.001;

  // Reconstruction NRIC -> vertex position
  int gaussNewtonIterations = 100;
  int initVertexIndex = 0; // Fixed initial vertex for frame-based reconstruction
  int initFrameIndex = -1; // Fixed initial triangle for frame-based reconstruction, -1 = automatic choice

  // Constriction
  RealType constrictionFactor = 0.8;
  // List of edges to constrict, each edge is given by the two indices of its vertices
  // This format ist for example used by OpenFlipper
  std::vector<int> constrictionList{ 0, 643, 642, 0, 643, 163, 162, 642, 163, 646, 646, 43, 42, 650, 650, 162, 43, 652,
                                     652, 166, 166, 655, 655, 13, 12, 665, 665, 170, 170, 667, 667, 42, 13, 673, 673,
                                     172, 172, 676, 676, 46, 46, 682, 682, 175, 175, 685, 685, 1, 1, 686, 686, 176, 176,
                                     688, 688, 47, 47, 691, 691, 178, 178, 694, 694, 14, 14, 703, 703, 181, 181, 706,
                                     706, 49, 49, 712, 712, 184, 184, 715, 715, 2, 2, 716, 716, 185, 185, 718, 718, 50,
                                     50, 721, 721, 187, 187, 724, 724, 12, 750, 1, 192, 750, 51, 756, 756, 192, 15, 769,
                                     769, 198, 198, 771, 771, 51, 0, 786, 786, 202, 202, 788, 788, 54, 54, 791, 791,
                                     204, 204, 794, 794, 16, 16, 803, 803, 207, 207, 806, 806, 56, 56, 812, 812, 210,
                                     210, 815, 815, 5, 5, 816, 816, 211, 211, 818, 818, 57, 57, 821, 821, 213, 213, 824,
                                     824, 15, 850, 0, 218, 850, 58, 856, 856, 218, 17, 869, 869, 224, 224, 871, 871, 58,
                                     2, 886, 886, 228, 228, 888, 888, 61, 61, 891, 891, 230, 230, 894, 894, 18, 18, 903,
                                     903, 233, 233, 906, 906, 63, 63, 912, 912, 236, 236, 915, 915, 3, 3, 916, 916, 237,
                                     237, 918, 918, 64, 64, 921, 921, 239, 239, 924, 924, 17, 950, 0, 244, 950, 65, 956,
                                     956, 244, 19, 969, 969, 250, 250, 971, 971, 65, 3, 986, 986, 254, 254, 988, 988,
                                     68, 68, 991, 991, 256, 256, 994, 994, 20, 20, 1003, 1003, 259, 259, 1006, 1006, 70,
                                     70, 1012, 1012, 262, 262, 1015, 1015, 4, 4, 1016, 1016, 263, 263, 1018, 1018, 71,
                                     71, 1021, 1021, 265, 265, 1024, 1024, 19, 4, 1082, 1082, 278, 278, 1084, 1084, 74,
                                     74, 1087, 1087, 280, 280, 1090, 1090, 21, 21, 1099, 1099, 283, 283, 1102, 1102, 76,
                                     76, 1108, 1108, 286, 286, 1111, 1111, 5, 1142, 1, 292, 1142, 77, 1148, 1148, 292,
                                     22, 1161, 1161, 298, 298, 1163, 1163, 77, 5, 1178, 1178, 302, 302, 1180, 1180, 80,
                                     80, 1183, 1183, 304, 304, 1186, 1186, 23, 23, 1195, 1195, 307, 307, 1198, 1198, 82,
                                     82, 1204, 1204, 310, 310, 1207, 1207, 10, 10, 1208, 1208, 311, 311, 1210, 1210, 83,
                                     83, 1213, 1213, 313, 313, 1216, 1216, 22, 1242, 2, 318, 1242, 84, 1248, 1248, 318,
                                     24, 1261, 1261, 324, 324, 1263, 1263, 84, 1, 1278, 1278, 328, 328, 1280, 1280, 87,
                                     87, 1283, 1283, 330, 330, 1286, 1286, 25, 25, 1295, 1295, 333, 333, 1298, 1298, 89,
                                     89, 1304, 1304, 336, 336, 1307, 1307, 6, 6, 1308, 1308, 337, 337, 1310, 1310, 90,
                                     90, 1313, 1313, 339, 339, 1316, 1316, 24, 1342, 3, 344, 1342, 91, 1348, 1348, 344,
                                     26, 1361, 1361, 350, 350, 1363, 1363, 91, 2, 1378, 1378, 354, 354, 1380, 1380, 94,
                                     94, 1383, 1383, 356, 356, 1386, 1386, 27, 27, 1395, 1395, 359, 359, 1398, 1398, 96,
                                     96, 1404, 1404, 362, 362, 1407, 1407, 7, 7, 1408, 1408, 363, 363, 1410, 1410, 97,
                                     97, 1413, 1413, 365, 365, 1416, 1416, 26, 1442, 4, 370, 1442, 98, 1448, 1448, 370,
                                     28, 1461, 1461, 376, 376, 1463, 1463, 98, 3, 1478, 1478, 380, 380, 1480, 1480, 101,
                                     101, 1483, 1483, 382, 382, 1486, 1486, 29, 29, 1495, 1495, 385, 385, 1498, 1498,
                                     103, 103, 1504, 1504, 388, 388, 1507, 1507, 8, 8, 1508, 1508, 389, 389, 1510, 1510,
                                     104, 104, 1513, 1513, 391, 391, 1516, 1516, 28, 1542, 5, 396, 1542, 105, 1548,
                                     1548, 396, 30, 1561, 1561, 402, 402, 1563, 1563, 105, 4, 1578, 1578, 406, 406,
                                     1580, 1580, 108, 108, 1583, 1583, 408, 408, 1586, 1586, 31, 31, 1595, 1595, 411,
                                     411, 1598, 1598, 110, 110, 1604, 1604, 414, 414, 1607, 1607, 9, 9, 1608, 1608, 415,
                                     415, 1610, 1610, 111, 111, 1613, 1613, 417, 417, 1616, 1616, 30, 10, 1674, 1674,
                                     430, 430, 1676, 1676, 114, 114, 1679, 1679, 432, 432, 1682, 1682, 32, 32, 1691,
                                     1691, 435, 435, 1694, 1694, 116, 116, 1700, 1700, 438, 438, 1703, 1703, 6, 6, 1766,
                                     1766, 452, 452, 1768, 1768, 119, 119, 1771, 1771, 454, 454, 1774, 1774, 33, 33,
                                     1783, 1783, 457, 457, 1786, 1786, 121, 121, 1792, 1792, 460, 460, 1795, 1795, 7, 7,
                                     1858, 1858, 474, 474, 1860, 1860, 124, 124, 1863, 1863, 476, 476, 1866, 1866, 34,
                                     34, 1875, 1875, 479, 479, 1878, 1878, 126, 126, 1884, 1884, 482, 482, 1887, 1887,
                                     8, 8, 1950, 1950, 496, 496, 1952, 1952, 129, 129, 1955, 1955, 498, 498, 1958, 1958,
                                     35, 35, 1967, 1967, 501, 501, 1970, 1970, 131, 131, 1976, 1976, 504, 504, 1979,
                                     1979, 9, 9, 2042, 2042, 518, 518, 2044, 2044, 134, 134, 2047, 2047, 520, 520, 2050,
                                     2050, 36, 36, 2059, 2059, 523, 523, 2062, 2062, 136, 136, 2068, 2068, 526, 526,
                                     2071, 2071, 10, 2102, 6, 532, 2102, 137, 2108, 2108, 532, 37, 2121, 2121, 538, 538,
                                     2123, 2123, 137, 10, 2138, 2138, 542, 542, 2140, 2140, 140, 140, 2143, 2143, 544,
                                     544, 2146, 2146, 38, 38, 2155, 2155, 547, 547, 2158, 2158, 142, 142, 2164, 2164,
                                     550, 550, 2167, 2167, 11, 11, 2168, 2168, 551, 551, 2170, 2170, 143, 143, 2173,
                                     2173, 553, 553, 2176, 2176, 37, 2202, 7, 558, 2202, 144, 2208, 2208, 558, 39, 2221,
                                     2221, 564, 564, 2223, 2223, 144, 11, 2260, 2260, 573, 573, 2262, 2262, 148, 148,
                                     2265, 2265, 575, 575, 2268, 2268, 39, 2294, 8, 580, 2294, 149, 2300, 2300, 580, 40,
                                     2313, 2313, 586, 586, 2315, 2315, 149, 11, 2352, 2352, 595, 595, 2354, 2354, 153,
                                     153, 2357, 2357, 597, 597, 2360, 2360, 40, 2386, 9, 602, 2386, 154, 2392, 2392,
                                     602, 41, 2405, 2405, 608, 608, 2407, 2407, 154, 11, 2444, 2444, 617, 617, 2446,
                                     2446, 158, 158, 2449, 2449, 619, 619, 2452, 2452, 41 };


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
  std::vector<int> strangulatedEdges;
  for ( int i = 0; i < constrictionList.size() / 2; i++ ) {
    int v_a = constrictionList[2 * i];
    int v_b = constrictionList[2 * i + 1];

    for ( int edgeIdx = 0; edgeIdx < numEdges; edgeIdx++ ) {
      int e1 = Topology.getAdjacentNodeOfEdge( edgeIdx, 0 );
      int e2 = Topology.getAdjacentNodeOfEdge( edgeIdx, 1 );

      if (( e1 == v_a && e2 == v_b ) || ( e1 == v_b && e2 == v_a ))
        strangulatedEdges.push_back( edgeIdx );
    }
  }

  std::cout << " - Number of strangulated edges: " << strangulatedEdges.size() << std::endl;

  VectorType constrictedGeometryNRIC = GeometryNRIC;
  for ( const int edgeIdx : strangulatedEdges )
    constrictedGeometryNRIC[edgeIdx] *= constrictionFactor;

  // Define functionals
  VectorType Alphas( 1 );
  Alphas.setConstant( 1. );

  // Abuse elastic mean functionals to get energy to a fixed mesh as simple objective
  ElasticMeanFunctional<DefaultConfigurator> EMF( W, GeometryNRIC, Alphas, 1 );
  ElasticMeanFunctionalGradient<DefaultConfigurator> EMG( W, GeometryNRIC, Alphas, 1 );
  ElasticMeanFunctionalHessian<DefaultConfigurator> EMH( W, GeometryNRIC, Alphas, 1 );


  AugmentedLagrangeMethod<DefaultConfigurator> constrainedSolverAL( EMF, EMG, EMH, intOp, intGrad, intHess,
                                                                    constrainedIterations, 1e-8, false );
  constrainedSolverAL.setParameter( "penalty_factor", 100. );
  constrainedSolverAL.setParameter( "constraint_decrease_exponent", 0.9 );
  constrainedSolverAL.setParameter( "constraint_tolerance", 5e-10 );
  constrainedSolverAL.setParameter( "optimality_tolerance", 1e-8 );
  constrainedSolverAL.setParameter( "inner_method", "LSN" );
  constrainedSolverAL.setParameter( "inner__tau_increase", 10. );
  constrainedSolverAL.setParameter( "inner__direction_beta", 1.e-3 );
  constrainedSolverAL.setParameter( "inner__minimal_stepsize", 1.e-15 );

  constrainedSolverAL.setBoundaryMask( strangulatedEdges );

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
  if ( !OpenMesh::IO::write_mesh( output, outputPrefix + "constricted_mesh.ply" ))
    throw std::runtime_error( "Failed to write file: " + outputPrefix + "constricted_mesh.ply" );
  std::cout << "Done." << std::endl;

}