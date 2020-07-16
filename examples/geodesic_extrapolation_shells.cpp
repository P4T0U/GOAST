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

typedef ShellDeformation<DefaultConfigurator, NonlinearMembraneDeformation<DefaultConfigurator>, SimpleBendingDeformation<DefaultConfigurator> > ShellDeformationType;

/** 
 * \brief Computation of discrete geodesic interpolation in shell space
 * \author Heeren
 * 
 * For two input shapes \f$ s_0 \f$ and \f$ s_1 \f$ compute sequence 
 * \f$ s_0, s_1, \ldots, s_K \f$, such that \f$ (s_0, \ldots, s_K) \f$ is 
 * a discrete geodesic of length \f$ K+1 \f$.
 * 
 * In detail, the sequence is constructed iteratively, by computing 3-point geodesics 
 * \f$ (s_{k-2}, s_{k-1}, s_k) \f$ for given \f$ s_{k-2} \f$ and \f$ s_{k-1} \f$ 
 * for \f$ k = 2,3, \ldots, K \f$.
 * 
 * \cite HeRuSc14
 * 
 */
int main(int argc, char *argv[])
{

try{

  TriMesh mesh1, mesh2;
  
  
  if (!OpenMesh::IO::read_mesh(mesh1, "../../data/finger/finger0.ply"))     
  {
    std::cerr << "read error\n";
    exit(1);
  }
  
  if (!OpenMesh::IO::read_mesh(mesh2, "../../data/finger/finger_var.ply")) 
  {
    std::cerr << "read error\n";
    exit(1);
  }

  double bendWeight = 0.001;
  int numShootingSteps = 3;
  OptimizationParameters<DefaultConfigurator> optPars;
  optPars.setNewtonIterations( 50 );
  optPars.setQuietMode( SHOW_ALL );
  bool forwardShooting = true;
 
  MeshTopologySaver Topol( mesh1 );
  typename DefaultConfigurator::VectorType positionGeom, variationGeom;
  getGeometry( mesh1, positionGeom );
  getGeometry( mesh2, variationGeom );
  
  ShellDeformationType W(Topol, bendWeight );
  GeodesicExtrapolation<DefaultConfigurator> extrapolationOp( Topol, positionGeom, variationGeom, W, optPars, false );
  
  // boundary mask
  std::vector<int> indices;
  Topol.getVertexNRingVertices( 10, 1302,  indices );
  std::cerr << "#bdry nodes = " << indices.size() << std::endl;
  extendBoundaryMask( Topol.getNumVertices(), indices );    
  extrapolationOp.setBoundaryMask( indices );
  
  
  auto t_start = std::chrono::high_resolution_clock::now();
  extrapolationOp.execute( numShootingSteps, forwardShooting, "fingerExtrapolation" );
  auto t_end = std::chrono::high_resolution_clock::now();
  std::cout << std::fixed << "Done in " << std::chrono::duration<double, std::ratio<1> >(t_end - t_start).count() << " seconds." << std::endl;

  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}