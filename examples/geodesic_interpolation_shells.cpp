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
 * For two input shapes \f$ s_A \f$ and \f$ s_B \f$ in dense correspondence and a prescribed length \f$ K+1 \f$,
 * the discrete geodesic \f$ (s_0, \ldots, s_K) \f$ connecting \f$ s_0 = s_A \f$ and \f$ s_K = s_B \f$,
 * is defined as the minimizer of the discrete path energy \f$ E[s_1, \ldots, s_{K-1}] = K \sum_{k=1}^K W[s_{k-1}, s_k] \f$. 
 * 
 * Here, the deformation energy W is given as an elastic shell energy (membrane plus bending term, controlled by bending weight parameter).
 * 
 * \cite HeRuWa12 and \cite HeRuSc14
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
  
  if (!OpenMesh::IO::read_mesh(mesh2, "../../data/finger/finger1.ply")) 
  {
    std::cerr << "read error\n";
    exit(1);
  }

  double bendWeight = 0.001;
  int length = 5;
  int initializationScheme = 2; // 2 = linear interpolation of nodal values
  OptimizationParameters<DefaultConfigurator> optPars;
  optPars.setGradientIterations( 100 );
  optPars.setBFGSIterations( 75 );
  optPars.setNewtonIterations( 35 );
  optPars.setQuietMode( SHOW_ALL );
  
  MeshTopologySaver Topol( mesh1 );
  typename DefaultConfigurator::VectorType startGeom, endGeom;
  getGeometry( mesh1, startGeom );
  getGeometry( mesh2, endGeom );
  
  ShellDeformationType W(Topol, bendWeight );
  GeodesicInterpolation<DefaultConfigurator> interpolationOp( Topol, startGeom, endGeom, W, length, optPars, false );
  
  
  std::vector<int> indices;
  Topol.getVertexNRingVertices( 10, 1302,  indices );
  std::cerr << "#bdry nodes = " << indices.size() << std::endl;
  extendBoundaryMask( Topol.getNumVertices(), indices );    
  interpolationOp.setBoundaryMask( indices );
  
  auto t_start = std::chrono::high_resolution_clock::now();
  interpolationOp.execute( "fingerInterpolation", initializationScheme );
  auto t_end = std::chrono::high_resolution_clock::now();
  std::cout << std::fixed << "Done in " << std::chrono::duration<double, std::ratio<1> >(t_end - t_start).count() << " seconds." << std::endl;

  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}