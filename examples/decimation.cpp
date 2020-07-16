// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include <iostream>
#include <string>

#include <goast/Core.h>
#include <goast/Multiresolution.h>


/**
 * \brief Example to compute simultaneous mesh decimation
 * \author Heeren
 * 
 * Load a sequence of meshes in dense correspondence,
 * specify percentage \f$ \theta > 0 \f$ of vertices as decimation target.
 * Then the algorithm performs a simultaneous mesh decimation by edge collapses
 * while keeping the dense correspondence.
 * 
 * By default, edge collapses are performed according to the quadric error metric
 * introduced in \cite GaHe97
 * 
 */
int main(int argc, char *argv[])
{

try{

  std::cerr << "=================================================================================" << std::endl;
  std::cerr << "EXAMPLE SIMULTANEOUS MESH DECIMATION" << std::endl;
  std::cerr << "=================================================================================" << std::endl << std::endl;
  
  // load meshes and copy to vector
  std::vector<TriMesh> meshes;
  int numMeshes = 6;
  for( int i = 0; i < numMeshes; i++ ){   
    TriMesh mesh;  
    std::ostringstream filename;
    filename << "../../data/cactus/cactus" << i << ".ply";
    OpenMesh::IO::read_mesh(mesh, filename.str() );  
    meshes.push_back(mesh);
  }
  std::cerr << "Original vertex number = " << meshes[0].n_vertices() << std::endl;
  
  // specify Theta (i.e. number of nodes n is reduced to Theta * n)
  double Theta = 0.15;
  // perform decimation by edge collapses
  std::cerr << "Simulatineous decimation of " << numMeshes << " to " << 100. * Theta << " percent of vertices." << std::endl;
  SimultaneousDecimater<DefaultConfigurator> simDec( meshes );  
  simDec.calcDecimation( Theta );
  std::cerr << "Decimated vertex number = " << simDec.getCoarseMesh( 0 ).n_vertices() << std::endl;       
  
  // saving
  for( int i = 0; i < numMeshes; i++ ){ 
    TriMesh coarseMesh = simDec.getCoarseMesh( i );
    std::ostringstream filename;
    filename << "cactus_15pc_" << i << ".ply";
    OpenMesh::IO::write_mesh(coarseMesh, filename.str() );
  }
  
  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}