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

//==============================================================================================================
typedef typename DefaultConfigurator::VectorType VectorType;

/** 
 * \brief Register a template shape to a given reference shape
 * \author Heeren
 * 
 */
int main(int argc, char *argv[])
{

try{

    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "EXAMPLE REGISTRATION" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;
    
    TriMesh referenceMesh, templateMesh;
    OpenMesh::IO::read_mesh(referenceMesh, "../../data/finger/finger0.ply");   
    OpenMesh::IO::read_mesh(templateMesh, "../../data/finger/finger1.ply"); 
    std::cerr << "Num of vertices = " << referenceMesh.n_vertices() << std::endl;
        
    VectorType referenceGeometry, templateGeometry;
    getGeometry( referenceMesh, referenceGeometry ); 
    getGeometry( templateMesh, templateGeometry ); 
    MeshTopologySaver Topol( referenceMesh );
    OpenMesh::IO::write_mesh( referenceMesh, "registeredFinger_reference.ply"); 
        
    std::cerr << "Perform momentum registration..." << std::endl;
    VectorType OptimalRBMdofs;
    MomentumRegistrationOperator<DefaultConfigurator>( Topol, referenceGeometry ).execute( templateGeometry, OptimalRBMdofs );
    setGeometry( templateMesh, templateGeometry ); 
    OpenMesh::IO::write_mesh(templateMesh, "registeredFinger_momentum.ply");  
    
    std::cerr << "Perform full least squares registration..." << std::endl;
    LeastSquaresRBMRegistrationOp<DefaultConfigurator>( referenceGeometry ).execute( templateGeometry, OptimalRBMdofs );
    setGeometry( templateMesh, templateGeometry ); 
    OpenMesh::IO::write_mesh(templateMesh, "registeredFinger_leastSquares.ply");  

    // least sqaures regristration is limited to specified vertex set
    std::cerr << "Perform partial least squares registration..." << std::endl;
    std::vector<int> indices;
    Topol.getVertexNRingVertices( 15, 1302,  indices );
    std::cerr << "#fixed nodes = " << indices.size() << std::endl;
    LeastSquaresRBMRegistrationOp<DefaultConfigurator>( referenceGeometry, indices, 10000, 1e-8 ).execute( templateGeometry, OptimalRBMdofs );
    setGeometry( templateMesh, templateGeometry ); 
    OpenMesh::IO::write_mesh(templateMesh, "registeredFinger_leastSquares_partial.ply"); 
  
  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}