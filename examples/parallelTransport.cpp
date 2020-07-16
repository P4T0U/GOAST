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

//==============================================================================================================
typedef ShellDeformation<DefaultConfigurator, NonlinearMembraneDeformation<DefaultConfigurator>, SimpleBendingDeformation<DefaultConfigurator> > ShellDeformationType;


/** 
 * \brief Discrete parallel transport via Schild's ladder in shell space 
 * \author Heeren
 * 
 * Corresponds to Fig. 1 in \cite HeRuSc14
 * 
 */
int main(int argc, char *argv[]){

  try{

    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "EXAMPLE PARALLEL TRANSPORT" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;
    
    //! first compute discrete path along which the variation is transported
    
    // load two end meshes
    TriMesh startMesh, endMesh;
    OpenMesh::IO::read_mesh(startMesh, "../../data/faces/neutral.ply"); 
    OpenMesh::IO::read_mesh(endMesh, "../../data/faces/smile.ply");
    
    //! get geometries
    typename DefaultConfigurator::VectorType startGeom, endGeom;
    getGeometry( startMesh, startGeom );
    getGeometry( endMesh, endGeom );
    
    MeshTopologySaver Topol( startMesh );
    typename DefaultConfigurator::RealType bendingWeight = 1.0;
    ShellDeformationType W(Topol, bendingWeight );
    int numLocalDofs = 3 * Topol.getNumVertices();
    
    
    int lengthOfPath = 5;
    int initializationScheme = 2; // 2 = linear interpolation of nodal values  
    OptimizationParameters<DefaultConfigurator> optPars;
    optPars.setGradientIterations( 50 );
    optPars.setQuietMode( SHOW_TERMINATION_INFO );    
    
    // boundary nodes
    std::vector<int> indices;  
    Topol.fillPartialBoundaryMask( indices, 6137 );
    std::cerr << "#bdry nodes = " << indices.size() << std::endl;
    extendBoundaryMask( Topol.getNumVertices(), indices );    
    
    // set up interpolation operator
    GeodesicInterpolation<DefaultConfigurator> interpolationOp( Topol, startGeom, endGeom, W, lengthOfPath, optPars, false );
    interpolationOp.setBoundaryMask( indices );
    
    auto t_start = std::chrono::high_resolution_clock::now();
    typename DefaultConfigurator::VectorType fullPath;
    interpolationOp.execute( fullPath, initializationScheme, "smoothedPathFaces" );  
    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout << std::fixed << "Done in " << std::chrono::duration<double, std::ratio<1> >(t_end - t_start).count() << " seconds." << std::endl;
    
    //! convert into more convenient form
    std::vector<typename DefaultConfigurator::VectorType> shapesAlongPath;
    shapesAlongPath.push_back( startGeom );
    for(int i = 0; i < lengthOfPath - 2; i++ )
        shapesAlongPath.push_back(  fullPath.segment( i*numLocalDofs, numLocalDofs ) );    
    shapesAlongPath.push_back( endGeom );
    
    //! now start actual parallel transport
    
    // load variation to be transported
    TriMesh varMesh;
    OpenMesh::IO::read_mesh(varMesh, "../../data/faces/disgust.ply"); 
    typename DefaultConfigurator::VectorType initialVarGeom, finalVarGeom;
    getGeometry( varMesh, initialVarGeom );
    
    optPars.setGradientIterations( 250 );
    optPars.setBFGSIterations( 0 );
    optPars.setNewtonIterations( 25 );
    optPars.setInitializationScheme( 2 );
    
    // Schild' ladder iteration
    for( int i = 0; i < lengthOfPath - 1; i ++ ){
        std::cerr << "\n\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=" << std::endl;
        std::cerr << "Step " << i+1 << " of " << lengthOfPath - 1 << " of Schild's ladder" << std::endl;
        std::cerr << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=" << std::endl;
        parTranspSingleStepShort<DefaultConfigurator>( Topol, shapesAlongPath[i], shapesAlongPath[i+1], initialVarGeom, indices, W, optPars, finalVarGeom, false );
        initialVarGeom = finalVarGeom;
    }
    
    // save transported variation
    setGeometry( startMesh, initialVarGeom );
    OpenMesh::IO::write_mesh(startMesh, "parTransportedDisgustFace.ply");
/*    
    // Schild' ladder iteration with deformation transfer
    getGeometry( varMesh, initialVarGeom );
    for( int i = 0; i < lengthOfPath - 1; i ++ ){
        std::cerr << "\n\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=" << std::endl;
        std::cerr << "Step " << i+1 << " of " << lengthOfPath - 1 << " of Schild's ladder with deformation transfer" << std::endl;
        std::cerr << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=" << std::endl;
        std::vector<int> mask;  
        Topol.fillPartialBoundaryMask( mask, 6137 );
        parTranspDeformationTransfer<DefaultConfigurator>( Topol, shapesAlongPath[i], shapesAlongPath[i+1], initialVarGeom, mask, finalVarGeom, false );
        initialVarGeom = finalVarGeom;
    }
    
    // save transported variation
    setGeometry( startMesh, initialVarGeom );
    OpenMesh::IO::write_mesh(startMesh, "parTransportedDisgustFace_deformationTransfer.ply");
*/    
  } 
  catch ( BasicException &el ){
        std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
  }
  
  return 0;
}